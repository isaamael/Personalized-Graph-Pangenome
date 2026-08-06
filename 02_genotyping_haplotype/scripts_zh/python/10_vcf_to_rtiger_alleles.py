#!/usr/bin/env python3
"""Convert VCF allele depths to RTIGER counts: REF to MM and ALT to TS."""
from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
from collections import Counter

csv.field_size_limit(min(sys.maxsize, 2**31 - 1))

try:
    import numpy as np
    from cyvcf2 import VCF
except ImportError:
    sys.exit("cyvcf2+numpy required; conda activate genome")


class PrepError(Exception):
    pass


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--tag", default="default")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--out-expdesign", required=True)
    ap.add_argument("--out-sites", default="", help="Optional marker-polarity table")
    ap.add_argument("--tmpdir", default="")
    ap.add_argument("--sort-parallel", type=int, default=4)
    return ap.parse_args()


def resolve_tmpdir(user_tmp: str, outdir: str) -> str:
    if user_tmp:
        return user_tmp
    return os.environ.get("RTIGER_TMPDIR", "") or outdir


def vcf_has_ad(header: str) -> bool:
    return any("##FORMAT" in ln and "ID=AD" in ln for ln in header.split("\n"))


def is_biallelic_marker(var) -> bool:
    alts = var.ALT or []
    return len(alts) == 1 and bool(alts[0])


def parse_ad(var, sample_i: int) -> tuple[int, int]:
    try:
        ad = var.format("AD")
    except (KeyError, ValueError, TypeError):
        return 0, 0
    if ad is None:
        return 0, 0
    row = ad[sample_i]
    if row is None or len(row) < 2:
        return 0, 0
    ref_ad = int(row[0]) if row[0] is not None and row[0] >= 0 else 0
    alt_ad = int(row[1]) if row[1] is not None and row[1] >= 0 else 0
    return ref_ad, alt_ad


def gt_called(gts: np.ndarray, i: int) -> bool:
    return gts[i, 0] >= 0 and gts[i, 1] >= 0


def write_variant_rows(
    w,
    var,
    gts: np.ndarray,
    samples: list[str],
    counters: Counter,
) -> tuple[int, int]:
    n_written = 0
    n_skip = 0
    chrom = var.CHROM
    pos = str(var.POS)
    ref_a = var.REF
    alt_a = (var.ALT or [""])[0]

    for i, sample in enumerate(samples):
        if not gt_called(gts, i):
            n_skip += 1
            counters["skip_gt"] += 1
            continue
        ref_c, alt_f = parse_ad(var, i)
        w.writerow([chrom, pos, ref_a, ref_c, alt_a, alt_f, sample])
        n_written += 1
        counters["sum_depth"] += ref_c + alt_f
        if ref_c + alt_f <= 0:
            counters["zero_total_rows"] += 1
    return n_written, n_skip


def split_long_table(long_path: str, outdir: str, samples: list[str]) -> dict[str, Counter]:
    os.makedirs(outdir, exist_ok=True)
    stats: dict[str, Counter] = {s: Counter() for s in samples}
    cur_sample = None
    cur_fh = None
    try:
        with open(long_path, newline="") as f:
            for row in csv.reader(f, delimiter="\t"):
                chrom, pos, mm_a, mm_c, ts_a, ts_c, sample = row
                if sample != cur_sample:
                    if cur_fh is not None:
                        cur_fh.close()
                    cur_sample = sample
                    p = os.path.join(outdir, f"{sample}.allele_count.tsv")
                    cur_fh = open(p, "w")
                cur_fh.write(f"{chrom}\t{pos}\t{mm_a}\t{mm_c}\t{ts_a}\t{ts_c}\n")
                stats[sample]["n_markers"] += 1
                stats[sample]["sum_depth"] += int(mm_c) + int(ts_c)
    finally:
        if cur_fh is not None:
            cur_fh.close()
    return stats


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(os.path.dirname(args.out_expdesign) or ".", exist_ok=True)
    tmp = resolve_tmpdir(args.tmpdir, args.outdir)
    os.makedirs(tmp, exist_ok=True)
    long_path = os.path.join(tmp, f"{args.dataset}.allele.long.tsv")
    counters: Counter = Counter()
    site_rows: list[list[str]] = []

    vcf = VCF(args.vcf)
    if not vcf_has_ad(vcf.raw_header):
        vcf.close()
        raise PrepError("VCF missing FORMAT:AD — cannot build RTIGER allele-counts from read depths")

    samples = list(vcf.samples)
    n_written = 0
    n_skip_gt = 0
    n_sites = 0

    with open(long_path, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        for var in vcf:
            if not is_biallelic_marker(var):
                counters["skip_multiallelic"] += 1
                continue
            n_sites += 1
            if args.out_sites:
                site_rows.append([var.CHROM, str(var.POS), var.REF, (var.ALT or [""])[0]])
            gts = var.genotype.array()[:, :2]
            nw, ns = write_variant_rows(w, var, gts, samples, counters)
            n_written += nw
            n_skip_gt += ns
    vcf.close()

    if n_written == 0:
        raise PrepError("no allele-count records written; check thinned VCF and AD fields")

    if args.out_sites and site_rows:
        os.makedirs(os.path.dirname(args.out_sites) or ".", exist_ok=True)
        with open(args.out_sites, "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["chrom", "pos", "vcf_ref", "vcf_alt", "refc_is_mm", "altf_is_ts"])
            for row in site_rows:
                w.writerow(row + ["1", "1"])

    sorted_path = long_path + ".sorted"
    sort_env = os.environ.copy()
    sort_env["LC_ALL"] = "C"
    subprocess.run([
        "sort", "-t", "\t", "-k7,7", "-k1,1V", "-k2,2n",
        f"--parallel={max(1, args.sort_parallel)}", "-T", tmp,
        long_path, "-o", sorted_path,
    ], check=True, env=sort_env)
    stats = split_long_table(sorted_path, args.outdir, samples)

    with open(args.out_expdesign, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["files", "name"])
        for s in sorted(samples):
            w.writerow([os.path.join(args.outdir, f"{s}.allele_count.tsv"), s])

    for p in (long_path, sorted_path):
        try:
            os.remove(p)
        except OSError:
            pass

    print(
        f"sites={n_sites} samples={len(samples)} records={n_written} skip_gt={n_skip_gt} "
        f"skip_zero_total={counters['zero_total_rows']} skip_multiallelic={counters['skip_multiallelic']}"
    )
    print(f"expDesign={args.out_expdesign}")


if __name__ == "__main__":
    try:
        main()
    except PrepError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)
