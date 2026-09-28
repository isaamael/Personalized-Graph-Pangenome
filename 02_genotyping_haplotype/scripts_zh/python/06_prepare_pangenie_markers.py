#!/usr/bin/env python3
"""筛选在标记面板中表现为亲本间固定差异的共线区 SyRI SNP。"""
from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import os
import sys
from typing import Dict, List, Set, Tuple


def load_collinear(bed: str) -> Dict[str, Tuple[List[int], List[int]]]:
    """将以 0 为起点的左闭右开区间按染色体整理为 (starts, ends)，并按起点排序。"""
    raw: Dict[str, List[Tuple[int, int]]] = {}
    with open(bed) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#") or line.startswith("track"):
                continue
            c, s, e = line.split()[:3]
            raw.setdefault(c, []).append((int(s), int(e)))
    out: Dict[str, Tuple[List[int], List[int]]] = {}
    for c, lst in raw.items():
        lst.sort()
        out[c] = ([a[0] for a in lst], [a[1] for a in lst])
    return out


def in_collinear(chr_iv: Tuple[List[int], List[int]], pos1: int) -> bool:
    starts, ends = chr_iv
    if not starts:
        return False
    x0 = pos1 - 1
    i = bisect.bisect_right(starts, x0) - 1
    if i < 0:
        return False
    return starts[i] <= x0 < ends[i]


def open_any(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--syri-snp", required=True)
    ap.add_argument("--collinear-bed", required=True)
    ap.add_argument("--panel-vcf", required=True)
    ap.add_argument("--out-markers", required=True)
    ap.add_argument("--out-targets", required=True, help="CHROM POS for bcftools -T")
    ap.add_argument("--out-stats", required=True)
    args = ap.parse_args()

    print("[step06] load collinear …", flush=True)
    col = load_collinear(args.collinear_bed)
    n_iv = sum(len(v[0]) for v in col.values())

    # SyRI 中 ref_allele=TS（第 3 列）、alt_allele=MM（第 4 列）；标记面板中 REF=MM、ALT=TS
    print("[step06] load SyRI SNP (full ~1.62M, not thin) …", flush=True)
    syri_keys: Dict[Tuple[str, int, str, str], str] = {}
    n_syri = n_syri_col = 0
    with open(args.syri_snp) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            n_syri += 1
            chrom = row["mm_chr"]
            pos = int(row["mm_pos"])
            chr_iv = col.get(chrom)
            if chr_iv is None or not in_collinear(chr_iv, pos):
                continue
            n_syri_col += 1
            ref_mm = row["alt_allele"]  # MM
            alt_ts = row["ref_allele"]  # TS
            if len(ref_mm) != 1 or len(alt_ts) != 1:
                continue
            syri_keys[(chrom, pos, ref_mm, alt_ts)] = row["id"]

    print(f"[step06] syri_all={n_syri} syri_collinear_snp={n_syri_col} keyed={len(syri_keys)}", flush=True)

    print("[step06] scan panel for MM=0|0 TS=1|1 biallelic SNP …", flush=True)
    n_panel = n_bi_snp = n_diag = n_hit = 0
    rows_out = []
    with open_any(args.panel_vcf) as fh:
        samples = []
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                try:
                    i_mm = samples.index("MM")
                    i_ts = samples.index("TS")
                except ValueError as e:
                    raise SystemExit(f"panel missing MM/TS samples: {samples[:5]}") from e
                continue
            n_panel += 1
            f = line.rstrip("\n").split("\t")
            chrom, pos_s, _id, ref, alt = f[0], f[1], f[2], f[3], f[4]
            if "," in alt or len(ref) != 1 or len(alt) != 1:
                continue
            n_bi_snp += 1
            fmt = f[8].split(":")
            try:
                gi = fmt.index("GT")
            except ValueError:
                gi = 0

            def gt_tok(cell: str) -> str:
                return cell.split(":")[gi]

            gt_mm = gt_tok(f[9 + i_mm])
            gt_ts = gt_tok(f[9 + i_ts])
            if gt_mm.replace("/", "|") not in ("0|0",) or gt_ts.replace("/", "|") not in ("1|1",):
                # 同时接受未定相基因型 0/0 和 1/1
                if gt_mm.replace("|", "/") != "0/0" or gt_ts.replace("|", "/") != "1/1":
                    continue
            n_diag += 1
            pos = int(pos_s)
            key = (chrom, pos, ref, alt)
            sid = syri_keys.get(key)
            if sid is None:
                continue
            n_hit += 1
            rows_out.append((chrom, pos, ref, alt, sid))

    rows_out.sort(key=lambda x: (x[0], x[1]))
    os.makedirs(os.path.dirname(args.out_markers) or ".", exist_ok=True)
    with open(args.out_markers, "w") as fo:
        fo.write("marker_id\tchrom\tpos\tref\talt\tsyri_id\n")
        for i, (chrom, pos, ref, alt, sid) in enumerate(rows_out, 1):
            fo.write(f"m{i}\t{chrom}\t{pos}\t{ref}\t{alt}\t{sid}\n")
    with open(args.out_targets, "w") as fo:
        for chrom, pos, *_ in rows_out:
            fo.write(f"{chrom}\t{pos}\n")

    with open(args.out_stats, "w") as fo:
        fo.write("metric\tcount\tdenom\trate\tnote\n")
        fo.write(f"syri_snp_all\t{n_syri}\t{n_syri}\t1.0000\tfull SyRI not thin markers_sim\n")
        fo.write(f"collinear_intervals\t{n_iv}\t{n_iv}\t1.0000\tMM collinear.bed\n")
        fo.write(f"syri_snp_collinear\t{n_syri_col}\t{n_syri}\t{n_syri_col/n_syri:.4f}\t\n")
        fo.write(f"panel_records\t{n_panel}\t{n_panel}\t1.0000\t\n")
        fo.write(f"panel_biallelic_snp\t{n_bi_snp}\t{n_panel}\t{n_bi_snp/max(n_panel,1):.4f}\t\n")
        fo.write(f"panel_MM00_TS11\t{n_diag}\t{n_bi_snp}\t{n_diag/max(n_bi_snp,1):.4f}\t\n")
        fo.write(
            f"ancestry_markers\t{n_hit}\t{n_syri_col}\t{n_hit/max(n_syri_col,1):.4f}\t"
            f"exact CHROM+POS+REF+ALT (REF=MM ALT=TS)\n"
        )

    print(f"[DONE] markers={n_hit} → {args.out_markers}", flush=True)
    if n_hit == 0:
        print("[FAIL] zero ancestry markers", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
