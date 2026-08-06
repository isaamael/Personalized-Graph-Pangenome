#!/usr/bin/env python3
"""Export Viterbi ancestry segments in GS and R/qtl-compatible formats.

States 0, 1 and 2 represent homozygous MM, heterozygous and homozygous TS.
One marker per genomic bin is selected by KC and distance to the bin centre.
"""
from __future__ import annotations

import argparse
import hashlib
import os
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from typing import Dict, List, Optional, Tuple

GT2ST = {0: "mat", 1: "het", 2: "pat"}
ST_CODE = {"NA": 0, "mat": 1, "het": 2, "pat": 3}
CODE_ST = ("NA", "mat", "het", "pat")


def chrom_sort_key(c: str):
    if c.startswith("chr") and c[3:].isdigit():
        return (0, int(c[3:]))
    return (1, c)


def sha256_first(path: str, n: int = 64 * 1024 * 1024) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        h.update(f.read(n))
    return h.hexdigest()


def hardlink_or_copy(src: str, dst: str) -> None:
    if os.path.lexists(dst):
        os.remove(dst)
    try:
        os.link(src, dst)
    except OSError:
        shutil.copy2(src, dst)


def load_samples(bcftools: str, vcf: str) -> List[str]:
    out = subprocess.check_output([bcftools, "query", "-l", vcf], text=True)
    samples = [ln.strip() for ln in out.splitlines() if ln.strip()]
    if not samples:
        raise SystemExit(f"[FAIL] no samples in {vcf}")
    return samples


def mean_kc_from_fields(fields: List[str]) -> float:
    s = 0.0
    n = 0
    for x in fields:
        if not x or x == ".":
            continue
        try:
            s += float(x)
            n += 1
        except ValueError:
            continue
    return s / n if n else -1.0


def select_gs_markers(
    bcftools: str, vcf: str, bin_bp: int
) -> Tuple[List[str], Dict[str, List[int]], List[Tuple[str, int, int, float, int]]]:
    """50kb-bin thin: max meanKC, then nearest bin center.

    Returns labels, by_chr_pos, site_rows(chrom,pos,bin,mean_kc,dist_center).
    bin_bp<=0 → keep all sites (KC still parsed for stats file).
    """
    cmd = [bcftools, "query", "-f", "%CHROM\t%POS[\t%KC]\n", vcf]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    assert proc.stdout is not None

    # best per (chrom, bin): (-kc, dist, pos, mean_kc)
    best: Dict[Tuple[str, int], Tuple[float, int, int, float]] = {}
    n_in = 0

    def consider(chrom: str, pos: int, kc: float) -> None:
        nonlocal n_in
        n_in += 1
        if bin_bp <= 0:
            bid = pos  # unique
            center = pos
        else:
            bid = (pos - 1) // bin_bp
            center = bid * bin_bp + (bin_bp + 1) // 2
        dist = abs(pos - center)
        key = (chrom, bid)
        # sort key for winner: higher kc, then smaller dist, then smaller pos
        score = (-kc, dist, pos, kc)
        cur = best.get(key)
        if cur is None or score < (-cur[0], cur[1], cur[2], cur[3]):
            # store as (kc, dist, pos, kc) with positive kc for later
            best[key] = (kc, dist, pos, kc)

    for line in proc.stdout:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 2:
            continue
        chrom, pos_s = parts[0], parts[1]
        pos = int(pos_s)
        kc = mean_kc_from_fields(parts[2:])
        consider(chrom, pos, kc)
    if proc.wait() != 0:
        raise SystemExit(f"[FAIL] bcftools query KC {vcf}")

    by_chr: Dict[str, List[int]] = defaultdict(list)
    site_rows: List[Tuple[str, int, int, float, int]] = []
    for (chrom, bid), (kc, dist, pos, _) in best.items():
        by_chr[chrom].append(pos)
        site_rows.append((chrom, pos, bid if bin_bp > 0 else 0, kc, dist))

    labels: List[str] = []
    for chrom in sorted(by_chr, key=chrom_sort_key):
        by_chr[chrom].sort()
        for p in by_chr[chrom]:
            labels.append(f"{chrom}:{p}")
    site_rows.sort(key=lambda r: (chrom_sort_key(r[0]), r[1]))
    print(f"[INFO] GS select: in={n_in} kept={len(labels)} bin_bp={bin_bp}", flush=True)
    return labels, dict(by_chr), site_rows


def load_segments(path: str) -> Dict[str, List[dict]]:
    by: Dict[str, List[dict]] = defaultdict(list)
    with open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(hdr)}
        need = ("sample", "chrom", "start", "end", "gt", "n_markers")
        for k in need:
            if k not in idx:
                raise SystemExit(f"[FAIL] segments missing col {k}: {hdr}")
        for line in f:
            p = line.rstrip("\n").split("\t")
            gt = int(p[idx["gt"]])
            if gt not in GT2ST:
                raise SystemExit(f"[FAIL] bad gt={gt} in {path}")
            by[p[idx["sample"]]].append(
                {
                    "chr": p[idx["chrom"]],
                    "start": int(p[idx["start"]]),
                    "end": int(p[idx["end"]]),
                    "gt": gt,
                    "state": GT2ST[gt],
                    "n_markers": int(p[idx["n_markers"]]),
                }
            )
    for sid in by:
        by[sid].sort(key=lambda x: (chrom_sort_key(x["chr"]), x["start"], x["end"]))
    return dict(by)


def load_transitions(path: str) -> List[dict]:
    rows = []
    with open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(hdr)}
        need = ("sample", "chrom", "left_pos", "right_pos", "gt_before", "gt_after")
        for k in need:
            if k not in idx:
                raise SystemExit(f"[FAIL] transitions missing col {k}: {hdr}")
        for line in f:
            p = line.rstrip("\n").split("\t")
            gb, ga = int(p[idx["gt_before"]]), int(p[idx["gt_after"]])
            if gb not in GT2ST or ga not in GT2ST:
                raise SystemExit(f"[FAIL] bad gt in transitions: {gb}/{ga}")
            rows.append(
                {
                    "sample": p[idx["sample"]],
                    "chr": p[idx["chrom"]],
                    "pos": int(p[idx["right_pos"]]),
                    "state_before": GT2ST[gb],
                    "state_after": GT2ST[ga],
                }
            )
    return rows


def load_seqlengths(fai: str) -> List[Tuple[str, int]]:
    rows = []
    with open(fai) as f:
        for line in f:
            chrom, length = line.split("\t")[:2]
            if chrom.startswith("chr") and chrom[3:].isdigit():
                n = int(chrom[3:])
                if 1 <= n <= 12:
                    rows.append((f"chr{n}", int(length)))
    rows.sort(key=lambda x: int(x[0][3:]))
    if len(rows) != 12:
        raise SystemExit(f"[FAIL] expect 12 chroms in fai, got {len(rows)}")
    return rows


def paint_sample(
    segs: List[dict], by_chr_pos: Dict[str, List[int]], n_markers: int, chrom_order: List[str]
) -> bytearray:
    """Paint GS markers; CO gaps split at midpoint (RTIGER-consistent, no false NA)."""
    out = bytearray(n_markers)
    offset = 0
    chrom_off: Dict[str, int] = {}
    for chrom in chrom_order:
        chrom_off[chrom] = offset
        offset += len(by_chr_pos[chrom])

    segs_by_chr: Dict[str, List[dict]] = defaultdict(list)
    for s in segs:
        segs_by_chr[s["chr"]].append(s)

    for chrom in chrom_order:
        positions = by_chr_pos.get(chrom, [])
        if not positions:
            continue
        sl = segs_by_chr.get(chrom, [])
        if not sl:
            continue
        off = chrom_off[chrom]
        nsl = len(sl)
        j = 0
        for i, pos in enumerate(positions):
            while j < nsl and sl[j]["end"] < pos:
                j += 1
            if j < nsl and sl[j]["start"] <= pos <= sl[j]["end"]:
                out[off + i] = ST_CODE[sl[j]["state"]]
                continue
            if j == 0:
                # before first segment → extend first state
                out[off + i] = ST_CODE[sl[0]["state"]]
                continue
            if j >= nsl:
                # after last segment → extend last state
                out[off + i] = ST_CODE[sl[-1]["state"]]
                continue
            # gap between sl[j-1] and sl[j]
            left = sl[j - 1]
            right = sl[j]
            mid = (left["end"] + right["start"]) // 2
            st = left["state"] if pos <= mid else right["state"]
            out[off + i] = ST_CODE[st]
    return out


def write_state_matrix(
    path: str,
    samples: List[str],
    by_seg: Dict[str, List[dict]],
    labels: List[str],
    by_chr_pos: Dict[str, List[int]],
    tag: str,
) -> Tuple[int, int, int, int, int]:
    n_markers = len(labels)
    chrom_order = sorted(by_chr_pos, key=chrom_sort_key)
    n_na = n_mat = n_het = n_pat = 0
    print(f"[INFO] write {tag}: n_markers={n_markers}", flush=True)
    with open(path, "w", buffering=1 << 20) as fo:
        fo.write("sample\t" + "\t".join(labels) + "\n")
        for si, sid in enumerate(samples):
            codes = paint_sample(by_seg[sid], by_chr_pos, n_markers, chrom_order)
            ctr = Counter(codes)
            n_na += ctr.get(0, 0)
            n_mat += ctr.get(1, 0)
            n_het += ctr.get(2, 0)
            n_pat += ctr.get(3, 0)
            fo.write(sid + "\t" + "\t".join(CODE_ST[c] for c in codes) + "\n")
            step = 10 if n_markers > 100000 else 40
            if (si + 1) % step == 0 or si + 1 == len(samples):
                print(f"[INFO] {tag} rows {si + 1}/{len(samples)}", flush=True)
    return n_markers, n_na, n_mat, n_het, n_pat


def promote_legacy_thin_to_dense(gs_dir: str) -> Optional[str]:
    """If legacy full-matrix viterbi_thin exists and dense missing, rename → dense."""
    thin = os.path.join(gs_dir, "viterbi_thin.tsv")
    dense = os.path.join(gs_dir, "viterbi_dense.tsv")
    if not os.path.isfile(thin) or os.path.getsize(thin) == 0:
        return None
    if os.path.isfile(dense) and os.path.getsize(dense) > 0:
        return dense
    with open(thin) as f:
        ncols = len(f.readline().rstrip("\n").split("\t")) - 1
    # full parentStable ~5e5; 50kb thin ~1e4
    if ncols < 100000:
        return None
    os.rename(thin, dense)
    print(f"[INFO] renamed legacy thin → dense (n_markers≈{ncols})", flush=True)
    return dense


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--eval-dir", required=True)
    ap.add_argument("--out-dir", default="")
    ap.add_argument("--segments", default="")
    ap.add_argument("--transitions", default="")
    ap.add_argument("--ancestry-vcf", default="")
    ap.add_argument("--fai", required=True)
    ap.add_argument("--bcftools", default="bcftools")
    ap.add_argument("--pgg-dest", default="")
    ap.add_argument("--pgg-prefix", default="pangenie_haplotypes")
    ap.add_argument("--skip-thin", action="store_true", help="skip GS matrices")
    ap.add_argument(
        "--gs-bin-bp",
        type=int,
        default=50000,
        help="GS thin bin size; 0=keep all parentStable sites as thin",
    )
    ap.add_argument(
        "--write-dense",
        action="store_true",
        help="also paint full parentStable → viterbi_dense.tsv (slow)",
    )
    ap.add_argument(
        "--reuse-dense",
        action="store_true",
        help="do not regenerate dense; rename legacy thin if needed",
    )
    args = ap.parse_args()
    if not args.write_dense:
        args.reuse_dense = True  # default: keep existing dense / rename

    eval_dir = os.path.abspath(args.eval_dir)
    out_dir = os.path.abspath(args.out_dir or os.path.join(eval_dir, "export_gs_rqtl"))
    seg_path = args.segments or os.path.join(eval_dir, "viterbi_segments_minspan.tsv")
    tr_path = args.transitions or os.path.join(eval_dir, "transitions.tsv")
    vcf = args.ancestry_vcf or os.path.join(eval_dir, "merge.ancestry_parentStable.vcf.gz")

    for p, name in ((seg_path, "segments"), (tr_path, "transitions"), (vcf, "ancestry-vcf"), (args.fai, "fai")):
        if not os.path.isfile(p) or os.path.getsize(p) == 0:
            raise SystemExit(f"[FAIL] missing {name}: {p}")

    print(f"[INFO] eval={eval_dir}", flush=True)
    print(f"[INFO] out={out_dir}", flush=True)
    print(f"[INFO] gs_bin_bp={args.gs_bin_bp}", flush=True)

    samples = load_samples(args.bcftools, vcf)
    print(f"[INFO] n_samples={len(samples)}", flush=True)
    by_seg = load_segments(seg_path)
    missing = [s for s in samples if s not in by_seg]
    if missing:
        raise SystemExit(f"[FAIL] {len(missing)} samples lack segments e.g. {missing[:5]}")
    extra = set(by_seg) - set(samples)
    if extra:
        print(f"[WARN] segments have {len(extra)} samples not in VCF (ignored)", flush=True)

    trans = load_transitions(tr_path)
    seqlen = load_seqlengths(args.fai)

    gs_dir = os.path.join(out_dir, "for_gs")
    rq_dir = os.path.join(out_dir, "for_rqtl")
    os.makedirs(gs_dir, exist_ok=True)
    os.makedirs(rq_dir, exist_ok=True)

    map_path = os.path.join(gs_dir, "viterbi_sample_map.tsv")
    with open(map_path, "w") as fo:
        fo.write("viterbi_name\tsample_id\n")
        for i, sid in enumerate(samples, 1):
            fo.write(f"Sample_{i}\t{sid}\n")
    map_rq = os.path.join(rq_dir, "viterbi_sample_map.tsv")
    hardlink_or_copy(map_path, map_rq)

    n_seg = 0
    raw_path = os.path.join(out_dir, "raw_viterbi_segments.tsv")
    rq_seg = os.path.join(rq_dir, "viterbi_segments.tsv")
    with open(raw_path, "w") as fr, open(rq_seg, "w") as fq:
        fr.write("sample\tchr\tfirst_marker\tlast_marker\tn_markers\tstate\n")
        fq.write("sample\tchr\tstart\tend\tstate\n")
        for sid in samples:
            for s in by_seg[sid]:
                fr.write(
                    f"{sid}\t{s['chr']}\t{s['start']}\t{s['end']}\t{s['n_markers']}\t{s['state']}\n"
                )
                fq.write(f"{sid}\t{s['chr']}\t{s['start']}\t{s['end']}\t{s['state']}\n")
                n_seg += 1

    rq_tr = os.path.join(rq_dir, "all_transitions.tsv")
    n_tr = 0
    with open(rq_tr, "w") as fo:
        fo.write("sample\tchr\tpos\tstate_before\tstate_after\n")
        for r in trans:
            if r["sample"] not in by_seg:
                continue
            fo.write(
                f"{r['sample']}\t{r['chr']}\t{r['pos']}\t{r['state_before']}\t{r['state_after']}\n"
            )
            n_tr += 1

    sl_path = os.path.join(out_dir, "seqlengths.tsv")
    with open(sl_path, "w") as fo:
        fo.write("chrom\tlength\n")
        for chrom, length in seqlen:
            fo.write(f"{chrom}\t{length}\n")

    thin_path = os.path.join(gs_dir, "viterbi_thin.tsv")
    dense_path = os.path.join(gs_dir, "viterbi_dense.tsv")
    sites_path = os.path.join(gs_dir, "gs_thin_sites.tsv")
    # GS-specific parent-stable markers are kept separate from the full marker set.
    ps50_path = os.path.join(gs_dir, "parent_stable_sites_50kb.tsv")
    thin_wrote = False
    n_markers = 0
    n_na = n_mat = n_het = n_pat = 0
    n_dense_markers = 0

    if not args.skip_thin:
        if args.reuse_dense:
            promote_legacy_thin_to_dense(gs_dir)

        print("[INFO] select GS thin markers (meanKC + bin center)…", flush=True)
        labels, by_chr_pos, site_rows = select_gs_markers(
            args.bcftools, vcf, args.gs_bin_bp
        )
        with open(sites_path, "w") as fo, open(ps50_path, "w") as fp:
            fo.write("chrom\tpos\tbin\tmean_kc\tdist_center\n")
            fp.write("chrom\tpos\tbin\tmean_kc\tdist_center\tset\n")
            for chrom, pos, bid, kc, dist in site_rows:
                fo.write(f"{chrom}\t{pos}\t{bid}\t{kc:.6f}\t{dist}\n")
                fp.write(
                    f"{chrom}\t{pos}\t{bid}\t{kc:.6f}\t{dist}\tparentStable_gs50kb\n"
                )
        # Export a two-column marker list for downstream tools.
        ancestry_dir = os.path.dirname(os.path.dirname(eval_dir))
        shared_dir = os.path.join(ancestry_dir, "shared")
        if os.path.isdir(shared_dir):
            shared_ps50 = os.path.join(shared_dir, "parent_stable_sites_gs50kb.tsv")
            with open(shared_ps50, "w") as fo:
                for chrom, pos, _, _, _ in site_rows:
                    fo.write(f"{chrom}\t{pos}\n")
            print(f"[INFO] wrote {shared_ps50}", flush=True)

        n_markers, n_na, n_mat, n_het, n_pat = write_state_matrix(
            thin_path, samples, by_seg, labels, by_chr_pos, "thin"
        )
        thin_wrote = True

        if args.write_dense and not args.reuse_dense:
            print("[INFO] write full dense matrix…", flush=True)
            d_labels, d_by, _ = select_gs_markers(args.bcftools, vcf, 0)
            n_dense_markers, _, _, _, _ = write_state_matrix(
                dense_path, samples, by_seg, d_labels, d_by, "dense"
            )
        elif os.path.isfile(dense_path):
            with open(dense_path) as f:
                n_dense_markers = len(f.readline().rstrip("\n").split("\t")) - 1
            print(f"[INFO] reuse dense n_markers={n_dense_markers}", flush=True)
    else:
        print("[INFO] skip GS matrices", flush=True)

    meta_path = os.path.join(out_dir, "EXPORT_META.tsv")
    prev_meta: Dict[str, str] = {}
    if args.skip_thin and os.path.isfile(meta_path):
        with open(meta_path) as fm:
            next(fm, None)
            for line in fm:
                if "\t" in line:
                    k, v = line.rstrip("\n").split("\t", 1)
                    prev_meta[k] = v

    n_cells = len(samples) * n_markers if n_markers else 0
    with open(meta_path, "w") as fo:
        fo.write("key\tvalue\n")
        fo.write(f"eval_dir\t{eval_dir}\n")
        fo.write(f"segments\t{seg_path}\n")
        fo.write(f"transitions\t{tr_path}\n")
        fo.write(f"ancestry_vcf\t{vcf}\n")
        fo.write(f"n_samples\t{len(samples)}\n")
        fo.write(f"n_segments\t{n_seg}\n")
        fo.write(f"n_transitions\t{n_tr}\n")
        fo.write(f"gs_bin_bp\t{args.gs_bin_bp}\n")
        fo.write("gs_select\tmax_meanKC_then_nearest_bin_center\n")
        fo.write("gs_paint\tclosed_segment_or_transition_midpoint_or_terminal_extend\n")
        if thin_wrote:
            fo.write(f"n_markers_thin\t{n_markers}\n")
            fo.write(f"n_markers\t{n_markers}\n")
            fo.write(f"n_cells\t{n_cells}\n")
            fo.write(f"n_mat\t{n_mat}\n")
            fo.write(f"n_het\t{n_het}\n")
            fo.write(f"n_pat\t{n_pat}\n")
            fo.write(f"n_NA\t{n_na}\n")
            fo.write(f"NA_rate\t{n_na / n_cells:.6f}\n")
            fo.write(f"NA_n\t{n_na}\n")
            fo.write(f"NA_d\t{n_cells}\n")
            if n_dense_markers:
                fo.write(f"n_markers_dense\t{n_dense_markers}\n")
        else:
            for k in (
                "n_markers",
                "n_markers_thin",
                "n_markers_dense",
                "n_cells",
                "n_mat",
                "n_het",
                "n_pat",
                "n_NA",
                "NA_rate",
                "NA_n",
                "NA_d",
                "gs_bin_bp",
            ):
                if k in prev_meta:
                    fo.write(f"{k}\t{prev_meta[k]}\n")
        fo.write(f"segments_sha256_first_64MiB\t{sha256_first(seg_path)}\n")
        fo.write(f"coord\tMM_chr1\n")
        fo.write("state_map\t0=mat(Hom-MM);1=het;2=pat(Hom-TS)\n")

    readme = os.path.join(out_dir, "README.md")
    with open(readme, "w") as fo:
        fo.write("# PanGenie Canon → GS / Rqtl export\n\n")
        fo.write("- Source: Canon `viterbi_segments_minspan` + parentStable ancestry VCF\n")
        fo.write("- Coord: MM `chr1`…`chr12` (same habit as GFA)\n")
        fo.write("- States: mat/het/pat (= Hom-MM / Het / Hom-TS)\n")
        fo.write(
            f"- `for_gs/viterbi_thin.tsv`: minspan painted on "
            f"**{args.gs_bin_bp} bp** bin thin (max meanKC → nearest center)\n"
        )
        fo.write(
            "- GS paint: in-segment closed; inter-segment gap at midpoint "
            "`(end_i+start_{i+1})//2` (≤mid→left); terminal extend first/last "
            "(RTIGER-consistent; no false NA in CO gaps)\n"
        )
        fo.write("- `for_gs/viterbi_dense.tsv`: **full** parentStable paint (optional/heavy)\n")
        fo.write(
            "- `for_gs/parent_stable_sites_50kb.tsv` (+ `gs_thin_sites.tsv`): "
            "**GS-only** parentStable subset; also `shared/parent_stable_sites_gs50kb.tsv`\n"
        )
        fo.write(
            "- Do **not** mix with Canon full set "
            "`shared/parent_stable_sites.tsv` / `merge.ancestry_parentStable.vcf.gz`\n"
        )
        fo.write(
            "- GS select uses FORMAT/KC only (mean); **not** offspring GT / missingness\n"
        )
        fo.write("- `for_rqtl/`: RTIGER step06 column schema (unchanged by GS thin)\n")
        fo.write("- Local only; see EXPORT_META.tsv\n")

    if args.pgg_dest:
        dest = os.path.abspath(args.pgg_dest)
        os.makedirs(dest, exist_ok=True)
        pref = args.pgg_prefix
        copies = [
            (raw_path, f"{pref}_raw_viterbi_segments.tsv"),
            (map_path, f"{pref}_viterbi_sample_map.tsv"),
            (rq_seg, f"{pref}_for_rqtl_viterbi_segments.tsv"),
            (rq_tr, f"{pref}_for_rqtl_all_transitions.tsv"),
        ]
        if not args.skip_thin and os.path.isfile(thin_path):
            copies.append((thin_path, f"{pref}_viterbi_thin.tsv"))
            copies.append((sites_path, f"{pref}_gs_thin_sites.tsv"))
            if os.path.isfile(ps50_path):
                copies.append((ps50_path, f"{pref}_parent_stable_sites_50kb.tsv"))
        if os.path.isfile(dense_path):
            copies.append((dense_path, f"{pref}_viterbi_dense.tsv"))
        for src, name in copies:
            dst = os.path.join(dest, name)
            hardlink_or_copy(src, dst)
            print(f"[INFO] PGG {dst}", flush=True)

    print(
        f"[DONE] n_seg={n_seg} n_tr={n_tr} n_thin={n_markers} n_dense={n_dense_markers} → {out_dir}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
