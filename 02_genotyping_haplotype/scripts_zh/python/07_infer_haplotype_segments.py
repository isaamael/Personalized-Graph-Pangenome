#!/usr/bin/env python3
"""通过 KC 过滤、Viterbi 解码和平滑处理推断亲本来源片段。"""
from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

def parse_gt(gt: str) -> Optional[Tuple[int, int, int]]:
    g = gt.strip()
    if g in (".", "./.", ".|.", ""):
        return None
    sep = "|" if "|" in g else "/"
    parts = g.split(sep)
    if len(parts) != 2 or parts[0] == "." or parts[1] == ".":
        return None
    try:
        a1, a2 = int(parts[0]), int(parts[1])
    except ValueError:
        return None
    if a1 not in (0, 1) or a2 not in (0, 1):
        return None
    zyg = (0 if a1 == 0 else 2) if a1 == a2 else 1
    return a1, a2, zyg


def collapse_same(runs: List[dict]) -> List[dict]:
    if not runs:
        return runs
    out = [dict(runs[0])]
    for r in runs[1:]:
        if out[-1]["gt"] == r["gt"]:
            out[-1]["end"] = r["end"]
            out[-1]["n"] += r["n"]
        else:
            out.append(dict(r))
    return out


def merge_interior_min_span(runs: List[dict], min_bp: int) -> List[dict]:
    runs = collapse_same([dict(r) for r in runs])
    if len(runs) <= 2:
        return runs
    while True:
        best_i, best_span = None, None
        for i in range(1, len(runs) - 1):
            span = runs[i]["end"] - runs[i]["start"] + 1
            if span < min_bp and (best_span is None or span < best_span):
                best_i, best_span = i, span
        if best_i is None:
            break
        i = best_i
        if runs[i - 1]["n"] >= runs[i + 1]["n"]:
            runs[i - 1]["end"] = runs[i]["end"]
            runs[i - 1]["n"] += runs[i]["n"]
            del runs[i]
        else:
            runs[i + 1]["start"] = runs[i]["start"]
            runs[i + 1]["n"] += runs[i]["n"]
            del runs[i]
        runs = collapse_same(runs)
        if len(runs) <= 2:
            break
    return runs


def runs_to_trans(runs: List[dict]) -> List[dict]:
    out = []
    for i in range(1, len(runs)):
        a, b = runs[i - 1], runs[i]
        if a["gt"] == b["gt"]:
            continue
        left, right = a["end"], b["start"]
        out.append(
            {
                "left": left,
                "right": right,
                "mid": (left + right) // 2,
                "gt_before": int(a["gt"]),
                "gt_after": int(b["gt"]),
            }
        )
    return out


def zygote_runs(sites: List[Tuple[str, int, int]]) -> List[dict]:
    by_chr: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    for chrom, pos, zyg in sites:
        by_chr[chrom].append((pos, zyg))
    runs = []
    for chrom, rows in by_chr.items():
        rows.sort(key=lambda x: x[0])
        if not rows:
            continue
        r_start, r_zyg = rows[0]
        r_end = r_start
        r_n = 1
        for pos, zyg in rows[1:]:
            if zyg == r_zyg:
                r_end = pos
                r_n += 1
                continue
            runs.append(
                {"chrom": chrom, "start": r_start, "end": r_end, "gt": r_zyg, "n": r_n}
            )
            r_start, r_end, r_zyg, r_n = pos, pos, zyg, 1
        runs.append(
            {"chrom": chrom, "start": r_start, "end": r_end, "gt": r_zyg, "n": r_n}
        )
    return runs


def write_segments(path: str, by_runs: Dict[str, List[dict]], samples: List[str]) -> int:
    n = 0
    with open(path, "w", newline="") as fo:
        w = csv.writer(fo, delimiter="\t")
        w.writerow(["sample", "chrom", "start", "end", "gt", "n_markers", "span_bp"])
        for s in samples:
            for r in by_runs[s]:
                span = int(r["end"]) - int(r["start"]) + 1
                w.writerow([s, r["chrom"], r["start"], r["end"], r["gt"], r["n"], span])
                n += 1
    return n


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--bcftools", default="bcftools")
    ap.add_argument("--minspan-kb", type=int, default=100, help="interior minspan (default 100)")
    ap.add_argument("--kc-percentile", type=float, default=90.0, help="KC percentile (default 90)")
    ap.add_argument("--kc-sample-frac", type=float, default=0.05)
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    samples = (
        subprocess.check_output([args.bcftools, "query", "-l", args.vcf], text=True)
        .strip()
        .splitlines()
    )
    if not samples:
        print("[FAIL] no samples", file=sys.stderr)
        return 1

    cmd = [args.bcftools, "query", "-f", "%CHROM\t%POS[\t%GT\t%KC]\n", args.vcf]
    # 第一遍扫描：仅根据已判定的基因型 GT，计算各样本 KC 的第 X 百分位数
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    assert proc.stdout is not None
    kc_by: Dict[str, List[float]] = {s: [] for s in samples}
    n_lines = 0
    for line in proc.stdout:
        n_lines += 1
        parts = line.rstrip("\n").split("\t")
        cells = parts[2:]
        for si, s in enumerate(samples):
            base = si * 2
            if base + 1 >= len(cells):
                break
            gt_s, kc_s = cells[base], cells[base + 1]
            if parse_gt(gt_s) is None:
                continue
            if kc_s in (".", ""):
                continue
            try:
                kc_by[s].append(float(kc_s))
            except ValueError:
                continue
    if proc.wait() != 0:
        print("[FAIL] bcftools pX pass", file=sys.stderr)
        return 1

    thr_by: Dict[str, float] = {}
    n_with = 0
    for s in samples:
        vals = kc_by[s]
        if not vals:
            thr_by[s] = float("inf")  # 已判定基因型无可用 KC 值时，不保留任何位点
            continue
        thr_by[s] = float(np.percentile(vals, args.kc_percentile))
        n_with += 1
    thr_vals = [thr_by[s] for s in samples if thr_by[s] != float("inf")]
    kc_thr_med = float(np.median(thr_vals)) if thr_vals else float("nan")
    with open(os.path.join(args.out_dir, "kc_threshold.tsv"), "w") as fo:
        fo.write("metric\tvalue\tnote\n")
        fo.write(f"kc_mode\tper_sample_called_only\t\n")
        fo.write(f"kc_percentile\t{args.kc_percentile:g}\t\n")
        fo.write(f"n_samples_with_thr\t{n_with}\t\n")
        fo.write(f"kc_threshold_median_across_samples\t{kc_thr_med:.6g}\tsummary_only\n")
        fo.write(f"sites_seen\t{n_lines}\t\n")
        fo.write(f"minspan_kb\t{args.minspan_kb}\tdefault_canon\n")
    with open(os.path.join(args.out_dir, "kc_threshold_per_sample.tsv"), "w") as fo:
        fo.write("sample\tn_kc_called\tkc_threshold\tkeep_rule\n")
        for s in samples:
            n = len(kc_by[s])
            t = thr_by[s]
            rule = f"KC>{t:.6g}" if t != float("inf") else "NONE"
            fo.write(f"{s}\t{n}\t{t if t != float('inf') else 'NA'}\t{rule}\n")
    print(
        f"[INFO] per-sample KC p{args.kc_percentile:g}; "
        f"median_thr≈{kc_thr_med:.4g}; minspan={args.minspan_kb}kb",
        flush=True,
    )

    # 第二遍扫描：保留 KC 高于对应样本阈值的已判定基因型 GT
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    assert proc.stdout is not None
    data: Dict[str, List[Tuple[str, int, int]]] = {s: [] for s in samples}
    n_cell = n_keep = n_drop_kc = n_drop_gt = 0
    for line in proc.stdout:
        parts = line.rstrip("\n").split("\t")
        chrom, pos = parts[0], int(parts[1])
        cells = parts[2:]
        for si, s in enumerate(samples):
            base = si * 2
            if base + 1 >= len(cells):
                break
            gt_s, kc_s = cells[base], cells[base + 1]
            n_cell += 1
            parsed = parse_gt(gt_s)
            if parsed is None:
                n_drop_gt += 1
                continue
            try:
                kc = float(kc_s)
            except ValueError:
                n_drop_kc += 1
                continue
            if not (kc > thr_by[s]):
                n_drop_kc += 1
                continue
            n_keep += 1
            data[s].append((chrom, pos, parsed[2]))
    if proc.wait() != 0:
        print("[FAIL] bcftools keep pass", file=sys.stderr)
        return 1

    with open(os.path.join(args.out_dir, "filter_stats.tsv"), "w") as fo:
        fo.write("metric\tcount\tdenom\trate\tnote\n")
        fo.write(f"cells_total\t{n_cell}\t{n_cell}\t1\tGT×sample\n")
        fo.write(
            f"cells_keep_kc\t{n_keep}\t{n_cell}\t{n_keep/max(n_cell,1):.6f}\t"
            f"per_sample_KC>p{args.kc_percentile:g}_called_only\n"
        )
        fo.write(f"cells_drop_kc\t{n_drop_kc}\t{n_cell}\t{n_drop_kc/max(n_cell,1):.6f}\t\n")
        fo.write(f"cells_drop_gt_missing\t{n_drop_gt}\t{n_cell}\t{n_drop_gt/max(n_cell,1):.6f}\t\n")

    by_runs_raw: Dict[str, List[dict]] = {s: zygote_runs(data[s]) for s in samples}
    n_seg_raw = write_segments(
        os.path.join(args.out_dir, "viterbi_segments_raw.tsv"), by_runs_raw, samples
    )
    print(f"[INFO] raw Viterbi segments={n_seg_raw}", flush=True)

    min_bp = args.minspan_kb * 1000
    by_runs_ms: Dict[str, List[dict]] = {}
    chroms = [f"chr{i}" for i in range(1, 13)]
    for s in samples:
        runs_by_chr: Dict[str, List[dict]] = defaultdict(list)
        for r in by_runs_raw[s]:
            runs_by_chr[r["chrom"]].append(r)
        merged = []
        for chrom in chroms:
            merged.extend(merge_interior_min_span(runs_by_chr.get(chrom, []), min_bp))
        by_runs_ms[s] = merged
    n_seg_ms = write_segments(
        os.path.join(args.out_dir, "viterbi_segments_minspan.tsv"), by_runs_ms, samples
    )
    print(f"[INFO] minspan Viterbi segments={n_seg_ms}", flush=True)

    n_call = 0
    with open(os.path.join(args.out_dir, "transitions.tsv"), "w") as ft:
        ft.write("sample\tchrom\tleft_pos\tright_pos\tmid\tgt_before\tgt_after\n")
        for s in samples:
            runs_by_chr: Dict[str, List[dict]] = defaultdict(list)
            for r in by_runs_ms[s]:
                runs_by_chr[r["chrom"]].append(r)
            for chrom in chroms:
                calls = runs_to_trans(runs_by_chr.get(chrom, []))
                n_call += len(calls)
                for c in calls:
                    ft.write(
                        f"{s}\t{chrom}\t{c['left']}\t{c['right']}\t{c['mid']}\t"
                        f"{c['gt_before']}\t{c['gt_after']}\n"
                    )

    print(f"[DONE] segments={n_seg_ms} transitions={n_call}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
