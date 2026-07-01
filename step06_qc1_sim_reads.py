#!/usr/bin/env python3
"""QC1：sim reads + truth CO。mode=fast（默认，<5min）| deep（全量 fq 扫描）。"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import random
import re
import subprocess
import sys
from typing import Dict, List, Tuple


def open_fq(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def stream_count_bases(path: str) -> Tuple[int, int]:
    cmd = (
        f'(pigz -dc "{path}" 2>/dev/null || gzip -dc "{path}") | '
        "awk 'NR%4==2{b+=length($0);n++} END{print n+0, b+0}'"
    )
    proc = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.strip() or f"count failed: {path}")
    parts = proc.stdout.strip().split()
    return int(parts[0]), int(parts[1])


def est_depth_from_size(r1: str, r2: str, expected_depth: float, ref_pair_bytes: float) -> float:
    total = os.path.getsize(r1) + os.path.getsize(r2)
    if ref_pair_bytes <= 0:
        return 0.0
    return expected_depth * total / ref_pair_bytes


def read_head_pairs(path: str, max_pairs: int) -> Tuple[List[str], List[int]]:
    ids: List[str] = []
    lens: List[int] = []
    with open_fq(path) as fh:
        while len(ids) < max_pairs:
            h = fh.readline()
            if not h:
                break
            seq = fh.readline().strip()
            fh.readline()
            fh.readline()
            ids.append(h[1:].strip().split()[0])
            lens.append(len(seq))
    return ids, lens


def read_random_ids(path: str, max_n: int, seed: int, rate: float = 0.0005) -> List[str]:
    cmd = (
        f'(pigz -dc "{path}" 2>/dev/null || gzip -dc "{path}") | '
        f"awk -v seed={seed} -v rate={rate} "
        "'BEGIN{srand(seed+0)} NR%4==1 {if(rand()<rate) print substr($0,2)}'"
    )
    proc = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.strip() or f"sample ids failed: {path}")
    ids = [x.split()[0] for x in proc.stdout.splitlines() if x.strip()]
    if len(ids) > max_n:
        ids = random.Random(seed).sample(ids, max_n)
    return ids


def norm_pair_id(rid: str) -> str:
    return re.sub(r"/[12]$", "", rid)


def check_pairing(r1_ids: List[str], r2_ids: List[str]) -> Tuple[int, int]:
    n = min(len(r1_ids), len(r2_ids))
    mis = sum(1 for i in range(n) if norm_pair_id(r1_ids[i]) != norm_pair_id(r2_ids[i]))
    return n, mis


def dup_rate(ids: List[str]) -> Tuple[int, int, float]:
    n = len(ids)
    u = len(set(ids))
    return n, u, (1.0 - u / n) if n else 0.0


def prefix_stats(ids: List[str]) -> Tuple[float, float, int]:
    mm = ts = bare = 0
    for rid in ids:
        if rid.startswith("MM_"):
            mm += 1
        elif rid.startswith("TS_"):
            ts += 1
        else:
            bare += 1
    n = len(ids) or 1
    return mm / n, ts / n, bare


def discover_samples(reads_dir: str, pf_depth: str, n_f2: int, n_pf: int) -> List[Tuple[str, str, str, str]]:
    out: List[Tuple[str, str, str, str]] = []
    for i in range(1, n_f2 + 1):
        sid = f"simF2_{i:03d}"
        out.append((sid, "simF2", os.path.join(reads_dir, "simF2", sid, f"{sid}.R1.fq.gz"),
                    os.path.join(reads_dir, "simF2", sid, f"{sid}.R2.fq.gz")))
    for i in range(1, n_pf + 1):
        sid = f"pseudoF1_{i:02d}"
        out.append((sid, "simF1",
                    os.path.join(reads_dir, "pseudoF1", f"{pf_depth}x", f"{sid}.R1.fq.gz"),
                    os.path.join(reads_dir, "pseudoF1", f"{pf_depth}x", f"{sid}.R2.fq.gz")))
    return out


def discover_spot(n_f2: int, n_pf: int, spot_f2: int, spot_f1: int, seed: int) -> set:
    rng = random.Random(seed)
    f2 = [f"simF2_{i:03d}" for i in range(1, n_f2 + 1)]
    f1 = [f"pseudoF1_{i:02d}" for i in range(1, n_pf + 1)]
    pick = set(f2 if spot_f2 >= n_f2 else rng.sample(f2, min(spot_f2, len(f2))))
    pick |= set(f1 if spot_f1 >= n_pf else rng.sample(f1, min(spot_f1, len(f1))))
    return pick


def load_cohort_stats(path: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    if not os.path.isfile(path):
        return out
    with open(path) as fh:
        for line in fh:
            if line.startswith("metric") or not line.strip():
                continue
            k, v = line.rstrip("\n").split("\t", 1)
            out[k] = v
    return out


def qc_truth(
    truth_dir: str, hap_ratio_min: float, hap_ratio_max: float,
    co_max: int, co_mean_lo: float, co_mean_hi: float,
) -> Tuple[List[dict], dict, bool]:
    summary_path = os.path.join(truth_dir, "truth_co_summary.tsv")
    cohort_path = os.path.join(truth_dir, "truth_co_cohort.tsv")
    if not os.path.isfile(summary_path):
        return [], {"status": "missing", "fail_reason": "truth_co_summary_missing"}, False

    cohort = load_cohort_stats(cohort_path)
    expected_co = float(cohort.get("expected_co_per_gamete", "8.0"))
    rows: List[dict] = []
    co_vals: List[int] = []
    with open(summary_path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            n_co = int(row["n_parent_co"])
            ratio = float(row["hap_len_ratio"])
            co_vals.append(n_co)
            len_pass = hap_ratio_min <= ratio <= hap_ratio_max
            co_pass = 0 <= n_co <= co_max
            rows.append({
                "sample": row["sample"], "gamete": row["gamete"],
                "n_parent_co": n_co, "n_frag": row["n_frag"],
                "hap_len_bp": row["hap_len_bp"], "hap_len_ratio": row["hap_len_ratio"],
                "expected_co_genome": f"{expected_co:.2f}",
                "co_pass": "PASS" if co_pass else "FAIL",
                "hap_len_pass": "PASS" if len_pass else "FAIL",
                "overall_pass": "PASS" if (co_pass and len_pass) else "FAIL",
            })

    mean_co = sum(co_vals) / len(co_vals) if co_vals else 0.0
    lo, hi = expected_co * co_mean_lo, expected_co * co_mean_hi
    cohort_co_pass = lo <= mean_co <= hi
    n_pass = sum(1 for r in rows if r["overall_pass"] == "PASS")
    summary = {
        "status": "ok", "n_gamete": len(rows), "n_pass": n_pass,
        "expected_co_per_gamete": f"{expected_co:.2f}",
        "obs_co_mean": f"{mean_co:.2f}",
        "obs_co_min": str(min(co_vals) if co_vals else 0),
        "obs_co_max": str(max(co_vals) if co_vals else 0),
        "cohort_co_pass": "PASS" if cohort_co_pass else "FAIL",
        "overall_pass": "PASS" if cohort_co_pass and n_pass == len(rows) else "FAIL",
    }
    return rows, summary, summary["overall_pass"] == "PASS"


def qc_one_fast(
    sample: str, cohort: str, r1: str, r2: str, expected_depth: float,
    ref_pair_bytes: float, max_pairs: int, depth_tol: float, dup_tol: float,
    prefix_min: float, do_spot: bool,
) -> Dict[str, object]:
    row: Dict[str, object] = {
        "sample": sample, "cohort": cohort, "mode": "fast", "r1": r1, "r2": r2,
        "expected_depth": f"{expected_depth:.2f}", "status": "missing", "overall_pass": "FAIL",
    }
    if not (os.path.isfile(r1) and os.path.isfile(r2)):
        row["fail_reason"] = "fq_missing"
        return row

    s1, s2 = os.path.getsize(r1), os.path.getsize(r2)
    size_ok = min(s1, s2) >= 100 * 1024 * 1024
    if not size_ok:
        row.update({"status": "ok", "fail_reason": "fq_too_small", "size_pass": "FAIL"})
        return row

    est = est_depth_from_size(r1, r2, expected_depth, ref_pair_bytes)
    depth_ratio = est / expected_depth if expected_depth > 0 else 0.0
    depth_pass = abs(depth_ratio - 1.0) <= depth_tol

    row.update({
        "status": "ok", "size_r1_mb": f"{s1/1024/1024:.1f}", "size_r2_mb": f"{s2/1024/1024:.1f}",
        "size_pass": "PASS", "est_depth": f"{est:.3f}", "depth_method": "file_size",
        "depth_ratio": f"{depth_ratio:.3f}", "depth_pass": "PASS" if depth_pass else "FAIL",
    })

    if not do_spot:
        row["pair_pass"] = "NA"
        row["prefix_pass"] = "NA" if cohort == "simF2" else "NA"
        row["overall_pass"] = "PASS" if depth_pass else "FAIL"
        return row

    r1_ids, _ = read_head_pairs(r1, max_pairs)
    r2_ids, _ = read_head_pairs(r2, max_pairs)
    pair_n, pair_mis = check_pairing(r1_ids, r2_ids)
    pair_pass = pair_mis == 0
    _, _, id_dup = dup_rate(r1_ids)
    id_dup_pass = id_dup <= dup_tol

    if cohort == "simF1":
        prefix_pass = "NA"
        mm_pct = ts_pct = bare_n = "NA"
    else:
        prefix_pass = True
        mm_pct = ts_pct = bare_n = "NA"

    overall = depth_pass and pair_pass and id_dup_pass and (prefix_pass in (True, "NA"))
    row.update({
        "pair_checked": pair_n, "pair_mismatch": pair_mis,
        "pair_pass": "PASS" if pair_pass else "FAIL",
        "id_dup_pass": "PASS" if id_dup_pass else "FAIL",
        "prefix_mm_pct": mm_pct if cohort == "simF1" else "NA",
        "prefix_ts_pct": ts_pct if cohort == "simF1" else "NA",
        "prefix_bare_n": bare_n if cohort == "simF1" else "NA",
        "prefix_pass": prefix_pass if cohort == "simF1" else "NA",
        "spot_checked": "1",
        "overall_pass": "PASS" if overall else "FAIL",
    })
    return row


def qc_one_deep(
    sample: str, cohort: str, r1: str, r2: str, expected_depth: float, genome_size: int,
    max_pairs: int, depth_tol: float, dup_tol: float, prefix_min: float,
) -> Dict[str, object]:
    row: Dict[str, object] = {
        "sample": sample, "cohort": cohort, "mode": "deep", "r1": r1, "r2": r2,
        "expected_depth": f"{expected_depth:.2f}", "status": "missing", "overall_pass": "FAIL",
    }
    if not (os.path.isfile(r1) and os.path.isfile(r2)):
        row["fail_reason"] = "fq_missing"
        return row
    try:
        n1, b1 = stream_count_bases(r1)
        n2, b2 = stream_count_bases(r2)
    except RuntimeError as e:
        row["fail_reason"] = str(e)
        return row
    if n1 != n2:
        row.update({"status": "ok", "pair_pass": "FAIL", "fail_reason": "read_count_mismatch", "overall_pass": "FAIL"})
        return row

    obs_depth = (b1 + b2) / genome_size if genome_size > 0 else 0.0
    depth_ratio = obs_depth / expected_depth if expected_depth > 0 else 0.0
    depth_pass = abs(depth_ratio - 1.0) <= depth_tol
    r1_ids, _ = read_head_pairs(r1, max_pairs)
    r2_ids, _ = read_head_pairs(r2, max_pairs)
    pair_n, pair_mis = check_pairing(r1_ids, r2_ids)
    pair_pass = pair_mis == 0
    id_n, id_u, id_dup = dup_rate(r1_ids)
    id_dup_pass = id_dup <= dup_tol

    if cohort == "simF1":
        try:
            pid_ids = read_random_ids(r1, max_pairs, hash(sample) % 100000)
            mm_pct, ts_pct, bare_n = prefix_stats(pid_ids)
            prefix_pass = len(pid_ids) >= 1000 and bare_n == 0 and mm_pct >= prefix_min and ts_pct >= prefix_min
        except RuntimeError:
            mm_pct, ts_pct, bare_n = prefix_stats(r1_ids)
            prefix_pass = bare_n == 0 and mm_pct >= prefix_min and ts_pct >= prefix_min
    else:
        mm_pct = ts_pct = bare_n = "NA"
        prefix_pass = True

    overall = depth_pass and pair_pass and id_dup_pass and prefix_pass
    row.update({
        "status": "ok", "n_reads": n1, "est_depth": f"{obs_depth:.3f}", "depth_method": "full_scan",
        "depth_ratio": f"{depth_ratio:.3f}", "depth_pass": "PASS" if depth_pass else "FAIL",
        "pair_checked": pair_n, "pair_mismatch": pair_mis, "pair_pass": "PASS" if pair_pass else "FAIL",
        "id_dup_pass": "PASS" if id_dup_pass else "FAIL",
        "prefix_mm_pct": f"{mm_pct:.3f}" if cohort == "simF1" else "NA",
        "prefix_ts_pct": f"{ts_pct:.3f}" if cohort == "simF1" else "NA",
        "prefix_pass": ("PASS" if prefix_pass else "FAIL") if cohort == "simF1" else "NA",
        "overall_pass": "PASS" if overall else "FAIL",
    })
    return row


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reads-dir", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--truth-dir", default="")
    ap.add_argument("--truth-out-tsv", default="")
    ap.add_argument("--genome-size", type=int, required=True)
    ap.add_argument("--mode", choices=("fast", "deep"), default="fast")
    ap.add_argument("--art-haplo-depth", type=float, default=5.0)
    ap.add_argument("--pf-depth", type=float, default=10.0)
    ap.add_argument("--pf-depth-label", default="10")
    ap.add_argument("--ref-pair-bytes-f2", type=float, default=5.0e9)
    ap.add_argument("--ref-pair-bytes-f1", type=float, default=5.8e9)
    ap.add_argument("--n-f2", type=int, default=50)
    ap.add_argument("--n-pf", type=int, default=10)
    ap.add_argument("--spot-f2", type=int, default=3)
    ap.add_argument("--spot-f1", type=int, default=3)
    ap.add_argument("--max-pairs", type=int, default=20000)
    ap.add_argument("--depth-tol", type=float, default=0.20)
    ap.add_argument("--dup-tol", type=float, default=0.001)
    ap.add_argument("--prefix-min", type=float, default=0.35)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--samples", default="")
    ap.add_argument("--reads-only", action="store_true")
    ap.add_argument("--truth-only", action="store_true")
    ap.add_argument("--hap-ratio-min", type=float, default=0.95)
    ap.add_argument("--hap-ratio-max", type=float, default=1.05)
    ap.add_argument("--co-max", type=int, default=25)
    ap.add_argument("--co-mean-lo", type=float, default=0.5)
    ap.add_argument("--co-mean-hi", type=float, default=1.5)
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out_tsv) or ".", exist_ok=True)
    truth_dir = args.truth_dir or os.path.dirname(args.out_tsv.replace("/qc/", "/truth/"))
    truth_out = args.truth_out_tsv or os.path.join(os.path.dirname(args.out_tsv), "QC1_truth_co.tsv")

    truth_ok = True
    if not args.reads_only:
        truth_rows, truth_summary, truth_ok = qc_truth(
            truth_dir, args.hap_ratio_min, args.hap_ratio_max,
            args.co_max, args.co_mean_lo, args.co_mean_hi,
        )
        truth_fields = [
            "sample", "gamete", "n_parent_co", "n_frag", "hap_len_bp", "hap_len_ratio",
            "expected_co_genome", "co_pass", "hap_len_pass", "overall_pass",
        ]
        with open(truth_out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=truth_fields, delimiter="\t")
            w.writeheader()
            w.writerows(truth_rows)
        with open(truth_out + ".summary.tsv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(truth_summary.keys()), delimiter="\t")
            w.writeheader()
            w.writerow(truth_summary)

    reads_ok = True
    if not args.truth_only:
        expected_f2 = 2.0 * args.art_haplo_depth
        samples = discover_samples(args.reads_dir, args.pf_depth_label, args.n_f2, args.n_pf)
        if args.samples.strip():
            want = {s.strip() for s in args.samples.split(",") if s.strip()}
            samples = [x for x in samples if x[0] in want]

        spot = discover_spot(args.n_f2, args.n_pf, args.spot_f2, args.spot_f1, args.seed)
        rows = []
        for sample, cohort, r1, r2 in samples:
            exp = expected_f2 if cohort == "simF2" else args.pf_depth
            ref_bytes = args.ref_pair_bytes_f2 if cohort == "simF2" else args.ref_pair_bytes_f1
            if args.mode == "fast":
                rows.append(qc_one_fast(
                    sample, cohort, r1, r2, exp, ref_bytes,
                    args.max_pairs, args.depth_tol, args.dup_tol, args.prefix_min,
                    do_spot=(sample in spot),
                ))
            else:
                rows.append(qc_one_deep(
                    sample, cohort, r1, r2, exp, args.genome_size,
                    args.max_pairs, args.depth_tol, args.dup_tol, args.prefix_min,
                ))

        fields = [
            "sample", "cohort", "mode", "status", "expected_depth", "est_depth", "depth_method",
            "depth_ratio", "depth_pass", "size_r1_mb", "size_r2_mb", "size_pass",
            "pair_checked", "pair_mismatch", "pair_pass", "id_dup_pass",
            "prefix_mm_pct", "prefix_ts_pct", "prefix_bare_n", "prefix_pass",
            "spot_checked", "overall_pass", "r1", "r2", "fail_reason",
        ]
        with open(args.out_tsv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)

        n_pass = sum(1 for r in rows if r.get("overall_pass") == "PASS")
        n_spot = sum(1 for r in rows if r.get("spot_checked") == "1")
        n_fail = len(rows) - n_pass
        reads_ok = n_fail == 0
        print(
            f"QC1 mode={args.mode} reads samples={len(rows)} pass={n_pass} fail={n_fail} spot={n_spot} out={args.out_tsv}",
            file=sys.stderr,
        )

    ok = (args.reads_only and reads_ok) or (args.truth_only and truth_ok) or (truth_ok and reads_ok)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
