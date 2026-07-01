#!/usr/bin/env python3
"""QC2：比对产物（BAM/GAM）抽查 — 比对率、深度、F1 混样比例、共线性 marker AF。"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import random
import re
import subprocess
import sys
from collections import Counter
from typing import Dict, List, Optional, Tuple

SAMTOOLS = os.environ.get("SAMTOOLS", "samtools")
VG_BIN = os.environ.get("VG", "vg")


def run_cmd(cmd: List[str], timeout: int = 3600) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)


def parse_flagstat(path: str) -> Dict[str, int]:
    out: Dict[str, int] = {}
    proc = run_cmd([SAMTOOLS, "flagstat", path])
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.strip() or f"flagstat failed: {path}")
    for line in proc.stdout.splitlines():
        m = re.match(r"^(\d+)\s+\+\s+\d+\s+(.+)$", line)
        if m:
            out[m.group(2).strip()] = int(m.group(1))
        else:
            m2 = re.match(r"^(\d+)\s+(.+)$", line)
            if m2:
                out[m2.group(2).strip()] = int(m2.group(1))
    return out


def parse_idxstats(path: str) -> Tuple[int, int, int]:
    """返回 (mapped, unmapped, total)。"""
    proc = run_cmd([SAMTOOLS, "idxstats", path])
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.strip() or f"idxstats failed: {path}")
    mapped = unmapped = 0
    for line in proc.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        n = int(parts[2])
        if parts[0] == "*":
            unmapped = n
        else:
            mapped += n
    return mapped, unmapped, mapped + unmapped


def parse_samtools_stats_depth(path: str, genome_size: int) -> float:
    proc = run_cmd([SAMTOOLS, "stats", path])
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.strip() or f"stats failed: {path}")
    mapped = avg_len = bases = 0
    for line in proc.stdout.splitlines():
        if not line.startswith("SN\t"):
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        key, val = parts[1], parts[2].split()[0]
        if key == "average length:":
            avg_len = float(val)
        elif key == "reads mapped:":
            mapped = int(val)
        elif key == "bases mapped (cigar):":
            bases = int(val)
    if bases > 0 and genome_size > 0:
        return bases / genome_size
    if mapped > 0 and avg_len > 0 and genome_size > 0:
        return mapped * avg_len / genome_size
    return 0.0


def bam_rates(flag: Dict[str, int]) -> Tuple[float, float, float]:
    total = flag.get("in total", 0)
    mapped = flag.get("mapped", 0)
    paired = flag.get("read1", 0) + flag.get("read2", 0)
    proper = flag.get("properly paired", 0)
    map_rate = mapped / total if total else 0.0
    pair_rate = proper / (paired / 2) if paired else 0.0
    return total, map_rate, pair_rate


def open_fq(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def read_head_pairs(path: str, max_pairs: int) -> Tuple[List[str], List[int]]:
    ids: List[str] = []
    lens: List[int] = []
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
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


def read_random_ids(path: str, max_n: int, seed: int, rate: float = 0.001) -> List[str]:
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
        rng = random.Random(seed)
        ids = rng.sample(ids, max_n)
    return ids


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


def sample_markers(panel: str, n: int, seed: int) -> List[dict]:
    rows: List[dict] = []
    with open(panel) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("breakpoint_ok", "1") != "1":
                continue
            if row.get("in_sv_interior", "0") == "1":
                continue
            rows.append(row)
    rng = random.Random(seed)
    if len(rows) <= n:
        return rows
    return rng.sample(rows, n)


def write_marker_bed(markers: List[dict], bed_path: str) -> None:
    with open(bed_path, "w") as fh:
        for m in markers:
            pos = int(m["mm_pos"])
            fh.write(f"{m['mm_chr']}\t{pos - 1}\t{pos}\t{m['marker_id']}\n")


def marker_pos_map(markers: List[dict]) -> Dict[Tuple[str, int], dict]:
    return {(m["mm_chr"], int(m["mm_pos"])): m for m in markers}


def clean_mpileup_bases(bases: str) -> str:
    # 去掉 read 起止标记 ^Q 与 $，以及 indel 块 +3AAA / -2NN
    bases = re.sub(r"\^.", "", bases)
    bases = bases.replace("$", "")
    bases = re.sub(r"[+-]\d+[A-Za-z]+", "", bases)
    return bases


def count_ref_alt(bases: str, alt: str) -> Tuple[int, int]:
    ref_n = bases.count(".") + bases.count(",")
    alt_u, alt_l = alt.upper(), alt.lower()
    alt_n = sum(1 for c in bases if c in (alt_u, alt_l))
    return ref_n, alt_n


def mpileup_af(bam: str, ref: str, bed: str, markers: List[dict]) -> Dict[str, Tuple[int, float]]:
    proc = run_cmd([
        SAMTOOLS, "mpileup", "-f", ref, "-l", bed, "-Q", "0", "--min-BQ", "0",
        "-d", "10000", bam,
    ])
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.strip() or "mpileup failed")
    pos_map = marker_pos_map(markers)
    out: Dict[str, Tuple[int, float]] = {}
    for line in proc.stdout.splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 5:
            continue
        chrom, pos_s, _ref, depth_s, bases_raw = fields[:5]
        pos = int(pos_s)
        m = pos_map.get((chrom, pos))
        if not m:
            continue
        depth = int(depth_s)
        bases = clean_mpileup_bases(bases_raw)
        ref_n, alt_n = count_ref_alt(bases, m["alt"])
        denom = ref_n + alt_n
        af = alt_n / denom if denom else 0.0
        out[m["marker_id"]] = (depth, af)
    return out


def classify_af(af: float, het_lo: float, het_hi: float) -> str:
    if af <= het_lo:
        return "hom_ref"
    if af >= het_hi:
        return "hom_alt"
    if het_lo < af < het_hi:
        return "het"
    return "unk"


def gam_stream_stats(gam: str, max_records: int) -> Dict[str, object]:
    proc = subprocess.Popen(
        [VG_BIN, "view", "-a", gam],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )
    n = mapped = proper = 0
    ident_sum = 0.0
    try:
        assert proc.stdout is not None
        for line in proc.stdout:
            if n >= max_records:
                break
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue
            n += 1
            mq = int(rec.get("mapping_quality", 0))
            if mq >= 1:
                mapped += 1
            ann = rec.get("annotation") or {}
            if ann.get("proper_pair"):
                proper += 1
            ident_sum += float(rec.get("identity", 0.0))
        proc.terminate()
        proc.wait(timeout=10)
    except Exception:
        proc.kill()
        raise
    map_rate = mapped / n if n else 0.0
    pair_rate = proper / n if n else 0.0
    mean_ident = ident_sum / n if n else 0.0
    return {
        "gam_records_sampled": n,
        "gam_map_rate": f"{map_rate:.4f}",
        "gam_proper_pair_rate": f"{pair_rate:.4f}",
        "gam_mean_identity": f"{mean_ident:.4f}",
    }


def discover_spot(n_f2: int, n_pf: int, spot_f2: int, spot_f1: int, seed: int) -> List[Tuple[str, str]]:
    rng = random.Random(seed)
    f2 = [f"simF2_{i:03d}" for i in range(1, n_f2 + 1)]
    f1 = [f"pseudoF1_{i:02d}" for i in range(1, n_pf + 1)]
    pick = [(s, "simF2") for s in (f2 if spot_f2 >= n_f2 else rng.sample(f2, min(spot_f2, len(f2))))]
    pick += [(s, "simF1") for s in (f1 if spot_f1 >= n_pf else rng.sample(f1, min(spot_f1, len(f1))))]
    return pick


def paths_for(cohort: str, sid: str, sim_f2: str, sim_f1: str, reads_dir: str, pf_depth: str) -> dict:
    root = sim_f2 if cohort == "simF2" else sim_f1
    p = {
        "bwa_bam": os.path.join(root, "bwa_dv", "01_bwa", f"{sid}.bam"),
        "gfa_bam": os.path.join(root, "gfa_dv", "01_bam", f"{sid}.bam"),
        "gam": os.path.join(root, "gfa_vgcall", "01_gam", f"{sid}.gam"),
    }
    if cohort == "simF2":
        p["r1"] = os.path.join(reads_dir, "simF2", sid, f"{sid}.R1.fq.gz")
        p["r2"] = os.path.join(reads_dir, "simF2", sid, f"{sid}.R2.fq.gz")
    else:
        p["r1"] = os.path.join(reads_dir, "pseudoF1", f"{pf_depth}x", f"{sid}.R1.fq.gz")
        p["r2"] = os.path.join(reads_dir, "pseudoF1", f"{pf_depth}x", f"{sid}.R2.fq.gz")
    return p


def depth_from_flagstat(flag: Dict[str, int], genome_size: int, read_len: int = 150) -> float:
    mapped = flag.get("mapped", 0)
    if mapped > 0 and genome_size > 0:
        return mapped * read_len / genome_size
    return 0.0


def cohort_inventory(
    cohort: str, n: int, fmt: str, sim_f2: str, sim_f1: str, min_bam_mb: int, min_gam_mb: int,
) -> List[dict]:
    root = sim_f2 if cohort == "simF2" else sim_f1
    rows = []
    for i in range(1, n + 1):
        sid = fmt % i
        bam = os.path.join(root, "bwa_dv", "01_bwa", f"{sid}.bam")
        gam = os.path.join(root, "gfa_vgcall", "01_gam", f"{sid}.gam")
        gfa_bam = os.path.join(root, "gfa_dv", "01_bam", f"{sid}.bam")
        bam_sz = os.path.getsize(bam) if os.path.isfile(bam) else 0
        gam_sz = os.path.getsize(gam) if os.path.isfile(gam) else 0
        gfa_sz = os.path.getsize(gfa_bam) if os.path.isfile(gfa_bam) else 0
        bam_ok = bam_sz > min_bam_mb * 1024 * 1024
        gam_ok = gam_sz > min_gam_mb * 1024 * 1024
        gfa_ok = gfa_sz > min_bam_mb * 1024 * 1024 if cohort == "simF1" else True
        rows.append({
            "sample": sid, "cohort": cohort,
            "bwa_bam": "PASS" if bam_ok else ("missing" if bam_sz == 0 else "FAIL"),
            "gam": "PASS" if gam_ok else ("missing" if gam_sz == 0 else "FAIL"),
            "gfa_bam": "PASS" if gfa_ok else ("missing" if gfa_sz == 0 else "NA"),
            "bwa_mb": f"{bam_sz/1024/1024:.1f}" if bam_sz else "0",
            "gam_mb": f"{gam_sz/1024/1024:.1f}" if gam_sz else "0",
        })
    return rows


def qc_one(
    sample: str, cohort: str, paths: dict, genome_size: int, expected_depth: float,
    ref_fa: str, markers: List[dict], marker_bed: str,
    max_fq_ids: int, max_gam_records: int,
    map_min: float, pair_min: float, depth_tol: float,
    prefix_tol: float, het_lo: float, het_hi: float,
    mode: str = "fast",
) -> Tuple[dict, List[dict]]:
    row: Dict[str, object] = {
        "sample": sample, "cohort": cohort, "mode": mode, "expected_depth": f"{expected_depth:.1f}",
        "status": "ok", "overall_pass": "FAIL",
    }
    af_rows: List[dict] = []
    do_deep = mode == "deep"

    # F1 混样：fast 跳过（QC1 已检）；deep 全文件随机抽
    if do_deep and cohort == "simF1" and os.path.isfile(paths["r1"]):
        try:
            ids = read_random_ids(paths["r1"], max_fq_ids, hash(sample) % 100000, rate=0.0005)
            mm_pct, ts_pct, bare = prefix_stats(ids)
            ratio = mm_pct / ts_pct if ts_pct > 0 else 0.0
            prefix_pass = len(ids) >= 1000 and bare == 0 and abs(ratio - 1.0) <= prefix_tol
            row.update({
                "prefix_n_sampled": len(ids),
                "prefix_mm_pct": f"{mm_pct:.3f}",
                "prefix_ts_pct": f"{ts_pct:.3f}",
                "prefix_mm_ts_ratio": f"{ratio:.3f}",
                "prefix_pass": "PASS" if prefix_pass else "FAIL",
            })
        except RuntimeError as e:
            row["prefix_error"] = str(e)
    elif cohort == "simF1":
        row["prefix_pass"] = "NA"

    # BWA BAM
    bam = paths["bwa_bam"]
    if os.path.isfile(bam) and os.path.getsize(bam) > 1000:
        try:
            if do_deep:
                fl = parse_flagstat(bam)
                total, map_r, pair_r = bam_rates(fl)
                depth = parse_samtools_stats_depth(bam, genome_size)
                row["bwa_depth_method"] = "stats"
            else:
                mapped, unmapped, total = parse_idxstats(bam)
                map_r = mapped / total if total else 0.0
                pair_r = map_r
                depth = mapped * 150 / genome_size if genome_size else 0.0
                row["bwa_depth_method"] = "idxstats_est"
            row.update({
                "bwa_total_reads": total if do_deep else mapped + unmapped,
                "bwa_map_rate": f"{map_r:.4f}",
                "bwa_pair_rate": f"{pair_r:.4f}",
                "bwa_depth": f"{depth:.3f}",
                "bwa_map_pass": "PASS" if map_r >= map_min else "FAIL",
                "bwa_pair_pass": "PASS" if pair_r >= pair_min else "FAIL",
                "bwa_depth_pass": "PASS" if abs(depth / expected_depth - 1.0) <= depth_tol else "FAIL",
            })
            if do_deep:
                af_bwa = mpileup_af(bam, ref_fa, marker_bed, markers)
                for m in markers:
                    mid = m["marker_id"]
                    dp, af = af_bwa.get(mid, (0, float("nan")))
                    gt = classify_af(af, het_lo, het_hi) if dp > 0 else "nodata"
                    af_rows.append({
                        "sample": sample, "cohort": cohort, "source": "bwa_bam",
                        "marker_id": mid, "chr": m["mm_chr"], "pos": m["mm_pos"],
                        "ref": m["ref"], "alt": m["alt"], "depth": dp,
                        "af": f"{af:.3f}" if dp > 0 else "NA",
                        "gt_class": gt,
                    })
        except RuntimeError as e:
            row["bwa_error"] = str(e)
    else:
        row["bwa_status"] = "missing"

    # GFA surject BAM（simF1）
    gfa_bam = paths["gfa_bam"]
    if do_deep and cohort == "simF1" and os.path.isfile(gfa_bam) and os.path.getsize(gfa_bam) > 1000:
        try:
            fl = parse_flagstat(gfa_bam)
            total, map_r, pair_r = bam_rates(fl)
            depth = parse_samtools_stats_depth(gfa_bam, genome_size)
            row.update({
                "gfa_bam_map_rate": f"{map_r:.4f}",
                "gfa_bam_pair_rate": f"{pair_r:.4f}",
                "gfa_bam_depth": f"{depth:.3f}",
            })
            af_gfa = mpileup_af(gfa_bam, ref_fa, marker_bed, markers)
            for m in markers:
                mid = m["marker_id"]
                dp, af = af_gfa.get(mid, (0, float("nan")))
                gt = classify_af(af, het_lo, het_hi) if dp > 0 else "nodata"
                af_rows.append({
                    "sample": sample, "cohort": cohort, "source": "gfa_surject_bam",
                    "marker_id": mid, "chr": m["mm_chr"], "pos": m["mm_pos"],
                    "ref": m["ref"], "alt": m["alt"], "depth": dp,
                    "af": f"{af:.3f}" if dp > 0 else "NA",
                    "gt_class": gt,
                })
        except RuntimeError as e:
            row["gfa_bam_error"] = str(e)
    elif cohort == "simF1" and os.path.isfile(gfa_bam) and os.path.getsize(gfa_bam) > 1000 and not do_deep:
        try:
            mapped, _, total = parse_idxstats(gfa_bam)
            map_r = mapped / total if total else 0.0
            row.update({
                "gfa_bam_map_rate": f"{map_r:.4f}",
                "gfa_bam_pair_rate": f"{map_r:.4f}",
                "gfa_bam_depth": f"{mapped * 150 / genome_size:.3f}" if genome_size else "0",
            })
        except RuntimeError as e:
            row["gfa_bam_error"] = str(e)

    # GAM：流式抽查（避免 vg stats -a OOM）
    gam = paths["gam"]
    if os.path.isfile(gam) and os.path.getsize(gam) > 1000:
        try:
            gs = gam_stream_stats(gam, max_gam_records)
            row.update(gs)
            row["gam_map_pass"] = "PASS" if float(gs["gam_map_rate"]) >= map_min else "FAIL"
            row["gam_pair_pass"] = "PASS" if float(gs["gam_proper_pair_rate"]) >= pair_min else "FAIL"
            est_depth = os.path.getsize(gam) / max(genome_size, 1) * 0.15
            row["gam_depth_est"] = f"{est_depth:.3f}"
        except Exception as e:
            row["gam_error"] = str(e)
    else:
        row["gam_status"] = "missing"

    # AF 分布汇总
    bwa_afs = [r for r in af_rows if r["source"] == "bwa_bam" and r["gt_class"] != "nodata"]
    if not do_deep:
        row["af_pattern_pass"] = "NA"
    elif bwa_afs:
        cls = Counter(r["gt_class"] for r in bwa_afs)
        row["af_hom_ref"] = cls.get("hom_ref", 0)
        row["af_het"] = cls.get("het", 0)
        row["af_hom_alt"] = cls.get("hom_alt", 0)
        if cohort == "simF1":
            het_frac = cls.get("het", 0) / len(bwa_afs)
            row["af_pattern_pass"] = "PASS" if het_frac >= 0.85 else "FAIL"
        else:
            has_three = sum(1 for k in ("hom_ref", "het", "hom_alt") if cls.get(k, 0) > 0)
            row["af_pattern_pass"] = "PASS" if has_three >= 2 else "FAIL"

    has_bam = os.path.isfile(bam) and os.path.getsize(bam) > 1000
    has_gam = os.path.isfile(gam) and os.path.getsize(gam) > 1000
    if not has_bam and not has_gam:
        row["align_status"] = "pending"
        if row.get("prefix_pass") == "PASS":
            row["overall_pass"] = "NA"
        elif "prefix_pass" in row:
            row["overall_pass"] = "FAIL"
        else:
            row["overall_pass"] = "NA"
        return row, af_rows

    passes = []
    for k in ("bwa_map_pass", "bwa_pair_pass", "bwa_depth_pass", "gam_map_pass", "gam_pair_pass",
              "prefix_pass", "af_pattern_pass"):
        if k in row and row[k] not in ("NA", ""):
            passes.append(row[k] == "PASS")
    row["overall_pass"] = "PASS" if passes and all(passes) else "FAIL"
    return row, af_rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", required=True)
    ap.add_argument("--ref-fa", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-af-tsv", required=True)
    ap.add_argument("--sim-f2-root", required=True)
    ap.add_argument("--sim-f1-root", required=True)
    ap.add_argument("--reads-dir", required=True)
    ap.add_argument("--genome-size", type=int, required=True)
    ap.add_argument("--art-haplo-depth", type=float, default=5.0)
    ap.add_argument("--pf-depth", type=float, default=10.0)
    ap.add_argument("--pf-depth-label", default="10")
    ap.add_argument("--n-f2", type=int, default=50)
    ap.add_argument("--n-pf", type=int, default=10)
    ap.add_argument("--spot-f2", type=int, default=3)
    ap.add_argument("--spot-f1", type=int, default=3)
    ap.add_argument("--n-markers", type=int, default=40)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--max-fq-ids", type=int, default=50000)
    ap.add_argument("--max-gam-records", type=int, default=30000)
    ap.add_argument("--map-min", type=float, default=0.95)
    ap.add_argument("--pair-min", type=float, default=0.90)
    ap.add_argument("--depth-tol", type=float, default=0.25)
    ap.add_argument("--prefix-tol", type=float, default=0.15)
    ap.add_argument("--het-lo", type=float, default=0.25)
    ap.add_argument("--het-hi", type=float, default=0.75)
    ap.add_argument("--mode", choices=("fast", "deep"), default="fast")
    ap.add_argument("--inventory-out", default="")
    ap.add_argument("--min-bam-mb", type=int, default=500)
    ap.add_argument("--min-gam-mb", type=int, default=500)
    ap.add_argument("--samples", default="")
    ap.add_argument("--tmp-dir", required=True)
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out_tsv) or ".", exist_ok=True)
    os.makedirs(args.tmp_dir, exist_ok=True)

    markers = sample_markers(args.panel, args.n_markers, args.seed) if args.mode == "deep" else []
    marker_bed = os.path.join(args.tmp_dir, "qc2_markers.bed")
    if markers:
        write_marker_bed(markers, marker_bed)

    inv_out = args.inventory_out or args.out_tsv.replace("align_summary", "align_inventory")
    inv_rows = cohort_inventory("simF2", args.n_f2, "simF2_%03d", args.sim_f2_root, args.sim_f1_root,
                                args.min_bam_mb, args.min_gam_mb)
    inv_rows += cohort_inventory("simF1", args.n_pf, "pseudoF1_%02d", args.sim_f2_root, args.sim_f1_root,
                                 args.min_bam_mb, args.min_gam_mb)
    with open(inv_out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(inv_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(inv_rows)

    if args.samples.strip():
        samples = []
        for s in args.samples.split(","):
            s = s.strip()
            if not s:
                continue
            cohort = "simF1" if s.startswith("pseudoF1") else "simF2"
            samples.append((s, cohort))
    else:
        samples = discover_spot(args.n_f2, args.n_pf, args.spot_f2, args.spot_f1, args.seed)

    expected_f2 = 2.0 * args.art_haplo_depth
    summary_rows: List[dict] = []
    all_af: List[dict] = []

    for sid, cohort in samples:
        paths = paths_for(cohort, sid, args.sim_f2_root, args.sim_f1_root,
                          args.reads_dir, args.pf_depth_label)
        exp = expected_f2 if cohort == "simF2" else args.pf_depth
        row, af_rows = qc_one(
            sid, cohort, paths, args.genome_size, exp, args.ref_fa, markers, marker_bed,
            args.max_fq_ids, args.max_gam_records,
            args.map_min, args.pair_min, args.depth_tol, args.prefix_tol,
            args.het_lo, args.het_hi, mode=args.mode,
        )
        summary_rows.append(row)
        all_af.extend(af_rows)

    sum_fields = [
        "sample", "cohort", "mode", "status", "expected_depth",
        "bwa_map_rate", "bwa_pair_rate", "bwa_depth", "bwa_depth_method",
        "bwa_map_pass", "bwa_pair_pass", "bwa_depth_pass",
        "gfa_bam_map_rate", "gfa_bam_pair_rate", "gfa_bam_depth",
        "gam_map_rate", "gam_proper_pair_rate", "gam_mean_identity", "gam_map_pass", "gam_pair_pass",
        "prefix_mm_pct", "prefix_ts_pct", "prefix_mm_ts_ratio", "prefix_n_sampled", "prefix_pass",
        "af_hom_ref", "af_het", "af_hom_alt", "af_pattern_pass", "overall_pass", "align_status",
        "bwa_status", "gam_status", "bwa_error", "gfa_bam_error", "gam_error",
    ]
    with open(args.out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=sum_fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(summary_rows)

    af_fields = ["sample", "cohort", "source", "marker_id", "chr", "pos", "ref", "alt", "depth", "af", "gt_class"]
    with open(args.out_af_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=af_fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(all_af)

    n_pass = sum(1 for r in summary_rows if r.get("overall_pass") == "PASS")
    n_na = sum(1 for r in summary_rows if r.get("overall_pass") == "NA")
    n_fail = sum(1 for r in summary_rows if r.get("overall_pass") == "FAIL")
    print(
        f"QC2 mode={args.mode} spot={len(summary_rows)} pass={n_pass} na={n_na} fail={n_fail} "
        f"inventory={len(inv_rows)} out={args.out_tsv} inv={inv_out}",
        file=sys.stderr,
    )
    sys.exit(0 if n_fail == 0 else 1)


if __name__ == "__main__":
    main()
