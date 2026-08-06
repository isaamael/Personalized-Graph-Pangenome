#!/usr/bin/env python3
"""Build simF1/simF2 haplotype-block truth sets and marker projections."""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import os
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Tuple


@dataclass
class Marker:
    marker_id: str
    mm_chr: str
    mm_pos: int
    ts_chr: str
    ts_pos: int


def load_markers(path: str) -> Dict[str, List[Marker]]:
    by_chr: Dict[str, List[Marker]] = defaultdict(list)
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["mm_chr"].startswith("chr0"):
                continue
            by_chr[row["mm_chr"]].append(Marker(
                marker_id=row["marker_id"],
                mm_chr=row["mm_chr"],
                mm_pos=int(row["mm_pos"]),
                ts_chr=row["ts_chr"],
                ts_pos=int(row["ts_pos"]),
            ))
    for c in by_chr:
        by_chr[c].sort(key=lambda m: m.mm_pos)
    return by_chr


def load_f1_samples(path: str) -> List[str]:
    with open(path) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        if "sample" in (r.fieldnames or []):
            return [row["sample"].strip() for row in r if row["sample"].strip()]
    return load_sample_list(path)


def discover_f2_samples(truth_dir: str) -> List[str]:
    sids = []
    for fn in sorted(os.listdir(truth_dir)):
        if fn.endswith(".gameteA.fragments.tsv"):
            sids.append(fn.replace(".gameteA.fragments.tsv", ""))
    return sids


def load_merged_blocks(path: str) -> Dict[str, List[dict]]:
    by_chr: Dict[str, List[dict]] = defaultdict(list)
    if not os.path.isfile(path):
        return by_chr
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    rows.sort(key=lambda r: (r["mm_chr"], int(r["frag_idx"])))
    for mm_chr, grp in _groupby_chr(rows):
        cur = None
        for r in grp:
            parent = r["parent"]
            mm_lo, mm_hi = int(r["mm_lo"]), int(r["mm_hi"])
            ts_lo, ts_hi = int(r["ts_lo"]), int(r["ts_hi"])
            if cur and cur["parent"] == parent:
                cur["mm_hi"] = max(cur["mm_hi"], mm_hi)
                cur["ts_hi"] = max(cur["ts_hi"], ts_hi)
                cur["marker_R"] = r["marker_R"]
            else:
                if cur:
                    by_chr[mm_chr].append(cur)
                cur = {
                    "parent": parent,
                    "mm_lo": mm_lo, "mm_hi": mm_hi,
                    "ts_lo": ts_lo, "ts_hi": ts_hi,
                    "marker_L": r["marker_L"], "marker_R": r["marker_R"],
                }
        if cur:
            by_chr[mm_chr].append(cur)
    return by_chr


def _groupby_chr(rows):
    cur_chr = None
    buf = []
    for r in rows:
        if r["mm_chr"] != cur_chr:
            if buf:
                yield cur_chr, buf
            cur_chr = r["mm_chr"]
            buf = [r]
        else:
            buf.append(r)
    if buf:
        yield cur_chr, buf


def parent_at(blocks: List[dict], pos: int) -> str:
    if not blocks:
        return "MM"
    lo = [b["mm_lo"] for b in blocks]
    i = bisect.bisect_right(lo, pos) - 1
    if i < 0:
        i = 0
    b = blocks[i]
    if pos > b["mm_hi"] and i + 1 < len(blocks):
        return blocks[i + 1]["parent"]
    return b["parent"]


def gt_code(hap_a: str, hap_b: str) -> int:
    a = 0 if hap_a == "MM" else 2
    b = 0 if hap_b == "MM" else 2
    return a if a == b else 1


def build_f2_blocks(
    sample: str,
    truth_dir: str,
    markers_by_chr: Dict[str, List[Marker]],
) -> Tuple[List[dict], List[dict]]:
    ba = load_merged_blocks(os.path.join(truth_dir, f"{sample}.gameteA.fragments.tsv"))
    bb = load_merged_blocks(os.path.join(truth_dir, f"{sample}.gameteB.fragments.tsv"))
    block_rows: List[dict] = []
    marker_rows: List[dict] = []
    for mm_chr in sorted(markers_by_chr.keys(), key=lambda c: int(c.replace("chr", ""))):
        mk = markers_by_chr[mm_chr]
        if not mk:
            continue
        cur = None
        for m in mk:
            pa = parent_at(ba.get(mm_chr, []), m.mm_pos)
            pb = parent_at(bb.get(mm_chr, []), m.mm_pos)
            gt = gt_code(pa, pb)
            marker_rows.append({
                "sample": sample, "marker_id": m.marker_id,
                "gt": gt, "hap_A": pa, "hap_B": pb,
                "mm_chr": mm_chr, "mm_pos": m.mm_pos,
            })
            if cur and cur["hap_A"] == pa and cur["hap_B"] == pb:
                cur["mm_end"] = m.mm_pos
                cur["ts_end"] = m.ts_pos
                cur["marker_R"] = m.marker_id
            else:
                if cur:
                    block_rows.append(cur)
                cur = {
                    "sample": sample, "chr": mm_chr,
                    "mm_start": m.mm_pos, "mm_end": m.mm_pos,
                    "ts_start": m.ts_pos, "ts_end": m.ts_pos,
                    "hap_A": pa, "hap_B": pb, "gt": gt,
                    "marker_L": m.marker_id, "marker_R": m.marker_id,
                }
        if cur:
            block_rows.append(cur)
    return block_rows, marker_rows


def build_f1_blocks(
    samples: List[str],
    markers_by_chr: Dict[str, List[Marker]],
) -> Tuple[List[dict], List[dict]]:
    block_rows: List[dict] = []
    marker_rows: List[dict] = []
    chrs = sorted(markers_by_chr.keys(), key=lambda c: int(c.replace("chr", "")))
    for sample in samples:
        for mm_chr in chrs:
            mk = markers_by_chr[mm_chr]
            if not mk:
                continue
            block_rows.append({
                "sample": sample, "chr": mm_chr,
                "mm_start": mk[0].mm_pos, "mm_end": mk[-1].mm_pos,
                "ts_start": mk[0].ts_pos, "ts_end": mk[-1].ts_pos,
                "hap_A": "TS", "hap_B": "MM", "gt": 1,
                "marker_L": mk[0].marker_id, "marker_R": mk[-1].marker_id,
            })
            for m in mk:
                marker_rows.append({
                    "sample": sample, "marker_id": m.marker_id,
                    "gt": 1, "hap_A": "TS", "hap_B": "MM",
                    "mm_chr": mm_chr, "mm_pos": m.mm_pos,
                })
    return block_rows, marker_rows


def write_blocks(path: str, rows: List[dict]) -> None:
    cols = [
        "sample", "chr", "mm_start", "mm_end", "ts_start", "ts_end",
        "hap_A", "hap_B", "gt", "marker_L", "marker_R",
    ]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def write_markers_by_chr(out_dir: str, rows: List[dict]) -> None:
    os.makedirs(out_dir, exist_ok=True)
    by_chr: Dict[str, List[dict]] = defaultdict(list)
    for r in rows:
        by_chr[r["mm_chr"]].append(r)
    cols = ["sample", "marker_id", "gt", "hap_A", "hap_B", "mm_chr", "mm_pos"]
    for mm_chr in sorted(by_chr.keys(), key=lambda c: int(c.replace("chr", ""))):
        path = os.path.join(out_dir, f"{mm_chr}.tsv.gz")
        with gzip.open(path, "wt", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
            w.writeheader()
            w.writerows(by_chr[mm_chr])


def write_summary(path: str, f1_n: int, f2_n: int, f2_blocks: List[dict], f1_blocks: List[dict]) -> None:
    f2_gt = defaultdict(int)
    for r in f2_blocks:
        f2_gt[r["gt"]] += 1
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["cohort", "n_samples", "n_blocks", "gt0_blocks", "gt1_blocks", "gt2_blocks"])
        w.writerow(["simF1", f1_n, len(f1_blocks), 0, len(f1_blocks), 0])
        w.writerow([
            "simF2", f2_n, len(f2_blocks),
            f2_gt.get(0, 0), f2_gt.get(1, 0), f2_gt.get(2, 0),
        ])


def link_if_missing(src: str, dst: str) -> None:
    if not os.path.isfile(src):
        return
    if os.path.lexists(dst):
        os.remove(dst)
    os.symlink(os.path.basename(src), dst)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth-dir", required=True)
    ap.add_argument("--panel", required=True)
    ap.add_argument("--f1-list", required=True)
    args = ap.parse_args()

    truth_dir = args.truth_dir
    markers_by_chr = load_markers(args.panel)
    f1_samples = load_f1_samples(args.f1_list)
    f2_samples = discover_f2_samples(truth_dir)
    if not f2_samples:
        raise SystemExit(f"no F2 fragments in {truth_dir}")

    f2_block_rows: List[dict] = []
    f2_marker_rows: List[dict] = []
    for sid in f2_samples:
        br, mr = build_f2_blocks(sid, truth_dir, markers_by_chr)
        f2_block_rows.extend(br)
        f2_marker_rows.extend(mr)

    f1_block_rows, f1_marker_rows = build_f1_blocks(f1_samples, markers_by_chr)

    f2_blocks_path = os.path.join(truth_dir, "simF2_hap_blocks.tsv")
    f1_blocks_path = os.path.join(truth_dir, "simF1_hap_blocks.tsv")
    write_blocks(f2_blocks_path, f2_block_rows)
    write_blocks(f1_blocks_path, f1_block_rows)

    marker_dir = os.path.join(truth_dir, "truth_markers_by_chr")
    write_markers_by_chr(marker_dir, f1_marker_rows + f2_marker_rows)

    summary_path = os.path.join(truth_dir, "sim_truth_summary.tsv")
    write_summary(summary_path, len(f1_samples), len(f2_samples), f2_block_rows, f1_block_rows)

    link_if_missing(
        os.path.join(truth_dir, "truth_co_detail.tsv"),
        os.path.join(truth_dir, "simF2_co_detail.tsv"),
    )
    link_if_missing(
        os.path.join(truth_dir, "truth_co_summary.tsv"),
        os.path.join(truth_dir, "simF2_co_summary.tsv"),
    )
    link_if_missing(
        os.path.join(truth_dir, "truth_breakpoints.bed"),
        os.path.join(truth_dir, "simF2_co_breakpoints_MM.bed"),
    )

    print(f"simF1 blocks: {f1_blocks_path} n={len(f1_block_rows)} samples={len(f1_samples)}")
    print(f"simF2 blocks: {f2_blocks_path} n={len(f2_block_rows)} samples={len(f2_samples)}")
    print(f"markers: {marker_dir} n={len(f1_marker_rows) + len(f2_marker_rows)}")
    print(f"summary: {summary_path}")


if __name__ == "__main__":
    main()
