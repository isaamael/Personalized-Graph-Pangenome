#!/usr/bin/env python3
"""step01: SyRI SNP marker + SV 块表 + 遗传图（SV 来源: TS_vs_MM/03_stats）。"""

from __future__ import annotations

import argparse
import bisect
import csv
import os
import sys
from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple

SKIP_TYPES = {"SNP", "SYN", "SYNAL"}


@dataclass
class SvBlock:
    sv_id: str
    typ: str
    mm_chr: str
    mm_start: int
    mm_end: int
    ts_chr: str
    ts_start: int
    ts_end: int
    parent: str
    length: int
    dest_mm_chr: str = "-"
    dest_mm_start: int = -1
    dest_mm_end: int = -1
    dest_ts_chr: str = "-"
    dest_ts_start: int = -1
    dest_ts_end: int = -1


@dataclass
class ChrMergedSv:
    starts: List[int]
    ends: List[int]


def parse_int(x: str) -> int:
    if x in ("", "-", "."):
        return -1
    return int(x)


def load_chr_lengths(fai: str) -> Dict[str, int]:
    out: Dict[str, int] = {}
    with open(fai) as fh:
        for line in fh:
            c, ln, *_ = line.split()
            out[c] = int(ln)
    return out


def load_sv_from_stats_tsv(path: str, min_sv_bp: int) -> Tuple[List[SvBlock], Dict[str, SvBlock]]:
    blocks: List[SvBlock] = []
    by_id: Dict[str, SvBlock] = {}
    trans_primary: Dict[str, SvBlock] = {}

    with open(path) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            typ = row["syri_type"]
            if typ in SKIP_TYPES:
                continue
            if typ.endswith("AL") and typ != "TRANSAL":
                continue
            ln = int(row["length_bp"])
            if ln < min_sv_bp:
                continue

            mm_chr = row["mm_chr"] if row["mm_chr"] not in ("", "-") else "-"
            mm_s, mm_e = parse_int(row["mm_start"]), parse_int(row["mm_end"])
            ts_chr = row["ts_chr"] if row["ts_chr"] not in ("", "-", ".") else "-"
            ts_s, ts_e = parse_int(row["ts_start"]), parse_int(row["ts_end"])
            sv_id = row["id"]
            parent = row["parent"] if row["parent"] not in ("", "-") else "-"

            blk = SvBlock(
                sv_id=sv_id, typ=typ,
                mm_chr=mm_chr, mm_start=mm_s, mm_end=mm_e,
                ts_chr=ts_chr, ts_start=ts_s, ts_end=ts_e,
                parent=parent, length=ln,
            )
            if typ == "TRANS":
                trans_primary[sv_id] = blk
            if typ == "TRANSAL" and parent != "-":
                if parent in trans_primary:
                    p = trans_primary[parent]
                    p.dest_mm_chr = mm_chr
                    p.dest_mm_start = mm_s
                    p.dest_mm_end = mm_e
                    p.dest_ts_chr = ts_chr
                    p.dest_ts_start = ts_s
                    p.dest_ts_end = ts_e
                    continue
            blocks.append(blk)
            by_id[sv_id] = blk

    return blocks, by_id


def merge_mm_intervals(blocks: List[SvBlock]) -> Dict[str, ChrMergedSv]:
    raw: Dict[str, List[Tuple[int, int]]] = {}
    for b in blocks:
        if b.mm_chr == "-" or b.mm_start < 0 or b.mm_end < b.mm_start:
            continue
        raw.setdefault(b.mm_chr, []).append((b.mm_start, b.mm_end))

    merged: Dict[str, ChrMergedSv] = {}
    for chrom, ivs in raw.items():
        ivs.sort()
        starts: List[int] = []
        ends: List[int] = []
        for lo, hi in ivs:
            if not starts or lo > ends[-1] + 1:
                starts.append(lo)
                ends.append(hi)
            else:
                ends[-1] = max(ends[-1], hi)
        merged[chrom] = ChrMergedSv(starts=starts, ends=ends)
    return merged


def inside_sv(mm_chr: str, mm_pos: int, merged: Dict[str, ChrMergedSv]) -> bool:
    iv = merged.get(mm_chr)
    if not iv or not iv.starts:
        return False
    idx = bisect.bisect_right(iv.starts, mm_pos) - 1
    if idx < 0:
        return False
    return iv.starts[idx] <= mm_pos <= iv.ends[idx]


def load_markers(snp_tsv: str, merged: Dict[str, ChrMergedSv], chrs: Set[str]) -> List[dict]:
    markers: List[dict] = []
    with open(snp_tsv) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            mm_chr = row["mm_chr"]
            if mm_chr not in chrs:
                continue
            mm_pos = int(row["mm_pos"])
            in_sv = inside_sv(mm_chr, mm_pos, merged)
            markers.append({
                "marker_id": row["id"],
                "mm_chr": mm_chr,
                "mm_pos": mm_pos,
                "ts_chr": row["ts_chr"],
                "ts_pos": int(row["ts_pos"]),
                "ref": row["ref_allele"],
                "alt": row["alt_allele"],
                "in_sv_interior": int(in_sv),
                "breakpoint_ok": int(not in_sv),
            })
    markers.sort(key=lambda m: (m["mm_chr"], m["mm_pos"]))
    return markers


def write_genetic_map(fai: str, out_path: str) -> None:
    lens = load_chr_lengths(fai)
    with open(out_path, "w") as fh:
        fh.write("chr\tposition\trate_cM_per_Mb\n")
        for c in sorted(lens, key=lambda x: int(x.replace("chr", ""))):
            fh.write(f"{c}\t{lens[c]}\t1.0\n")


def thin_markers_for_sim(markers: List[dict], bin_size: int, keep: int) -> List[dict]:
    from collections import defaultdict

    bins: dict = defaultdict(list)
    for m in markers:
        if not m["breakpoint_ok"]:
            continue
        bins[(m["mm_chr"], m["mm_pos"] // bin_size)].append(m)

    result = []
    for (_, _), ms in sorted(bins.items()):
        center = (ms[0]["mm_pos"] // bin_size) * bin_size + bin_size // 2
        ms_sorted = sorted(ms, key=lambda m: abs(m["mm_pos"] - center))
        result.extend(ms_sorted[:keep])
    return sorted(result, key=lambda m: (m["mm_chr"], m["mm_pos"]))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--snp-tsv", required=True)
    ap.add_argument("--syri-sv-tsv", required=True)
    ap.add_argument("--mm-fai", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--min-sv-bp", type=int, default=1000)
    ap.add_argument("--pilot-chr", default="")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    chrs = {f"chr{i}" for i in range(1, 13)}
    if args.pilot_chr:
        chrs = {args.pilot_chr}

    blocks, _ = load_sv_from_stats_tsv(args.syri_sv_tsv, args.min_sv_bp)
    blocks = [b for b in blocks if b.mm_chr in chrs or b.ts_chr in chrs]
    merged = merge_mm_intervals(blocks)
    n_merged = sum(len(v.starts) for v in merged.values())

    markers = load_markers(args.snp_tsv, merged, chrs)

    mpath = os.path.join(args.out_dir, "markers.tsv")
    with open(mpath, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=[
            "marker_id", "mm_chr", "mm_pos", "ts_chr", "ts_pos",
            "ref", "alt", "in_sv_interior", "breakpoint_ok",
        ], delimiter="\t")
        w.writeheader()
        w.writerows(markers)

    sim_bin = int(os.environ.get("SIM_BIN_SIZE", 100_000))
    sim_keep = int(os.environ.get("SIM_KEEP_PER_BIN", 2))
    sim_markers = thin_markers_for_sim(markers, sim_bin, sim_keep)
    sim_path = os.path.join(args.out_dir, "markers_sim.tsv")
    with open(sim_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=[
            "marker_id", "mm_chr", "mm_pos", "ts_chr", "ts_pos",
            "ref", "alt", "in_sv_interior", "breakpoint_ok",
        ], delimiter="\t")
        w.writeheader()
        w.writerows(sim_markers)

    spath = os.path.join(args.out_dir, "sv_blocks.tsv")
    with open(spath, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([
            "sv_id", "type", "mm_chr", "mm_start", "mm_end",
            "ts_chr", "ts_start", "ts_end", "parent", "length",
            "dest_mm_chr", "dest_mm_start", "dest_mm_end",
            "dest_ts_chr", "dest_ts_start", "dest_ts_end",
        ])
        for b in blocks:
            w.writerow([
                b.sv_id, b.typ, b.mm_chr, b.mm_start, b.mm_end,
                b.ts_chr, b.ts_start, b.ts_end, b.parent, b.length,
                b.dest_mm_chr, b.dest_mm_start, b.dest_mm_end,
                b.dest_ts_chr, b.dest_ts_start, b.dest_ts_end,
            ])

    mpath_iv = os.path.join(args.out_dir, "sv_merged_intervals.tsv")
    with open(mpath_iv, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["mm_chr", "mm_start", "mm_end"])
        for chrom in sorted(merged, key=lambda c: int(c.replace("chr", ""))):
            iv = merged[chrom]
            for s, e in zip(iv.starts, iv.ends):
                w.writerow([chrom, s, e])

    write_genetic_map(args.mm_fai, os.path.join(args.out_dir, "genetic_map.txt"))

    bp_ok = sum(1 for m in markers if m["breakpoint_ok"])
    print(
        f"markers_full={len(markers)} breakpoint_ok={bp_ok} "
        f"markers_sim={len(sim_markers)} sv_blocks={len(blocks)} "
        f"sv_merged_intervals={n_merged} min_sv_bp={args.min_sv_bp}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
