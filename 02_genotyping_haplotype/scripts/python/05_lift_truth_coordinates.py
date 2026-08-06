#!/usr/bin/env python3
"""Lift MM-coordinate truth sets and marker panels to SL6 with a PAF alignment."""
import argparse
import csv
import gzip
import os
import sys
from collections import defaultdict


def norm_sl6_chr(chrom):
    c = str(chrom).strip()
    if c.startswith("chr") and not c.startswith("chr0"):
        n = c.replace("chr", "")
        if n.isdigit():
            return f"chr{n.zfill(2)}"
    if c.startswith("chr0"):
        return c
    return c


def lift_point(lo, mm_chr, mm_pos):
    c, p = lo.convert(mm_chr, mm_pos, "mm", "sl6")
    if c is None or p is None:
        return None, None
    return norm_sl6_chr(c), int(p)


def lift_interval(lo, mm_chr, start, end):
    c1, p1 = lift_point(lo, mm_chr, start)
    c2, p2 = lift_point(lo, mm_chr, end)
    if c1 is None or c2 is None or c1 != c2:
        return None
    if p1 > p2:
        p1, p2 = p2, p1
    if p2 <= p1:
        return None
    return c1, p1, p2


def lift_markers_panel(lo, panel_path, out_path):
    rows_out = []
    n_fail = 0
    with open(panel_path) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        cols = r.fieldnames
        for row in r:
            mc, mp = row["mm_chr"], int(row["mm_pos"])
            sc, sp = lift_point(lo, mc, mp)
            if sc is None:
                n_fail += 1
                continue
            rows_out.append({
                "marker_id": row["marker_id"],
                "mm_chr": mc,
                "mm_pos": mp,
                "sl6_chr": sc,
                "sl6_pos": sp,
            })
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["marker_id", "mm_chr", "mm_pos", "sl6_chr", "sl6_pos"], delimiter="\t")
        w.writeheader()
        w.writerows(rows_out)
    return len(rows_out), n_fail


def lift_truth_markers(lo, mm_dir, sl6_dir):
    os.makedirs(sl6_dir, exist_ok=True)
    stats = {}
    mid_map = {}
    for fn in sorted(os.listdir(mm_dir)):
        if not fn.endswith(".tsv.gz"):
            continue
        mm_chr = fn.replace(".tsv.gz", "")
        in_path = os.path.join(mm_dir, fn)
        out_path = os.path.join(sl6_dir, fn)
        rows_out = []
        n_in = n_fail = 0
        with gzip.open(in_path, "rt") as fh:
            r = csv.DictReader(fh, delimiter="\t")
            for row in r:
                n_in += 1
                mc = row["mm_chr"]
                mp = int(row["mm_pos"])
                sc, sp = lift_point(lo, mc, mp)
                if sc is None:
                    n_fail += 1
                    continue
                rows_out.append({
                    "sample": row["sample"],
                    "marker_id": row["marker_id"],
                    "gt": row["gt"],
                    "hap_A": row["hap_A"],
                    "hap_B": row["hap_B"],
                    "sl6_chr": sc,
                    "sl6_pos": sp,
                })
                mid_map[row["marker_id"]] = (mc, mp, sc, sp)
        with gzip.open(out_path, "wt", newline="") as fh:
            w = csv.DictWriter(
                fh,
                fieldnames=["sample", "marker_id", "gt", "hap_A", "hap_B", "sl6_chr", "sl6_pos"],
                delimiter="\t",
            )
            w.writeheader()
            w.writerows(rows_out)
        stats[mm_chr] = (n_in, len(rows_out), n_fail)
    return stats, mid_map


def lift_hap_blocks(lo, in_path, out_path):
    rows_out = []
    n_in = n_fail = 0
    with open(in_path) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            n_in += 1
            hit = lift_interval(lo, row["chr"], int(row["mm_start"]), int(row["mm_end"]))
            if hit is None:
                n_fail += 1
                continue
            sc, ss, se = hit
            rows_out.append({
                "sample": row["sample"],
                "chr": sc,
                "sl6_start": ss,
                "sl6_end": se,
                "ts_start": row["ts_start"],
                "ts_end": row["ts_end"],
                "hap_A": row["hap_A"],
                "hap_B": row["hap_B"],
                "gt": row["gt"],
                "marker_L": row["marker_L"],
                "marker_R": row["marker_R"],
            })
    cols = [
        "sample", "chr", "sl6_start", "sl6_end", "ts_start", "ts_end",
        "hap_A", "hap_B", "gt", "marker_L", "marker_R",
    ]
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows_out)
    return n_in, len(rows_out), n_fail


def lift_breakpoints(lo, in_path, out_path):
    rows_out = []
    n_in = n_fail = 0
    with open(in_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            p = line.strip().split("\t")
            if len(p) < 3:
                continue
            n_in += 1
            mc, s, e = p[0], int(p[1]), int(p[2])
            hit = lift_interval(lo, mc, s, e) if e > s else None
            if hit is None:
                sc, sp = lift_point(lo, mc, s)
                if sc is None:
                    n_fail += 1
                    continue
                hit = (sc, sp, sp + 1)
            sc, ss, se = hit
            rest = p[3:] if len(p) > 3 else []
            rows_out.append("\t".join([sc, str(ss), str(se)] + rest))
    with open(out_path, "w") as fh:
        fh.write("# ref=SL6\n")
        fh.writelines(x + "\n" for x in rows_out)
    return n_in, len(rows_out), n_fail


def link_if_missing(src, dst):
    if not os.path.isfile(src):
        return
    if os.path.lexists(dst):
        os.remove(dst)
    os.symlink(os.path.basename(src), dst)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth-dir", required=True)
    ap.add_argument("--panel-dir", required=True)
    ap.add_argument("--liftover-dir", required=True)
    args = ap.parse_args()

    sys.path.insert(0, args.liftover_dir)
    from liftover_utils import Liftover

    lo = Liftover()
    truth = args.truth_dir
    panel = args.panel_dir

    panel_n, panel_fail = lift_markers_panel(
        lo, os.path.join(panel, "markers_sim.tsv"), os.path.join(panel, "markers_lifted.tsv")
    )
    print(f"markers_lifted: {panel_n} ok, {panel_fail} fail")

    mm_mk = os.path.join(truth, "truth_markers_by_chr")
    sl6_mk = os.path.join(truth, "truth_markers_by_chr_SL6")
    mk_stats, _ = lift_truth_markers(lo, mm_mk, sl6_mk)
    for c, (ni, no, nf) in mk_stats.items():
        print(f"  truth_markers {c}: in={ni} out={no} fail={nf}")

    block_stats = {}
    for tag in ("simF1", "simF2"):
        src = os.path.join(truth, f"{tag}_hap_blocks.tsv")
        dst = os.path.join(truth, f"{tag}_hap_blocks_SL6.tsv")
        ni, no, nf = lift_hap_blocks(lo, src, dst)
        block_stats[tag] = (ni, no, nf)
        print(f"{tag}_hap_blocks_SL6: in={ni} out={no} fail={nf}")

    bp_src = os.path.join(truth, "truth_breakpoints.bed")
    bp_dst = os.path.join(truth, "truth_breakpoints_SL6.bed")
    bp_ni, bp_no, bp_nf = lift_breakpoints(lo, bp_src, bp_dst)
    print(f"truth_breakpoints_SL6: in={bp_ni} out={bp_no} fail={bp_nf}")

    link_if_missing(bp_dst, os.path.join(truth, "simF2_co_breakpoints_SL6.bed"))

    summary = os.path.join(truth, "truth_SL6_lift_summary.tsv")
    with open(summary, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["artifact", "n_in", "n_out", "n_fail"])
        w.writerow(["markers_lifted", panel_n + panel_fail, panel_n, panel_fail])
        for c, (ni, no, nf) in mk_stats.items():
            w.writerow([f"truth_markers_{c}", ni, no, nf])
        for tag, (ni, no, nf) in block_stats.items():
            w.writerow([f"{tag}_hap_blocks_SL6", ni, no, nf])
        w.writerow(["truth_breakpoints_SL6", bp_ni, bp_no, bp_nf])
    print(f"summary: {summary}")


if __name__ == "__main__":
    main()
