#!/usr/bin/env python3
"""step02: PedigreeSim F2 重组 → 片段表 + parent-CO 真值（MM 坐标）。"""

from __future__ import annotations

import argparse
import bisect
import csv
import os
import subprocess
import sys
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
    breakpoint_ok: bool


@dataclass
class ChrMarkers:
    all_m: List[Marker]
    bp: List[Marker]
    positions: List[int]
    bp_positions: List[int]

    @classmethod
    def from_markers(cls, markers: List[Marker]) -> "ChrMarkers":
        bp = [m for m in markers if m.breakpoint_ok]
        return cls(
            all_m=markers,
            bp=bp,
            positions=[m.mm_pos for m in markers],
            bp_positions=[m.mm_pos for m in bp],
        )

    def at_or_before(self, mm_pos: int) -> Marker:
        idx = bisect.bisect_right(self.positions, mm_pos) - 1
        if idx < 0:
            return self.all_m[0]
        return self.all_m[idx]


def load_chr_lengths(fai: str) -> Dict[str, int]:
    out: Dict[str, int] = {}
    with open(fai) as fh:
        for line in fh:
            c, ln, *_ = line.split()
            out[c] = int(ln)
    return out


def load_markers(path: str, pilot_chr: str) -> Dict[str, List[Marker]]:
    by_chr: Dict[str, List[Marker]] = defaultdict(list)
    with open(path) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            if pilot_chr and row["mm_chr"] != pilot_chr:
                continue
            by_chr[row["mm_chr"]].append(Marker(
                marker_id=row["marker_id"],
                mm_chr=row["mm_chr"],
                mm_pos=int(row["mm_pos"]),
                ts_chr=row["ts_chr"],
                ts_pos=int(row["ts_pos"]),
                breakpoint_ok=bool(int(row["breakpoint_ok"])),
            ))
    for c in by_chr:
        by_chr[c].sort(key=lambda m: m.mm_pos)
    return by_chr


def cm_from_mm(mm_pos: int, rate: float = 1.0) -> float:
    return mm_pos / 1_000_000.0 * rate


def mm_from_cm(cm: float, chr_len: int) -> int:
    pos = int(round(cm * 1_000_000.0))
    return max(1, min(pos, chr_len))


def snap_to_bp(mm_pos: int, bp_positions: List[int], bp: List[Marker]) -> int:
    if not bp:
        return mm_pos
    idx = bisect.bisect_left(bp_positions, mm_pos)
    if idx <= 0:
        return bp[0].mm_pos
    if idx >= len(bp):
        return bp[-1].mm_pos
    if abs(bp[idx - 1].mm_pos - mm_pos) <= abs(bp[idx].mm_pos - mm_pos):
        return bp[idx - 1].mm_pos
    return bp[idx].mm_pos


def founder_parent(fa: int) -> str:
    return "TS" if fa in (0, 1) else "MM"


def parse_hsb_recomb(hbw: List[str]) -> List[float]:
    out: List[float] = []
    for x in hbw:
        if x == "NA":
            break
        out.append(float(x))
    return out


def frag_len(f: dict) -> int:
    if f["parent"] == "TS":
        return f["ts_hi"] - f["ts_lo"] + 1
    return f["mm_hi"] - f["mm_lo"] + 1


def write_pedigreesim_inputs(
    work: str,
    markers_by_chr: Dict[str, List[Marker]],
    chr_lens: Dict[str, int],
    seed: int,
    n_f2: int,
) -> None:
    os.makedirs(work, exist_ok=True)
    chrom_path = os.path.join(work, "chromosomes.txt")
    with open(chrom_path, "w") as fh:
        fh.write("CHROMOSOME\tLENGTH\tCENTROMERE\n")
        for mm_chr in sorted(markers_by_chr.keys(), key=lambda c: int(c.replace("chr", ""))):
            clen = chr_lens.get(mm_chr, markers_by_chr[mm_chr][-1].mm_pos)
            cM = cm_from_mm(clen)
            fh.write(f"{mm_chr}\t{cM:.6f}\t{cM / 2.0:.6f}\n")

    par_path = os.path.join(work, "PedigreeSim.par")
    with open(par_path, "w") as fh:
        fh.write("PLOIDY=2\n")
        fh.write("MAPFUNCTION=HALDANE\n")
        fh.write(f"SEED={seed}\n")
        fh.write("CHROMFILE=chromosomes.txt\n")
        fh.write("OUTPUT=sim\n")
        fh.write("POPTYPE=F2\n")
        fh.write(f"POPSIZE={n_f2}\n")


def run_pedigreesim(java: str, jar: str, ps_dir: str, work: str) -> None:
    jar_link = os.path.join(work, "PedigreeSim.jar")
    lib_link = os.path.join(work, "lib")
    if not os.path.lexists(jar_link):
        os.symlink(jar, jar_link)
    if not os.path.lexists(lib_link):
        os.symlink(os.path.join(ps_dir, "lib"), lib_link)

    cmd = [java, "-Xmx4g", "-jar", "PedigreeSim.jar", "PedigreeSim.par"]
    proc = subprocess.run(cmd, cwd=work, capture_output=True, text=True)
    if proc.returncode != 0:
        sys.stderr.write(proc.stdout)
        sys.stderr.write(proc.stderr)
        raise RuntimeError(f"PedigreeSim failed rc={proc.returncode}")
    if not os.path.isfile(os.path.join(work, "sim.hsa")):
        sys.stderr.write(proc.stdout)
        sys.stderr.write(proc.stderr)
        raise RuntimeError("PedigreeSim did not produce sim.hsa")


def parse_hsa_hsb(
    hsa_path: str,
    hsb_path: str,
    chrom_order: List[str],
) -> Dict[str, Dict[str, List[Tuple[List[int], List[float]]]]]:
    with open(hsa_path) as ha, open(hsb_path) as hb:
        hsa_lines = [ln.strip() for ln in ha if ln.strip()]
        hsb_lines = [ln.strip() for ln in hb if ln.strip()]

    n_chr = len(chrom_order)
    ploidy = 2
    out: Dict[str, Dict[str, List[Tuple[List[int], List[float]]]]] = {}
    i = 0
    while i < len(hsa_lines):
        words = hsa_lines[i].split()
        ind = words[0]
        haplos: Dict[str, List[Tuple[List[int], List[float]]]] = {}
        for k in range(n_chr * ploidy):
            w = hsa_lines[i + k].split()
            hbw = hsb_lines[i + k].split()
            chrom = w[1]
            founders = [int(x) for x in w[3:] if x != "NA"]
            recomb_cm = parse_hsb_recomb(hbw)
            haplos.setdefault(chrom, []).append((founders, recomb_cm))
        out[ind] = haplos
        i += n_chr * ploidy
    return out


def haplo_to_frags(
    founders: List[int],
    recomb_cm: List[float],
    mm_chr: str,
    mm_chr_len: int,
    ts_chr_lens: Dict[str, int],
    mm_chr_lens: Dict[str, int],
    chr_mk: ChrMarkers,
) -> Tuple[List[dict], List[dict]]:
    if len(founders) == 0:
        return [], []
    if len(recomb_cm) != len(founders) - 1:
        raise RuntimeError(
            f"{mm_chr}: hsb cM count {len(recomb_cm)} != founder_segments-1 {len(founders) - 1}"
        )

    bp = chr_mk.bp
    bounds_cm = [0.0] + list(recomb_cm) + [cm_from_mm(mm_chr_len)]
    mm_bounds = [1]
    for j in range(1, len(bounds_cm) - 1):
        mm_bounds.append(snap_to_bp(mm_from_cm(bounds_cm[j], mm_chr_len), chr_mk.bp_positions, bp))
    mm_bounds.append(mm_chr_len)
    for j in range(1, len(mm_bounds)):
        if mm_bounds[j] < mm_bounds[j - 1]:
            raise RuntimeError(
                f"{mm_chr}: non-monotonic mm snap {mm_bounds[j - 1]} -> {mm_bounds[j]}"
            )

    merged_bounds: List[Tuple[str, int, int]] = []
    for si, fa in enumerate(founders):
        parent = founder_parent(fa)
        mm_lo, mm_hi = mm_bounds[si], mm_bounds[si + 1]
        if merged_bounds and merged_bounds[-1][0] == parent:
            merged_bounds[-1] = (parent, merged_bounds[-1][1], mm_hi)
        else:
            merged_bounds.append((parent, mm_lo, mm_hi))

    frags: List[dict] = []
    co_rows: List[dict] = []
    for i, (parent, mm_lo, mm_hi) in enumerate(merged_bounds):
        is_first = i == 0
        is_last = i == len(merged_bounds) - 1
        mL = chr_mk.at_or_before(mm_lo)
        mR = chr_mk.at_or_before(mm_hi) if not is_last else chr_mk.at_or_before(mm_chr_len)
        mL_id = "START" if is_first else mL.marker_id
        mR_id = "END" if is_last else mR.marker_id
        if parent == "TS":
            ts_chr = mL.ts_chr
            ts_lo = mL.ts_pos if is_first else mL.ts_pos + 1
            ts_hi = ts_chr_lens.get(ts_chr, mm_hi) if is_last else mR.ts_pos
        else:
            ts_chr = mm_chr
            ts_lo = mm_lo if is_first else mm_lo + 1
            ts_hi = mm_chr_len if is_last else mm_hi
        if is_last:
            mm_hi = mm_chr_lens.get(mm_chr, mm_chr_len)
        frags.append({
            "parent": parent,
            "mm_chr": mm_chr,
            "mm_lo": mm_lo, "mm_hi": mm_hi,
            "ts_chr": ts_chr, "ts_lo": max(1, ts_lo), "ts_hi": ts_hi,
            "marker_L": mL_id, "marker_R": mR_id,
        })
        if i < len(merged_bounds) - 1:
            nxt = merged_bounds[i + 1][0]
            co_rows.append({
                "mm_chr": mm_chr,
                "mm_pos": mm_hi,
                "parent_L": parent,
                "parent_R": nxt,
                "marker": mR_id,
            })
    return frags, co_rows


def write_fragments(path: str, sample: str, gamete: str, frags: List[dict]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([
            "sample", "gamete", "frag_idx", "parent", "mm_chr",
            "mm_lo", "mm_hi", "ts_chr", "ts_lo", "ts_hi",
            "marker_L", "marker_R", "trans_sv_id",
        ])
        for i, f in enumerate(frags):
            w.writerow([
                sample, gamete, i, f["parent"], f["mm_chr"],
                f["mm_lo"], f["mm_hi"], f["ts_chr"], f["ts_lo"], f["ts_hi"],
                f["marker_L"], f["marker_R"], f.get("trans_sv_id", "-"),
            ])


def f2_ps_id(sid: int, n_f2: int) -> str:
    width = len(str(n_f2))
    return f"F2_{sid:0{width}d}"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--mm-fai", required=True)
    ap.add_argument("--ts-fai", required=True)
    ap.add_argument("--markers-file", default="markers_sim.tsv")
    ap.add_argument("--n-f2", type=int, default=50)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--pilot-chr", default="")
    ap.add_argument("--java", default="java")
    ap.add_argument("--pedigreesim-jar", required=True)
    ap.add_argument("--pedigreesim-dir", required=True)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    markers_path = os.path.join(args.panel_dir, args.markers_file)
    if not os.path.isfile(markers_path):
        raise FileNotFoundError(f"markers file not found: {markers_path}")
    markers_by_chr = load_markers(markers_path, args.pilot_chr)
    chr_markers = {c: ChrMarkers.from_markers(ms) for c, ms in markers_by_chr.items()}
    mm_chr_lens = load_chr_lengths(args.mm_fai)
    ts_chr_lens = load_chr_lengths(args.ts_fai)
    chrom_order = sorted(markers_by_chr.keys(), key=lambda c: int(c.replace("chr", "")))
    genome_len = sum(mm_chr_lens.get(c, 0) for c in chrom_order)
    genome_cM = sum(cm_from_mm(mm_chr_lens.get(c, 0)) for c in chrom_order)

    work = os.path.join(args.out_dir, "pedigreesim_work")
    write_pedigreesim_inputs(work, markers_by_chr, mm_chr_lens, args.seed, args.n_f2)
    run_pedigreesim(args.java, args.pedigreesim_jar, args.pedigreesim_dir, work)

    haplos = parse_hsa_hsb(
        os.path.join(work, "sim.hsa"),
        os.path.join(work, "sim.hsb"),
        chrom_order,
    )

    co_summary_path = os.path.join(args.out_dir, "truth_co_summary.tsv")
    co_detail_path = os.path.join(args.out_dir, "truth_co_detail.tsv")
    bp_path = os.path.join(args.out_dir, "truth_breakpoints.bed")
    co_summary_rows: List[dict] = []

    with open(bp_path, "w") as bp_fh, \
         open(co_detail_path, "w", newline="") as co_fh:
        bp_fh.write("# ref=MM\n")
        co_w = csv.DictWriter(
            co_fh,
            fieldnames=["sample", "gamete", "mm_chr", "mm_pos", "parent_L", "parent_R", "marker"],
            delimiter="\t",
        )
        co_w.writeheader()

        for sid in range(1, args.n_f2 + 1):
            ps_id = f2_ps_id(sid, args.n_f2)
            sample = f"simF2_{sid:03d}"
            if ps_id not in haplos:
                raise RuntimeError(f"PedigreeSim output missing {ps_id}")

            for gamete, hap_idx in (("A", 0), ("B", 1)):
                all_frags: List[dict] = []
                n_co = 0
                for mm_chr in chrom_order:
                    chr_mk = chr_markers[mm_chr]
                    clen = mm_chr_lens.get(mm_chr, chr_mk.all_m[-1].mm_pos)
                    founders, recomb_cm = haplos[ps_id][mm_chr][hap_idx]
                    frags, co_rows = haplo_to_frags(
                        founders, recomb_cm, mm_chr, clen,
                        ts_chr_lens, mm_chr_lens, chr_mk,
                    )
                    all_frags.extend(frags)
                    n_co += len(co_rows)
                    for co in co_rows:
                        co_w.writerow({
                            "sample": sample, "gamete": gamete, **co,
                        })
                        pos = co["mm_pos"]
                        bp_fh.write(f"{mm_chr}\t{pos}\t{pos + 1}\t{sample}\t0\t+\n")

                hap_len = sum(frag_len(f) for f in all_frags)
                write_fragments(
                    os.path.join(args.out_dir, f"{sample}.gamete{gamete}.fragments.tsv"),
                    sample, gamete, all_frags,
                )
                co_summary_rows.append({
                    "sample": sample,
                    "gamete": gamete,
                    "n_parent_co": n_co,
                    "n_frag": len(all_frags),
                    "hap_len_bp": hap_len,
                    "hap_len_ratio": f"{hap_len / genome_len:.4f}" if genome_len else "0",
                })

    with open(co_summary_path, "w", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=["sample", "gamete", "n_parent_co", "n_frag", "hap_len_bp", "hap_len_ratio"],
            delimiter="\t",
        )
        w.writeheader()
        w.writerows(co_summary_rows)

    co_a = [r["n_parent_co"] for r in co_summary_rows if r["gamete"] == "A"]
    co_b = [r["n_parent_co"] for r in co_summary_rows if r["gamete"] == "B"]
    ratio_a = [float(r["hap_len_ratio"]) for r in co_summary_rows if r["gamete"] == "A"]
    stats_path = os.path.join(args.out_dir, "truth_co_cohort.tsv")
    with open(stats_path, "w") as fh:
        fh.write("metric\tvalue\n")
        fh.write(f"genome_len_bp\t{genome_len}\n")
        fh.write(f"genome_cM\t{genome_cM:.3f}\n")
        fh.write(f"expected_co_per_gamete\t{genome_cM / 100:.2f}\n")
        fh.write(f"n_f2\t{args.n_f2}\n")
        fh.write(f"co_gameteA_min\t{min(co_a)}\n")
        fh.write(f"co_gameteA_max\t{max(co_a)}\n")
        fh.write(f"co_gameteA_mean\t{sum(co_a)/len(co_a):.2f}\n")
        fh.write(f"co_gameteB_min\t{min(co_b)}\n")
        fh.write(f"co_gameteB_max\t{max(co_b)}\n")
        fh.write(f"co_gameteB_mean\t{sum(co_b)/len(co_b):.2f}\n")
        fh.write(f"hap_len_ratioA_min\t{min(ratio_a):.4f}\n")
        fh.write(f"hap_len_ratioA_max\t{max(ratio_a):.4f}\n")
        fh.write(f"hap_len_ratioA_mean\t{sum(ratio_a)/len(ratio_a):.4f}\n")

    meta = os.path.join(args.out_dir, "sim_meta.tsv")
    with open(meta, "w") as fh:
        fh.write(f"seed\t{args.seed}\n")
        fh.write(f"n_f2\t{args.n_f2}\n")
        fh.write(f"pilot_chr\t{args.pilot_chr or 'all'}\n")
        fh.write("ref_coord\tMM\n")
        fh.write("engine\tPedigreeSim\n")
        fh.write(f"markers_file\t{args.markers_file}\n")
        fh.write(f"co_gameteA_mean\t{sum(co_a)/len(co_a):.2f}\n")

    print(
        f"Done n_f2={args.n_f2} seed={args.seed} "
        f"coA_mean={sum(co_a)/len(co_a):.1f} coB_mean={sum(co_b)/len(co_b):.1f} "
        f"hap_ratioA_mean={sum(ratio_a)/len(ratio_a):.3f}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
