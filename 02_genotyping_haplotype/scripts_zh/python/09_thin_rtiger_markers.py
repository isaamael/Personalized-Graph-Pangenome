#!/usr/bin/env python3
"""按基因组窗口稀疏化 VCF，同时完整保留 FORMAT 列和样本列。"""
from __future__ import annotations

import argparse
import gzip
import subprocess
import sys
from dataclasses import dataclass, field
from typing import BinaryIO, Optional

try:
    import numpy as np
    from cyvcf2 import VCF
except ImportError:
    sys.exit("cyvcf2+numpy required; conda activate genome")

KEEP_METHODS = ("missing", "dp", "dp_gq", "balanced", "het5", "het1", "het0", "first")
HET_METHODS = frozenset({"het5", "het1", "het0"})
GQ_METHODS = frozenset({"dp_gq", "balanced"})
BIN_PCTL_METHODS = frozenset({"dp_gq", "balanced"})


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True)
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--bin-bp", type=int, default=50000)
    ap.add_argument(
        "--keep-method",
        choices=list(KEEP_METHODS) + ["least_missing"],
        default="balanced",
    )
    ap.add_argument(
        "--bin-pctl",
        type=float,
        default=-1.0,
        help="Within-bin quantile threshold for balanced/dp_gq selection",
    )
    ap.add_argument("--soft-dp-k", type=float, default=10.0)
    ap.add_argument("--fill-indel-sv", type=int, choices=[0, 1], default=1)
    ap.add_argument("--sv-min-len", type=int, default=50)
    ap.add_argument("--stats-tsv", default="")
    return ap.parse_args()


def normalize_method(method: str) -> str:
    return "missing" if method == "least_missing" else method


def default_bin_pctl(method: str, user_pctl: float) -> float:
    if user_pctl >= 0:
        return user_pctl
    return 0.75 if method in BIN_PCTL_METHODS else 0.0


def vcf_has_gq(header: str) -> bool:
    return any("##FORMAT" in ln and "ID=GQ" in ln for ln in header.split("\n"))


def bin_id(pos: int, bin_bp: int) -> int:
    return (pos - 1) // bin_bp


def miss_rate(gts: np.ndarray) -> float:
    n = gts.shape[0]
    if n == 0:
        return 1.0
    return float(np.sum((gts[:, 0] < 0) | (gts[:, 1] < 0))) / n


def het_rate(gts: np.ndarray) -> float:
    called = (gts[:, 0] >= 0) & (gts[:, 1] >= 0)
    n = int(called.sum())
    if n == 0:
        return 0.0
    return float(np.sum(gts[called, 0] != gts[called, 1])) / n


def mean_dp(var) -> float:
    d = var.gt_depths
    if d is None or len(d) == 0:
        return 0.0
    valid = d[d >= 0]
    return float(np.mean(valid)) if len(valid) else 0.0


def mean_gq(var, need: bool) -> float:
    if not need:
        return -1.0
    try:
        gq = np.asarray(var.format("GQ"), dtype=np.float64)
    except (KeyError, ValueError, TypeError):
        return -1.0
    valid = gq[gq >= 0]
    return float(np.mean(valid)) if len(valid) else -1.0


def has_svtype(var) -> bool:
    info = var.INFO
    if isinstance(info, dict):
        return "SVTYPE" in info
    if info:
        s = str(info)
        return "SVTYPE=" in s or s.startswith("SVTYPE=")
    return False


def variant_class(var, sv_min_len: int) -> Optional[str]:
    alts = var.ALT or []
    if len(alts) != 1:
        return None
    ref, alt = var.REF, alts[0]
    if not ref or not alt:
        return None
    if has_svtype(var) or max(len(ref), len(alt)) >= sv_min_len:
        return "sv"
    if len(ref) == 1 and len(alt) == 1 and var.is_snp:
        return "snp"
    return "indel"


@dataclass
class Cand:
    chrom: str
    pos: int
    ref: str
    alt: str
    gts: np.ndarray
    miss: float
    mean_dp: float
    mean_gq: float
    het: float
    vclass: str
    bin_idx: int = 0


def make_cand(var, vclass: str, need_het: bool, need_gq: bool, bid: int) -> Cand:
    gts = var.genotype.array()[:, :2]
    alts = var.ALT or []
    return Cand(
        var.CHROM,
        int(var.POS),
        var.REF,
        alts[0],
        gts.copy(),
        miss_rate(gts),
        mean_dp(var),
        mean_gq(var, need_gq),
        het_rate(gts) if need_het else -1.0,
        vclass,
        bid,
    )


def soft_score(cand: Cand, has_gq: bool, dp_k: float) -> float:
    miss_w = max(0.0, 1.0 - cand.miss)
    dp_w = cand.mean_dp / (cand.mean_dp + dp_k) if cand.mean_dp > 0 else 0.0
    gq_w = min(cand.mean_gq, 99.0) / 99.0 if has_gq and cand.mean_gq >= 0 else 1.0
    return gq_w * dp_w * miss_w


def score_key(cand: Cand, method: str, has_gq: bool, dp_k: float) -> tuple:
    if method == "first":
        return (cand.pos, cand.miss)
    if method == "missing":
        return (cand.miss, -cand.mean_dp, cand.pos)
    if method == "dp":
        return (-cand.mean_dp, cand.miss, cand.pos)
    if method == "dp_gq":
        if has_gq and cand.mean_gq >= 0:
            return (-cand.mean_gq, -cand.mean_dp, cand.miss, cand.pos)
        return (-cand.mean_dp, cand.miss, cand.pos)
    if method == "balanced":
        return (-soft_score(cand, has_gq, dp_k), cand.miss, cand.pos)
    if method == "het5":
        return (abs(cand.het - 0.5), cand.miss, cand.pos)
    if method == "het1":
        return (-cand.het, cand.miss, cand.pos)
    if method == "het0":
        return (cand.het, cand.miss, cand.pos)
    return (cand.miss, cand.pos)


def filter_by_bin_pctl(pool: list[Cand], has_gq: bool, bin_pctl: float) -> tuple[list[Cand], dict]:
    meta = {"n_pool": len(pool), "p_dp": "", "p_gq": "", "p_miss": "", "fallback": "N"}
    if bin_pctl <= 0 or len(pool) <= 1:
        return pool, meta
    pct = bin_pctl * 100.0
    p_dp = float(np.percentile([c.mean_dp for c in pool], pct))
    p_miss = float(np.percentile([c.miss for c in pool], 100.0 - pct))
    meta["p_dp"] = f"{p_dp:.2f}"
    meta["p_miss"] = f"{p_miss:.4f}"
    eligible = [c for c in pool if c.mean_dp >= p_dp and c.miss <= p_miss]
    gq_vals = [c.mean_gq for c in pool if c.mean_gq >= 0]
    if has_gq and gq_vals:
        p_gq = float(np.percentile(gq_vals, pct))
        meta["p_gq"] = f"{p_gq:.2f}"
        eligible = [c for c in eligible if c.mean_gq >= 0 and c.mean_gq >= p_gq]
    if not eligible:
        eligible = pool
        meta["fallback"] = "Y"
    return eligible, meta


def pick_from_pool(
    pool: list[Cand],
    method: str,
    has_gq: bool,
    dp_k: float,
    bin_pctl: float,
) -> tuple[Optional[Cand], dict]:
    if not pool:
        return None, {}
    eligible, meta = filter_by_bin_pctl(pool, has_gq, bin_pctl)
    pick = min(eligible, key=lambda c: score_key(c, method, has_gq, dp_k))
    meta["n_eligible"] = len(eligible)
    return pick, meta


def select_sites(
    inp: str,
    bin_bp: int,
    keep_method: str,
    sv_min_len: int,
    soft_dp_k: float,
    bin_pctl: float,
    fill_indel_sv: bool,
    stats_tsv: str,
) -> tuple[set[tuple[str, int, str, str]], int, dict]:
    method = normalize_method(keep_method)
    bin_pctl = default_bin_pctl(method, bin_pctl)
    need_het = method in HET_METHODS
    need_gq = method in GQ_METHODS or bin_pctl > 0
    vcf = VCF(inp)
    has_gq = vcf_has_gq(vcf.raw_header)
    # 同一位置存在多条记录时，保留选定的等位基因。
    selected: set[tuple[str, int, str, str]] = set()
    n_in = 0
    n_snp = n_indel = n_sv = n_fallback = n_skip_empty = 0
    cur_chr: Optional[str] = None
    cur_bin: Optional[int] = None
    snp_pool: list[Cand] = []
    indel_sv_pool: list[Cand] = []
    stat_rows: list[str] = []

    def settle() -> None:
        nonlocal n_snp, n_indel, n_sv, n_fallback, n_skip_empty
        nonlocal snp_pool, indel_sv_pool
        pick, meta = pick_from_pool(snp_pool, method, has_gq, soft_dp_k, bin_pctl)
        from_snp = pick is not None
        if pick is None and fill_indel_sv:
            pick, meta = pick_from_pool(indel_sv_pool, method, has_gq, soft_dp_k, bin_pctl)
        if pick is None:
            if not snp_pool and (indel_sv_pool or not fill_indel_sv):
                n_skip_empty += 1
        else:
            if meta.get("fallback") == "Y":
                n_fallback += 1
            selected.add((pick.chrom, pick.pos, pick.ref, pick.alt))
            if from_snp:
                n_snp += 1
            elif pick.vclass == "indel":
                n_indel += 1
            else:
                n_sv += 1
            if stats_tsv:
                ss = soft_score(pick, has_gq, soft_dp_k)
                stat_rows.append(
                    f"{pick.chrom}\t{pick.bin_idx}\t{pick.pos}\t{pick.vclass}\t"
                    f"{pick.mean_dp:.2f}\t{pick.mean_gq:.2f}\t{pick.miss:.4f}\t{ss:.6f}\t"
                    f"{meta.get('n_pool','')}\t{meta.get('n_eligible','')}\t"
                    f"{meta.get('p_dp','')}\t{meta.get('p_gq','')}\t{meta.get('p_miss','')}\t"
                    f"{meta.get('fallback','')}"
                )
        snp_pool.clear()
        indel_sv_pool.clear()

    for var in vcf:
        n_in += 1
        vclass = variant_class(var, sv_min_len)
        if vclass is None:
            continue
        chrom = var.CHROM
        pos = int(var.POS)
        bid = bin_id(pos, bin_bp)
        if cur_chr is None:
            cur_chr, cur_bin = chrom, bid
        elif chrom != cur_chr:
            settle()
            cur_chr, cur_bin = chrom, bid
        elif bid != cur_bin:
            settle()
            cur_bin = bid
        cand = make_cand(var, vclass, need_het, need_gq, bid)
        if vclass == "snp":
            snp_pool.append(cand)
        elif fill_indel_sv:
            indel_sv_pool.append(cand)
    settle()
    vcf.close()

    if stats_tsv and stat_rows:
        with open(stats_tsv, "w") as f:
            f.write(
                "chrom\tbin\tpos\tvclass\tmean_dp\tmean_gq\tmiss\tsoft_score\t"
                "n_pool\tn_eligible\tp_dp\tp_gq\tp_miss\tfallback\n"
            )
            f.write("\n".join(stat_rows) + "\n")

    meta = {
        "n_in": n_in,
        "n_out": len(selected),
        "n_snp": n_snp,
        "n_indel": n_indel,
        "n_sv": n_sv,
        "n_skip_empty": n_skip_empty,
        "n_fallback": n_fallback,
        "method": method,
        "bin_pctl": bin_pctl,
        "has_gq": has_gq,
    }
    return selected, n_in, meta


def write_thinned_vcf(inp: str, out: str, selected: set[tuple[str, int, str, str]]) -> int:
    n_out = 0
    opener = gzip.open if inp.endswith(".gz") else open
    with opener(inp, "rt") as fin, open(out, "wb") as fout:
        proc = subprocess.Popen(["bgzip", "-c"], stdin=subprocess.PIPE, stdout=fout, text=False)
        assert proc.stdin is not None
        bgzip_in: BinaryIO = proc.stdin
        for line in fin:
            if line.startswith("#"):
                bgzip_in.write(line.encode("utf-8"))
                continue
            cols = line.split("\t", 5)
            if len(cols) < 5:
                continue
            try:
                pos = int(cols[1])
            except ValueError:
                continue
            # 精确匹配 REF 和 ALT，避免同一位置出现重复记录。
            if (cols[0], pos, cols[3], cols[4]) in selected:
                if not line.endswith("\n"):
                    line += "\n"
                bgzip_in.write(line.encode("utf-8"))
                n_out += 1
        bgzip_in.close()
        if proc.wait() != 0:
            sys.exit("bgzip failed")
    subprocess.run(["tabix", "-p", "vcf", out], check=True)
    return n_out


def run_thin(
    inp: str,
    out: str,
    bin_bp: int,
    keep_method: str,
    sv_min_len: int,
    soft_dp_k: float,
    bin_pctl: float,
    fill_indel_sv: int = 1,
    stats_tsv: str = "",
) -> int:
    selected, n_in, meta = select_sites(
        inp, bin_bp, keep_method, sv_min_len, soft_dp_k, bin_pctl, bool(fill_indel_sv), stats_tsv
    )
    n_out = write_thinned_vcf(inp, out, selected)
    if n_out != len(selected):
        sys.exit(f"thin write mismatch: selected={len(selected)} written={n_out}")
    print(
        f"thin: in={n_in} kept={n_out} snp={meta['n_snp']} indel_fill={meta['n_indel']} "
        f"sv_fill={meta['n_sv']} skip_empty_bin={meta['n_skip_empty']} "
        f"fill_indel_sv={int(fill_indel_sv)} bin_bp={bin_bp} method={meta['method']} "
        f"bin_pctl={meta['bin_pctl']} gq={'on' if meta['has_gq'] else 'off'} "
        f"fallback_bins={meta['n_fallback']} fmt=preserve"
        + (f" soft_dp_k={soft_dp_k}" if meta["method"] == "balanced" else "")
    )
    return n_out


def main():
    args = parse_args()
    if args.bin_bp <= 0:
        sys.exit("bin_bp must be > 0")
    if args.soft_dp_k <= 0:
        sys.exit("soft_dp_k must be > 0")
    if args.bin_pctl > 1:
        sys.exit("bin_pctl must be in [0,1]")
    run_thin(
        args.input,
        args.output,
        args.bin_bp,
        args.keep_method,
        args.sv_min_len,
        args.soft_dp_k,
        args.bin_pctl,
        args.fill_indel_sv,
        args.stats_tsv,
    )


if __name__ == "__main__":
    main()
