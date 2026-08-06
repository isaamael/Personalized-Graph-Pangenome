#!/usr/bin/env python3
"""PGGB–SyRI concordance — recommended final matching rules.

Class          Size              Match
SNP            1 bp SNV          exact CHROM,POS,REF,ALT (upper)
Small indel    1–19 bp           exact REF/ALT at POS; OR same type, |dPOS|≤3,
                                 len_sim≥0.7 (size need not be equal)
Midsize indel  20–49 bp          same type, |dPOS|≤10, len_sim≥0.7
Large INS/DEL  ≥50 bp            same type, |dPOS|≤200, len_sim≥0.7
Complex        MNP/INV/DUP/…     excluded from P/R/F1

seq_sim not enforced (allele representation / left-norm mismatch).
"""
from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path

OUTDIR = WORK0 = BEDS = Path()
GENOME_BP = 0
CHRS = [f"chr{i}" for i in range(1, 13)]
CLASSES = ["SNP", "Small", "Midsize", "Large"]
POS_TOL = {"SNP": 0, "Small": 3, "Midsize": 10, "Large": 200}
LEN_SIM_MIN = {"Small": 0.7, "Midsize": 0.7, "Large": 0.7}
SEQ_SIM_MIN = {"Small": 0.7, "Midsize": 0.7, "Large": 0.7}


def open_vcf(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def parse_info(info: str) -> dict:
    d = {}
    for x in info.split(";"):
        if "=" in x:
            k, v = x.split("=", 1)
            d[k] = v
    return d


def indel_allele_seq(ref: str, alt: str, typ: str) -> str | None:
    """Inserted or deleted sequence for seq_sim; None if unusable."""
    if not ref or not alt or alt.startswith("<"):
        return None
    if any(c not in "ACGTNacgtn" for c in ref + alt):
        return None
    ref, alt = ref.upper(), alt.upper()
    if typ == "INS":
        # left-anchored: REF shorter
        if len(alt) <= len(ref):
            return None
        # common: REF = anchor, ALT = anchor+ins
        if alt.startswith(ref):
            return alt[len(ref):]
        if ref.startswith(alt):
            return None
        return alt  # fallback raw
    if typ == "DEL":
        if len(ref) <= len(alt):
            return None
        if ref.startswith(alt):
            return ref[len(alt):]
        return ref
    return None


def seq_sim(a: str | None, b: str | None) -> float | None:
    """None = cannot evaluate (skip seq filter)."""
    if a is None or b is None or a == "" or b == "":
        return None
    if a == b:
        return 1.0
    # simple identity / max-len (no alignment needed for short; for long use overlap ratio)
    la, lb = len(a), len(b)
    if la == 0 or lb == 0:
        return 0.0
    # character identity on min prefix + length penalty via min/max already in len_sim
    m = min(la, lb)
    match = sum(1 for i in range(m) if a[i] == b[i])
    return match / max(la, lb)


def len_sim(a: int, b: int) -> float:
    if a <= 0 or b <= 0:
        return 0.0
    return min(a, b) / max(a, b)


def classify_record(ref: str, alt: str, info: dict):
    if "," in alt:
        return None
    if alt.startswith("<"):
        typ = alt.strip("<>").upper()
        if typ not in ("INS", "DEL"):
            return None
        end = int(float(info.get("END", 0) or 0))
        svlen = abs(int(float(info.get("SVLEN", 0) or 0)))
        pos = int(info.get("_POS", 0))
        if svlen <= 0 and end > pos:
            svlen = end - pos
        if svlen <= 0:
            return None
        if svlen <= 19:
            cls = "Small"
        elif svlen <= 49:
            cls = "Midsize"
        else:
            cls = "Large"
        return cls, typ, svlen, None  # no usable seq

    if any(c not in "ACGTNacgtn" for c in ref + alt):
        return None
    ref_u, alt_u = ref.upper(), alt.upper()
    lr, la = len(ref_u), len(alt_u)
    if lr == 1 and la == 1:
        return "SNP", "SNP", 1, None
    if lr == la:
        return None  # MNP / complex
    sz = abs(lr - la)
    typ = "DEL" if lr > la else "INS"
    if sz <= 19:
        cls = "Small"
    elif sz <= 49:
        cls = "Midsize"
    else:
        cls = "Large"
    seq = indel_allele_seq(ref_u, alt_u, typ)
    return cls, typ, sz, seq


def load_bed(bed: Path | None):
    if bed is None:
        return None
    d = defaultdict(list)
    with open(bed) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            a = line.split()
            d[a[0]].append((int(a[1]), int(a[2])))
    for c in d:
        d[c].sort()
    return d


def in_mask(mask, chrom, pos1):
    if mask is None:
        return True
    ivs = mask.get(chrom)
    if not ivs:
        return False
    x = pos1 - 1
    lo, hi = 0, len(ivs)
    while lo < hi:
        mid = (lo + hi) // 2
        if ivs[mid][1] <= x:
            lo = mid + 1
        else:
            hi = mid
    return lo < len(ivs) and ivs[lo][0] <= x < ivs[lo][1]


def load_variants(paths, mask, require_samechr_orig: bool):
    data = {c: {k: [] for k in CLASSES} for c in CHRS}
    n_keep = n_cx = n_mask = n_cross = 0
    for path in paths:
        with open_vcf(path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                a = line.rstrip("\n").split("\t")
                chrom, pos_s, _id, ref, alt, _q, _f, info = a[:8]
                if chrom not in CHRS:
                    continue
                pos = int(pos_s)
                if not in_mask(mask, chrom, pos):
                    n_mask += 1
                    continue
                if require_samechr_orig:
                    m = re.search(r"ORIG_TS=([^;\t]+)", info)
                    if m and m.group(1).split(":")[0] != chrom:
                        n_cross += 1
                        continue
                inf = parse_info(info)
                inf["_POS"] = pos
                hit = classify_record(ref, alt, inf)
                if hit is None:
                    n_cx += 1
                    continue
                cls, typ, sz, seq = hit
                data[chrom][cls].append({
                    "pos": pos, "size": sz, "type": typ,
                    "ref": ref.upper(), "alt": alt.upper(),
                    "seq": seq, "matched": False,
                })
                n_keep += 1
    for c in CHRS:
        for cls in CLASSES:
            data[c][cls].sort(key=lambda x: (x["pos"], x["size"]))
    print(f"keep={n_keep} complex={n_cx} mask={n_mask} cross={n_cross}", flush=True)
    return data


def match_snp(P, S):
    idx = {}
    for i, s in enumerate(S):
        idx.setdefault((s["pos"], s["ref"], s["alt"]), []).append(i)
    tp = 0
    used = set()
    for p in P:
        for i in idx.get((p["pos"], p["ref"], p["alt"]), []):
            if i not in used:
                used.add(i)
                tp += 1
                break
    return tp, len(P) - tp, len(S) - tp


def pair_ok(p, s, cls: str) -> bool:
    if p["type"] != s["type"]:
        return False
    if abs(p["pos"] - s["pos"]) > POS_TOL[cls]:
        return False
    if len_sim(p["size"], s["size"]) < LEN_SIM_MIN[cls]:
        return False
    # Allele-seq similarity optional: SyRI remapped vs PGGB often differ in
    # left-normalization / representation → do NOT require seq_sim (user OK).
    return True


def match_small(P, S):
    """Exact allele first; else same type + Δpos≤3 + len_sim≥0.7."""
    exact = {}
    for i, s in enumerate(S):
        exact.setdefault((s["pos"], s["ref"], s["alt"]), []).append(i)
    used = [False] * len(S)
    matched_p = [False] * len(P)
    tp = 0
    for pi, p in enumerate(P):
        for i in exact.get((p["pos"], p["ref"], p["alt"]), []):
            if not used[i]:
                used[i] = True
                matched_p[pi] = True
                tp += 1
                break
    j0 = 0
    for pi, p in enumerate(P):
        if matched_p[pi]:
            continue
        while j0 < len(S) and S[j0]["pos"] < p["pos"] - POS_TOL["Small"]:
            j0 += 1
        best = -1
        best_key = None
        j = j0
        while j < len(S) and S[j]["pos"] <= p["pos"] + POS_TOL["Small"]:
            if not used[j] and pair_ok(p, S[j], "Small"):
                d = abs(S[j]["pos"] - p["pos"])
                dd = abs(S[j]["size"] - p["size"])
                key = (d, dd)
                if best_key is None or key < best_key:
                    best_key = key
                    best = j
            j += 1
        if best >= 0:
            used[best] = True
            tp += 1
    return tp, len(P) - tp, len(S) - tp


def match_fuzzy(P, S, cls: str):
    used = [False] * len(S)
    tol = POS_TOL[cls]
    j0 = 0
    tp = 0
    for p in P:
        while j0 < len(S) and S[j0]["pos"] < p["pos"] - tol:
            j0 += 1
        best = -1
        best_d = None
        j = j0
        while j < len(S) and S[j]["pos"] <= p["pos"] + tol:
            if not used[j] and pair_ok(p, S[j], cls):
                d = abs(S[j]["pos"] - p["pos"])
                dd = abs(S[j]["size"] - p["size"])
                if best_d is None or d < best_d[0] or (d == best_d[0] and dd < best_d[1]):
                    best_d = (d, dd)
                    best = j
            j += 1
        if best >= 0:
            used[best] = True
            tp += 1
    return tp, len(P) - tp, len(S) - tp


def prf(tp, fp, fn):
    p = tp / (tp + fp) if (tp + fp) else 0.0
    r = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = (2 * p * r / (p + r)) if (p + r) else 0.0
    return p, r, f1


def bed_bp(path: Path) -> int:
    return sum(int(a[2]) - int(a[1]) for a in (l.split() for l in open(path)))


def eval_stratum(name, mask_bed, mask_bp, syri_paths, pggb_paths):
    mask = load_bed(mask_bed)
    print(f"\n=== {name} ===", flush=True)
    print("SyRI", flush=True)
    syri = load_variants(syri_paths, mask, True)
    print("PGGB", flush=True)
    pggb = load_variants(pggb_paths, mask, False)
    rows = []
    for cls in CLASSES:
        tp = fp = fn = np_ = ns_ = 0
        for c in CHRS:
            P, S = pggb[c][cls], syri[c][cls]
            np_ += len(P)
            ns_ += len(S)
            if cls == "SNP":
                t, f, n = match_snp(P, S)
            elif cls == "Small":
                t, f, n = match_small(P, S)
            else:
                t, f, n = match_fuzzy(P, S, cls)
            tp += t
            fp += f
            fn += n
        p, r, f1 = prf(tp, fp, fn)
        size_lab = {"SNP": "1bp SNV", "Small": "1-19", "Midsize": "20-49", "Large": ">=50"}[cls]
        rows.append({
            "Stratum": name,
            "Mask_Mb": mask_bp / 1e6,
            "Genome_fraction": mask_bp / GENOME_BP,
            "Class": cls,
            "Size_bp": size_lab,
            "Pos_tol_bp": POS_TOL[cls],
            "PGGB": np_, "SyRI": ns_,
            "TP": tp, "FP": fp, "FN": fn,
            "Precision": p, "Recall": r, "F1": f1,
        })
        print(f"{cls}: P={p:.4f} R={r:.4f} F1={f1:.4f} TP={tp} FP={fp} FN={fn} nP={np_} nS={ns_}", flush=True)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--concordance-dir", required=True, type=Path)
    parser.add_argument("--genome-bp", required=True, type=int)
    args = parser.parse_args()
    global OUTDIR, WORK0, BEDS, GENOME_BP
    OUTDIR = args.concordance_dir
    WORK0 = OUTDIR / "prepared_variants"
    BEDS = OUTDIR / "beds"
    GENOME_BP = args.genome_bp
    OUTDIR.mkdir(parents=True, exist_ok=True)
    syri_all = [
        WORK0 / "syri_mm/syri_mm.snp.vcf.gz",
        WORK0 / "syri_mm/syri_mm.indel_lt50.vcf.gz",
        WORK0 / "syri_mm/syri_mm.sv_ge50.vcf.gz",
    ]
    pggb_all = [WORK0 / "pggb/pggb.chr.vcf.gz"]
    strata = [
        ("A_genomewide", None, GENOME_BP),
        ("B_collinear_1to1", BEDS / "mask_B.bed", bed_bp(BEDS / "mask_B.bed")),
        ("C_collinear_nonrep", BEDS / "mask_C.bed", bed_bp(BEDS / "mask_C.bed")),
    ]
    all_rows = []
    for name, bed, bp in strata:
        all_rows.extend(eval_stratum(name, bed, bp, syri_all, pggb_all))

    summary = OUTDIR / "pggb_syri_concordance_final_rules_ABC_summary_table.tsv"
    detail = OUTDIR / "pggb_syri_concordance_final_rules_SNP_Small_Midsize_Large_ABC.tsv"
    tcols = ["Stratum", "Mask_Mb", "Genome_fraction", "Class", "PGGB", "SyRI",
             "TP", "FP", "FN", "Precision", "Recall", "F1"]
    dcols = tcols[:4] + ["Size_bp", "Pos_tol_bp"] + tcols[4:]

    with open(summary, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=tcols, delimiter="\t")
        w.writeheader()
        for r in all_rows:
            w.writerow({
                "Stratum": r["Stratum"],
                "Mask_Mb": f"{r['Mask_Mb']:.3f}",
                "Genome_fraction": f"{r['Genome_fraction']:.4f}",
                "Class": r["Class"],
                "PGGB": r["PGGB"], "SyRI": r["SyRI"],
                "TP": r["TP"], "FP": r["FP"], "FN": r["FN"],
                "Precision": f"{r['Precision']:.4f}",
                "Recall": f"{r['Recall']:.4f}",
                "F1": f"{r['F1']:.4f}",
            })

    with open(detail, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=dcols, delimiter="\t")
        w.writeheader()
        for r in all_rows:
            w.writerow({
                "Stratum": r["Stratum"],
                "Mask_Mb": f"{r['Mask_Mb']:.3f}",
                "Genome_fraction": f"{r['Genome_fraction']:.4f}",
                "Class": r["Class"],
                "Size_bp": r["Size_bp"],
                "Pos_tol_bp": r["Pos_tol_bp"],
                "PGGB": r["PGGB"], "SyRI": r["SyRI"],
                "TP": r["TP"], "FP": r["FP"], "FN": r["FN"],
                "Precision": f"{r['Precision']:.6f}",
                "Recall": f"{r['Recall']:.6f}",
                "F1": f"{r['F1']:.6f}",
            })

    print("wrote", summary)
    print("wrote", detail)


if __name__ == "__main__":
    main()
