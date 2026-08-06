#!/usr/bin/env python3
"""Project reference centromere intervals to assembly coords via minimap2 PAF."""
import argparse
from pathlib import Path


def norm_chr(name):
    s = name.strip()
    if s.lower().startswith("chr"):
        n = s[3:].lstrip("0") or "0"
        return f"chr{int(n)}"
    return s


def load_paf(path):
    blocks = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 10:
                continue
            blocks.append({
                "qname": c[0], "qlen": int(c[1]), "qs": int(c[2]), "qe": int(c[3]),
                "strand": c[4],
                "tname": c[5], "tlen": int(c[6]), "ts": int(c[7]), "te": int(c[8]),
            })
    return blocks


def overlap(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0))


def lift_interval(blocks, ref_chr, rs, re):
    """PAF from minimap2 REF QRY: col1-3=ref, col6-8=qry."""
    ref_chr = norm_chr(ref_chr)
    rs, re = int(rs), int(re)
    best = None
    for b in blocks:
        if norm_chr(b["qname"]) != ref_chr:
            continue
        ov = overlap(rs, re, b["qs"], b["qe"])
        if ov <= 0:
            continue
        if best is None or ov > best["ov"]:
            best = {**b, "ov": ov}
    if best is None:
        return None

    b = best
    span = b["qe"] - b["qs"]
    if span <= 0:
        return None
    frac0 = (max(rs, b["qs"]) - b["qs"]) / span
    frac1 = (min(re, b["qe"]) - b["qs"]) / span
    tspan = b["te"] - b["ts"]
    q0 = int(b["ts"] + frac0 * tspan)
    q1 = int(b["ts"] + frac1 * tspan)
    if b["strand"] == "-":
        q0, q1 = b["te"] - int(frac1 * tspan), b["te"] - int(frac0 * tspan)
    if q0 > q1:
        q0, q1 = q1, q0
    return {
        "qry_chr": norm_chr(b["tname"]),
        "qry_start": q0,
        "qry_end": q1,
        "ref_chr": ref_chr,
        "ref_start": rs,
        "ref_end": re,
        "overlap_bp": b["ov"],
        "strand": b["strand"],
    }


def parse_cen_pos(path, species, accession):
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("Clades"):
                continue
            p = [x.strip() for x in line.split("\t")]
            if len(p) < 7:
                continue
            if p[1] != species or p[2] != accession:
                continue
            chr_name = norm_chr(p[3])
            start = int(p[4].replace(",", ""))
            end = int(p[5].replace(",", ""))
            rows.append((chr_name, start, end, p[6].replace(",", "")))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cen-pos", required=True)
    ap.add_argument("--species", required=True)
    ap.add_argument("--accession", required=True)
    ap.add_argument("--paf", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    blocks = load_paf(args.paf)
    cens = parse_cen_pos(args.cen_pos, args.species, args.accession)

    bed_path = outdir / f"{args.sample}_centromeres_final.bed"
    tsv_path = outdir / f"{args.sample}.centromere_projection.tsv"
    unmapped_path = outdir / f"{args.sample}.centromere_unmapped.tsv"

    with open(tsv_path, "w") as tf, open(bed_path, "w") as bf, open(unmapped_path, "w") as uf:
        tf.write("sample\tref_chr\tref_start\tref_end\tqry_chr\tqry_start\tqry_end\toverlap_bp\tstrand\tstatus\n")
        uf.write("sample\tref_chr\tref_start\tref_end\tstatus\n")
        ok = 0
        for chr_name, rs, re, length in cens:
            hit = lift_interval(blocks, chr_name, rs, re)
            if hit is None:
                uf.write(f"{args.sample}\t{chr_name}\t{rs}\t{re}\tUNMAPPED\n")
                tf.write(f"{args.sample}\t{chr_name}\t{rs}\t{re}\t.\t.\t.\t0\t.\tUNMAPPED\n")
                continue
            ok += 1
            bf.write(f"{hit['qry_chr']}\t{hit['qry_start']}\t{hit['qry_end']}\t{args.sample}_{hit['qry_chr']}_cen\n")
            tf.write(
                f"{args.sample}\t{chr_name}\t{rs}\t{re}\t{hit['qry_chr']}\t{hit['qry_start']}\t{hit['qry_end']}"
                f"\t{hit['overlap_bp']}\t{hit['strand']}\tOK\n"
            )

    link = outdir / f"{args.sample}_projected_centromeres.bed"
    if link.exists() or link.is_symlink():
        link.unlink()
    link.symlink_to(bed_path.name)

    print(f"{args.sample}: mapped {ok}/{len(cens)} centromeres -> {bed_path}")
    if ok < len(cens):
        raise SystemExit(f"UNMAPPED centromeres: {len(cens) - ok} (see {unmapped_path})")


if __name__ == "__main__":
    main()
