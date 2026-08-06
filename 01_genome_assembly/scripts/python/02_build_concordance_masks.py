#!/usr/bin/env python3
"""Stratified PGGB–SyRI concordance: A=genomewide (existing), B=1-1 collinear, C=B∩non-repetitive.

B mask (MM coords):
  include SYN (query) ∪ SYNAL (query)
  exclude INV*/TRANS*/DUP*/HDR/NOTAL (query)
  exclude chromosome ends (END_MARGIN_BP)
  exclude coord-failure proxy: remapped ShV with ORIG_TS chrom != MM chrom

C = B minus:
  TE (EDTA TEanno)
  tandem (SyRI TDM + EDTA tandem/simple if present)
  satellite / centromere (MM_centromeres_final.bed)
  segmental duplication (SyRI DUP*/INVDP* on MM)
  assembly gap ± GAP_FLANK_BP
"""
from __future__ import annotations

import argparse
import re
import subprocess
from pathlib import Path

SYRI_OUT = MM_FAI = MM_FA = EDTA_GFF = CEN_BED = Path()


END_MARGIN_BP = 100_000
GAP_FLANK_BP = 1_000
GAP_MIN_N = 10
CHRS = [f"chr{i}" for i in range(1, 13)]

EXCLUDE_STRUCT = {
    "INV", "INVAL", "INVDP", "INVDPAL", "INVTR", "INVTRAL",
    "TRANS", "TRANSAL",
    "DUP", "DUPAL",
    "HDR", "NOTAL",
}
DUP_SD = {"DUP", "DUPAL", "INVDP", "INVDPAL", "CPG", "CPL"}
TE_SKIP_TYPES = {"target_site_duplication"}  # tiny TSD footprints; keep as TE-adjacent via parent TE


def run(cmd, check=True):
    print("+", " ".join(map(str, cmd)), flush=True)
    return subprocess.run(cmd, check=check)


def load_fai(path: Path) -> dict[str, int]:
    d = {}
    with open(path) as fh:
        for line in fh:
            c, ln = line.split()[:2]
            d[c] = int(ln)
    return d


def syri_mm_intervals(syri_out: Path, types: set[str]) -> list[tuple[str, int, int, str]]:
    """Return 0-based half-open MM (query) intervals for given SyRI types."""
    rows = []
    with open(syri_out) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            a = line.rstrip("\n").split("\t")
            if len(a) < 11:
                continue
            typ = a[10]
            if typ not in types:
                continue
            qc, qs, qe = a[5], a[6], a[7]
            if qc in ("-", ".") or qs in ("-", ".") or qe in ("-", "."):
                continue
            if qc not in CHRS:
                continue
            s, e = int(qs), int(qe)
            if e < s:
                s, e = e, s
            # SyRI coords 1-based inclusive → BED 0-based half-open
            rows.append((qc, s - 1, e, typ))
    return rows


def write_bed(path: Path, rows: list[tuple], merge: bool = False):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".raw")
    with open(tmp, "w") as fh:
        for r in rows:
            fh.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3] if len(r) > 3 else '.'}\n")
    if merge:
        run(["bash", "-lc",
             f"sort -k1,1 -k2,2n {tmp} | bedtools merge -i - > {path}"])
        tmp.unlink(missing_ok=True)
    else:
        run(["bash", "-lc", f"sort -k1,1 -k2,2n {tmp} -o {path}"])
        tmp.unlink(missing_ok=True)


def bed_bp(path: Path) -> int:
    n = 0
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            a = line.split()
            n += int(a[2]) - int(a[1])
    return n


def build_gap_bed(fa: Path, out: Path, min_n: int = GAP_MIN_N):
    if out.exists() and out.stat().st_size > 0:
        print(f"reuse {out}")
        return
    chrom = None
    seq: list[str] = []
    rows = 0
    total = 0

    def flush(ch, parts, wo):
        nonlocal rows, total
        s = "".join(parts).upper()
        i = 0
        L = len(s)
        while i < L:
            if s[i] == "N":
                j = i
                while j < L and s[j] == "N":
                    j += 1
                if j - i >= min_n:
                    wo.write(f"{ch}\t{i}\t{j}\n")
                    rows += 1
                    total += j - i
                i = j
            else:
                i += 1

    with open(fa) as fh, open(out, "w") as wo:
        for line in fh:
            if line.startswith(">"):
                if chrom is not None:
                    flush(chrom, seq, wo)
                chrom = line[1:].strip().split()[0]
                seq = []
            else:
                seq.append(line.strip())
        if chrom is not None:
            flush(chrom, seq, wo)
    print(f"wrote {out} gaps={rows} bp={total}")


def edta_te_beds(gff: Path, te_bed: Path, tandem_bed: Path):
    te_rows = []
    tandem_rows = []
    tandem_re = re.compile(r"tandem|simple_repeat|low_complexity|satellite|Satellite", re.I)
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            a = line.rstrip("\n").split("\t")
            if len(a) < 9:
                continue
            chrom, typ, start, end, attrs = a[0], a[2], int(a[3]), int(a[4]), a[8]
            if chrom not in CHRS:
                continue
            if typ in TE_SKIP_TYPES:
                continue
            bed = (chrom, start - 1, end, typ)
            if tandem_re.search(typ) or tandem_re.search(attrs):
                tandem_rows.append(bed)
            else:
                # all other EDTA features as TE/repeat
                te_rows.append(bed)
    write_bed(te_bed, te_rows, merge=True)
    write_bed(tandem_bed, tandem_rows, merge=True)
    print(f"TE intervals raw={len(te_rows)} merged_bp={bed_bp(te_bed)}")
    print(f"tandem/satellite-like raw={len(tandem_rows)} merged_bp={bed_bp(tandem_bed)}")


def build_masks(beddir: Path, lens: dict[str, int]):
    beddir.mkdir(parents=True, exist_ok=True)

    syn = syri_mm_intervals(SYRI_OUT, {"SYN", "SYNAL"})
    write_bed(beddir / "syn_synal.mm.bed", syn, merge=True)

    excl = syri_mm_intervals(SYRI_OUT, EXCLUDE_STRUCT)
    write_bed(beddir / "exclude_struct.mm.bed", excl, merge=True)

    # chromosome ends
    ends = []
    for c, L in lens.items():
        if c not in CHRS:
            continue
        ends.append((c, 0, min(END_MARGIN_BP, L), "tel5"))
        if L > END_MARGIN_BP:
            ends.append((c, max(0, L - END_MARGIN_BP), L, "tel3"))
    write_bed(beddir / "chr_ends.bed", ends, merge=True)

    # B = syn − struct − ends
    run(["bash", "-lc",
         f"bedtools subtract -a {beddir/'syn_synal.mm.bed'} -b {beddir/'exclude_struct.mm.bed'} "
         f"| bedtools subtract -a - -b {beddir/'chr_ends.bed'} "
         f"| bedtools sort -i - | bedtools merge -i - > {beddir/'mask_B.bed'}"])

    # C extras
    build_gap_bed(MM_FA, beddir / "MM_gaps_ge10.bed", GAP_MIN_N)
    run(["bash", "-lc",
         f"bedtools slop -i {beddir/'MM_gaps_ge10.bed'} -g {MM_FAI} -b {GAP_FLANK_BP} "
         f"| bedtools sort -i - | bedtools merge -i - > {beddir/'gaps_flank.bed'}"])

    edta_te_beds(EDTA_GFF, beddir / "te.edta.bed", beddir / "tandem_sat_like.edta.bed")

    tdm = syri_mm_intervals(SYRI_OUT, {"TDM"})
    write_bed(beddir / "tdm.syri.bed", tdm, merge=True)

    sd = syri_mm_intervals(SYRI_OUT, DUP_SD)
    write_bed(beddir / "segdup.syri.bed", sd, merge=True)

    # centromere / satellite proxy
    run(["bash", "-lc",
         f"sort -k1,1 -k2,2n {CEN_BED} | bedtools merge -i - > {beddir/'centromere.bed'}"])

    # union of C excludes
    run(["bash", "-lc",
         f"cat {beddir/'te.edta.bed'} {beddir/'tandem_sat_like.edta.bed'} "
         f"{beddir/'tdm.syri.bed'} {beddir/'centromere.bed'} "
         f"{beddir/'segdup.syri.bed'} {beddir/'gaps_flank.bed'} "
         f"| bedtools sort -i - | bedtools merge -i - > {beddir/'exclude_C.union.bed'}"])

    run(["bash", "-lc",
         f"bedtools subtract -a {beddir/'mask_B.bed'} -b {beddir/'exclude_C.union.bed'} "
         f"| bedtools sort -i - | bedtools merge -i - > {beddir/'mask_C.bed'}"])

    for name in ("syn_synal.mm.bed", "mask_B.bed", "mask_C.bed",
                 "exclude_struct.mm.bed", "exclude_C.union.bed"):
        p = beddir / name
        print(f"{name}\tintervals\t{sum(1 for _ in open(p))}\tbp\t{bed_bp(p)}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--syri-out", required=True, type=Path)
    parser.add_argument("--mm-fasta", required=True, type=Path)
    parser.add_argument("--edta-gff", required=True, type=Path)
    parser.add_argument("--centromere-bed", required=True, type=Path)
    parser.add_argument("--concordance-dir", required=True, type=Path)
    args = parser.parse_args()

    global SYRI_OUT, MM_FAI, MM_FA, EDTA_GFF, CEN_BED
    SYRI_OUT = args.syri_out
    MM_FA = args.mm_fasta
    MM_FAI = Path(str(MM_FA) + ".fai")
    EDTA_GFF = args.edta_gff
    CEN_BED = args.centromere_bed

    beds = args.concordance_dir / "beds"
    build_masks(beds, load_fai(MM_FAI))


if __name__ == "__main__":
    main()
