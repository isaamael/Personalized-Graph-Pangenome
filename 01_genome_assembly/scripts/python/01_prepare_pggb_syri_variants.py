#!/usr/bin/env python3
"""Prepare comparable MM-reference VCFs from SyRI and PGGB.

SyRI VCF is TS-referenced. Sequence ShV (SNP/indel) are remapped to MM via
INFO ChrB/StartB with allele swap (TS↔MM) to match PGGB REF=MM ALT=TS.
Symbolic SVs are written as INS/DEL with END on MM when possible for fuzzy
position/length matching.
"""
from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

CHRS = [f"chr{i}" for i in range(1, 13)]


def run(cmd, check=True):
    print("+", " ".join(map(str, cmd)), flush=True)
    return subprocess.run(cmd, check=check)


def parse_info(info: str) -> dict:
    d = {}
    for x in info.split(";"):
        if "=" in x:
            k, v = x.split("=", 1)
            d[k] = v
        else:
            d[x] = True
    return d


def load_fai(fai: Path) -> dict:
    lens = {}
    with open(fai) as fh:
        for line in fh:
            c, ln = line.split()[:2]
            lens[c] = int(ln)
    return lens


def syri_to_mm_vcfs(syri_vcf: Path, out_dir: Path, chr_lens: dict | None = None):
    """Write snp / indel<50 / indel+sv>=50 MM-ref VCFs from SyRI VCF."""
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "snp": out_dir / "syri_mm.snp.vcf",
        "indel": out_dir / "syri_mm.indel_lt50.vcf",
        "sv": out_dir / "syri_mm.sv_ge50.vcf",
    }
    headers = []
    contig_seen = set()
    rows = {k: [] for k in paths}
    chr_lens = chr_lens or {}

    with open(syri_vcf) as fh:
        for line in fh:
            if line.startswith("##"):
                if line.startswith("##contig="):
                    # keep only chr1-12 style; rewrite later from MM fai if needed
                    headers.append(line)
                    if "ID=" in line:
                        contig_seen.add(line.split("ID=")[1].split(",")[0].split(">")[0])
                elif not line.startswith("##fileDate") and not line.startswith("##source"):
                    headers.append(line)
                continue
            if line.startswith("#CHROM"):
                headers.append('##INFO=<ID=ORIG_TS,Number=1,Type=String,Description="Original TS-ref POS/REF/ALT">\n')
                headers.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTS\n")
                continue
            a = line.rstrip("\n").split("\t")
            if len(a) < 8:
                continue
            chrom, pos, vid, ref, alt, qual, filt, info = a[:8]
            inf = parse_info(info)
            qb = inf.get("ChrB", ".")
            sb = inf.get("StartB", ".")
            eb = inf.get("EndB", ".")
            if qb not in CHRS or sb in (".", "-") or sb == "":
                continue
            try:
                mm_pos = int(sb)
            except ValueError:
                continue

            gt = a[9] if len(a) > 9 else "1"
            orig = f"ORIG_TS={chrom}:{pos}:{ref}>{alt}"

            # Sequence alleles: SNP / indel
            if not alt.startswith("<") and all(x in "ACGTNacgtn" for x in ref + alt.replace(",", "")):
                # allele swap to MM-ref
                mm_ref, mm_alt = alt, ref
                # For multi-base, SyRI INS (insertion in query) has longer ALT on TS-ref;
                # after swap: longer REF on MM = DEL relative to MM when MM has the insertion...
                # SyRI REF=TS: INS means insertion in QRY(MM) → MM has extra seq → on MM-ref this is
                # not a simple allele swap for indels spanning different coordinates.
                # For ShV SNPs (len1/len1) swap is exact. For indels use length on MM StartB
                # and swapped alleles when both are sequence.
                sz = abs(len(mm_ref) - len(mm_alt))
                rec = f"{qb}\t{mm_pos}\t{vid}\t{mm_ref}\t{mm_alt}\t.\tPASS\t{orig}\tGT\t{gt}\n"
                if len(mm_ref) == 1 and len(mm_alt) == 1:
                    rows["snp"].append((qb, mm_pos, rec))
                elif sz < 50:
                    rows["indel"].append((qb, mm_pos, rec))
                else:
                    rows["sv"].append((qb, mm_pos, rec))
                continue

            # Symbolic INS/DEL on MM via ChrB; type flip TS-ref → MM-ref
            # SyRI INS (ins in MM vs TS) ↔ DEL on MM-ref perspective is subtle;
            # Keep symbolic ALT but place it on MM coordinates for fuzzy matching.
            if alt in ("<INS>", "<DEL>"):
                try:
                    end = int(eb) if eb not in (".", "-") else mm_pos
                except ValueError:
                    end = mm_pos
                if end < mm_pos:
                    mm_pos, end = end, mm_pos
                length = end - mm_pos + 1
                if length < 50 and alt.startswith("<"):
                    continue
                # Flip label for MM-ref orientation (same as concordance script)
                new_alt = "<DEL>" if alt == "<INS>" else "<INS>"
                rec = (
                    f"{qb}\t{mm_pos}\t{vid}\tN\t{new_alt}\t.\tPASS\t"
                    f"END={end};SVTYPE={new_alt[1:-1]};SVLEN={length};{orig}\tGT\t{gt}\n"
                )
                rows["sv"].append((qb, mm_pos, rec))

    for key, path in paths.items():
        with open(path, "w") as out:
            out.write("##fileformat=VCFv4.2\n")
            out.write("##source=syri_to_mm_for_pggb_comparison\n")
            for c in CHRS:
                if c in chr_lens:
                    out.write(f"##contig=<ID={c},length={chr_lens[c]}>\n")
                else:
                    out.write(f"##contig=<ID={c}>\n")
            out.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End">\n')
            out.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">\n')
            out.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">\n')
            out.write('##INFO=<ID=ORIG_TS,Number=1,Type=String,Description="Original TS-ref">\n')
            out.write('##ALT=<ID=INS,Description="Insertion">\n')
            out.write('##ALT=<ID=DEL,Description="Deletion">\n')
            out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTS\n")
            for _, _, rec in sorted(rows[key], key=lambda x: (CHRS.index(x[0]) if x[0] in CHRS else 99, x[1])):
                out.write(rec)
        print(f"wrote {path} n={len(rows[key])}")
    return paths


def bgzip_index(vcf: Path):
    gz = Path(str(vcf) + ".gz")
    if gz.exists():
        gz.unlink()
    run(["bgzip", "-f", str(vcf)])
    gz = Path(str(vcf) + ".gz")
    run(["tabix", "-f", "-p", "vcf", str(gz)])
    return gz


def prep_pggb(pggb: Path, out_dir: Path):
    """Rename MM#0#chr* → chr* and split by size class."""
    out_dir.mkdir(parents=True, exist_ok=True)
    rename = out_dir / "chr_rename.txt"
    with open(rename, "w") as fh:
        for c in CHRS:
            fh.write(f"MM#0#{c}\t{c}\n")
    renamed = out_dir / "pggb.chr.vcf.gz"
    run([
        "bcftools", "annotate", "--rename-chrs", str(rename),
        "-Oz", "-o", str(renamed), str(pggb),
    ])
    run(["bcftools", "index", "-f", "-t", str(renamed)])

    def extract(expr, outname):
        outp = out_dir / outname
        run([
            "bcftools", "view", "-i", expr, "-Oz", "-o", str(outp), str(renamed),
        ])
        run(["bcftools", "index", "-f", "-t", str(outp)])
        return outp

    # SNP: both alleles length 1
    snp = extract('strlen(REF)=1 && strlen(ALT)=1', "pggb.snp.vcf.gz")
    # indel <50
    indel = extract(
        '((strlen(REF)!=strlen(ALT)) && abs(strlen(REF)-strlen(ALT))<50 && strlen(REF)>0 && strlen(ALT)>0 && (strlen(REF)>1 || strlen(ALT)>1))',
        "pggb.indel_lt50.vcf.gz",
    )
    # SV >=50 sequence
    sv = extract(
        '(abs(strlen(REF)-strlen(ALT))>=50)',
        "pggb.sv_ge50.vcf.gz",
    )
    return {"snp": snp, "indel": indel, "sv": sv, "all": renamed}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--syri-vcf", required=True, type=Path)
    parser.add_argument("--pggb", required=True, type=Path)
    parser.add_argument("--mm-fa", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()

    work = args.outdir / "prepared_variants"
    work.mkdir(parents=True, exist_ok=True)
    fai = Path(str(args.mm_fa) + ".fai")
    if not fai.exists():
        run(["samtools", "faidx", str(args.mm_fa)])

    syri_vcfs = syri_to_mm_vcfs(args.syri_vcf, work / "syri_mm", load_fai(fai))
    for vcf in syri_vcfs.values():
        bgzip_index(vcf)
    prep_pggb(args.pggb, work / "pggb")


if __name__ == "__main__":
    main()
