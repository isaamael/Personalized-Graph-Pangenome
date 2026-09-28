#!/usr/bin/env python3
"""将 SyRI 和 PGGB 变异整理为可比较的、以 MM 为参考的 VCF 文件。

SyRI VCF 以 TS 为参考。通过 INFO 中的 ChrB/StartB 将序列型短变异（ShV，SNP/插入缺失）
重新映射到 MM，并交换等位基因（TS↔MM），以匹配 PGGB 的 REF=MM、ALT=TS。
对于使用符号表示的结构变异，尽可能转换为 MM 坐标下带 END 字段的 INS/DEL，
用于按位置和长度进行模糊匹配。
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
    """将 SyRI VCF 按 SNP、<50 bp 插入缺失及 ≥50 bp 插入缺失/结构变异输出为以 MM 为参考的 VCF。"""
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
                    # 仅保留 chr1–chr12 格式；如有需要，随后根据 MM 的 FAI 索引重写
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

            # 序列型等位基因：SNP / 插入缺失
            if not alt.startswith("<") and all(x in "ACGTNacgtn" for x in ref + alt.replace(",", "")):
                # 交换等位基因，使 MM 成为参考
                mm_ref, mm_alt = alt, ref
                # 对于多碱基变异，SyRI INS（查询序列中的插入）在以 TS 为参考时具有较长的 ALT；
                # 交换后，MM 中较长的 REF 表示相对于 MM 的缺失（DEL），因为插入序列位于 MM 中。
                # 当 SyRI 的 REF=TS 时，INS 表示查询序列 MM 相对 TS 多出一段序列；
                # 若插入缺失跨越不同坐标，仅交换等位基因并不能完成以 MM 为参考的转换。
                # 对于两个等位基因均为单碱基的 ShV SNP，交换可准确完成转换。对于插入缺失，使用 MM 的 StartB 和变异长度，
                # 并在两个等位基因均以序列表示时交换等位基因。
                sz = abs(len(mm_ref) - len(mm_alt))
                rec = f"{qb}\t{mm_pos}\t{vid}\t{mm_ref}\t{mm_alt}\t.\tPASS\t{orig}\tGT\t{gt}\n"
                if len(mm_ref) == 1 and len(mm_alt) == 1:
                    rows["snp"].append((qb, mm_pos, rec))
                elif sz < 50:
                    rows["indel"].append((qb, mm_pos, rec))
                else:
                    rows["sv"].append((qb, mm_pos, rec))
                continue

            # 通过 ChrB 将符号型 INS/DEL 定位到 MM；由 TS 参考转换为 MM 参考时交换类型
            # SyRI 的 INS（MM 相对 TS 的插入）与 MM 参考下的 DEL 之间的对应关系需注意参考方向；
            # 保留 ALT 的符号表示，并将其定位到 MM 坐标，用于模糊匹配。
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
                # 按 MM 参考方向交换标签（与一致性计算脚本保持一致）
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
    """将 MM#0#chr* 重命名为 chr*，并按变异大小分组。"""
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

    # SNP：两个等位基因的长度均为 1
    snp = extract('strlen(REF)=1 && strlen(ALT)=1', "pggb.snp.vcf.gz")
    # 长度 <50 bp 的插入缺失
    indel = extract(
        '((strlen(REF)!=strlen(ALT)) && abs(strlen(REF)-strlen(ALT))<50 && strlen(REF)>0 && strlen(ALT)>0 && (strlen(REF)>1 || strlen(ALT)>1))',
        "pggb.indel_lt50.vcf.gz",
    )
    # 以序列表示、长度 ≥50 bp 的结构变异
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
