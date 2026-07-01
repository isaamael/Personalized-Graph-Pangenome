#!/bin/bash
# step04：单样本 simF2 — pysam 拼单倍型 + 按染色体并发 ART + pigz 流式压缩
# 用法: bash step04_art_f2.sh simF2_NNN
set -euo pipefail
source /public/home/xuruiqiang/work/xinpian/genotype/env.sh

# 本步参数：N_F2=${N_F2} ART_HAPLO_DEPTH=${ART_HAPLO_DEPTH}x ART_CHR_JOBS=${ART_CHR_JOBS} PIGZ_THREADS=${PIGZ_THREADS}

SID="${1:?usage: step04_art_f2.sh simF2_NNN}"
[[ "${SID}" =~ ^simF2_([0-9]+)$ ]] || { echo "invalid sample: ${SID}" >&2; exit 1; }
TASK_IDX=$((10#${BASH_REMATCH[1]} - 1))
[[ "${TASK_IDX}" -lt "${N_F2}" ]] || { echo "sample index >= N_F2=${N_F2}" >&2; exit 1; }

OUT_R1="${READS_DIR}/simF2/${SID}/${SID}.R1.fq.gz"
OUT_R2="${READS_DIR}/simF2/${SID}/${SID}.R2.fq.gz"
if fq_ok "${OUT_R1}" && fq_ok "${OUT_R2}"; then
    workflow_log SKIP sim "step04_art_f2 sample=${SID}"
    exit 0
fi

# 预检：truth、panel、亲本 FASTA、工具
require_file "${TRUTH_DIR}/${SID}.gameteA.fragments.tsv"
require_file "${PANEL_DIR}/sv_blocks.tsv"
require_file "${MM_FA}"
require_file "${TS_FA}"
require_exec "${PYTHON}"
require_exec "${ART_ILLUMINA}"
require_exec "${PIGZ_BIN}"
require_exec "${SAMTOOLS}"

workflow_log INFO sim "step04_art_f2 start sample=${SID}"
# 执行：pysam 提取单倍型 FASTA（workers=${EXTRACT_WORKERS}）
"${PYTHON}" "${GENOTYPE_00_SCRIPT}/step03_extract_native.py" \
    --truth-dir "${TRUTH_DIR}" --panel-dir "${PANEL_DIR}" --out-dir "${HAPLO_DIR}" \
    --mm-fa "${MM_FA}" --ts-fa "${TS_FA}" --sample "${SID}" --workers 1

FA_A="${HAPLO_DIR}/${SID}.gameteA.fa"
FA_B="${HAPLO_DIR}/${SID}.gameteB.fa"
require_file "${FA_A}" 1000
require_file "${FA_B}" 1000

OUT_DIR="${READS_DIR}/simF2/${SID}"
ART_WORK="${OUT_DIR}/art_work"
SEED_A=$((SEED + TASK_IDX))
SEED_B=$((SEED + TASK_IDX + ART_SEED_OFFSET))
require_dir "${OUT_DIR}"
rm -rf "${ART_WORK}"
mkdir -p "${ART_WORK}"

grep '^>' "${FA_A}" | sed 's/>//' > "${ART_WORK}/chr_list.txt"

# 执行：单线程预建 .fai，避免 xargs 并发 faidx 时索引竞争损坏
"${SAMTOOLS}" faidx "${FA_A}"
"${SAMTOOLS}" faidx "${FA_B}"

# 执行：按染色体并发 ART（-P ${ART_CHR_JOBS}），每染色体 hap_depth=${ART_HAPLO_DEPTH}
run_art_chr() {
    local chr="$1"
    "${SAMTOOLS}" faidx "${FA_A}" "${chr}" > "${ART_WORK}/${SID}.A.${chr}.fa"
    "${SAMTOOLS}" faidx "${FA_B}" "${chr}" > "${ART_WORK}/${SID}.B.${chr}.fa"
    "${ART_ILLUMINA}" -ss "${ART_PLATFORM}" -l 150 -p -f "${ART_HAPLO_DEPTH}" -m 500 -s 50 -rs "${SEED_A}" \
        -i "${ART_WORK}/${SID}.A.${chr}.fa" -o "${ART_WORK}/${SID}.A.${chr}" > "${ART_WORK}/${SID}.A.${chr}.log" 2>&1
    "${ART_ILLUMINA}" -ss "${ART_PLATFORM}" -l 150 -p -f "${ART_HAPLO_DEPTH}" -m 500 -s 50 -rs "${SEED_B}" \
        -i "${ART_WORK}/${SID}.B.${chr}.fa" -o "${ART_WORK}/${SID}.B.${chr}" > "${ART_WORK}/${SID}.B.${chr}.log" 2>&1
}
export -f run_art_chr
export SID FA_A FA_B ART_WORK ART_ILLUMINA ART_PLATFORM ART_HAPLO_DEPTH SEED_A SEED_B SAMTOOLS

xargs -P "${ART_CHR_JOBS}" -I {} bash -c 'run_art_chr "$1"' _ {} < "${ART_WORK}/chr_list.txt"

# 执行：按染色体顺序合并 R1/R2（ART 输出为 ${prefix}1.fq / ${prefix}2.fq），pigz 流式压缩
{
    while read -r chr; do
        [[ -f "${ART_WORK}/${SID}.A.${chr}1.fq" ]] && cat "${ART_WORK}/${SID}.A.${chr}1.fq"
        [[ -f "${ART_WORK}/${SID}.B.${chr}1.fq" ]] && cat "${ART_WORK}/${SID}.B.${chr}1.fq"
    done < "${ART_WORK}/chr_list.txt"
} | "${PIGZ_BIN}" -p "${PIGZ_THREADS}" -c > "${OUT_R1}"

{
    while read -r chr; do
        [[ -f "${ART_WORK}/${SID}.A.${chr}2.fq" ]] && cat "${ART_WORK}/${SID}.A.${chr}2.fq"
        [[ -f "${ART_WORK}/${SID}.B.${chr}2.fq" ]] && cat "${ART_WORK}/${SID}.B.${chr}2.fq"
    done < "${ART_WORK}/chr_list.txt"
} | "${PIGZ_BIN}" -p "${PIGZ_THREADS}" -c > "${OUT_R2}"

rm -rf "${ART_WORK}" "${FA_A}" "${FA_B}" "${FA_A}.fai" "${FA_B}.fai"

fq_ok "${OUT_R1}" && fq_ok "${OUT_R2}" \
    || { workflow_log FAIL sim "step04_art_f2 bad output sample=${SID}"; exit 1; }
workflow_log DONE sim "step04_art_f2 sample=${SID} seed_a=${SEED_A} seed_b=${SEED_B}"
