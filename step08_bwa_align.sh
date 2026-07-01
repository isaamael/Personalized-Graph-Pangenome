#!/bin/bash
# step08：BWA-MEM 比对 → sorted BAM（SL6 线性参考）
# 用法: bash step08_bwa_align.sh <simF2|simF1> <sample_id>
set -euo pipefail
source /public/home/xuruiqiang/work/xinpian/genotype/env.sh

# 本步参数：cohort BWA_REF=${BWA_REF_FA} SKIP_EXISTING=${SKIP_EXISTING}

COHORT="${1:?usage: step08_bwa_align.sh simF2|simF1 sample_id}"
SID="${2:?usage: step08_bwa_align.sh simF2|simF1 sample_id}"
ROOT="$(cohort_root "${COHORT}")"
mapfile -t FQ < <(sample_fq_paths "${COHORT}" "${SID}")
FQ1="${FQ[0]}" FQ2="${FQ[1]}"
BAM="${ROOT}/bwa_dv/01_bwa/${SID}.bam"
TMP="${ROOT}/bwa_dv/tmp_${SID}"

if [[ "${SKIP_EXISTING}" == "1" ]] && bam_ok "${BAM}"; then
    workflow_log SKIP align "step08_bwa cohort=${COHORT} sample=${SID}"
    exit 0
fi

require_exec "${BWA}"
require_exec "${SAMTOOLS}"
require_file "${BWA_REF_FA}"
require_file "${BWA_REF_FA}.sa"
require_file "${FQ1}"
require_file "${FQ2}"
require_dir "${ROOT}/bwa_dv/01_bwa"

# 预检：清空并注册 temp 清理 trap
prep_work_tmp "${TMP}"
setup_work_tmp_trap "${TMP}"
export TMPDIR="${TMP}/samtools"
mkdir -p "${TMPDIR}"

workflow_log INFO align "step08_bwa start cohort=${COHORT} sample=${SID} tmp=${TMP}"
# 执行：bwa mem | samtools view | sort → temp BAM，成功后 mv 到最终路径
RG="@RG\\tID:${SID}\\tSM:${SID}\\tPL:ILLUMINA"
BAM_PART="${TMP}/${SID}.bam.part"
"${BWA}" mem -t "${BWA_MEM_THREADS}" -M -R "${RG}" "${BWA_REF_FA}" "${FQ1}" "${FQ2}" \
    | "${SAMTOOLS}" view -bS -@"${BWA_VIEW_THREADS}" - \
    | "${SAMTOOLS}" sort -T "${TMPDIR}/sort" -@"${BWA_SORT_THREADS}" -m 2G -o "${BAM_PART}" -
mv -f "${BAM_PART}" "${BAM}"
"${SAMTOOLS}" index -@"${BWA_INDEX_THREADS}" "${BAM}"

bam_ok "${BAM}" || { workflow_log FAIL align "step08_bwa bad bam cohort=${COHORT} sample=${SID}"; exit 1; }
workflow_log DONE align "step08_bwa cohort=${COHORT} sample=${SID} bam=${BAM}"
