#!/bin/bash
# step09：vg giraffe 比对
#   simF2 → 仅 GAM（gfa_vgcall 下游）
#   simF1 → GAM + surject BAM（gfa_dv 下游）
# 用法: bash step09_giraffe_align.sh <simF2|simF1> <sample_id>
set -euo pipefail
source /public/home/xuruiqiang/work/xinpian/genotype/env.sh

# 本步参数：GIRAFFE_GBZ SKIP_EXISTING=${SKIP_EXISTING} SLURM_CPUS=${SLURM_GFA_CPUS}

COHORT="${1:?usage: step09_giraffe_align.sh simF2|simF1 sample_id}"
SID="${2:?usage: step09_giraffe_align.sh simF2|simF1 sample_id}"
ROOT="$(cohort_root "${COHORT}")"
mapfile -t FQ < <(sample_fq_paths "${COHORT}" "${SID}")
FQ1="${FQ[0]}" FQ2="${FQ[1]}"
GAM="${ROOT}/gfa_vgcall/01_gam/${SID}.gam"
BAM="${ROOT}/gfa_dv/01_bam/${SID}.bam"
TMP="${ROOT}/gfa_align/tmp_${SID}"
NCPU="${SLURM_CPUS_PER_TASK:-${SLURM_GFA_CPUS}}"
ST_THREADS="$(samtools_threads_cap "${NCPU}" "${GFA_SAMTOOLS_THREADS}")"
RG_LINE="@RG\\tID:${SID}\\tSM:${SID}\\tPL:ILLUMINA"

# simF2 不做 surject；simF1 需要 surject BAM 供 gfa_dv
DO_SURJECT_BAM=0
[[ "${COHORT}" == "simF1" ]] && DO_SURJECT_BAM=1

need_gam=1 need_bam="${DO_SURJECT_BAM}"
[[ "${SKIP_EXISTING}" == "1" ]] && gam_ok "${GAM}" && need_gam=0
[[ "${SKIP_EXISTING}" == "1" && "${DO_SURJECT_BAM}" == "1" ]] && bam_ok "${BAM}" && [[ -f "${BAM}.bai" ]] && need_bam=0
if [[ "${need_gam}" == "0" && "${need_bam}" == "0" ]]; then
    workflow_log SKIP align "step09_giraffe cohort=${COHORT} sample=${SID}"
    exit 0
fi

require_exec "${VG}"
require_file "${GIRAFFE_GBZ}"
require_file "${GIRAFFE_DIST}"
require_file "${GIRAFFE_MIN}"
require_file "${FQ1}"
require_file "${FQ2}"
require_dir "${ROOT}/gfa_vgcall/01_gam"
if [[ "${DO_SURJECT_BAM}" == "1" ]]; then
    require_exec "${SAMTOOLS}"
    require_file "${REF_GBZ}"
    require_file "${RENAME_SED}"
    require_dir "${ROOT}/gfa_dv/01_bam"
fi

# 预检：清空 temp、注册 trap、指定 vg/samtools 临时目录
prep_work_tmp "${TMP}"
setup_work_tmp_trap "${TMP}"
export TMPDIR="${TMP}/vg"
mkdir -p "${TMPDIR}"
rm -f "${TMP}"/*.part

lock="${GAM}.lock"
exec 200>"${lock}"
flock -n 200 || { echo "ERROR: ${GAM} locked" >&2; exit 3; }

workflow_log INFO align "step09_giraffe start cohort=${COHORT} sample=${SID} cpus=${NCPU} surject=${DO_SURJECT_BAM}"
if [[ "${need_gam}" == "1" ]]; then
    # 执行：vg giraffe → .part，成功后 mv 到 GAM
    GAM_PART="${TMP}/${SID}.gam.part"
    "${VG}" giraffe -Z "${GIRAFFE_GBZ}" -d "${GIRAFFE_DIST}" -m "${GIRAFFE_MIN}" \
        -f "${FQ1}" -f "${FQ2}" -t "${NCPU}" --sample "${SID}" \
        > "${GAM_PART}" 2> "${TMP}/giraffe.log"
    mv -f "${GAM_PART}" "${GAM}"
fi

if [[ "${need_bam}" == "1" ]]; then
    # 执行：surject → rename → sort → addreplacerg → temp BAM（仅 simF1）
    BAM_SORT="${TMP}/sorted.bam"
    BAM_RG="${TMP}/${SID}.bam.part"
    "${VG}" surject -x "${REF_GBZ}" -b -t "${NCPU}" "${GAM}" 2> "${TMP}/surject.log" \
        | "${SAMTOOLS}" view -h - \
        | sed -f "${RENAME_SED}" \
        | "${SAMTOOLS}" sort -T "${TMPDIR}/sort" -@"${ST_THREADS}" -m 2G -O BAM -o "${BAM_SORT}" -
    "${SAMTOOLS}" addreplacerg -r "${RG_LINE}" \
        -@"${ST_THREADS}" -O BAM -o "${BAM_RG}" "${BAM_SORT}"
    mv -f "${BAM_RG}" "${BAM}"
    "${SAMTOOLS}" index -@"${ST_THREADS}" "${BAM}"
fi

if [[ "${DO_SURJECT_BAM}" == "1" ]]; then
    gam_ok "${GAM}" && bam_ok "${BAM}" \
        || { workflow_log FAIL align "step09_giraffe bad output cohort=${COHORT} sample=${SID}"; exit 1; }
else
    gam_ok "${GAM}" \
        || { workflow_log FAIL align "step09_giraffe bad gam cohort=${COHORT} sample=${SID}"; exit 1; }
fi
workflow_log DONE align "step09_giraffe cohort=${COHORT} sample=${SID} gam=${GAM} surject_bam=${DO_SURJECT_BAM}"
