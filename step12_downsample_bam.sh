#!/bin/bash
# step12：从全深度 BAM 按 read-name 下采样（samtools view -N）
# 用法: bash step12_downsample_bam.sh <bwa|gfa_dv> <simF2|simF1> <sample_id>
set -euo pipefail
source /public/home/xuruiqiang/work/xinpian/genotype/env.sh

STRATEGY="${1:?usage: step12_downsample_bam.sh bwa|gfa_dv simF2|simF1 sample_id}"
COHORT="${2:?usage: step12_downsample_bam.sh bwa|gfa_dv simF2|simF1 sample_id}"
SID="${3:?usage: step12_downsample_bam.sh bwa|gfa_dv simF2|simF1 sample_id}"

if [[ "${STRATEGY}" == "gfa_dv" && "${COHORT}" == "simF2" ]]; then
    echo "ERROR: simF2 无 surject BAM，gfa_dv 仅 simF1" >&2
    exit 2
fi

case "${STRATEGY}" in
    bwa) DEPTHS="${DSAMPLE_DEPTHS_BWA}" ;;
    gfa_dv) DEPTHS="${DSAMPLE_DEPTHS_GFA}" ;;
    *) echo "unknown strategy: ${STRATEGY}" >&2; exit 1 ;;
esac

NDIR="$(cohort_names_dir "${COHORT}")"
SRC="$(downsample_bam_src "${STRATEGY}" "${COHORT}" "${SID}")"
TMP="$(cohort_root "${COHORT}")/downsample/tmp_${STRATEGY}_${SID}"
ST_SORT="${DSAMPLE_SORT_THREADS}"
ST_VIEW="${DSAMPLE_VIEW_THREADS}"

require_exec "${SAMTOOLS}"
require_file "${SRC}"
require_file "${SRC}.bai"
bam_ok "${SRC}" || { echo "ERROR: source BAM missing: ${SRC}" >&2; exit 2; }
names_ready "${NDIR}" "${SID}" "${DEPTHS}" || {
    echo "ERROR: 先跑 step11_gen_downsample_names.sh ${COHORT} ${SID}" >&2
    exit 2
}

prep_work_tmp "${TMP}"
setup_work_tmp_trap "${TMP}"
export TMPDIR="${TMP}/samtools"
mkdir -p "${TMPDIR}"

workflow_log INFO downsample "step12_bam start strategy=${STRATEGY} cohort=${COHORT} sample=${SID}"
for depth in ${DEPTHS}; do
    dstr="$(depth_label "${depth}")"
    OUT="$(downsample_bam_out "${STRATEGY}" "${COHORT}" "${SID}" "${depth}")"
    require_dir "$(dirname "${OUT}")"

    if [[ "${SKIP_EXISTING}" == "1" ]] && ds_bam_ok "${OUT}" && [[ -f "${OUT}.bai" ]]; then
        workflow_log SKIP downsample "step12_bam exists strategy=${STRATEGY} sample=${SID} depth=${dstr}x"
        continue
    fi

    NAMES="$(names_file_for "${NDIR}" "${SID}" "${depth}")"
    NAMES_EXP="${TMP}/names_${dstr}x.txt"
    PART="${TMP}/sorted_${dstr}x.bam.part"

    expand_read_names_bam "${NAMES}" > "${NAMES_EXP}"
    "${SAMTOOLS}" view -@"${ST_VIEW}" -b -N "${NAMES_EXP}" "${SRC}" \
        | "${SAMTOOLS}" sort -T "${TMPDIR}/sort_${dstr}" -@"${ST_SORT}" -m 1G \
            -o "${PART}" -
    rm -f "${NAMES_EXP}"
    mv -f "${PART}" "${OUT}"
    "${SAMTOOLS}" index -@"${ST_SORT}" "${OUT}"

    ds_bam_ok "${OUT}" || { workflow_log FAIL downsample "step12_bam bad out=${OUT}"; exit 1; }
done
workflow_log DONE downsample "step12_bam strategy=${STRATEGY} cohort=${COHORT} sample=${SID}"
