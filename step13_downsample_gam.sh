#!/bin/bash
# step13：GAM 下采样（vg filter -N -e；name 须与 GAM 完全一致）
# GAM read name 带 /1 /2 → expand_read_names 将 base ID 扩展为两行
# 用法: bash step13_downsample_gam.sh <simF2|simF1> <sample_id>
set -euo pipefail
source /public/home/xuruiqiang/work/xinpian/genotype/env.sh

COHORT="${1:?usage: step13_downsample_gam.sh simF2|simF1 sample_id}"
SID="${2:?usage: step13_downsample_gam.sh simF2|simF1 sample_id}"
NDIR="$(cohort_names_dir "${COHORT}")"
SRC="$(downsample_gam_src "${COHORT}" "${SID}")"
TMP="$(cohort_root "${COHORT}")/downsample/tmp_gam_${SID}"
NCPU="${SLURM_CPUS_PER_TASK:-${SLURM_DSAMPLE_CPUS}}"

require_exec "${VG}"
require_file "${SRC}"
gam_ok "${SRC}" || { echo "ERROR: source GAM missing: ${SRC}" >&2; exit 2; }
names_ready "${NDIR}" "${SID}" "${DSAMPLE_DEPTHS_GFA}" || {
    echo "ERROR: 先跑 step11_gen_downsample_names.sh ${COHORT} ${SID}" >&2
    exit 2
}

prep_work_tmp "${TMP}"
setup_work_tmp_trap "${TMP}"
export TMPDIR="${TMP}/vg"
mkdir -p "${TMPDIR}"

workflow_log INFO downsample "step13_gam start cohort=${COHORT} sample=${SID}"
for depth in ${DSAMPLE_DEPTHS_GFA}; do
    dstr="$(depth_label "${depth}")"
    OUT="$(downsample_gam_out "${COHORT}" "${SID}" "${depth}")"
    require_dir "$(dirname "${OUT}")"

    if [[ "${SKIP_EXISTING}" == "1" ]] && ds_gam_ok "${OUT}"; then
        workflow_log SKIP downsample "step13_gam exists sample=${SID} depth=${dstr}x"
        continue
    fi

    NAMES="$(names_file_for "${NDIR}" "${SID}" "${depth}")"
    NAMES_EXP="${TMP}/names_${dstr}x.txt"
    PART="${TMP}/${dstr}x.gam.part"
    SRC_BYTES="$(stat -c%s "${SRC}")"

    expand_read_names_gam "${NAMES}" > "${NAMES_EXP}"
    N_NAMES="$(wc -l < "${NAMES_EXP}")"
    [[ "${N_NAMES}" -gt 0 ]] || { echo "ERROR: empty names list ${NAMES_EXP}" >&2; exit 2; }

    "${VG}" filter -N "${NAMES_EXP}" -e -t "${NCPU}" "${SRC}" > "${PART}" 2> "${TMP}/filter_${dstr}.log"
    rm -f "${NAMES_EXP}"
    mv -f "${PART}" "${OUT}"

    OUT_BYTES="$(stat -c%s "${OUT}" 2>/dev/null || echo 0)"
    if [[ "${OUT_BYTES}" -lt $((DS_GAM_MIN_MB * 1024 * 1024)) ]]; then
        workflow_log FAIL downsample "step13_gam empty out=${OUT} names=${N_NAMES} src_mb=$((SRC_BYTES/1024/1024))"
        echo "ERROR: GAM filter output too small; check name suffix match (vg view -a ... | head)" >&2
        exit 1
    fi
    ds_gam_ok "${OUT}" || { workflow_log FAIL downsample "step13_gam bad out=${OUT}"; exit 1; }
done
workflow_log DONE downsample "step13_gam cohort=${COHORT} sample=${SID}"
