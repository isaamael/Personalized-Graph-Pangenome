#!/bin/bash
# step10 QC2：比对产物质检（mode=fast|deep，默认 fast）
# fast：全 cohort 文件清单 + 抽检 flagstat/GAM 流式；deep 含 mpileup marker AF
# 用法: bash step10_qc2_align.sh [sample1,sample2,...]
set -euo pipefail
source /public/home/xuruiqiang/work/xinpian/genotype/env.sh

QC2_MODE="${QC2_MODE:-fast}"
QC2_OUT="${QC2_OUT:-${COMMON_ROOT}/qc/QC2_align_summary.tsv}"
QC2_AF_OUT="${QC2_AF_OUT:-${COMMON_ROOT}/qc/QC2_marker_af.tsv}"
QC2_INV_OUT="${QC2_INV_OUT:-${COMMON_ROOT}/qc/QC2_align_inventory.tsv}"
[[ "${QC2_MODE}" == "deep" ]] && QC2_OUT="${QC2_OUT%.tsv}.deep.tsv" && QC2_AF_OUT="${QC2_AF_OUT%.tsv}.deep.tsv"
QC2_TMP="${QC2_TMP:-${GENOTYPE_TMP_ROOT}/qc2}"
QC2_SPOT_F2="${QC2_SPOT_F2:-3}"
QC2_SPOT_F1="${QC2_SPOT_F1:-3}"
QC2_N_MARKERS="${QC2_N_MARKERS:-40}"
QC2_MAX_GAM="${QC2_MAX_GAM_RECORDS:-30000}"
EXTRA=(--mode "${QC2_MODE}" --inventory-out "${QC2_INV_OUT}")
[[ "${QC2_MODE}" == "fast" ]] && EXTRA+=(--map-min 0.90)
[[ $# -ge 1 ]] && EXTRA+=(--samples "$1")
[[ "${QC2_MODE}" == "deep" ]] && EXTRA+=(--n-markers "${QC2_N_MARKERS}")

require_exec "${PYTHON}"
require_exec "${SAMTOOLS}"
require_exec "${VG}"
require_file "${SL6_REF_FA}"
require_dir "${COMMON_ROOT}/qc"
[[ "${QC2_MODE}" == "deep" ]] && require_file "${PANEL_DIR}/markers_sim.tsv"
prep_work_tmp "${QC2_TMP}"

workflow_log INFO align "step10_qc2 start mode=${QC2_MODE} spot_f2=${QC2_SPOT_F2} spot_f1=${QC2_SPOT_F1}"
"${PYTHON}" "${GENOTYPE_00_SCRIPT}/step10_qc2_align.py" \
    --panel "${PANEL_DIR}/markers_sim.tsv" \
    --ref-fa "${SL6_REF_FA}" \
    --out-tsv "${QC2_OUT}" \
    --out-af-tsv "${QC2_AF_OUT}" \
    --sim-f2-root "${SIM_F2_ROOT}" \
    --sim-f1-root "${SIM_F1_ROOT}" \
    --reads-dir "${READS_DIR}" \
    --genome-size "${GENOME_SIZE_MM}" \
    --art-haplo-depth "${ART_HAPLO_DEPTH}" \
    --pf-depth "${PF_DEPTH}" \
    --pf-depth-label "${PF_DEPTH}" \
    --n-f2 "${N_F2}" --n-pf "${N_PF}" \
    --spot-f2 "${QC2_SPOT_F2}" --spot-f1 "${QC2_SPOT_F1}" \
    --seed "${SEED}" \
    --max-gam-records "${QC2_MAX_GAM}" \
    --tmp-dir "${QC2_TMP}" \
    "${EXTRA[@]}"

require_file "${QC2_OUT}" 50
require_file "${QC2_INV_OUT}" 50
workflow_log DONE align "step10_qc2 mode=${QC2_MODE} out=${QC2_OUT} inv=${QC2_INV_OUT}"
