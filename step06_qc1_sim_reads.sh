#!/bin/bash
# step06 QC1：truth CO + reads 质检（mode=fast|deep，默认 fast）
set -euo pipefail
source /public/home/xuruiqiang/work/xinpian/genotype/env.sh

QC1_MODE="${QC1_MODE:-fast}"
QC1_OUT="${QC1_OUT:-${COMMON_ROOT}/qc/QC1_sim_reads.tsv}"
[[ "${QC1_MODE}" == "deep" ]] && QC1_OUT="${QC1_OUT%.tsv}.deep.tsv"
QC1_TRUTH_OUT="${QC1_TRUTH_OUT:-${COMMON_ROOT}/qc/QC1_truth_co.tsv}"
QC1_MAX_PAIRS="${QC1_MAX_PAIRS:-20000}"
QC1_DEPTH_TOL="${QC1_DEPTH_TOL:-0.20}"
QC1_SPOT_F2="${QC1_SPOT_F2:-3}"
QC1_SPOT_F1="${QC1_SPOT_F1:-3}"
EXTRA=(--mode "${QC1_MODE}" --spot-f2 "${QC1_SPOT_F2}" --spot-f1 "${QC1_SPOT_F1}" --seed "${SEED}")
[[ $# -ge 1 ]] && EXTRA+=(--samples "$1")
[[ -n "${QC1_TRUTH_ONLY:-}" ]] && EXTRA+=(--truth-only)
[[ -n "${QC1_READS_ONLY:-}" ]] && EXTRA+=(--reads-only)

require_exec "${PYTHON}"
require_dir "${COMMON_ROOT}/qc"
require_file "${TRUTH_DIR}/truth_co_summary.tsv"

workflow_log INFO sim "step06_qc1 start mode=${QC1_MODE} truth=${TRUTH_DIR} reads=${READS_DIR}"
"${PYTHON}" "${GENOTYPE_00_SCRIPT}/step06_qc1_sim_reads.py" \
    --reads-dir "${READS_DIR}" --out-tsv "${QC1_OUT}" \
    --truth-dir "${TRUTH_DIR}" --truth-out-tsv "${QC1_TRUTH_OUT}" \
    --genome-size "${GENOME_SIZE_MM}" \
    --art-haplo-depth "${ART_HAPLO_DEPTH}" --pf-depth "${PF_DEPTH}" \
    --pf-depth-label "${PF_DEPTH}" --n-f2 "${N_F2}" --n-pf "${N_PF}" \
    --max-pairs "${QC1_MAX_PAIRS}" --depth-tol "${QC1_DEPTH_TOL}" \
    "${EXTRA[@]}"

require_file "${QC1_TRUTH_OUT}" 100
workflow_log DONE sim "step06_qc1 mode=${QC1_MODE} truth=${QC1_TRUTH_OUT} reads=${QC1_OUT}"
