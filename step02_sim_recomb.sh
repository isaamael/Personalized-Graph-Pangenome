#!/bin/bash
# step02：PedigreeSim 模拟 N_F2 个 F2 重组，输出 fragments + breakpoints
set -euo pipefail
source /public/home/xuruiqiang/work/xinpian/genotype/env.sh

# 本步参数：N_F2=${N_F2} SEED=${SEED} PILOT_CHR=${PILOT_CHR:-all}

EXTRA=()
[[ -n "${PILOT_CHR}" ]] && EXTRA+=(--pilot-chr "${PILOT_CHR}")

# 预检：panel、PedigreeSim、参考索引
require_file "${PEDIGREESIM_JAR}"
require_file "${PANEL_DIR}/markers_sim.tsv"
require_file "${MM_FAI}"
require_file "${TS_FAI}"
require_exec "${PYTHON}"
require_exec "${JAVA}"
require_dir "${TRUTH_DIR}"

workflow_log INFO sim "step02_sim_recomb start n_f2=${N_F2} seed=${SEED}"
# 执行：PedigreeSim → fragments + parent-CO 真值（summary/detail/breakpoints）
"${PYTHON}" "${GENOTYPE_00_SCRIPT}/step02_sim_recomb.py" \
    --panel-dir "${PANEL_DIR}" --out-dir "${TRUTH_DIR}" \
    --mm-fai "${MM_FAI}" --ts-fai "${TS_FAI}" \
    --n-f2 "${N_F2}" --seed "${SEED}" --java "${JAVA}" \
    --pedigreesim-jar "${PEDIGREESIM_JAR}" --pedigreesim-dir "${PEDIGREESIM_DIR}" \
    "${EXTRA[@]}"

require_file "${TRUTH_DIR}/simF2_001.gameteA.fragments.tsv"
require_file "${TRUTH_DIR}/truth_co_summary.tsv"
require_file "${TRUTH_DIR}/truth_co_cohort.tsv"
workflow_log DONE sim "step02_sim_recomb n_f2=${N_F2} out=${TRUTH_DIR}"
