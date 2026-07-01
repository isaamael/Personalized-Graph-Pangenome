#!/bin/bash
# step01：SyRI SNP+SV 构建 marker panel（SV≥MIN_SV_BP，合并区间判 SNP 是否在 SV 内）
set -euo pipefail
source /public/home/xuruiqiang/work/xinpian/genotype/env.sh

# 本步参数：MIN_SV_BP=${MIN_SV_BP} PILOT_CHR=${PILOT_CHR:-all} SIM_BIN=${SIM_BIN_SIZE} keep=${SIM_KEEP_PER_BIN}

EXTRA=()
[[ -n "${PILOT_CHR}" ]] && EXTRA+=(--pilot-chr "${PILOT_CHR}")

# 预检：SyRI 统计表、参考索引
require_file "${SYRI_SNP_MM}"
require_file "${SYRI_SV_MM}"
require_file "${MM_FAI}"
require_exec "${PYTHON}"
require_dir "${PANEL_DIR}"

workflow_log INFO sim "step01_build_panel start min_sv_bp=${MIN_SV_BP}"
# 执行：Python 读 syri_sv_positions_MM.tsv + SNP 表，二分判断 breakpoint_ok
"${PYTHON}" "${GENOTYPE_00_SCRIPT}/step01_build_panel.py" \
    --snp-tsv "${SYRI_SNP_MM}" --syri-sv-tsv "${SYRI_SV_MM}" \
    --mm-fai "${MM_FAI}" --out-dir "${PANEL_DIR}" --min-sv-bp "${MIN_SV_BP}" \
    "${EXTRA[@]}"

require_file "${PANEL_DIR}/markers_sim.tsv"
workflow_log DONE sim "step01_build_panel out=${PANEL_DIR}"
