#!/bin/bash
# step11：双端配对 + shuffle 嵌套 name 列表（登录节点）
# 用法: bash step11_gen_downsample_names.sh <simF2|simF1> <sample_id>
# 前置：QC1 pair_pass=PASS（step06 已检 R1/R2 计数与 head 配对）
set -euo pipefail
source /public/home/xuruiqiang/work/xinpian/genotype/env.sh

COHORT="${1:?usage: step11_gen_downsample_names.sh simF2|simF1 sample_id}"
SID="${2:?usage: step11_gen_downsample_names.sh simF2|simF1 sample_id}"
NDIR="$(cohort_names_dir "${COHORT}")"
mapfile -t FQ < <(sample_fq_paths "${COHORT}" "${SID}")
EXP_DEPTH="$(cohort_expected_depth "${COHORT}")"
DEPTHS="$(echo ${DSAMPLE_DEPTHS_BWA} ${DSAMPLE_DEPTHS_GFA} | tr ' ' '\n' | sort -nu | tr '\n' ' ')"

require_exec "${PYTHON}"
require_file "${FQ[0]}"
require_file "${FQ[1]}"
require_dir "${NDIR}"

if [[ "${SKIP_EXISTING}" == "1" ]] && names_ready "${NDIR}" "${SID}" "${DEPTHS}"; then
    workflow_log SKIP downsample "step11_names cohort=${COHORT} sample=${SID}"
    exit 0
fi

workflow_log INFO downsample "step11_names start cohort=${COHORT} sample=${SID}"
"${PYTHON}" "${GENOTYPE_00_SCRIPT}/step11_gen_downsample_names.py" \
    --sample "${SID}" --r1 "${FQ[0]}" --r2 "${FQ[1]}" \
    --names-dir "${NDIR}" --expected-depth "${EXP_DEPTH}" \
    --depths "${DEPTHS}" --seed "${SEED}"

names_ready "${NDIR}" "${SID}" "${DEPTHS}" \
    || { workflow_log FAIL downsample "step11_names incomplete cohort=${COHORT} sample=${SID}"; exit 1; }
workflow_log DONE downsample "step11_names cohort=${COHORT} sample=${SID}"
