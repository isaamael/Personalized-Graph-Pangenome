#!/bin/bash
# step05：pseudoF1 — 同 seed 比例抽样混亲本 reads，read 名加 MM_/TS_ 前缀避免 ID 冲突
# 用法: bash step05_pseudo_f1.sh [pseudoF1_NN]
set -euo pipefail
source /public/home/xuruiqiang/work/xinpian/genotype/env.sh

# 本步参数：N_PF=${N_PF} PF_DEPTH=${PF_DEPTH}x GENOME_SIZE_MM=${GENOME_SIZE_MM} SEQKIT_JOBS=${SEQKIT_JOBS}

run_one() {
    local sample="$1" p1r1="$2" p1r2="$3" p2r1="$4" p2r2="$5" seed="$6"
    local pf_out pf_work tmp out_r1 out_r2
    pf_out="${READS_DIR}/pseudoF1"
    pf_work="${pf_out}/work"
    tmp="${pf_work}/${sample}"
    out_r1="${pf_out}/${PF_DEPTH}x/${sample}.R1.fq.gz"
    out_r2="${pf_out}/${PF_DEPTH}x/${sample}.R2.fq.gz"

    require_dir "${tmp}"
    require_dir "${pf_out}/${PF_DEPTH}x"
    for f in "${p1r1}" "${p1r2}" "${p2r1}" "${p2r2}"; do require_file "${f}"; done

    if fq_ok "${out_r1}" && fq_ok "${out_r2}"; then
        workflow_log SKIP sim "step05_pseudo_f1 exists sample=${sample}"
        return 0
    fi

    local n1 n2 frac1 frac2
    n1=$("${SEQKIT}" stat -j "${SEQKIT_JOBS}" -T "${p1r1}" | awk 'NR>1 {sum+=$4} END {print sum+0}')
    n2=$("${SEQKIT}" stat -j "${SEQKIT_JOBS}" -T "${p2r1}" | awk 'NR>1 {sum+=$4} END {print sum+0}')
    frac1=$(awk -v a="${n1}" -v b="${n2}" 'BEGIN{n=(a<b)?a:b; print (a>0)?n/a:0}')
    frac2=$(awk -v a="${n1}" -v b="${n2}" 'BEGIN{n=(a<b)?a:b; print (b>0)?n/b:0}')

    # 执行：同 seed 比例抽样 + 行首加 MM_/TS_ 前缀（seqkit replace -p '^'），R1/R2 配对一致
    "${SEQKIT}" sample -p "${frac1}" -s "${seed}" -j "${SEQKIT_JOBS}" "${p1r1}" \
        | "${SEQKIT}" replace -p '^' -r 'MM_' | "${PIGZ_BIN}" -p "${PIGZ_THREADS}" > "${tmp}/p1_R1.fq.gz"
    "${SEQKIT}" sample -p "${frac1}" -s "${seed}" -j "${SEQKIT_JOBS}" "${p1r2}" \
        | "${SEQKIT}" replace -p '^' -r 'MM_' | "${PIGZ_BIN}" -p "${PIGZ_THREADS}" > "${tmp}/p1_R2.fq.gz"
    "${SEQKIT}" sample -p "${frac2}" -s "${seed}" -j "${SEQKIT_JOBS}" "${p2r1}" \
        | "${SEQKIT}" replace -p '^' -r 'TS_' | "${PIGZ_BIN}" -p "${PIGZ_THREADS}" > "${tmp}/p2_R1.fq.gz"
    "${SEQKIT}" sample -p "${frac2}" -s "${seed}" -j "${SEQKIT_JOBS}" "${p2r2}" \
        | "${SEQKIT}" replace -p '^' -r 'TS_' | "${PIGZ_BIN}" -p "${PIGZ_THREADS}" > "${tmp}/p2_R2.fq.gz"

    # 执行：gzip 物理拼接（cat），避免 pigz 解压再压缩
    cat "${tmp}/p1_R1.fq.gz" "${tmp}/p2_R1.fq.gz" > "${tmp}/mix_R1.fq.gz"
    cat "${tmp}/p1_R2.fq.gz" "${tmp}/p2_R2.fq.gz" > "${tmp}/mix_R2.fq.gz"

    local total_bases fraction
    total_bases=$("${SEQKIT}" stat -j "${SEQKIT_JOBS}" -T "${tmp}/mix_R1.fq.gz" "${tmp}/mix_R2.fq.gz" \
        | awk 'NR>1 {sum+=$5} END {print sum+0}')
    fraction=$(awk -v d="${PF_DEPTH}" -v g="${GENOME_SIZE_MM}" -v t="${total_bases}" 'BEGIN{print (t>0)?d*g/t:1}')

    if awk -v f="${fraction}" 'BEGIN{exit !(f>=1)}'; then
        cp -f "${tmp}/mix_R1.fq.gz" "${out_r1}"
        cp -f "${tmp}/mix_R2.fq.gz" "${out_r2}"
    else
        # 执行：混合后按目标深度 ${PF_DEPTH}x 再抽样（seqkit 直读 .gz，同 seed 保配对）
        "${SEQKIT}" sample -p "${fraction}" -s "${seed}" -j "${SEQKIT_JOBS}" "${tmp}/mix_R1.fq.gz" \
            -o "${out_r1}"
        "${SEQKIT}" sample -p "${fraction}" -s "${seed}" -j "${SEQKIT_JOBS}" "${tmp}/mix_R2.fq.gz" \
            -o "${out_r2}"
    fi
    rm -rf "${tmp}"

    fq_ok "${out_r1}" && fq_ok "${out_r2}" \
        || { workflow_log FAIL sim "step05_pseudo_f1 bad output sample=${sample}"; return 1; }
    workflow_log DONE sim "step05_pseudo_f1 sample=${sample} depth=${PF_DEPTH}x seed=${seed}"
}

require_file "${PF_PAIRS_TSV}"
require_exec "${SEQKIT}"
require_exec "${PIGZ_BIN}"

if [[ $# -ge 1 ]]; then
    SID="$1"
    line=$(awk -F'\t' -v s="${SID}" 'NR>1 && $1==s {print; exit}' "${PF_PAIRS_TSV}")
    [[ -n "${line}" ]] || { echo "not in pairs: ${SID}" >&2; exit 1; }
    IFS=$'\t' read -r sample p1r1 p1r2 p2r1 p2r2 seed <<< "${line}"
    run_one "${sample}" "${p1r1}" "${p1r2}" "${p2r1}" "${p2r2}" "${seed}"
else
    workflow_log INFO sim "step05_pseudo_f1 batch n_pf=${N_PF} depth=${PF_DEPTH}x"
    tail -n +2 "${PF_PAIRS_TSV}" | while IFS=$'\t' read -r sample p1r1 p1r2 p2r1 p2r2 seed; do
        run_one "${sample}" "${p1r1}" "${p1r2}" "${p2r1}" "${p2r2}" "${seed}"
    done
    workflow_log DONE sim "step05_pseudo_f1 batch depth=${PF_DEPTH}x"
fi
