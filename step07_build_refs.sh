#!/bin/bash
# step07：BWA 索引 + PGGB/giraffe GBZ 构建（记录复用命令；产物已存在则仅预检）
# 用法: bash step07_build_refs.sh          # 默认只检查 + 打印命令
#       RUN=1 bash step07_build_refs.sh    # 用户 review 后显式重跑（慎用）
set -euo pipefail
source /public/home/xuruiqiang/work/xinpian/genotype/env.sh

# 本步参数：RUN=${RUN:-0} REF_ROOT=${REF_ROOT} PGGB_DIR=${PGGB_DIR}

RUN="${RUN:-0}"
PGGB_DIR="${PGGB_DIR:-${REF_ROOT}/10_pggb}"
PGGB_SCRIPT="/public/home/xuruiqiang/work/xinpian/genome_sl6/00_script"

workflow_log INFO refs "step07_build_refs start run=${RUN}"

# 预检：当前流程依赖的索引是否就绪
require_file "${BWA_REF_FA}"
require_file "${BWA_REF_FA}.fai"
require_file "${BWA_REF_FA}.sa"
require_file "${GIRAFFE_GBZ}"
require_file "${REF_GBZ}"
require_file "${GIRAFFE_DIST}"
require_file "${GIRAFFE_MIN}"
require_file "${RENAME_SED}"

cat <<EOF

========== step07 参考构建：复用命令（来自完成记录）==========
完成时间: 2026-06-05 BWA index; 2026-06-18 PGGB/giraffe (Slurm 2321765–2321781)
日志: ${PGGB_DIR}/00_logs/slurm_vg_2321766.out
      ${PGGB_DIR}/00_logs/slurm_vg_post_2321767.out

--- A. SL6 线性参考 BWA 索引 ---
# 参考: ${BWA_REF_FA}
${BWA} index -a bwis "${BWA_REF_FA}"
samtools faidx "${BWA_REF_FA}"

--- B. PGGB 分染色体建图（array 1–12）---
# prep: ${PGGB_SCRIPT}/11.pggb_prep.sh
# 投递整条链: ${PGGB_SCRIPT}/11_submit_pggb_pipeline.sh
# 单染色体命令（chr1 示例，高变 chr4/5/9/11/12 用 -p 95 其余 -p 97）:
#   pggb -i ${PGGB_DIR}/chr_inputs/chr1.fa.gz -o ${PGGB_DIR}/chr1/out \\
#     -s 10000 -l 50000 -p 97 -n 2 -k 47 -K 19 -F 0.001 -f 0 -B 10000000 \\
#     -j 0 -e 0 -G 700,900,1100 -P 1,19,39,3,81,1 -O 0.001 -d 100 -Q Consensus_ \\
#     -V 'MM:100000' -t 40 -T 40

--- C. 合并 GFA + deconstruct VCF ---
# ${PGGB_SCRIPT}/12_vg_combine.sh
${VG} combine \$(awk -F'\\t' '{print \$2}' ${PGGB_DIR}/pggb_run.manifest) > ${PGGB_DIR}/pggb.gfa
${VG} paths -x ${PGGB_DIR}/pggb.gfa -L > ${PGGB_DIR}/path.list
${VG} deconstruct --path-prefix MM --all-snarls --threads 40 ${PGGB_DIR}/pggb.gfa > ${PGGB_DIR}/pggb.vcf

--- D. giraffe 索引 + surject 用 ref GBZ ---
# ${PGGB_SCRIPT}/13_vg_postprocess.sh
${VG} autoindex --workflow giraffe -g ${PGGB_DIR}/pggb.gfa -p ${PGGB_DIR}/pggb -t 40 -T ${PGGB_DIR}/temp -M 100G
${VG} gbwt -Z --set-reference TS --gbz-format -g ${PGGB_DIR}/pggb.ref.gbz ${PGGB_DIR}/pggb.giraffe.gbz
${VG} gbwt -Z --set-reference MM --gbz-format -g ${PGGB_DIR}/pggb.ref_mm.gbz ${PGGB_DIR}/pggb.giraffe.gbz

--- E. surject 染色体 rename（step09 使用）---
# ${PGGB_DIR}/rename_mm.sed  (MM#0#chrN → chr0#chrN)

产物:
  ${BWA_REF_FA}.sa
  ${PGGB_DIR}/pggb.giraffe.gbz  ${PGGB_DIR}/pggb.dist  ${PGGB_DIR}/pggb.min
  ${PGGB_DIR}/pggb.ref_mm.gbz   ${RENAME_SED}
============================================================

EOF

if [[ "${RUN}" != "1" ]]; then
    workflow_log DONE refs "step07_build_refs skip run=0 refs_ok=1"
    exit 0
fi

workflow_log INFO refs "step07_build_refs RUN=1 executing"
require_exec "${BWA}"
require_exec "${SAMTOOLS}"
require_exec "${VG}"

if [[ ! -s "${BWA_REF_FA}.sa" ]]; then
    "${BWA}" index -a bwis "${BWA_REF_FA}"
    "${SAMTOOLS}" faidx "${BWA_REF_FA}"
fi

echo "[step07] PGGB 全链路请用 Slurm 脚本（勿在登录节点跑）:"
echo "  bash ${PGGB_SCRIPT}/11_submit_pggb_pipeline.sh"
workflow_log DONE refs "step07_build_refs RUN=1 bwa_only pggb=manual_slurm"
