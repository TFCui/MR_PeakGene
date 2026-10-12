#!/bin/bash
#SBATCH --job-name=MR_AllPairs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=360G
#SBATCH --time=04:00:00
#SBATCH --array=0-2
#SBATCH --output=/scratch/user/tfcui2025/CIMA/Claude_project/Scripts/Model/logs/%x_%A_%a.out

# 全量 peak-gene pairs 的 MR / SMR / 边际相关, 一个 array task 一个细胞类型 (流程见 Scripts/Pipeline.md 第 2 节)
#
# 提交前先建目录 (#SBATCH --output 在投递时就要存在):
#   cd /scratch/user/tfcui2025/CIMA/Claude_project && mkdir -p Scripts/Model/logs Results/MR tmp
#
# 冒烟测试 (只跑 Bn 前 200 个 pair): sbatch --array=0 --export=ALL,TEST_N=200 --time=01:00:00 Scripts/Model/MR_AllPairs_slurm.sh
# 正式运行 (3 个 CT):                sbatch Scripts/Model/MR_AllPairs_slurm.sh
# 中断后续跑: 重投同一命令即可 (已完成的 chunk 会跳过)

set -euo pipefail
module purge
module load GCC/13.3.0 OpenMPI/5.0.3 ArchR/1.0.3-R-4.4.2 SeuratDisk/20231104-R-4.4.2 R_tamu/4.4.2

PROJ=/scratch/user/tfcui2025/CIMA/Claude_project
CELLS=(Bn_TCL1A CD4_Tn_CCR7 CD8_Tn_CCR7)
CT=${CELLS[$SLURM_ARRAY_TASK_ID]}
export TEST_N=${TEST_N:-0}

echo "cell_type = $CT | TEST_N = $TEST_N | cores = $SLURM_CPUS_PER_TASK | job = ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
Rscript $PROJ/Scripts/Model/MR_AllPairs.R $CT
