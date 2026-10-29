#!/bin/bash
#SBATCH --job-name=QC2_FlipScan
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=200G
#SBATCH --time=04:00:00
#SBATCH --array=0-2
#SBATCH --output=/scratch/user/tfcui2025/CIMA/Claude_project/Scripts/QC/logs/%x_%A_%a.out
# QC2 第 1 步: 全基因组 MR-Egger flip 前后对比, 一个 array task 一个细胞类型 (流程见 Scripts/Pipeline.md 第 4.2 节)
# 冒烟测试: sbatch --array=0 --export=ALL,TEST_N=200 --cpus-per-task=8 --mem=64G --time=01:00:00 Scripts/QC/QC2_1_FlipScan_slurm.sh
# 正式运行: sbatch Scripts/QC/QC2_1_FlipScan_slurm.sh   (中断后重投同一命令即可续跑)

set -euo pipefail
module purge
module load GCC/13.3.0 OpenMPI/5.0.3 ArchR/1.0.3-R-4.4.2 SeuratDisk/20231104-R-4.4.2 R_tamu/4.4.2

PROJ=/scratch/user/tfcui2025/CIMA/Claude_project
CELLS=(Bn_TCL1A CD4_Tn_CCR7 CD8_Tn_CCR7)
CT=${CELLS[$SLURM_ARRAY_TASK_ID]}
export TEST_N=${TEST_N:-0}
echo "cell_type = $CT | TEST_N = $TEST_N | cores = $SLURM_CPUS_PER_TASK | job = ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
Rscript $PROJ/Scripts/QC/QC2_1_FlipScan.R $CT
