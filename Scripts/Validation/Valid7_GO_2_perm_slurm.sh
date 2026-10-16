#!/bin/bash
#SBATCH --job-name=Valid7_GO_2_perm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --array=0-3
#SBATCH --output=/scratch/user/tfcui2025/CIMA/Claude_project/Scripts/Validation/logs/%x_%A_%a.out

# Valid7_GO 第 2 步: 随机抽样 GSEA, array 0-3 = Bn / CD4 / CD8 / Combined, 每个 task N_PERM 次 (默认 5000)
# 依赖第 1 步 (Valid7_GO_1_gsea); 完成后: sbatch --dependency=afterok:<jobid> Valid7_GO_3_merge_slurm.sh
# 冒烟测试: sbatch --array=0,3 --export=ALL,N_PERM=96 --time=00:30:00 Valid7_GO_2_perm_slurm.sh (测试后删 tmp/Valid7_GO/perm/*_perm_00001-00096.rds)
# 超时或失败: 直接重投, 已完成的 batch 会跳过

set -euo pipefail
module purge
module load GCC/13.3.0 OpenMPI/5.0.3 ArchR/1.0.3-R-4.4.2 SeuratDisk/20231104-R-4.4.2 R_tamu/4.4.2

SETS=(Bn_TCL1A CD4_Tn_CCR7 CD8_Tn_CCR7 Combined)
S=${SETS[$SLURM_ARRAY_TASK_ID]}
DIR=/scratch/user/tfcui2025/CIMA/Claude_project/Scripts/Validation
echo "job = ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} | set = $S | cores = $SLURM_CPUS_PER_TASK | N_PERM = ${N_PERM:-5000} | node = $(hostname)"
Rscript $DIR/Valid7_GO_2_perm.R $S
