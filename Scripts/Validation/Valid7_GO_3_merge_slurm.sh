#!/bin/bash
#SBATCH --job-name=Valid7_GO_3_merge
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=00:30:00
#SBATCH --output=/scratch/user/tfcui2025/CIMA/Claude_project/Scripts/Validation/logs/%x_%j.out

# Valid7_GO 第 3 步: 合并观测与随机结果, 算经验 P, 覆盖 Results/Validation/Valid7_GO/GO_GSEA_MROnly.tsv.gz
# 提交: sbatch --dependency=afterok:<Valid7_GO_2_perm 的 jobid> Valid7_GO_3_merge_slurm.sh

set -euo pipefail
module purge
module load GCC/13.3.0 OpenMPI/5.0.3 ArchR/1.0.3-R-4.4.2 SeuratDisk/20231104-R-4.4.2 R_tamu/4.4.2

DIR=/scratch/user/tfcui2025/CIMA/Claude_project/Scripts/Validation
echo "job = $SLURM_JOB_ID | N_PERM = ${N_PERM:-5000} | node = $(hostname)"
Rscript $DIR/Valid7_GO_3_merge.R
