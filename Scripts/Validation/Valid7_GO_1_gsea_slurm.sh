#!/bin/bash
#SBATCH --job-name=Valid7_GO_1_gsea
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/user/tfcui2025/CIMA/Claude_project/Scripts/Validation/logs/%x_%j.out

# Valid7_GO 第 1 步: MR-Only gene 的观测 GSEA (4 组, exponent = 1); 之后: Valid7_GO_2_perm (随机) -> Valid7_GO_3_merge
# 提交前: mkdir -p /scratch/user/tfcui2025/CIMA/Claude_project/Scripts/Validation/logs

set -euo pipefail
module purge
module load GCC/13.3.0 OpenMPI/5.0.3 ArchR/1.0.3-R-4.4.2 SeuratDisk/20231104-R-4.4.2 R_tamu/4.4.2

DIR=/scratch/user/tfcui2025/CIMA/Claude_project/Scripts/Validation
echo "job = $SLURM_JOB_ID | cores = $SLURM_CPUS_PER_TASK | node = $(hostname)"
Rscript $DIR/Valid7_GO_1_gsea.R
