#!/bin/bash
#SBATCH --job-name=Valid2_SLDSC_3_merge
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/scratch/user/tfcui2025/CIMA/Claude_project/Scripts/Validation/logs/%x_%j.out

# Valid2_SLDSC_3_merge (流程见 Scripts/Pipeline.md 第 3 节; 提交顺序见 3.0.3)
# 提交前: mkdir -p /scratch/user/tfcui2025/CIMA/Claude_project/Scripts/Validation/logs

set -euo pipefail
module purge
module load GCC/13.3.0 OpenMPI/5.0.3 ArchR/1.0.3-R-4.4.2 SeuratDisk/20231104-R-4.4.2 R_tamu/4.4.2

DIR=/scratch/user/tfcui2025/CIMA/Claude_project/Scripts/Validation
echo "job = $SLURM_JOB_ID | cores = $SLURM_CPUS_PER_TASK | node = $(hostname)"
Rscript $DIR/Valid2_SLDSC_3_merge.R
