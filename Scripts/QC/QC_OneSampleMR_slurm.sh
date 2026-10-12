#!/bin/bash
#SBATCH --job-name=QC_OneSampleMR
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --array=0-2
#SBATCH --output=/scratch/user/tfcui2025/CIMA/Claude_project/Scripts/QC/logs/%x_%A_%a.out

# QC_OneSampleMR 第 2 步: null SNP 上 caQTL 与 eQTL z-score 的相关性 (流程见 Scripts/Pipeline.md 第 4.1 节)
# 提交前: mkdir -p /scratch/user/tfcui2025/CIMA/Claude_project/Scripts/QC/logs

set -euo pipefail
module purge
module load GCC/13.3.0 OpenMPI/5.0.3 ArchR/1.0.3-R-4.4.2 SeuratDisk/20231104-R-4.4.2 R_tamu/4.4.2

DIR=/scratch/user/tfcui2025/CIMA/Claude_project/Scripts/QC
CELLS=(Bn_TCL1A CD4_Tn_CCR7 CD8_Tn_CCR7)
CT=${CELLS[$SLURM_ARRAY_TASK_ID]}
echo "job = ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} | cell_type = $CT | cores = $SLURM_CPUS_PER_TASK | node = $(hostname)"
Rscript $DIR/QC_OneSampleMR.R $CT
