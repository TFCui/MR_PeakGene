#!/bin/bash
#SBATCH --job-name=Valid6_GWAS_2_annot
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/user/tfcui2025/CIMA/Claude_project/Results/Validation/Valid6_GWAS/logs/%x_%j.out

# Valid6_GWAS_2_annot: Fig 4B 候选 example 的注释表 (见脚本头部); 1 个 job, 3 个 CT 一起
# 提交前: mkdir -p /scratch/user/tfcui2025/CIMA/Claude_project/Results/Validation/Valid6_GWAS/logs

set -euo pipefail
module purge
module load GCC/13.3.0 OpenMPI/5.0.3 ArchR/1.0.3-R-4.4.2 SeuratDisk/20231104-R-4.4.2 R_tamu/4.4.2

# R 的 tempdir 放在项目内
export TMPDIR=/scratch/user/tfcui2025/CIMA/Claude_project/tmp/Valid6_GWAS/Rtmp_$SLURM_JOB_ID
mkdir -p $TMPDIR

DIR=/scratch/user/tfcui2025/CIMA/Claude_project/Scripts/Validation
echo "job = $SLURM_JOB_ID | cores = $SLURM_CPUS_PER_TASK | node = $(hostname)"
Rscript $DIR/Valid6_GWAS_2_annot.R
rm -rf $TMPDIR
