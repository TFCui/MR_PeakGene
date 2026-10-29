#!/bin/bash
#SBATCH --job-name=Plot10_Combined
#SBATCH --output=/scratch/user/tfcui2025/CIMA/Claude_project/Scripts/Visualization/logs/%x_%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=02:00:00
# Figure 4 (原 Fig10): S-LDSC Combined + MR-only link + GWAS pair; Figure_4 = 4 个 pair (1 行), Figure_S4 = 8 个 pair (2 行, 单页)
# 输出: Results/Figures/Figures_20261027/Figure_4.{pdf,svg} + Figure_S4.{pdf,svg}
# 提交: cd /scratch/user/tfcui2025/CIMA/Claude_project && mkdir -p Scripts/Visualization/logs && sbatch Scripts/Visualization/Plot10_Combined_slurm.sh

PROJ=/scratch/user/tfcui2025/CIMA/Claude_project
module purge
module load GCC/13.3.0 OpenMPI/5.0.3 ArchR/1.0.3-R-4.4.2 SeuratDisk/20231104-R-4.4.2 R_tamu/4.4.2

# R 的 tempdir 也放在项目内
export TMPDIR=$PROJ/tmp/plink_fig10/Rtmp_$SLURM_JOB_ID
mkdir -p $TMPDIR
cd $PROJ

echo "===== Plot10 $(date)"
Rscript Scripts/Visualization/Plot10_Combined.R && echo "===== Plot10 OK" || echo "!!!!! FAILED: Plot10"
rm -rf $TMPDIR
echo "===== DONE $(date)"
