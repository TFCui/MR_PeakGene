#!/bin/bash
#SBATCH --job-name=Plot8_Link
#SBATCH --output=/scratch/user/tfcui2025/CIMA/Claude_project/Scripts/Visualization/logs/%x_%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=02:00:00
# Fig8: Valid6_GWAS MR-only pairs 的 link 图 (一页 = 一个 pair, 只画这一条 link; pairs 见 Results/Validation/Valid6_GWAS/GWAS_MRonly_overlap.tsv)
# 输出: Results/Figures/Figures_20261019/Fig8_Link_GWAS.pdf
# 提交: cd /scratch/user/tfcui2025/CIMA/Claude_project && mkdir -p Scripts/Visualization/logs && sbatch Scripts/Visualization/Plot8_Link_slurm.sh

PROJ=/scratch/user/tfcui2025/CIMA/Claude_project
module purge
module load GCC/13.3.0 OpenMPI/5.0.3 ArchR/1.0.3-R-4.4.2 SeuratDisk/20231104-R-4.4.2 R_tamu/4.4.2

# R 的 tempdir 也放在项目内
export TMPDIR=$PROJ/tmp/plink_fig8/Rtmp_$SLURM_JOB_ID
mkdir -p $TMPDIR
cd $PROJ

echo "===== Plot8 $(date)"
Rscript Scripts/Visualization/Plot8_Link.R && echo "===== Plot8 OK" || echo "!!!!! FAILED: Plot8"
rm -rf $TMPDIR
echo "===== DONE $(date)"
