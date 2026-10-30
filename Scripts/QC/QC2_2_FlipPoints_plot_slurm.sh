#!/bin/bash
#SBATCH --job-name=QC2_FlipPoints
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/user/tfcui2025/CIMA/Claude_project/Scripts/QC/logs/%x_%j.out

# QC2_2: MR flip 前后对比图 (6 个 pair, 单页 3 行 x 4 panel), 输出 Results/QC/2_FlipPoints/FigureS2_estimates.tsv + Figures_20261102/Figure_S2.pdf (只出 PDF) (流程见 Scripts/Pipeline.md 第 4.2 节)
# 提交: cd /scratch/user/tfcui2025/CIMA/Claude_project && mkdir -p Scripts/QC/logs && sbatch Scripts/QC/QC2_2_FlipPoints_plot_slurm.sh

PROJ=/scratch/user/tfcui2025/CIMA/Claude_project

module purge
module load GCC/13.3.0 OpenMPI/5.0.3 ArchR/1.0.3-R-4.4.2 SeuratDisk/20231104-R-4.4.2 R_tamu/4.4.2

# R 的 tempdir 也放在项目内
export TMPDIR=$PROJ/tmp/QC2_FlipPoints/Rtmp_${SLURM_JOB_ID}
mkdir -p $TMPDIR
cd $PROJ

echo "===== QC2_FlipPoints $(date)"
Rscript Scripts/QC/QC2_2_FlipPoints_plot.R && echo "===== OK" || echo "!!!!! FAILED"

rm -rf $TMPDIR
echo "===== DONE $(date)"
