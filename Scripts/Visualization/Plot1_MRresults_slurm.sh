#!/bin/bash
#SBATCH --job-name=Plot1_MRresults
#SBATCH --output=/scratch/user/tfcui2025/CIMA/Claude_project/Scripts/Visualization/logs/%x_%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=01:00:00
# Fig1B/C + S3: 重写 Figure_1B_C (PLOT1_WRITE_1BC=1) 和 Figure_S3 到 Figures_20261027/
# 提交: cd /scratch/user/tfcui2025/CIMA/Claude_project && mkdir -p Scripts/Visualization/logs && sbatch Scripts/Visualization/Plot1_MRresults_slurm.sh

PROJ=/scratch/user/tfcui2025/CIMA/Claude_project

module purge
module load GCC/13.3.0 OpenMPI/5.0.3 ArchR/1.0.3-R-4.4.2 SeuratDisk/20231104-R-4.4.2 R_tamu/4.4.2

export TMPDIR=$PROJ/tmp/Rtmp_plot1_${SLURM_JOB_ID}
mkdir -p $TMPDIR
cd $PROJ

echo "===== Plot1_MRresults $(date)"
PLOT1_WRITE_1BC=1 Rscript Scripts/Visualization/Plot1_MRresults.R && echo "===== OK" || echo "!!!!! FAILED"

rm -rf $TMPDIR
echo "===== DONE $(date)"
