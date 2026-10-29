#!/bin/bash
#SBATCH --job-name=Table4_DataOverview
#SBATCH --output=/scratch/user/tfcui2025/CIMA/Claude_project/Scripts/Visualization/logs/%x_%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
# Table 4 (论文 S4 Table): 三个细胞类型的数据基本信息; 输出 Results/Figures/Figures_20261027/Table_4.xlsx + Table_4.tsv
# 提交: cd /scratch/user/tfcui2025/CIMA/Claude_project && sbatch Scripts/Visualization/Table4_DataOverview_slurm.sh

PROJ=/scratch/user/tfcui2025/CIMA/Claude_project
module purge
module load GCC/13.3.0 OpenMPI/5.0.3 ArchR/1.0.3-R-4.4.2 SeuratDisk/20231104-R-4.4.2 R_tamu/4.4.2
export TMPDIR=$PROJ/tmp/Table4_Rtmp_$SLURM_JOB_ID
mkdir -p $TMPDIR
cd $PROJ
echo "===== Table4 $(date)"
Rscript Scripts/Visualization/Table4_DataOverview.R && echo "===== OK" || echo "!!!!! FAILED"
rm -rf $TMPDIR
echo "===== DONE $(date)"
