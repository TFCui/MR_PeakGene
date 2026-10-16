#!/bin/bash
#SBATCH --job-name=Plot2_examples
#SBATCH --output=/scratch/user/tfcui2025/CIMA/Claude_project/Scripts/Visualization/logs/%x_%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=03:00:00
# Fig2 example pairs: 一个作业画全部选定的 pair (Example_pairs.tsv 中 pair_set = Fig2_selected, 跨三个细胞类型), 输出 Figures_20261019/Fig2_examples.pdf
# 提交: cd /scratch/user/tfcui2025/CIMA/Claude_project && mkdir -p Scripts/Visualization/logs && sbatch Scripts/Visualization/Plot2_examples_slurm.sh

PROJ=/scratch/user/tfcui2025/CIMA/Claude_project

module purge
module load GCC/13.3.0 OpenMPI/5.0.3 ArchR/1.0.3-R-4.4.2 SeuratDisk/20231104-R-4.4.2 R_tamu/4.4.2

# R 的 tempdir 也放在项目内
export TMPDIR=$PROJ/tmp/plink_fig2/Rtmp_${SLURM_JOB_ID}
mkdir -p $TMPDIR
cd $PROJ

echo "===== Plot2_examples_v2 $(date)"
Rscript Scripts/Visualization/Plot2_examples_v2.R && echo "===== OK" || echo "!!!!! FAILED"

rm -rf $TMPDIR
echo "===== DONE $(date)"
