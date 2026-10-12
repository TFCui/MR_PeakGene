#!/bin/bash
#SBATCH --job-name=Plot2_examples
#SBATCH --output=/scratch/user/tfcui2025/CIMA/Claude_project/Scripts/Visualization/logs/%x_%A_%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=03:00:00
#SBATCH --array=0-2
# Fig2 example pairs: 每个 array task 一个细胞类型, 依次画 v1 和 v2 (全部 example, --n 0)
# 提交: cd /scratch/user/tfcui2025/CIMA/Claude_project && mkdir -p Scripts/Visualization/logs && sbatch Scripts/Visualization/Plot2_examples_slurm.sh

CELLS=(Bn_TCL1A CD4_Tn_CCR7 CD8_Tn_CCR7)
CT=${CELLS[$SLURM_ARRAY_TASK_ID]}
PROJ=/scratch/user/tfcui2025/CIMA/Claude_project

module purge
module load GCC/13.3.0 OpenMPI/5.0.3 ArchR/1.0.3-R-4.4.2 SeuratDisk/20231104-R-4.4.2 R_tamu/4.4.2

# R 的 tempdir 也放在项目内
export TMPDIR=$PROJ/tmp/plink_fig2/Rtmp_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p $TMPDIR
cd $PROJ

echo "===== $CT v1 $(date)"
Rscript Scripts/Visualization/Plot2_examples_v1.R --celltype $CT --n 0 && echo "===== $CT v1 OK" || echo "!!!!! FAILED: $CT v1"
echo "===== $CT v2 $(date)"
Rscript Scripts/Visualization/Plot2_examples_v2.R --celltype $CT --n 0 && echo "===== $CT v2 OK" || echo "!!!!! FAILED: $CT v2"

rm -rf $TMPDIR
echo "===== DONE $(date)"
