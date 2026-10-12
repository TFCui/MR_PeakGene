#!/bin/bash
#SBATCH --job-name=Plot8_Link
#SBATCH --output=/scratch/user/tfcui2025/CIMA/Claude_project/Scripts/Visualization/logs/%x_%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=02:00:00
# Fig8: Valid6_GWAS pairs 的位点 link 图 (GWAS 模式, 一页 = 一个 link x 一个 trait), 3 个 CT 画在同一个 PDF
# 依赖: Results/Validation/Valid6_GWAS/GWAS_MRonly_overlap.tsv
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
