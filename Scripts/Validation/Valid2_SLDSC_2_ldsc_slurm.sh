#!/bin/bash
#SBATCH --job-name=Valid2_SLDSC_ldsc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=22
#SBATCH --mem=150G
#SBATCH --time=12:00:00
#SBATCH --array=0-20%20
#SBATCH --output=/scratch/user/tfcui2025/CIMA/Claude_project/Scripts/Validation/logs/%x_%A_%a.out

# Valid2_SLDSC 步骤 2 (流程见 Scripts/Pipeline.md 3.2); 移植自 Scripts_reference/D2_sldsc_slurm.sh
# 一个 array task = 一个 (CT, 集合); task id = 7 * CT_index + set_index
#   (1) ldsc --l2: 22 条染色体并行 (22 个 ldscore 都在就跳过)
#   (2) ldsc --h2: 9 个 GWAS 并行, 模型 = 该集合 (第一位) + baselineLD (.results 已存在就跳过)
# 依赖步骤 1: sbatch --dependency=afterok:<prep jobid> Valid2_SLDSC_2_ldsc_slurm.sh

module purge
module load Anaconda3/2025.12-2
eval "$(conda shell.bash hook)"
conda activate /scratch/user/tfcui2025/CIMA/R_package/S_LDSC/ldsc/ldsc_env
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

LDSC_DIR=/scratch/user/tfcui2025/CIMA/R_package/S_LDSC/ldsc
REF_DIR=/scratch/user/tfcui2025/CIMA/R_package/S_LDSC/reference_genome
PLINK=$REF_DIR/1000G_EUR_Phase3_plink/1000G.EUR.QC.
HM3=$REF_DIR/listHM3.txt
BASELINE=$REF_DIR/baselineLD_v2.2/baselineLD.
WEIGHTS=$REF_DIR/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.
FRQ=$REF_DIR/1000G_Phase3_frq/1000G.EUR.QC.
GWAS_DIR=/scratch/user/tfcui2025/CIMA/R_package/S_LDSC/GWAS_data/selected_sumstats

PROJ=/scratch/user/tfcui2025/CIMA/Claude_project
CELLS=(Bn_TCL1A CD4_Tn_CCR7 CD8_Tn_CCR7)
SETS=(ivw egger smr cct pseudo wnn pool)

# 9 个 GWAS: 标签和 sumstats 文件前缀一一对应
TRAITS=(AS CD IBD RA SLE MS T1D Asthma PBC)
TFILES=(CUSTOM.AS
        CUSTOM.CD
        PASS.Inflammatory_Bowel_Disease.deLange2017
        PASS.Rheumatoid_Arthritis.Saevarsdottir2022
        PASS.Lupus.Bentham2015
        PASS_Multiple_sclerosis
        PASS.Type_1_Diabetes.Chiou2021
        GBMI.Asthma
        PASS.Primary_Biliary_Cirrhosis.Cordell2015)

ID=$SLURM_ARRAY_TASK_ID
CT=${CELLS[$((ID / 7))]}
S=${SETS[$((ID % 7))]}
NAME=${CT}_${S}
PRE=$PROJ/tmp/Valid2_SLDSC/annot_ld/$CT/$NAME
OUT3=$PROJ/Results/Validation/Valid2_SLDSC/h2/$CT
mkdir -p $OUT3
echo "[$(date)] task $ID: $NAME"

# ---- (1) LD score ----
N_LD=$(ls ${PRE}.*.l2.ldscore.gz 2>/dev/null | wc -l)
if [ $N_LD -eq 22 ]; then
  echo "   skip l2"
else
  for CHR in {1..22}; do
    python $LDSC_DIR/ldsc.py --l2 \
      --bfile ${PLINK}${CHR} \
      --ld-wind-cm 1 \
      --annot ${PRE}.${CHR}.annot.gz \
      --thin-annot \
      --print-snps $HM3 \
      --out ${PRE}.${CHR} \
      > /dev/null 2>&1 &
  done
  wait
  N_LD=$(ls ${PRE}.*.l2.ldscore.gz 2>/dev/null | wc -l)
  echo "[$(date +%T)] ldscore: $N_LD / 22"
  if [ $N_LD -ne 22 ]; then echo "LD score 不完整, 看 ${PRE}.*.log"; exit 1; fi
  rm -f ${PRE}.*.log                                      # 成功才删 l2 的 log
fi

# ---- (2) S-LDSC h2: 9 个 GWAS 并行 ----
run_h2 () {                                               # $1 = trait 标签, $2 = sumstats 文件前缀
  local F=$OUT3/${NAME}_$1
  if [ -s $F.results ]; then echo "   skip h2 $1"; return; fi
  python $LDSC_DIR/ldsc.py --h2 $GWAS_DIR/$2.sumstats.gz \
    --ref-ld-chr ${PRE}.,${BASELINE} \
    --w-ld-chr $WEIGHTS \
    --overlap-annot \
    --frqfile-chr $FRQ \
    --print-coefficients \
    --out $F > /dev/null 2>&1
  [ -s $F.results ] && echo "   h2 $1: OK" || echo "   h2 $1: FAILED (看 $F.log)"
}

for i in "${!TRAITS[@]}"; do
  run_h2 ${TRAITS[$i]} ${TFILES[$i]} &
done
wait

N_RES=$(ls $OUT3/${NAME}_*.results 2>/dev/null | wc -l)
echo "[$(date)] task $ID done: $N_RES / 9 results"
[ $N_RES -eq 9 ]
