#!/bin/bash -l
#$ -cwd
#$ -l h_data=16G,h_rt=3:30:00,highp
#$ -j y
#$ -o ./job_out
#$ -t 1-22

# for trait in cholesterol ldl_direct; do 
# qsub s04_plink_assoc.sh ${trait}
# done

. /u/local/Modules/default/init/modules.sh
module load plink

i_chr=${SGE_TASK_ID}
trait=$1

pheno_file=01_plink/${trait}.pheno
covar_file=01_plink/covar.txt

num_cols=$(head -1 $covar_file | awk '{print NF}')
num_cols=$((num_cols-2))

mkdir -p 04_assoc/plink_assoc_imputed

  plink \
    --allow-no-sex \
    --vcf 02_impute/processed/chr${i_chr}/chr${i_chr}.sample.imputed.vcf.gz \
    --ci 0.95 \
    --covar ${covar_file} \
    --covar-number 1-${num_cols} \
    --linear hide-covar \
    --pheno ${pheno_file} \
    --out 04_assoc/plink_assoc_imputed/${trait}.${i_chr}
