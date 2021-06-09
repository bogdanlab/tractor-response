#!/bin/bash -l
#$ -cwd
#$ -l h_data=12G,h_rt=3:30:00,highp
#$ -j y
#$ -o ./job_out

# for trait in cholesterol ldl_direct; do 
# qsub s04_plink_assoc.sh ${trait}
# done

. /u/local/Modules/default/init/modules.sh
module load plink

trait=$1

pheno_file=01_plink/${trait}.pheno
covar_file=01_plink/covar.txt

num_cols=$(head -1 $covar_file | awk '{print NF}')
num_cols=$((num_cols-2))

mkdir -p 04_assoc/plink_assoc

for i_chr in $(seq 1 22); do
  plink \
    --allow-no-sex \
    --vcf 02_impute/processed/chr${i_chr}/chr${i_chr}.sample.typed.vcf.gz \
    --ci 0.95 \
    --covar ${covar_file} \
    --covar-number 1-${num_cols} \
    --linear hide-covar \
    --pheno ${pheno_file} \
    --out 04_assoc/plink_assoc/${trait}.${i_chr}
done
