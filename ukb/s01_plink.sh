#!/bin/bash -l
#$ -cwd
#$ -l h_data=8G,h_rt=0:60:00
#$ -j y
#$ -o ./job_out
#$ -t 1-22
. /u/local/Modules/default/init/modules.sh

module load plink
module load htslib

bfile_dir=/u/project/sriram/ukbiobank/data/geno/cal/
out_dir=01_plink/

mkdir -p ${out_dir}

i_chr=${SGE_TASK_ID}
plink --bfile ${bfile_dir}/${i_chr} \
    --keep-fam 00_data/admix_indiv.txt \
    --keep-allele-order \
    --make-bed \
    --hwe 1e-6 \
    --maf 0.01 \
    --geno 0.05 \
    --out ${out_dir}/chr${i_chr}

plink --bfile ${out_dir}/chr${i_chr} \
    --recode vcf \
    --out ${out_dir}/chr${i_chr}

bgzip -c ${out_dir}/chr${i_chr}.vcf > ${out_dir}/chr${i_chr}.vcf.gz
rm ${out_dir}/chr${i_chr}.vcf
