#!/bin/bash -l
#$ -cwd
#$ -l h_data=20G,h_rt=5:00:00,highp
#$ -j y
#$ -o ./job_out
#$ -t 1-22


i_chr=${SGE_TASK_ID}

dir_data=./02_impute/processed/chr${i_chr}/

mkdir -p ./03_rfmix/lanc/
./03_rfmix/rfmix/rfmix \
    -f ${dir_data}/chr${i_chr}.sample.typed.nochr.vcf.gz \
    -r ${dir_data}/chr${i_chr}.ref.typed.nochr.vcf.gz \
    --chromosome=${i_chr} \
    -m ./03_rfmix/metadata/sample_map.tsv \
    -g ./03_rfmix/metadata/genetic_map/chr${i_chr}.tsv \
    -e 1 -n 5 \
    -o ./03_rfmix/lanc/chr${i_chr}

