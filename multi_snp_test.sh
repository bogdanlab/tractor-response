#!/bin/bash -l
#$ -cwd
#$ -l h_data=4G,h_rt=0:30:00,highp
#$ -j y
#$ -o ./job_out

. /u/local/Modules/default/init/modules.sh
export PATH=~/project-pasaniuc/software/anaconda3/bin:$PATH
export PYTHONNOUSERSITE=True
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

region_i=${SGE_TASK_ID}

odds_ratio=1.2
anc_effect=0.0
case_prevalence=0.1
control_ratio=1.0

out_dir=out/multi_snp_test/cp_${case_prevalence}_cr_${control_ratio}/or_${odds_ratio}_ae_${anc_effect}

mkdir -p ${out_dir}

python experiment.py multi_snp_test_cli \
        --data_dir data/geno/finemap_300 \
        --odds_ratio ${odds_ratio} \
        --anc_effect ${anc_effect} \
        --region_size 41 \
	--region_i ${region_i} \
        --case_prevalence ${case_prevalence} \
        --control_ratio ${control_ratio} \
        --n_sim 50 \
        --out ${out_dir}/region_${region_i}.csv.gz

