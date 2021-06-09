#!/bin/bash -l
#$ -cwd
#$ -l h_data=16G,h_rt=8:00:00
#$ -j y
#$ -o ./job_out
#$ -t 1-22

. /u/local/Modules/default/init/modules.sh
export PATH=~/project-pasaniuc/software/miniconda3/bin:$PATH
export PYTHONNOUSERSITE=True

# for trait in cholesterol ldl_direct; do
# 
#
i_chr=${SGE_TASK_ID}
trait=$1
mkdir -p 04_assoc/assoc/

python s04_assoc.py \
    --ds 04_assoc/dataset/chr${i_chr}.zarr \
    --trait ${trait} \
    --out 04_assoc/assoc/${trait}.chr${i_chr}.csv 
