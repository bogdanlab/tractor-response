#!/bin/bash -l
#$ -cwd
#$ -l h_data=12G,h_rt=16:00:00
#$ -j y
#$ -o ./job_out

export PATH=~/project-pasaniuc/software/anaconda3/bin:$PATH
export PYTHONNOUSERSITE=True

. /u/local/Modules/default/init/modules.sh

snakemake \
    --snakefile Snakefile \
    --jobscript cluster_jobscript.sh \
    --cluster-sync "qsub -l h_data={resources.mem_gb}G,h_rt=00:{resources.time_min}:00 -o job_out" \
    --jobs 100 \
    --printshellcmds \
    --max-jobs-per-second 0.3 \
    --restart-times 5 \
    --latency-wait 20 \
    --default-resources mem_gb=16 time_min=10
