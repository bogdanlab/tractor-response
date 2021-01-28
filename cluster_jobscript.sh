#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -sync y

export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

{exec_job}
