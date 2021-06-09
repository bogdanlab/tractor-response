#!/bin/bash -l
#$ -cwd
#$ -l h_data=25G,h_rt=2:00:00
#$ -j y
#$ -o ./job_out

. /u/local/Modules/default/init/modules.sh && module load plink

cd 01_plink
rm -f merge_list.txt

for i_chr in $(seq 2 22); do
  echo chr${i_chr} >> merge_list.txt
done

plink --bfile chr1  --merge-list merge_list.txt --make-bed --out merged
plink --bfile merged --pca --out pca
