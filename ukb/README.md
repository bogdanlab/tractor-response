# Real data analysis
We are not able to share individual-level data so we describe procedures to select 
individuals and pre-process the dataset.


## Admixed individual selection
We run ![SCOPE](https://www.biorxiv.org/content/10.1101/2021.05.11.443705v1.full), 
a scalable program to infer admixture proportion in biobank-scale data. We used four 
population EUR, AFR, EAS, SAS, as ancestral population. We select individuals based on
(EUR > 0.05) & (AFR > 0.05) & (EAS < 0.05) & (SAS < 0.05). This yields XX individuals 
as admixed population in this study. See 00_data.ipynb notebook for details.

## Genotype processing
With the subset of admixed individuals, we filter for hwe < 1e-6, MAF > 0.01, 
genotype missing rate < 0.05 to select SNPs, using the following comamand.

```bash
plink --bfile ${bfile} \
    --keep-fam ${admix_id} \
    --keep-allele-order \
    --make-bed \
    --hwe 1e-6 \
    --maf 0.01 \
    --geno 0.05 \
    --out ${out}
```

Then we perform phasing and imputation using 
![TOPMed Imputation Server](https://imputation.biodatacatalyst.nhlbi.nih.gov/#!).

We perform post-imputation QC  
```bash
bcftools filter -i 'INFO/R2>0.8 & INFO/MAF > 0.005' ${vcf_input} -Oz -o ${vcf_output}
tabix -p vcf ${vcf_output}
```

## Local ancestry inference
We follow Atkinson et al.: we use AFR and EUR in 1000G reference panel and RFmix to infer
local ancestry. See `s03_rfmix.sh` for details. 

Alternatively, we also infer local ancestry with ![LAMP-LD](https://github.com/bogdanlab/lamp-ld)
with window size 200 SNPs, and 10 protypical states. Our inference results are robust to the
parameters. Results from LAMP-LD and RFmix are consistent (data not shown).

## Association testing
See `s04_assoc.sh`