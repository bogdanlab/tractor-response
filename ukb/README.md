# Real data analysis
We are not able to share individual-level data so we describe procedures to select 
individuals and pre-process the dataset.


## Admixed individual selection
We run [SCOPE](https://www.biorxiv.org/content/10.1101/2021.05.11.443705v1.full), 
a scalable program to infer admixture proportion in biobank-scale data. We used four 
population EUR, AFR, EAS, SAS, as ancestral population. We select individuals based on
(EUR > 0.05) & (AFR > 0.05) & (EAS < 0.05) & (SAS < 0.05). This yields 4327 individuals 
as admixed population in this study. See `00_data.ipynb` notebook for details.

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
[TOPMed Imputation Server](https://imputation.biodatacatalyst.nhlbi.nih.gov/#!).
We perform post-imputation QC to filter for imputation R2 > 0.8 and MAF > 0.5%.   
```bash
bcftools filter -i 'INFO/R2>0.8 & INFO/MAF > 0.005' ${vcf_input} -Oz -o ${vcf_output}
tabix -p vcf ${vcf_output}
```

## Local ancestry inference
We follow Tractor paper: we use AFR and EUR in 1000G reference panel and RFmix to infer
local ancestry. See `s03_rfmix.sh` for details. The inferred local ancestry will be used as input to SNP1 and Tractor.


## Association testing on known risk regions to lipid traits
We compared ATT, SNP1, Tractor on four well-known risk regions (LDLR, APOE, PCSK9, SORT1) of lipid traits (HDL, TC). For all three methods, we include age, sex, dilution factor, and top 10 PCs as covariates.

No inflation of p-values are observed for the three methods. We take +- 50kb window around the transcribed regions of each gene, and compare the association strength for each of the method.

See `s04_assoc.ipynb` for more details. See `results/risk_regions.xlsx` for numerical results.