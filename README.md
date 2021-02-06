# Response to "Tractor uses local ancestry to enable the inclusion of admixed individuals in GWAS and to boost power"

- See [power.ipynb](power.ipynb) for full results.

## Steps to replicate our results
1. Prepare simulated admixed genotype. Download the simulated genotype with 1000G [here](https://github.com/bogdanlab/tractor-response/releases/download/v1.0/raw_geno.zip) and unzip it to `raw_geno` folder.
2. Run the snakemake pipeline with `Snakefile`. One may need to run the experiments in parallel to speed up.


## Required softwares
- Download MIXSCORE from [here](https://cdn1.sph.harvard.edu/wp-content/uploads/sites/181/2013/02/mixscore-1.3.tar.gz) to `software/` folder, so the `mixscore_wrapper` function in [admix.py](admix.py) can find the mixscore executable program.
- As an alternative to run the experiments, use our genotype simulation [pipeline](https://github.com/bogdanlab/admixed_genotype_simulation) to simulate admixed individuals with real genotypes from 1000 Genomes project.
