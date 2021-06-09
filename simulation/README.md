# Simulation
## Results
See [results.ipynb](results.ipynb), including
1. main figure.
2. qq plot and false positive ratio in the null simulation.
3. extended results of power.
4. association testing when the causal SNP is not typed

## Steps to replicate our results
1. Clone this repository and install dependency with `pip install requirements.txt`. Also download MIXSCORE from [here](https://cdn1.sph.harvard.edu/wp-content/uploads/sites/181/2013/02/mixscore-1.3.tar.gz) to `software/` folder, so the `mixscore_wrapper` function in [admix.py](admix.py) can find the mixscore executable program.
1. Prepare simulated admixed genotype. Download the simulated genotype with 1000G [here](https://github.com/bogdanlab/tractor-response/releases/download/v1.0/raw_geno.zip) and unzip it to `raw_geno` folder.
2. Run the snakemake pipeline with `Snakefile`. One may need to run the experiments in parallel to speed up.

## Also see
- As an alternative to run the experiments, use our genotype simulation [pipeline](https://github.com/bogdanlab/admixed_genotype_simulation) to simulate admixed individuals with real genotypes from 1000 Genomes project.
