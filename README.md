# Response to "Tractor uses local ancestry to enable the inclusion of admixed individuals in GWAS and to boost power"

- See [power.ipynb](power.ipynb) for the power simulation
- See [Snakefile](Snakefile) for replicating the pipeline

## Required softwares
- Download MIXSCORE from [here](https://cdn1.sph.harvard.edu/wp-content/uploads/sites/181/2013/02/mixscore-1.3.tar.gz) to `software/` folder, so the `mixscore_wrapper` function in [admix.py](admix.py) can find the mixscore executable program.
- Use our genotype simulation [pipeline](https://github.com/bogdanlab/admixed_genotype_simulation) to simulate admixed individuals with real genotypes from 1000 Genomes project.
