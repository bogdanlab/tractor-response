# Response to "Tractor uses local ancestry to enable the inclusion of admixed individuals in GWAS and to boost power"

See our manuscript [On powerful GWAS in admixed populations, Nat. Genet. 2021](https://www.nature.com/articles/s41588-021-00953-5).

> This repository is mainly for reproducing results of the manuscript. Go to [admix-kit software package](https://github.com/KangchengHou/admix-kit) if you want to perform similar analysis. 

See [simulation study](
https://nbviewer.jupyter.org/github/bogdanlab/tractor-response/blob/main/simulation/results.ipynb) for
1. main figure.
2. qq plot and false positive ratio in the null simulation.
3. extended results of power.
4. association testing when the causal SNP is not typed

See [real data study](https://nbviewer.jupyter.org/github/bogdanlab/tractor-response/blob/main/ukb/04_assoc.ipynb) 
for analyses of known risk reigons (LDLR, APOE, PCSK9, SORT1) for lipids traits (LDL, TC) in 4,327 individuals
with European and African ancestries.


See directory [simulation](simulation), [ukb](ukb) for details of simulation and real data analysis.
