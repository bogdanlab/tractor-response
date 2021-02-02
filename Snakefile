from os.path import join
import os
import yaml
import experiment
import admix
import itertools
import numpy as np
import pandas as pd
config = {
    "N_TEST_CHUNK": 100
}

RAW_DATA_DIR = "/u/project/pasaniuc/pasaniucdata/admixture/kangcheng/genotype_simulation/out/kg_3k"

GENO_PREFIX_LIST = [
    "EUR_0.2_AFR_0.8_7_80000"
]

PARAM_LIST = [param for param in itertools.product([1.0, 1.15, 1.2, 1.25], [0.0])]

rule all:
    input:
        expand("out/single_snp_test_tractor/{geno_prefix}/{case_control_param}/{sim_param}/chunk_{chunk_i}.csv",
            geno_prefix=GENO_PREFIX_LIST,
            case_control_param=["0.1_1.0", "0.1_2.5"],
            sim_param=['_'.join([str(i) for i in p]) for p in PARAM_LIST],
            chunk_i=np.arange(config["N_TEST_CHUNK"])),
        expand("out/single_snp_test_tractor/{geno_prefix}/{case_control_param}/{sim_param}/summary.csv",
            geno_prefix=GENO_PREFIX_LIST,
            case_control_param=["0.1_1.0", "0.1_2.5"],
            sim_param=['_'.join([str(i) for i in p]) for p in PARAM_LIST])


rule format_data:
    resources:
        mem_gb=8,
        time_min=20
    input:
        legend = join(RAW_DATA_DIR, "legend.txt"),
        anc = join(RAW_DATA_DIR, "{sim_prefix}/admix.hanc"),
        phgeno = join(RAW_DATA_DIR, "{sim_prefix}/admix.phgeno"),
    output:
        anc = "data/geno/{sim_prefix}/anc.npy",
        phgeno = "data/geno/{sim_prefix}/phgeno.npy",
        legend = "data/geno/{sim_prefix}/legend.csv"
    run:
        import numpy as np
        import pandas as pd

        legend = pd.read_csv(input.legend, delim_whitespace=True)
        ancestry = admix.read_int_mat(input.anc)
        phgeno = admix.read_int_mat(input.phgeno).T

        np.save(output.anc, ancestry)
        np.save(output.phgeno, phgeno)
        legend.to_csv(output.legend, index=False)


rule single_snp_test_tractor:
    resources:
        mem_gb=5,
        time_min=180
    input:
        anc = "data/geno/{geno_prefix}/anc.npy",
        phgeno = "data/geno/{geno_prefix}/phgeno.npy",
        legend = "data/geno/{geno_prefix}/legend.csv"
    output:
        "out/single_snp_test_tractor/{geno_prefix}/{case_prevalence}_{control_ratio}/{odds_ratio}_{anc_effect}/chunk_{chunk_i}.csv"
    run:
        import numpy as np
        # read phgeno, anc
        anc = np.load(input.anc)
        phgeno = np.load(input.phgeno)
        # calculating global ancestry
        global_ancestry = anc.reshape((anc.shape[0] // 2, anc.shape[1] * 2)).mean(axis=1)
        print(global_ancestry)
        print(np.mean(global_ancestry), np.std(global_ancestry))
        assert np.all(anc.shape == phgeno.shape)
        
        n_snp = anc.shape[1]
        chunk_index = np.array_split(np.arange(n_snp), config["N_TEST_CHUNK"])[int(wildcards.chunk_i)]
        print("chunk_index", chunk_index)
        phgeno = phgeno[:, chunk_index]
        anc = anc[:, chunk_index]
        print(f"case_prevalence: {wildcards.case_prevalence}, control_ratio: {wildcards.control_ratio}")
        print(f"odds_ratio: {wildcards.odds_ratio}, anc_effect: {wildcards.anc_effect}")
        rls = experiment.single_snp_test_tractor(phgeno=phgeno, 
                                                 anc=anc,
                                                 theta=global_ancestry, 
                                                 odds_ratio=float(wildcards.odds_ratio),
                                                 anc_effect=float(wildcards.anc_effect),
                                                 case_prevalence=float(wildcards.case_prevalence),
                                                 control_ratio=float(wildcards.control_ratio),
                                                 seed=int(wildcards.chunk_i), n_sim=50)

        # save
        rls.to_csv(output[0], index=False)

rule single_snp_test_summary:
    resources:
        mem_gb=8,
        time_min=10
    input:
        expand("out/single_snp_test_tractor/{{prefix}}/chunk_{chunk_i}.csv",
                chunk_i=np.arange(config["N_TEST_CHUNK"])),
    output:
        "out/single_snp_test_tractor/{prefix}/summary.csv"
    run:
        import numpy as np
        
        rls_list = []
        total_n_snp = 0
        for chunk_i in range(config["N_TEST_CHUNK"]):
            chunk_rls = pd.read_csv(input[chunk_i])
            chunk_n_snp = len(chunk_rls["SNP_I"].unique())
            chunk_rls["SNP_I"] += total_n_snp
            total_n_snp += chunk_n_snp
            rls_list.append(chunk_rls)
        print(f"total_n_snp: {total_n_snp}")
        pd.concat(rls_list).to_csv(output[0], index=False)



# rule finemap:
#     resources:
#         mem_gb=8,
#         time_min=30
#     input:
#         anc = "data/geno/{geno_prefix}/anc.npy",
#         phgeno = "data/geno/{geno_prefix}/phgeno.npy",
#         legend = "data/geno/{geno_prefix}/legend.csv"
#     output:
#         score= "out/finemap/{geno_prefix}/{var_g}_{cov}/chunk_{chunk_i}.csv",
#         beta = "out/finemap/{geno_prefix}/{var_g}_{cov}/chunk_{chunk_i}.beta.npy",
#     run:
#         # read phgeno, anc
#         anc = np.load(input.anc)
#         phgeno = np.load(input.phgeno)

#         # calculating global ancestry
#         global_ancestry = anc.reshape((anc.shape[0] // 2, anc.shape[1] * 2)).mean(axis=1)
#         assert np.all(anc.shape == phgeno.shape)

#         # seperate chunks
#         n_snp = anc.shape[1]
#         chunk_index = np.array_split(np.arange(n_snp), config["N_FINEMAP_CHUNK"])[int(wildcards.chunk_i)]
#         print(chunk_index)
#         phgeno = phgeno[:, chunk_index]
#         anc = anc[:, chunk_index]
#         print(phgeno.shape, anc.shape)
#         print(f"var_g: {wildcards.var_g}, cov: {wildcards.cov}")
#         score_df, beta = experiment.finemap(phgeno=phgeno, 
#                                             anc=anc, 
#                                             theta=global_ancestry, 
#                                             var_g=float(wildcards.var_g), cov=float(wildcards.cov), 
#                                             seed=int(wildcards.chunk_i), n_sim=10)

#         # save
#         score_df.to_csv(output.score, index=False)
#         np.save(output.beta, beta)



