import pandas as pd
import numpy as np
from os.path import join, dirname
import os
import matplotlib.pyplot as plt
import statsmodels.api as sm
import scipy.stats as stats
import seaborn as sns
import itertools
import pickle
import fire
import admix


def single_snp_test(phgeno, anc, theta, var_g, cov, n_sim=50, seed=1234, case_prevalence=0.2):
    """
    Run one experiment, both simulation and testing
    
    Args
    -----
    phgeno: phased genotype (2 * n_indiv, n_snp), odd and even rows represent the 
        two haplotypes per individual
    anc: local ancestry (2 * n_indiv, n_snp), same as `phgeno`
    theta: global ancestry proportion
    var_g: variance component parameter
    cov: covariance parameter
    
    Returns
    -----
    score_df: DataFrame
    """
    np.random.seed(seed)
    # form individuals
    n_indiv = anc.shape[0] // 2
    n_snp = anc.shape[1]
    anc = anc.reshape((n_indiv, n_snp * 2))
    phgeno = phgeno.reshape((n_indiv, n_snp * 2))
    
    sim_list = admix.simulate_phenotype_cc_1snp(phgeno, anc, var_g=var_g, cov=cov, n_sim=n_sim, case_prevalence=case_prevalence)
    score_df_list = []
    for snp_i, snp_sim in enumerate(sim_list):
        print(snp_i)
        for sim_i in range(n_sim):
            mix_pheno = snp_sim[2][:, sim_i]
            mix_anc = admix.add_up_haplotype(anc[:, [snp_i, n_snp + snp_i]])
            mix_geno = admix.add_up_haplotype(phgeno[:, [snp_i, n_snp + snp_i]])
            study_index = admix.sample_case_control(mix_pheno)
            score_df = admix.mixscore_wrapper(pheno=mix_pheno[study_index], 
                                              anc=mix_anc[study_index, :], 
                                              geno=mix_geno[study_index, :], 
                                              theta=theta[study_index])
            tractor_geno = admix.convert_anc_count(anc=anc[:, [snp_i, n_snp + snp_i]], phgeno=phgeno[:, [snp_i, n_snp + snp_i]])
            score_df["TRACTOR"] = admix.tractor(pheno=mix_pheno[study_index],
                                                anc=mix_anc[study_index, :], 
                                                geno=tractor_geno[study_index, :], 
                                                theta=theta[study_index])
            score_df["SIM_I"] = sim_i
            score_df["SNP_I"] = snp_i
            score_df_list.append(score_df)
    return pd.concat(score_df_list)

def finemap(phgeno, anc, theta, var_g, cov, n_sim=50, seed=1234, case_prevalence=0.2):
    """
    Run one experiment, both simulation and testing
    
    Args
    -----
    phgeno: phased genotype (2 * n_indiv, n_snp), odd and even rows represent the 
        two haplotypes per individual
    anc: local ancestry (2 * n_indiv, n_snp), same as `phgeno`
    theta: global ancestry proportion
    var_g: variance component parameter
    cov: covariance parameter
    
    Returns
    -----
    score_df: DataFrame
    beta: ground truth effect sizes vector
    """
    np.random.seed(seed)
    
    # form individuals
    n_indiv = anc.shape[0] // 2
    n_snp = anc.shape[1]
    anc = anc.reshape((n_indiv, n_snp * 2))
    phgeno = phgeno.reshape((n_indiv, n_snp * 2))
    
    beta, phe_g, phe = admix.simulate_phenotype_cc(phgeno, anc, var_g=var_g, cov=cov, n_sim=n_sim, n_causal=1, case_prevalence=0.2)
    score_df_list = []
    
    mix_anc = admix.add_up_haplotype(anc)
    mix_geno = admix.add_up_haplotype(phgeno)
    for sim_i in range(n_sim):
        print(sim_i)
        sim_phe = phe[:, sim_i]
        # sample case control
        study_index = admix.sample_case_control(sim_phe)
        score_df = admix.mixscore_wrapper(pheno=sim_phe[study_index], 
                                          anc=mix_anc[study_index, :], 
                                          geno=mix_geno[study_index, :], 
                                          theta=theta[study_index])
        print(score_df)
        score_df["TRACTOR"] = admix.tractor_multi_snp(pheno=sim_phe[study_index], 
                                                      phgeno=phgeno[study_index, :], 
                                                      anc=anc[study_index, :], 
                                                      theta=theta[study_index])
        score_df["SNP_I"] = np.arange(n_snp)
        score_df["SIM_I"] = sim_i
        score_df_list.append(score_df)
    score_df = pd.concat(score_df_list)
    return score_df, beta