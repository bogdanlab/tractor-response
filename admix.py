import numpy as np
import scipy
import warnings
import functools
from scipy.special import logit, expit
import tempfile
import subprocess
from os.path import join
import pandas as pd
import matplotlib.pyplot as plt
import numbers
from scipy.optimize import fsolve
import statsmodels.api as sm
from scipy import stats

def deprecated(func):
    """This is a decorator which can be used to mark functions
    as deprecated. It will result in a warning being emitted
    when the function is used."""
    @functools.wraps(func)
    def new_func(*args, **kwargs):
        warnings.simplefilter('always', DeprecationWarning)  # turn off filter
        warnings.warn("Call to deprecated function {}.".format(func.__name__),
                      category=DeprecationWarning,
                      stacklevel=2)
        warnings.simplefilter('default', DeprecationWarning)  # reset filter
        return func(*args, **kwargs)
    return new_func


def tractor(pheno, anc, geno, theta):
    """
    A reimplementation for Tractor
        (base) admixture only;
        (tractor) admixture + the number of copies of the risk allele on a EUR background + the number of copies on an AFR background.
    Args
    -----
    pheno: (n_indiv, ) phenotypes
    anc: (n_indiv, ) number of local ancestry
    geno: (n_indiv, n_anc) number of risk alleles for each ancestry
    theta: (n_indiv, ) global ancestry
    
    Returns
    -----
    p-value
    """
    base_design = np.hstack([sm.add_constant(anc), theta[:, np.newaxis]])
    tractor_design = np.hstack([base_design, geno])
    
#     base_model = sm.Logit(pheno, base_design).fit(disp=0, method="bfgs", maxiter=100)
#     tractor_model = sm.Logit(pheno, tractor_design).fit(disp=0, method="bfgs", maxiter=100)
    base_model = sm.Logit(pheno, base_design).fit(disp=0, method="bfgs", maxiter=200)
    tractor_model = sm.Logit(pheno, tractor_design).fit(disp=0, method="bfgs", maxiter=200)
    lr_stat = -2 * (base_model.llf - tractor_model.llf)
    return lr_stat 


def sample_case_control(pheno, ratio=1):
    """
    Sample case control from the population with a desired ratio
    Args
    -----
    pheno: (n_indiv, ) binary vector
    ratio: the ratio between control and case
    
    Returns
    -----
    index: (n_indiv, ) vector indicating whether i-th individual is sampled
    """
    case_index = np.where(pheno == 1)[0]
    control_index = np.random.choice(np.where(pheno == 0)[0], size=int(len(case_index) * ratio), replace=False)
    study_index = np.sort(np.concatenate([case_index, control_index]))
    return study_index

def add_up_haplotype(haplo):
    """
    Adding up the values from two haplotypes
    
    Args
    -----
    haplo: (n_indiv, 2 * n_snp) matrix
    
    Returns
    -----
    (n_indiv, n_snp) matrix with added up haplotypes
    """
    assert haplo.shape[1] % 2 == 0
    n_snp = haplo.shape[1] // 2
    return haplo[:, np.arange(n_snp)] + haplo[:, np.arange(n_snp) + n_snp]
    

def zsc2pval(zsc):
    return 1 - scipy.stats.norm.cdf(zsc)

def pval2zsc(pval):
    return -scipy.stats.norm.ppf(pval)

def chi2_to_logpval(chi2, dof=1):
    return stats.chi2.logsf(chi2, dof)

def read_int_mat(path):
    """
    Read a matrix of integer with [0-9], and with no delimiter.
    
    Args
    ----
    
    """
    with open(path) as f:
        mat = np.array([np.array([int(c) for c in line.strip()]) for line in f.readlines()], dtype=np.int8)
    return mat

def write_int_mat(path, mat):
    """
    Read a matrix of integer with [0-9], and with no delimiter.
    
    Args
    ----
    
    """
    np.savetxt(path, mat, fmt="%d", delimiter='')


def seperate_ld_blocks(anc, phgeno, legend, ld_blocks):
    assert len(legend) == anc.shape[1]
    assert len(legend) == phgeno.shape[1]
    
    rls_list = []
    for block_i, block in ld_blocks.iterrows():
        block_index = np.where((block.START <= legend.position) & (legend.position < block.STOP))[0]
        block_legend = legend.loc[block_index]
        block_anc = anc[:, block_index]
        block_phgeno = phgeno[:, block_index]
        rls_list.append((block_anc, block_phgeno, block_legend))
    return rls_list

def convert_anc_count(phgeno, anc):
    """
    Convert from ancestry and phased genotype to number of minor alles for each ancestry

    Args
    ----
    phgeno: n_indiv x 2n_snp, the first half columns contain the first haplotype,
        the second half columns contain the second haplotype
    anc: n_indiv x 2n_snp, match `phgeno`

    Returns
    ----
    geno: n_indiv x 2n_snp, the first half columns stores the number of minor alleles 
        from the first ancestry, the second half columns stores the number of minor
        alleles from the second ancestry
    """
    n_indiv = anc.shape[0]
    n_snp = anc.shape[1] // 2
    phgeno = phgeno.reshape((n_indiv * 2, n_snp))
    anc = anc.reshape((n_indiv * 2, n_snp))
    
    geno = np.zeros((n_indiv, n_snp * 2), dtype=np.int8)
    for indiv_i in range(n_indiv):
        for haplo_i in range(2 * indiv_i, 2 * indiv_i + 2):
            for anc_i in range(2):
                anc_snp_index = np.where(anc[haplo_i, :] == anc_i)[0]
                geno[indiv_i, anc_snp_index + anc_i * n_snp] += phgeno[haplo_i, anc_snp_index]
    return geno

def simulate_phenotype_cc(phgeno, anc, var_g, case_prevalence=0.1, n_causal=1, cov = 0., n_sim=30, seed=1234):

    """
    Simulate case control phenotype for admixture population 

    Args
    ----
    anc: #indiv x 2#snp, the first half columns contain the ancestry for the first haplotype,
                             the second half columns contain the ancestry for the second haplotype
    haplo: #indiv x 2#snp, matches `anc` 
    h2g: heritability
    n_causal: number of causal variants
    n_sim: number of simulation
    
    Returns
    ----
    beta: (2 * n_snp, n_sim) matrix
    phe: (n_indiv, n_sim) matrix
    """
    np.random.seed(seed)

    # phgeno
    # [0 1 0 | 1 1 0]
    # anc
    # [0 0 0 | 0 0 0]
    # geno
    # [1 2 0 | 0 0 0]
    geno = convert_anc_count(phgeno, anc)
    n_snp = geno.shape[1] // 2

    beta1 = np.zeros((n_snp, n_sim))
    beta2 = np.zeros((n_snp, n_sim))
    for i_sim in range(n_sim):
        cau = sorted(np.random.choice(np.arange(n_snp), size=n_causal, replace=False))
        beta = np.random.multivariate_normal(mean = [0., 0.], cov=[[1.0, cov],
                                                                    [cov, 1.0]], size=n_causal)
        beta1[cau, i_sim] = beta[:, 0]
        beta2[cau, i_sim] = beta[:, 1]

    # TODO: see popcorn paper for simulating beta1, beta2
    beta = np.vstack([beta1, beta2])
    phe_g = np.dot(geno, beta)
    phe = np.zeros_like(phe_g, dtype=np.int8)

    for sim_i in range(n_sim):
        scale = np.sqrt(var_g / np.var(phe_g[:, sim_i]))
        phe_g[:, sim_i] = phe_g[:, sim_i] * scale
        beta[:, sim_i] *= scale
        # find an intercept, such that the expectation is case_prevalence.
        func = lambda b : np.mean(expit(b + phe_g[:, sim_i])) - case_prevalence
        intercept = fsolve(func, logit(case_prevalence))
        phe[:, sim_i] = np.random.binomial(1, expit(intercept + phe_g[:, sim_i]))
        
    return beta, phe_g, phe

def simulate_phenotype_cc_1snp(phgeno, anc, var_g, cov=0.0, n_sim=10, case_prevalence=0.1):
    """
    Simulate case control phenotypes from phased genotype and ancestry (one SNP at a time)
    
    Args
    -----
    phgeno: phased genotype (n_indiv, 2 * n_snp)
    anc: local ancestry (n_indiv, 2 * n_snp)
    var_g: variance of genetic component
    cov: covariance between the allelic effects
    
    Returns
    -----
    n_snp list of tuple (beta, phe_g, phe)
    """
    n_indiv = phgeno.shape[0]
    n_snp = phgeno.shape[1] // 2
    
    return_list = []
    for snp_i in range(n_snp):
        snp_phgeno = phgeno[:, [snp_i, snp_i + n_snp]]
        snp_anc = anc[:, [snp_i, snp_i + n_snp]]
        beta, phe_g, phe = simulate_phenotype_cc(snp_phgeno, snp_anc, 
                                                 var_g=var_g, n_causal=1, cov=cov, n_sim=n_sim, case_prevalence=case_prevalence)
        return_list.append((beta, phe_g, phe))
    return return_list

def simulate_admix_phenotype_cont(phgeno, anc, h2g, n_causal, cov = 0., n_sim=30, seed=1234):
    """
    Simulate phenotype for admixture population [continuous]

    Args
    ----
    anc: #indiv x 2#snp, the first half columns contain the ancestry for the first haplotype,
                             the second half columns contain the ancestry for the second haplotype
    haplo: #indiv x 2#snp, matches `anc` 
    h2g: heritability
    n_causal: number of causal variants
    n_sim: number of simulation
    
    Returns
    ----
    beta: (2 * n_snp, n_sim) matrix
    phe: (n_indiv, n_sim) matrix
    """

    np.random.seed(seed)
    geno = convert_anc_count(phgeno, anc)
    n_snp = geno.shape[1] // 2

    beta1 = np.zeros((n_snp, n_sim))
    beta2 = np.zeros((n_snp, n_sim))
    for i_sim in range(n_sim):
        cau = sorted(np.random.choice(np.arange(n_snp), size=n_causal, replace=False))
        beta = np.random.multivariate_normal(mean = [0., 0.], cov=[[1.0, cov],
                                                                    [cov, 1.0]], size=n_causal)
        beta1[cau, i_sim] = beta[:, 0]
        beta2[cau, i_sim] = beta[:, 1]
        
    beta = np.vstack([beta1, beta2])

    phe_g = np.dot(geno, beta)
    phe_e = np.zeros_like(phe_g)

    for sim_i in range(n_sim):
        var_g = np.var(phe_g[:, sim_i])
        var_e = var_g * ( (1. / h2g) - 1)
        phe_e[:, sim_i] = np.random.normal(loc=0.0, scale=np.sqrt(var_e), size=n_indiv)
        
    phe = phe_g + phe_e
    return beta, phe_g, phe

@deprecated
def deprecated_simulate_phenotype_admix(anc, haplo, h2g, n_causal, cov = 0., n_sim=30, seed=1234):
    """
    Simulate phenotype for admixture population
    # Args
    anc: #indiv x 2#snp
        the first half columns contain the ancestry for the first haplotype,
        the second half columns contain the ancestry for the second haplotype
        index starting from 1
                            
    haplo: #indiv x 2#snp, matches `anc` 
    h2g: heritability
    n_causal: number of causal variants
    n_sim: number of simulation
    # Returns
    beta: (2 * n_snp, n_sim) matrix
    phe: (n_indiv, n_sim) matrix
    """
    np.random.seed(seed)
    # input checking
    assert anc.shape == haplo.shape
    n_indiv = anc.shape[0]
    n_snp = int(anc.shape[1] / 2)
    beta1 = np.zeros((n_snp, n_sim))
    beta2 = np.zeros((n_snp, n_sim))
    for i_sim in range(n_sim):
        cau = sorted(np.random.choice(np.arange(n_snp), size=n_causal, replace=False))
        beta = np.random.multivariate_normal(mean = [0., 0.], cov=[[1.0, cov],
                                                                    [cov, 1.0]], size=n_causal)
        beta1[cau, i_sim] = beta[:, 0]
        beta2[cau, i_sim] = beta[:, 1]

    beta = np.vstack([beta1, beta2])
    print(beta1.shape)
    phe_g = np.zeros((n_indiv, n_sim))
    phe_e = np.zeros((n_indiv, n_sim))
    
    haplo_index1 = np.arange(n_snp)
    haplo_index2 = np.arange(n_snp) + n_snp
    
    for i_sim in range(n_sim):
        sim_beta = np.zeros((n_indiv, n_snp * 2))
        for i_indiv in range(n_indiv):
            # copy effect sizes based on ancestry
            # first haplotype
            indiv_beta1 = np.zeros(n_snp)
            indiv_anc1 = anc[i_indiv, haplo_index1]
            
            indiv_beta1[indiv_anc1 == 0] = beta1[indiv_anc1 == 0, i_sim]
            indiv_beta1[indiv_anc1 == 1] = beta2[indiv_anc1 == 1, i_sim]
            
            # second haplotype
            indiv_beta2 = np.zeros(n_snp)
            indiv_anc2 = anc[i_indiv, haplo_index2]
            indiv_beta2[indiv_anc2 == 0] = beta1[indiv_anc2 == 0, i_sim]
            indiv_beta2[indiv_anc2 == 1] = beta2[indiv_anc2 == 1, i_sim]
            
            sim_beta[i_indiv, :] = np.concatenate([indiv_beta1, indiv_beta2])

        phe_g[:, i_sim] = np.sum(haplo * sim_beta, axis=1)
        var_g = np.var(phe_g[:, i_sim])
        var_e = var_g * ( (1. / h2g) - 1)
        phe_e[:, i_sim] = np.random.normal(loc=0.0, scale=np.sqrt(var_e), size=n_indiv)
    
    phe = phe_g + phe_e
    return beta, phe_g, phe


def mixscore_wrapper(pheno, anc, geno, theta, 
                    scores=["ADM", "ATT", "MIX", "SNP1", "SUM"],
                    mixscore_path="/u/project/pasaniuc/kangchen/tractor/software/mixscore-1.3/bin/mixscore",
                    verbose=False):
    """
    A python wrapper for mixscore
    
    Args
    ----
    pheno: phenotypes
    anc: ancestry
    geno: genotype
    theta: global ancestry component
    """
    
    tmp = tempfile.TemporaryDirectory()
    tmp_dir = tmp.name

    n_sample = len(pheno)
    n_snp = anc.shape[1]

    write_int_mat(join(tmp_dir, "pheno"), pheno.reshape((1, -1)))
    write_int_mat(join(tmp_dir, "anc"), anc.T)
    write_int_mat(join(tmp_dir, "geno"), geno.T)
    np.savetxt(join(tmp_dir, "theta"), theta, fmt='%.6f')

    param = {"nsamples": str(n_sample),
              "nsnps": str(n_snp),
              "phenofile": join(tmp_dir, "pheno"),
              "ancfile": join(tmp_dir, "anc"),
              "genofile": join(tmp_dir, "geno"),
              "thetafile": join(tmp_dir, "theta"),
              "outfile": join(tmp_dir, "out")}

    with open(join(tmp_dir, "param"), 'w') as f:
        f.writelines([k + ':' + param[k] + '\n' for k in param])
    
    rls_dict = {}
    for name in scores:
        if verbose:
            print(f"Calculating {name}...")
        
        cmd = ' '.join([mixscore_path, name, f"{tmp_dir}/param"])
        subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
        with open(param["outfile"]) as f:
            out = [line.strip() for line in f.readlines()]
        rls_dict[name] = out
    tmp.cleanup()
    score_df = pd.DataFrame(rls_dict).apply(pd.to_numeric, errors='coerce')
    return score_df