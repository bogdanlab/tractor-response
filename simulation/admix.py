import numpy as np
import scipy
from scipy.special import logit, expit
import tempfile
import subprocess
from os.path import join
import pandas as pd
from scipy.optimize import fsolve
import statsmodels.api as sm
from scipy import stats


def tractor(pheno, anc, geno, theta):
    """
    A reimplementation for Tractor
        (base) admixture only;
        (tractor) admixture + the number of copies of the risk allele on a EUR background + the number of copies on an AFR background.

        TEST1: y ~logit(localanc)
        TEST2: y ~logit (genotypes)
        TEST3: y ~logit (localanc + genotypes) , comparing to y~logit(localanc)
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
    # local ancestry
    m1_design = np.hstack([sm.add_constant(anc), theta[:, np.newaxis]])
    m1_model = sm.Logit(pheno, m1_design).fit(disp=0, method="bfgs", maxiter=200)

    # local ancestry + genotype (regardless of ancestry)
    m2_design = np.hstack([m1_design, geno.mean(axis=1)[:, np.newaxis]])
    m2_model = sm.Logit(pheno, m2_design).fit(
        disp=0,
        method="bfgs",
        maxiter=200,
        start_params=np.concatenate([m1_model.params, [0.0]]),
    )

    # local ancestry + genotype (ancestry aware)
    m3_design = np.hstack([m1_design, geno])
    m3_model = sm.Logit(pheno, m3_design).fit(
        disp=0,
        method="bfgs",
        maxiter=200,
        start_params=np.concatenate([m1_model.params, [0.0, 0.0]]),
    )

    # genotype (regardless of ancestry)
    att_design = np.hstack(
        [sm.add_constant(geno.mean(axis=1)[:, np.newaxis]), theta[:, np.newaxis]]
    )
    att_model = sm.Logit(pheno, att_design).fit(disp=0, method="bfgs", maxiter=200)

    rls_dict = {
        "ADM_LOGISTIC": m1_model.pvalues[1],
        "ATT_LOGISTIC": att_model.pvalues[1],
        "SNP1_LOGISTIC": stats.chi2.sf(-2 * (m1_model.llf - m2_model.llf), 1),
        "TRACTOR": stats.chi2.sf(-2 * (m1_model.llf - m3_model.llf), 2),
    }

    return rls_dict


def tractor_multi_snp(pheno, phgeno, anc, theta):
    """
    A convenient function for fitting multiple SNPs with the tractor model
    Args
    -----
    pheno: (n_indiv, ) phenotypes
    phgeno: (n_indiv, 2xn_snp) phased genotype
    anc: (n_indiv, 2xn_snp) local ancestry
    theta: (n_indiv, ) global ancestry

    Returns
    -----
    (n_snp, ) p-value
    """
    n_indiv = len(pheno)
    n_snp = phgeno.shape[1] // 2
    t_geno = convert_anc_count(anc=anc, phgeno=phgeno)
    t_anc = add_up_haplotype(anc)
    df = {
        "ADM_LOGISTIC": [],
        "ATT_LOGISTIC": [],
        "SNP1_LOGISTIC": [],
        "TRACTOR": [],
    }
    for snp_i in range(n_snp):
        pval_dict = tractor(
            pheno=pheno,
            anc=t_anc[:, snp_i],
            geno=t_geno[:, [snp_i, n_snp + snp_i]],
            theta=theta,
        )
        for name in df:
            df[name].append(pval_dict[name])
    return pd.DataFrame(df)


def sample_case_control(pheno, control_ratio=1):
    """
    Sample case control from the population with a desired ratio
    Args
    -----
    pheno: (n_indiv, ) binary vector
    ratio: the ratio of control / case

    Returns
    -----
    index: (n_indiv, ) vector indicating whether i-th individual is sampled
    """
    case_index = np.where(pheno == 1)[0]
    control_index = np.random.choice(
        np.where(pheno == 0)[0],
        size=int(len(case_index) * control_ratio),
        replace=False,
    )
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
        mat = np.array(
            [np.array([int(c) for c in line.strip()]) for line in f.readlines()],
            dtype=np.int8,
        )
    return mat


def write_int_mat(path, mat):
    """
    Read a matrix of integer with [0-9], and with no delimiter.

    Args
    ----

    """
    np.savetxt(path, mat, fmt="%d", delimiter="")


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
                geno[indiv_i, anc_snp_index + anc_i * n_snp] += phgeno[
                    haplo_i, anc_snp_index
                ]
    return geno


def convert_anc_count2(phgeno, anc):
    """
    Convert from ancestry and phased genotype to number of minor alles for each ancestry
    version 2, it should lead to exact the same results as `convert_anc_count`

    anc = np.random.randint(0, 2, size=(10, 6))
    phgeno = np.random.randint(0, 2, size=(10, 6))
    count1 = admix.convert_anc_count(phgeno=phgeno, anc=anc)
    count2 = convert_anc_count2(phgeno = phgeno, anc=anc)
    assert np.all(count1 == count2)

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
    n_anc = 2
    geno = np.zeros_like(phgeno)
    for haplo_i in range(2):
        haplo_slice = slice(haplo_i * n_snp, (haplo_i + 1) * n_snp)
        haplo_phgeno = phgeno[:, haplo_slice]
        haplo_anc = anc[:, haplo_slice]
        for anc_i in range(n_anc):
            geno[:, (anc_i * n_snp) : ((anc_i + 1) * n_snp)][
                haplo_anc == anc_i
            ] += haplo_phgeno[haplo_anc == anc_i]

    return geno


def simulate_phenotype_cc_1snp_tractor(
    phgeno, anc, odds_ratio, theta, anc_effect, n_sim=10, case_prevalence=0.1
):
    """
    Simulate case control phenotypes from phased genotype and ancestry (one SNP at a time)
    This function is specifically to mimic how Tractor paper simulate, as follows:
    In this model, the probability of disease was set to
    −2.19 + log[allelic risk effect size] × the number of copies of the minor allele coming from an AFR ancestral background +
        0.5 × AFR admixture proportion.
    Args
    -----
    phgeno: phased genotype (n_indiv, 2 * n_snp)
    anc: local ancestry (n_indiv, 2 * n_snp)
    odds_ratio: odds ratio for the effect sizes, assumed to be the same across population
    anc_effect: the effect sizes associated with global ancestry theta * anc_effect
    Returns
    -----
    n_snp list of tuple (beta, phe_g, phe)
    """
    n_indiv = phgeno.shape[0]
    n_snp = phgeno.shape[1] // 2
    assert len(theta) == n_indiv
    return_list = []
    # simulate snp by snp
    for snp_i in range(n_snp):
        snp_phgeno = phgeno[:, [snp_i, snp_i + n_snp]]
        snp_anc = anc[:, [snp_i, snp_i + n_snp]]
        snp_geno = convert_anc_count(snp_phgeno, snp_anc)
        # snp_phe_g: (n_indiv, n_sim)
        # allelic risk effect size x number of minor alleles
        snp_phe_g = np.dot(snp_geno, np.log(odds_ratio) * np.ones((2, n_sim)))
        # anc_effect x global ancestry
        snp_phe_g += np.dot(theta[:, np.newaxis], anc_effect * np.ones((1, n_sim)))
        snp_phe = np.zeros_like(snp_phe_g, dtype=np.int8)

        for sim_i in range(n_sim):
            # find an intercept, such that the expectation is case_prevalence.
            func = lambda b: np.mean(expit(b + snp_phe_g[:, sim_i])) - case_prevalence
            intercept = fsolve(func, logit(case_prevalence))
            snp_phe[:, sim_i] = np.random.binomial(
                1, expit(intercept + snp_phe_g[:, sim_i])
            )
        return_list.append((snp_phe_g, snp_phe))
    return return_list


def mixscore_wrapper(
    pheno,
    anc,
    geno,
    theta,
    scores=["ADM", "ATT", "MIX", "SNP1", "SUM"],
    mixscore_path="./software/mixscore-1.3/bin/mixscore",
    verbose=False,
):
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
    np.savetxt(join(tmp_dir, "theta"), theta, fmt="%.6f")

    param = {
        "nsamples": str(n_sample),
        "nsnps": str(n_snp),
        "phenofile": join(tmp_dir, "pheno"),
        "ancfile": join(tmp_dir, "anc"),
        "genofile": join(tmp_dir, "geno"),
        "thetafile": join(tmp_dir, "theta"),
        "outfile": join(tmp_dir, "out"),
    }

    with open(join(tmp_dir, "param"), "w") as f:
        f.writelines([k + ":" + param[k] + "\n" for k in param])

    rls_dict = {}
    for name in scores:
        if verbose:
            print(f"Calculating {name}...")

        cmd = " ".join([mixscore_path, name, f"{tmp_dir}/param"])
        subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
        with open(param["outfile"]) as f:
            out = [line.strip() for line in f.readlines()]
        rls_dict[name] = out
    tmp.cleanup()
    score_df = pd.DataFrame(rls_dict).apply(pd.to_numeric, errors="coerce")
    # convert to p-value
    for name in score_df.columns:
        score_df[name] = stats.chi2.sf(score_df[name], df=(2 if name == "SUM" else 1))
    return score_df