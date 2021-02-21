import pandas as pd
import numpy as np
from os.path import join
import fire
import admix


def single_snp_test_tractor(
    phgeno,
    anc,
    theta,
    odds_ratio,
    anc_effect,
    n_sim=50,
    seed=1234,
    case_prevalence=0.1,
    control_ratio=2.5,
):
    """
    Run one experiment, both simulation and testing

    Args
    -----
    phgeno: phased genotype (2 * n_indiv, n_snp), odd and even rows represent the
        two haplotypes per individual
    anc: local ancestry (2 * n_indiv, n_snp), same as `phgeno`
    theta: global ancestry proportion

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

    sim_list = admix.simulate_phenotype_cc_1snp_tractor(
        phgeno=phgeno,
        anc=anc,
        odds_ratio=odds_ratio,
        theta=theta,
        anc_effect=anc_effect,
        n_sim=n_sim,
        case_prevalence=case_prevalence,
    )
    score_df_list = []
    for snp_i, snp_sim in enumerate(sim_list):
        print(snp_i)
        for sim_i in range(n_sim):
            mix_pheno = snp_sim[1][:, sim_i]
            mix_anc = admix.add_up_haplotype(anc[:, [snp_i, n_snp + snp_i]])
            mix_geno = admix.add_up_haplotype(phgeno[:, [snp_i, n_snp + snp_i]])
            study_index = admix.sample_case_control(
                mix_pheno, control_ratio=control_ratio
            )
            print(
                f"Number of cases: {np.sum(mix_pheno[study_index] == 1)}, number of controls: {np.sum(mix_pheno[study_index] == 0)}"
            )
            score_df = admix.mixscore_wrapper(
                pheno=mix_pheno[study_index],
                anc=mix_anc[study_index, :],
                geno=mix_geno[study_index, :],
                theta=theta[study_index],
            )
            # logistic regression results
            tractor_geno = admix.convert_anc_count(
                anc=anc[:, [snp_i, n_snp + snp_i]],
                phgeno=phgeno[:, [snp_i, n_snp + snp_i]],
            )

            logistic_rls = admix.tractor(
                pheno=mix_pheno[study_index],
                anc=mix_anc[study_index, :],
                geno=tractor_geno[study_index, :],
                theta=theta[study_index],
            )
            for name in logistic_rls:
                score_df[name] = logistic_rls[name]

            score_df["SIM_I"] = sim_i
            score_df["SNP_I"] = snp_i
            print(score_df)
            score_df_list.append(score_df)
    return pd.concat(score_df_list)


def multi_snp_test_cli(
    data_dir,
    odds_ratio,
    anc_effect,
    region_size,
    region_i,
    out,
    n_sim=50,
    case_prevalence=0.1,
    control_ratio=2.5,
):
    # read phgeno, anc
    phgeno = np.load(join(data_dir, "phgeno.npy"))
    anc = np.load(join(data_dir, "anc.npy"))
    legend = pd.read_csv(join(data_dir, "legend.csv"))
    # calculating global ancestry
    theta = anc.reshape((anc.shape[0] // 2, anc.shape[1] * 2)).mean(axis=1)
    print(theta)
    print(np.mean(theta), np.std(theta))
    assert np.all(anc.shape == phgeno.shape)

    n_snp = anc.shape[1]
    region_index = np.arange(region_size * (region_i - 1), region_size * region_i)
    print(region_index)
    print(legend.iloc[region_index, :])
    phgeno = phgeno[:, region_index]
    anc = anc[:, region_index]
    print(
        f"case_prevalence: {case_prevalence}, control_ratio: {control_ratio}, odds_ratio: {odds_ratio}, anc_effect: {anc_effect}"
    )
    rls = multi_snp_test(
        phgeno=phgeno,
        anc=anc,
        theta=theta,
        odds_ratio=odds_ratio,
        anc_effect=anc_effect,
        case_prevalence=case_prevalence,
        control_ratio=control_ratio,
        seed=region_i,
        n_sim=n_sim,
    )

    # save
    rls.to_csv(out, index=False)


def multi_snp_test(
    phgeno,
    anc,
    theta,
    odds_ratio,
    anc_effect,
    n_sim=50,
    seed=1234,
    case_prevalence=0.1,
    control_ratio=2.5,
):
    """
    With a regional input, simulate single causal SNP at the center, and generate association statistics for all SNPs
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

    center_snp_index = n_snp // 2

    sim_list = admix.simulate_phenotype_cc_1snp_tractor(
        phgeno=phgeno[:, [center_snp_index, center_snp_index + n_snp]],
        anc=anc[:, [center_snp_index, center_snp_index + n_snp]],
        odds_ratio=odds_ratio,
        theta=theta,
        anc_effect=anc_effect,
        n_sim=n_sim,
        case_prevalence=case_prevalence,
    )
    assert len(sim_list) == 1
    snp_sim = sim_list[0]
    score_df_list = []
    mix_anc = admix.add_up_haplotype(anc)
    mix_geno = admix.add_up_haplotype(phgeno)

    for sim_i in range(n_sim):
        mix_pheno = snp_sim[1][:, sim_i]
        study_index = admix.sample_case_control(mix_pheno, control_ratio=control_ratio)
        print(
            f"Number of cases: {np.sum(mix_pheno[study_index] == 1)}, number of controls: {np.sum(mix_pheno[study_index] == 0)}"
        )

        score_df = admix.mixscore_wrapper(
            pheno=mix_pheno[study_index],
            anc=mix_anc[study_index, :],
            geno=mix_geno[study_index, :],
            theta=theta[study_index],
        )
        logistic_rls = admix.tractor_multi_snp(
            pheno=mix_pheno[study_index],
            phgeno=phgeno[study_index, :],
            anc=anc[study_index, :],
            theta=theta[study_index],
        )
        for name in logistic_rls:
            score_df[name] = logistic_rls[name]
        score_df["SNP_I"] = np.arange(n_snp)
        score_df["SIM_I"] = sim_i
        score_df_list.append(score_df)
        print(score_df)

    return pd.concat(score_df_list)

# Example:
"""
root_dir=/u/project/pasaniuc/pasaniucdata/admixture/kangcheng/genotype_simulation/out/kg_3k/
prefix=EUR_0.2_AFR_0.8_7_20000
python format_data.py format_data \
    --raw_dir ${root_dir}/${prefix} \
    --legend ${root_dir}/legend.txt \
    --out_dir data/geno/${prefix}
"""

def format_data(legend, raw_dir, out_dir):
    os.makedirs(out_dir)
    legend = pd.read_csv(legend, delim_whitespace=True)
    ancestry = admix.read_int_mat(join(raw_dir, "admix.hanc"))
    phgeno = admix.read_int_mat(join(raw_dir, "admix.phgeno")).T
    np.save(join(out_dir, "anc.npy"), ancestry)
    np.save(join(out_dir, "phgeno.npy"), phgeno)
    legend.to_csv(join(out_dir, "legend.csv"), index=False)


if __name__ == "__main__":
    fire.Fire()
