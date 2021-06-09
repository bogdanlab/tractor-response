import numpy as np
import xarray as xr
from tqdm import tqdm
from typing import List
import statsmodels.api as sm
import pandas as pd
import admix
from scipy import stats

def fill_nan(mat):
    """
    fill nan with average of the mean
    """
    ret = mat.copy()
    col_mean = np.nanmean(mat, axis=0)
    inds = np.where(np.isnan(mat))
    ret[inds] = np.take(col_mean, inds[1])
    return ret

def marginal(ds: xr.Dataset, col_pheno: str, col_cov: List[str] = None, method: str = None, family: str = "gaussian"):
    """

    Parameters
    ----------
    ds
    pheno
    method
    cov
    family

    Returns
    -------
    pval
        Association p-values for each SNP being tested
    """
    
    assert family == "gaussian"
    if method is None:
        method = "ATT"
        
    assert method in ["ATT", "TRACTOR", "ADM", "SNP1"]
    
    pheno = ds[col_pheno].data
    cov = np.vstack([ds[col].data for col in col_cov]).T
    n_snp = len(ds['snp'])
    mask_indiv = ~np.isnan(pheno)
    
    # TODO: deal with missing genotype (current it is fine because of imputation)
    if method == "ATT":
        geno = np.sum(ds["geno"].data, axis=2)
        pvalues = []
        for i_snp in tqdm(range(n_snp)):
            design = np.hstack(
                [sm.add_constant(geno[:, i_snp][:, np.newaxis]), cov]
            )
            model = sm.GLM(pheno[mask_indiv], design[mask_indiv, :], family=sm.families.Gaussian()).fit()
            pvalues.append(model.pvalues[1])
        pvalues = np.array(pvalues)
        
    elif method == "SNP1":
        geno = np.sum(ds["geno"].data, axis=2)
        lanc = np.sum(ds["lanc"].data, axis=2)
        
        pvalues = []
        for i_snp in tqdm(range(n_snp)):

            design = np.hstack(
                [sm.add_constant(geno[:, i_snp][:, np.newaxis]), lanc[:, i_snp][:, np.newaxis], cov]
            )
            model = sm.GLM(pheno[mask_indiv], design[mask_indiv, :], family=sm.families.Gaussian()).fit()
            pvalues.append(model.pvalues[1])
        pvalues = np.array(pvalues)
    
    elif method == "TRACTOR":
        lanc = np.sum(ds["lanc"].data, axis=2)
        allele_per_anc = admix.data.compute_allele_per_anc(ds).compute()
        
        pvalues = []
        for i_snp in tqdm(range(n_snp)):
            
            design_null = np.hstack(
                [sm.add_constant(lanc[:, i_snp][:, np.newaxis]), cov]
            )
            model_null = sm.GLM(pheno[mask_indiv], design_null[mask_indiv, :], family=sm.families.Gaussian()).fit()
            design_alt = np.hstack([design_null, allele_per_anc[:, i_snp, :]])
            model_alt = sm.GLM(pheno[mask_indiv], design_alt[mask_indiv, :], family=sm.families.Gaussian()).fit(
                start_params=np.concatenate([model_null.params, [0.0, 0.0]])
            )
            pvalues.append(stats.chi2.sf(-2 * (model_null.llf - model_alt.llf), 2))
        pvalues = np.array(pvalues)
    
    elif method == "ADM":
        lanc = np.sum(ds["lanc"].data, axis=2)
        pvalues = []
        for i_snp in tqdm(range(n_snp)):
            design = np.hstack(
                [sm.add_constant(lanc[:, i_snp][:, np.newaxis]), cov]
            )
            model = sm.GLM(pheno[mask_indiv], design[mask_indiv, :], family=sm.families.Gaussian()).fit()
            pvalues.append(model.pvalues[1])
        pvalues = np.array(pvalues)
        
    else:
        raise NotImplementedError
        
    return pd.DataFrame({"SNP": ds.snp.values, "P": pvalues}).set_index("SNP")

def impute_dataset(vcf, ds):
    gt = vcf["calldata/GT"]
    assert (gt == -1).sum() == 0

    ds_imputed = xr.Dataset(
        data_vars={
            "geno": (("indiv", "snp", "haploid"), np.swapaxes(gt, 0, 1)),
        },
        coords={"snp": vcf["variants/ID"],
                "indiv": vcf["samples"],
                "chrom": ("snp", vcf["variants/CHROM"]),
                "position": ("snp", vcf["variants/POS"]),
                "snp_a0": ("snp", vcf["variants/REF"]),
                "snp_a1": ("snp", vcf["variants/ALT"][:, 0])},

        attrs={"n_anc": 2}
    )
    
    # fill in individual information
    dict_indiv = {}
    for col in ds["indiv"].coords:
        if col == "indiv":
            assert np.all(ds["indiv"] == ds_imputed["indiv"])
        else:
            dict_indiv[col] = ("indiv", ds[col])
    ds_imputed = ds_imputed.assign_coords(dict_indiv)
    
    # impute local ancestry
    
    # relevant typed region
    typed_start = np.where(ds['position'] < vcf['variants/POS'][0])[0][-1]
    typed_stop = np.where(ds['position'] > vcf['variants/POS'][-1])[0][0]
    ds_typed_subset = ds.isel(snp=slice(typed_start, typed_stop + 1))
    ds_typed_margin = ds_typed_subset.isel(snp=[0, -1])

    imputed_lanc = []
    for i_hap in range(2):
        df_typed_margin = pd.DataFrame(ds_typed_margin.lanc[:, i_hap].values.T, columns=ds_typed_margin.indiv.values, index=ds_typed_margin.snp.values)
        df_imputed = pd.DataFrame({
            "snp": ds_imputed.snp["snp"],
        }).set_index("snp")
        df_imputed = pd.concat([df_imputed, pd.DataFrame(columns=ds_imputed["indiv"].values, dtype=float)])
        # fill margin
        df_imputed = pd.concat([df_typed_margin.iloc[[0], :], df_imputed, df_typed_margin.iloc[[-1], :]], axis=0)
        df_imputed.index.name = "snp"
        # fill inside
        df_imputed.loc[ds_typed_subset.snp.values, ds_typed_subset.indiv.values] = ds_typed_subset["lanc"][:, :, i_hap].values.T
        # interpolate
        df_imputed = df_imputed.reset_index().interpolate(method="nearest").set_index("snp")
        
        imputed_lanc.append(df_imputed.loc[ds_imputed["snp"].values, ds_imputed["indiv"].values].values.astype(np.int8).T)
    
    ds_imputed = ds_imputed.assign(lanc=(("indiv", "snp", "haploid"), np.dstack(imputed_lanc)))
    return ds_imputed
