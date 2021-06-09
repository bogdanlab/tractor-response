import numpy as np
import allel
from os.path import join
import xarray as xr
import pandas as pd
from tqdm import tqdm
from functools import reduce


dir_plink = "01_plink/"
col_cov = ["AGE", "SEX", "dilution_factor"] + [f"PC{i + 1}" for i in range(10)]
col_trait = ["cholesterol", "ldl_direct"]

df_trait = [
    pd.read_csv(join(dir_plink, f"{trait}.pheno"), sep="\t") for trait in col_trait
]
df_trait = reduce(
    lambda left, right: pd.merge(left, right, on=["FID", "IID"]), df_trait
)
df_covar = pd.read_csv(join(dir_plink, f"covar.txt"), sep="\t")
df_pheno = pd.merge(df_trait, df_covar, on=["FID", "IID"])
df_pheno.index = df_pheno.FID.astype(str) + "_" + df_pheno.IID.astype(str)


for i_chr in tqdm(range(1, 23)):

    path_vcf = f"02_impute/processed/chr{i_chr}/chr{i_chr}.sample.typed.nochr.vcf.gz"
    vcf = allel.read_vcf(path_vcf)

    gt = vcf["calldata/GT"]
    assert (gt == -1).sum() == 0

    ds = xr.Dataset(
        data_vars={
            "geno": (("indiv", "snp", "haploid"), np.swapaxes(gt, 0, 1)),
        },
        coords={
            "snp": vcf["variants/ID"],
            "indiv": vcf["samples"],
            "chrom": ("snp", vcf["variants/CHROM"]),
            "position": ("snp", vcf["variants/POS"]),
            "snp_a0": ("snp", vcf["variants/REF"]),
            "snp_a1": ("snp", vcf["variants/ALT"][:, 0]),
        },
        attrs={"n_anc": 2},
    )

    # assign covariates
    dict_df = {
        col: ("indiv", df_pheno[col].reindex(ds.indiv)) for col in col_trait + col_cov
    }

    ds = ds.assign_coords(dict_df)

    # assign local ancestry
    rfmix = pd.read_csv(f"03_rfmix/lanc/chr{i_chr}.msp.tsv", sep="\t", skiprows=1)
    lanc_full = np.full((ds.dims["indiv"], ds.dims["snp"], ds.dims["haploid"]), np.nan)
    lanc0 = rfmix.loc[:, rfmix.columns.str.endswith(".0")].rename(
        columns=lambda x: x[:-2]
    )
    lanc1 = rfmix.loc[:, rfmix.columns.str.endswith(".1")].rename(
        columns=lambda x: x[:-2]
    )

    assert np.all(ds.indiv == lanc0.columns)
    for i_row, row in rfmix.iterrows():
        mask_row = (
            (row.spos <= ds.snp.position) & (ds.snp.position <= row.epos)
        ).values
        lanc_full[:, mask_row, 0] = lanc0.loc[i_row, :].values[:, np.newaxis]
        lanc_full[:, mask_row, 1] = lanc1.loc[i_row, :].values[:, np.newaxis]
    lanc_full = lanc_full.astype(np.int8)
    # EUR <-> AFR coding adjustment
    lanc_full = 1 - lanc_full
    ds = ds.assign(lanc=(("indiv", "snp", "haploid"), lanc_full))
    ds.to_zarr(f"04_assoc/dataset/chr{i_chr}.zarr")
