from utils import marginal
import xarray as xr
import fire
import statsmodels.api as sm
import pandas as pd
def assoc(ds, trait, out):
    ds = xr.open_zarr(ds, chunks=None)
    col_cov = ["AGE", "SEX", "dilution_factor"] + [f"PC{i + 1}" for i in range(10)]
    
    list_method = ["ATT", "SNP1", "TRACTOR", "ADM"]
    list_df = []
    for method in list_method:
        list_df.append(marginal(ds, col_pheno=trait, col_cov=col_cov, method=method))
    df_assoc = pd.concat(list_df, axis=1)
    df_assoc.columns = list_method
    df_assoc.to_csv(out)
    
if __name__ == "__main__":
    fire.Fire(assoc)
