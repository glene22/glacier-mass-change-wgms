"""
Global Glacier Mass Change Tables
------------------------------------------
Author: Emily Glen
Date: 2025-08-16
------------------------------------------

Description:
    This script processes the WGMS global gridded annual glacier 
    mass change dataset (.nc4 format). It generates:

      1. Global annual time series of glacier mass change (Gt),
         glacier area (km²), area-weighted mass balance (m w.e.),
         and cumulative mass change.
      2. Hemispheric (North/South) splits for the same metrics.
      3. A ranking of the top 10 years with the most negative 
         global mass change (highest losses).
      4. CSV outputs (rounded) for use in reports or plots.

User settings:
    - NcFile: Path to the dataset (.nc4 file)

Outputs:
    - global_annual_mass_change.csv
    - hemisphere_annual_mass_change.csv
    - top10_losses.csv
"""

# ========= Imports =========
import xarray as xr
import pandas as pd
import numpy as np
from pathlib import Path

# ========= User settings =========
# Repo-friendly default; will be overridden by main(data_dir, fig_dir) if present
NcFile = "data/wgms-amce-2025-02b/global-gridded/global-gridded-annual-glacier-mass-change.nc4"

# ========= NetCDF variable names =========
mwe_var        = "glacier_mass_change_mwe"   # Mass change in m water equivalent
uncert_mwe_var = "uncertainty_mwe"           # Uncertainty in m w.e.
gt_var         = "glacier_mass_change_gt"    # Mass change in gigatonnes
uncert_gt_var  = "uncertainty_gt"            # Uncertainty in gigatonnes
area_var       = "glacier_area_km2"          # Glacier area in km²

def main(data_dir: Path, fig_dir: Path) -> None:
    """Generate CSV tables; writes to fig_dir / 'tables' (or next to NcFile if you change one line)."""
    global NcFile

    # Prefer canonical path under data_dir; fall back to current NcFile
    candidate = data_dir / "wgms-amce-2025-02b" / "global-gridded" / "global-gridded-annual-glacier-mass-change.nc4"
    if candidate.exists():
        NcFile = str(candidate)

    # Open NetCDF
    ds = xr.open_dataset(NcFile)

    # Extract variables (filling NaNs only where safe)
    gt   = ds[gt_var].fillna(0)
    mwe  = ds[mwe_var]                 # keep NaNs for weighting
    area = ds[area_var].fillna(0)
    unc  = ds[uncert_gt_var].fillna(0) # using Gt uncertainty here   

    # Convert time axis to datetime, then extract integer years for index
    time = pd.to_datetime(ds["time"].values)
    years = pd.Index(pd.to_datetime(time).year, name="year")

    # Helper to sum over space
    Space_Dims = [d for d in gt.dims if d not in ["time"]]

    # ========= Global annual table =========
    global_gt = gt.sum(dim=Space_Dims)
    global_unc_rss = np.sqrt((unc ** 2).sum(dim=Space_Dims))
    global_area = area.sum(dim=Space_Dims)

    weighted_num = (mwe.fillna(0) * area).sum(dim=Space_Dims)
    weighted_den = global_area.where(global_area != 0)  # avoid divide-by-zero
    global_mwe_aw = (weighted_num / weighted_den)

    global_cum_gt = global_gt.cumsum(dim="time")

    df_global = pd.DataFrame({
        "mass_change_gt": global_gt.values,
        "uncertainty_gt_rss": global_unc_rss.values,
        "glacier_area_km2": global_area.values,
        "area_weighted_mwe": global_mwe_aw.values,
        "cumulative_mass_change_gt": global_cum_gt.values
    }, index=years)

    df_global_round = df_global.copy()
    df_global_round["mass_change_gt"] = df_global_round["mass_change_gt"].round(2)
    df_global_round["uncertainty_gt_rss"] = df_global_round["uncertainty_gt_rss"].round(2)
    df_global_round["glacier_area_km2"] = df_global_round["glacier_area_km2"].round(0)
    df_global_round["area_weighted_mwe"] = df_global_round["area_weighted_mwe"].round(4)
    df_global_round["cumulative_mass_change_gt"] = df_global_round["cumulative_mass_change_gt"].round(2)

    # ========= Hemispheric split (N/S) =========
    lat = ds["lat"]
    mask_N = xr.DataArray(lat >= 0, dims=("lat",))  # Northern Hemisphere
    mask_S = xr.DataArray(lat < 0, dims=("lat",))   # Southern Hemisphere

    def sum_where(mask, da):
        return da.where(mask).sum(dim=Space_Dims)

    def aw_mwe_where(mask):
        num = (mwe.fillna(0) * area).where(mask).sum(dim=Space_Dims)
        den = area.where(mask).sum(dim=Space_Dims)
        return (num / den.where(den != 0))

    hemi_gt_N  = sum_where(mask_N, gt)
    hemi_unc_N = np.sqrt(sum_where(mask_N, unc**2))
    hemi_area_N = sum_where(mask_N, area)
    hemi_mwe_N = aw_mwe_where(mask_N)

    hemi_gt_S  = sum_where(mask_S, gt)
    hemi_unc_S = np.sqrt(sum_where(mask_S, unc**2))
    hemi_area_S = sum_where(mask_S, area)
    hemi_mwe_S = aw_mwe_where(mask_S)

    df_hemi = pd.DataFrame({
        "mass_change_gt_N": hemi_gt_N.values,
        "uncertainty_gt_rss_N": hemi_unc_N.values,
        "glacier_area_km2_N": hemi_area_N.values,
        "area_weighted_mwe_N": hemi_mwe_N.values,
        "mass_change_gt_S": hemi_gt_S.values,
        "uncertainty_gt_rss_S": hemi_unc_S.values,
        "glacier_area_km2_S": hemi_area_S.values,
        "area_weighted_mwe_S": hemi_mwe_S.values,
    }, index=years)

    df_hemi_round = df_hemi.copy()
    for col in df_hemi_round.columns:
        if "mwe" in col:
            df_hemi_round[col] = df_hemi_round[col].round(4)
        elif "area" in col:
            df_hemi_round[col] = df_hemi_round[col].round(0)
        else:
            df_hemi_round[col] = df_hemi_round[col].round(2)

    # ========= Top 10 most negative years =========
    df_top10_round = (
        df_global.sort_values("mass_change_gt")
                 .head(10)
                 .round({
                     "mass_change_gt": 2,
                     "uncertainty_gt_rss": 2,
                     "glacier_area_km2": 0,
                     "area_weighted_mwe": 4,
                     "cumulative_mass_change_gt": 2
                 })
    )

    # ========= Write CSV =========
    outdir = (Path(fig_dir) / "tables") if fig_dir is not None else Path(NcFile).parent
    outdir.mkdir(parents=True, exist_ok=True)

    df_global_round.to_csv(outdir / "global_annual_mass_change.csv", index=True)
    df_hemi_round.to_csv(outdir / "hemisphere_annual_mass_change.csv", index=True)
    df_top10_round.to_csv(outdir / "top10_losses.csv", index=True)
    print(f"Saved CSVs in: {outdir}")


