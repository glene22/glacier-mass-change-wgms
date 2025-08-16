"""
Global Glacier Mass Change Exploration
------------------------------------------
Author: Emily Glen
Date: 2025-08-16
------------------------------------------

Description:
    This script explores the WGMS global gridded annual glacier 
    mass change dataset (WGMS-AMCE-2025-02b). It processes data 
    for 1976-2024 and generates:

      1. Global time series of annual and cumulative glacier 
         mass change (Gt) with uncertainties.
      2. Regional (RGI regions) area-weighted means of annual 
         mass balance (m w.e./yr) and cumulative mass change.
      3. Combined figure with global + regional panels 
         (cumulative + annual bar plots).

User settings:
    - NcFile : Path to the WGMS dataset (.nc4 file)
    - Outdir : Directory for saving outputs

Outputs:
    - global_plus_regions_timeseries.png : Figure of global 
      and regional glacier mass change time series
"""

# ========= Imports =========
import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import MaxNLocator
from pathlib import Path
from RGI_Regions import RGI_REGIONS

# ========= User settings =========
NcFile = "data/wgms-amce-2025-02b/global-gridded/global-gridded-annual-glacier-mass-change.nc4"
Outdir = "figures"
os.makedirs(Outdir, exist_ok=True)

# ---- Single source of truth for the time window ----
START_YEAR = "1976"
END_YEAR = "2024"

# Ensure default output folder exists (caller may override in main)
os.makedirs(Outdir, exist_ok=True)

# ========= Styles =========
# Centralized plot style for consistent visuals.
bar_alpha = 0.5
bar_pos_color = "powderblue"
bar_neg_color = "indianred"
fill_alpha = 0.5
line_zero_style = dict(color="grey", lw=0.6, linestyle="--")
err_style = dict(ecolor="grey", lw=0.5)
small_fontsize = 7

# ========= NetCDF variable names =========
mwe_var = "glacier_mass_change_mwe"      # Mass change (m water equivalent)
uncert_mwe_var = "uncertainty_mwe"       # Uncertainty (m w.e.)
gt_var = "glacier_mass_change_gt"        # Mass change (gigatonnes)
uncert_gt_var = "uncertainty_gt"         # Uncertainty (gigatonnes)
area_var = "glacier_area_km2"            # Glacier area (km²)

# ========= Helpers =========
def load_dataset(ncfile: str, start: str = START_YEAR, end: str = END_YEAR) -> xr.Dataset:
    """Open the NetCDF dataset and subset to the given time range [start, end]."""
    return xr.open_dataset(ncfile).sel(time=slice(start, end))

def area_weighted_mean(da, area):
    # Compute area-weighted mean over lat/lon.
    return (da * area).sum(("lat", "lon")) / area.sum(("lat", "lon"))

def compute_global_series(ds):
    # Build global annual and cumulative series.
    area = ds[area_var]
    # Area-weighted annual mass balance and its uncertainty (m w.e. yr⁻¹)
    mean_mwe = area_weighted_mean(ds[mwe_var], area)
    mean_unc_mwe = area_weighted_mean(ds[uncert_mwe_var], area)
    # Aggregate annual Gt and uncertainty across space
    annual_gt = ds[gt_var].sum(("lat", "lon"))
    # Combine cell uncertainties via RSS (assume independent errors)
    annual_unc_gt = np.sqrt((ds[uncert_gt_var]**2).sum(("lat","lon")))
    # Sum over years
    cumulative_gt = annual_gt.cumsum("time")
    return mean_mwe, mean_unc_mwe, annual_gt, annual_unc_gt, cumulative_gt

def compute_region_series(ds: xr.Dataset, bounds):
    """
    Build regional annual/cumulative series for a bounding box.
    """
    # Unpack [lat_min, lat_max, lon_min, lon_max] from the region spec
    lat_min, lat_max, lon_min, lon_max = bounds

    # Subset to the region box.
    reg = ds.sel(lat=slice(lat_max, lat_min), lon=slice(lon_min, lon_max))

    # Fast exit: if the region has no non-NaN values for m w.e., there's nothing to plot.
    if reg[mwe_var].count() == 0:
        return None

    # Area-weighted annual mass balance (m w.e. yr-1) and its uncertainty.
    area = reg[area_var]
    mean_mwe = area_weighted_mean(reg[mwe_var], area)
    mean_unc = area_weighted_mean(reg[uncert_mwe_var], area)

    # Annual regional mass change in Gt: spatial sum over lat/lon.
    annual_gt = reg[gt_var].sum(("lat", "lon"))

    # Uncertainty propagation in Gt across grid cells:
    # Assume cell errors are independent / combine by root-sum-of-squares (RSS).
    annual_unc_gt = np.sqrt((reg[uncert_gt_var] ** 2).sum(("lat", "lon")))

    # Cumulative mass change (Gt): running sum over the time axis.
    cumulative_gt = annual_gt.cumsum("time")

    return {
        "mean_mwe": mean_mwe,
        "mean_unc": mean_unc,
        "annual_gt": annual_gt,
        "annual_unc_gt": annual_unc_gt,
        "cumulative_gt": cumulative_gt,
        "reg": reg,
    }

# ========= Plotting =========
def plot_global(ax, ds, bar_axes):
    """
    Plot global cumulative Gt (line + uncertainty band) and annual m w.e. bars.

    - Primary axis (ax): cumulative mass change in Gt
    - Secondary axis (ax2): annual mass balance (m w.e. yr⁻¹)
    """
    mean_mwe, mean_unc_mwe, annual_gt, annual_unc_gt, cumulative_gt = compute_global_series(ds)
    # Line: cumulative Gt; shaded band: cumulative uncertainty (RSS over years)
    ax.plot(ds.time.dt.year, cumulative_gt.values, color="black", lw=1)
    ax.fill_between(
        ds.time.dt.year,
        (cumulative_gt - annual_unc_gt.cumsum("time")).values,
        (cumulative_gt + annual_unc_gt.cumsum("time")).values,
        color="lightgrey", alpha=fill_alpha
    )
    ax.set_ylabel("cumulative mass change (Gt)")
    ax.set_xlabel("year")
    ax.set_title(f"global glacier mass change ({START_YEAR}–{END_YEAR})", fontsize=12)
    # Bars: annual mass balance (m w.e. yr⁻¹) with error bars (uncertainty m w.e.)
    ax2 = ax.twinx()
    colors = [bar_pos_color if float(v) > 0 else bar_neg_color for v in mean_mwe.values]
    ax2.bar(
        ds.time.dt.year, mean_mwe.values,
        yerr=mean_unc_mwe.values,
        color=colors, alpha=bar_alpha, edgecolor="grey",
        linewidth=0.5, capsize=2, error_kw=err_style
    )
    ax2.axhline(0, **line_zero_style)
    ax2.set_ylabel("annual mass change (m w.e./yr)")
    bar_axes.append(ax2)

    # Manual decade ticks
    years = ds.time.dt.year
    xticks = np.arange((years.min() // 10) * 10, years.max() + 1, 10)
    ax.set_xticks(xticks)
    ax.set_xticklabels([str(x) for x in xticks], fontsize=9)

    # Summary box: mean annual Gt and final cumulative Gt
    ax.text(
        0.02, 0.02,
        f"{annual_gt.mean().item():+.1f} Gt/yr\n{cumulative_gt[-1].item():+.0f} Gt",
        transform=ax.transAxes,
        fontsize=small_fontsize,
        va="bottom", ha="left",
        bbox=dict(facecolor="white", alpha=0.6, edgecolor="none", pad=1)
    )
    bar_axes.append(ax2)

def plot_region(ax, name, series, bar_axes):
    """
    Plot regional cumulative Gt (line + uncertainty band) and annual m w.e. bars.

    - Primary axis (ax): cumulative mass change (Gt)
    - Secondary axis (ax2): annual mass balance (m w.e. yr⁻¹)
    """
    reg = series["reg"]
    mean_mwe = series["mean_mwe"]
    mean_unc = series["mean_unc"]
    annual_gt = series["annual_gt"]
    annual_unc_gt = series["annual_unc_gt"]
    cumulative_gt = series["cumulative_gt"]

    # Line: cumulative Gt; shaded band: cumulative uncertainty (RSS of annual σ)
    ax.plot(reg.time.dt.year, cumulative_gt.values, color="black", lw=1)
    ax.fill_between(
        reg.time.dt.year,
        (cumulative_gt - annual_unc_gt.cumsum("time")).values,
        (cumulative_gt + annual_unc_gt.cumsum("time")).values,
        color="lightgrey", alpha=fill_alpha
    )

    # Bars: annual mass balance (m w.e. yr⁻¹) with error bars (uncertainty m w.e.)
    ax2 = ax.twinx()
    colors = [bar_pos_color if float(v) > 0 else bar_neg_color for v in mean_mwe.values]
    ax2.bar(
        reg.time.dt.year, mean_mwe.values,
        yerr=mean_unc.values,
        color=colors, alpha=bar_alpha, edgecolor="grey",
        linewidth=0.5, capsize=2, error_kw=err_style
    )
    ax2.axhline(0, **line_zero_style)

    # Titles and decade ticks (every 10 years) for the regional panel
    ax.set_title(name, fontsize=9)
    years = reg.time.dt.year  # regional years (not the full dataset)
    xticks = np.arange((years.min() // 10) * 10, years.max() + 1, 10)
    ax.set_xticks(xticks)
    ax.set_xticklabels([str(x) for x in xticks], fontsize=9)

    # Axis label font sizes
    ax.tick_params(axis="y", labelsize=small_fontsize)
    ax2.tick_params(axis="y", labelsize=small_fontsize)

    # Summary box: mean annual Gt and final cumulative Gt
    ax.text(
        0.02, 0.02,
        f"{annual_gt.mean().item():+.1f} Gt/yr\n{cumulative_gt[-1].item():+.0f} Gt",
        transform=ax.transAxes,
        fontsize=small_fontsize,
        va="bottom", ha="left",
        bbox=dict(facecolor="white", alpha=0.6, edgecolor="none", pad=1)
    )
    bar_axes.append(ax2)

# ========= Main =========
def main(data_dir: Path, fig_dir: Path) -> None:
    # Always use NcFile as configured at the top
    Outdir = str(fig_dir)
    os.makedirs(Outdir, exist_ok=True)

    ds = load_dataset(NcFile)

    # Layout: 2 rows for global (spanning 3 cols) + 6 rows for 18 regions (3 per row)
    fig = plt.figure(figsize=(12, 17))
    gs = gridspec.GridSpec(8, 3, figure=fig)
    bar_axes = []

    # Global panel (top two rows, all columns)
    ax_g = fig.add_subplot(gs[0:2, 0:3])
    plot_global(ax_g, ds, bar_axes)

    # Regional panels: fixed layout (3×6), draw up to 18 regions — no messages
    for idx, (region, bounds) in enumerate(list(RGI_REGIONS.items())[:18]):
        row, col = idx // 3 + 2, idx % 3
        ax = fig.add_subplot(gs[row, col])
        series = compute_region_series(ds, bounds)
        plot_region(ax, region, series, bar_axes)

    # Harmonize secondary y-limits (annual m w.e.) across panels for comparability
    bar_axes[0].set_ylim(-1.5, 0.5)  # global
    for ax in bar_axes[1:]:
        ax.set_ylim(-2.5, 2.5)

    plt.tight_layout()
    plt.savefig(os.path.join(Outdir, "global_plus_regions_timeseries.png"), dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {os.path.join(Outdir, 'global_plus_regions_timeseries.png')}")
