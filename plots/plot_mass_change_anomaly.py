"""
Glacier Mass Change Anomaly Mapping
-----------------------------------
Author: Emily Glen
Date: 2025-08-15
------------------------------------------

Description:
    This script loads annual global glacier mass change data from a NetCDF file,
    subsets it to a specific RGI glacier region, and computes the anomaly for a
    given year relative to a baseline period. The anomaly is then plotted and
    saved as a PNG map.

User settings:
    Set these to control the run:
    - NC_FILE: Path to the dataset (.nc4 file)
    - REGION: Glacier region name (must be in RGI_REGIONS). 
    - Examples:
        [1] Alaska
        [2] Western Canada & US
        [3] Arctic Canada North
        [4] Arctic Canada South
        [5] Greenland
        [6] Iceland
        [7] Svalbard
        [8] Scandinavia
        [9] Russian Arctic
        [10] Central Europe
        [11] Caucasus & Middle East
        [12] Central Asia
        [13] South Asia West
        [14] South Asia East
        [15] Low Latitudes
        [16] Southern Andes
        [17] New Zealand
        [18] Antarctic & Subantarctic
    - YEAR: Single year to analyse (1976-2024)
    - BASELINE_START / BASELINE_END: Baseline period for comparison
    - OUTDIR: Directory where the map will be saved

Outputs:
    - PNG file showing the spatial anomaly for the selected year vs baseline.
"""

# ========= Imports =========
import os
import argparse
from pathlib import Path
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from RGI_Regions import RGI_REGIONS  # script with bounding boxes of RGI regions

# ========= User settings =========
NcFile = "data/wgms-amce-2025-02b/global-gridded/global-gridded-annual-glacier-mass-change.nc4"
Outdir = "figures"
Region = "Greenland"
Year = 2023
BaselineStart = 1991
BaselineEnd = 2020

mwe_var = "glacier_mass_change_mwe"  # Mass change in m water equivalent

# ========= Helpers =========
def subset_region(ds, bounds):
    lat_min, lat_max, lon_min, lon_max = bounds
    ds = ds.sortby("lat").sortby("lon")  # make sure coords are ascending
    a, b = sorted([lat_min, lat_max])    # lat bounds in order
    c, d = sorted([lon_min, lon_max])    # lon bounds in order
    return ds.sel(lat=slice(a, b), lon=slice(c, d))  # slice to bounding box

def mean_of_year(da, year):
    # Select all time steps in the given year
    sel = da.sel(time=str(year))
    # Compute the mean over time, ignoring NaNs
    return sel.mean("time", skipna=True)


def mean_of_baseline(da, start, end):
    # Select data between start and end years, then return the time-mean ignoring NaNs
    return da.sel(time=slice(f"{start}-01-01", f"{end}-12-31")).mean("time", skipna=True)


def compute_anomaly(ds, year, bstart, bend):
    # Return anomaly: mean(year) - mean(baseline)
    da = ds[mwe_var]
    return mean_of_year(da, year) - mean_of_baseline(da, bstart, bend)

# ========= Plotting =========
def plot_anomaly(anom, region, year, bstart, bend, out_png):
    """Plot glacier mass change anomaly with balanced colour scale."""
    lon, lat = anom.lon, anom.lat
    # Set up figure with a PlateCarree projection (identity for lon/lat in degrees, 
    # avoids reprojection and aligns coastlines/features correctly)
    fig = plt.figure(figsize=(8, 6))
    ax = plt.axes(projection=ccrs.PlateCarree()) 
    ax.set_extent([lon.min(), lon.max(), lat.min(), lat.max()])

    # Compute symmetric colour scale centred at zero (so that gain/loss is visually balanced)
    vmin = float(np.nanmin(anom))
    vmax = float(np.nanmax(anom))
    vlim = max(abs(vmin), abs(vmax))

    # Plot anomaly data
    im = ax.pcolormesh(
        lon, lat, anom,
        transform=ccrs.PlateCarree(),
        cmap="coolwarm_r",     
        shading="auto",
        vmin=-vlim, vmax=vlim # enforce symmetric range
    )

    # Add land and coastlines for reference
    ax.add_feature(cfeature.LAND, facecolor="lightgray", alpha=0.5)
    ax.coastlines(linewidth=0.6)

    # Add labelled gridlines
    gl = ax.gridlines(draw_labels=True, linestyle=":", linewidth=0.5, alpha=0.6)
    gl.right_labels = False
    gl.top_labels = False

    # Title
    ax.set_title(
        f"{region} - {year} glacier mass change anomaly\n"
        f"relative to {bstart}-{bend} mean",
        fontsize=11
    )

    # Colourbar
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Mass change anomaly (m w.e.)", fontsize=10)

    # Save and close
    plt.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_png}")


# ========= Main Function =========
def main(data_dir: Path, fig_dir: Path) -> None:
    """Run using data in data_dir and write outputs into fig_dir."""
    global NcFile, Outdir, Region, Year, BaselineStart, BaselineEnd

    # Prefer canonical dataset under data_dir; fall back to NcFile
    candidate = data_dir / "wgms-amce-2025-02b" / "global-gridded" / "global-gridded-annual-glacier-mass-change.nc4"
    if candidate.exists():
        NcFile = str(candidate)

    Outdir = str(fig_dir) if fig_dir is not None else Outdir
    os.makedirs(Outdir, exist_ok=True)

    # Validate region
    bounds = RGI_REGIONS.get(Region)
    if bounds is None:
        raise KeyError(f"Region '{Region}' not found. Options: {list(RGI_REGIONS.keys())}")

    print(f"Opening dataset and subsetting region: {Region}…")
    ds_reg = subset_region(xr.open_dataset(NcFile), bounds)

    print(f"Computing anomaly for {Year} vs {BaselineStart}-{BaselineEnd}…")
    anom = compute_anomaly(ds_reg, Year, BaselineStart, BaselineEnd)

    outdir = Path(Outdir); outdir.mkdir(parents=True, exist_ok=True)
    safe_region = Region.replace(" ", "_").replace("/", "-")
    out_png = outdir / f"anomaly_{safe_region}_{Year}_vs_{BaselineStart}-{BaselineEnd}.png"

    plot_anomaly(anom, Region, Year, BaselineStart, BaselineEnd, str(out_png))


