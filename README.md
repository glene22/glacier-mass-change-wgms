# WGMS Global Glacier Mass Change
**Emily Glen — 2025**

---
This repo downloads the **WGMS global gridded annual glacier mass change** dataset (version `WGMS-AMCE-2025-02b`) and generates:
- A **global + regional time-series** figure (annual + cumulative).
- A **spatial anomaly** map for a chosen region/year vs a baseline.
- **CSV tables** (global + hemispheres + top-10 loss years).

This repository uses the WGMS Annual Mass-Change Estimates (AMCE) global gridded product v2025-02b ([10.5904/wgms-amce-2025-02b](https://doi.org/10.5904/wgms-amce-2025-02b)) as described by Dussaillant et al. (2025),  ([10.5194/essd-17-1977-2025](https://doi.org/10.5194/essd-17-1977-2025)).

---

## 0. Get the code

Option A — SSH (recommended)  

```bash
git clone git@github.com:glene22/glacier-mass-change-wgms.git
cd glacier-mass-change-wgms
```

 Option B — HTTPS

```bash
git clone https://github.com/glene22/glacier-mass-change-wgms.git
cd glacier-mass-change-wgms
```

---

## 1. Contents

```
├── run_all.py                      # one-command download, unzip, and run all plots/exports
├── plots/
│   ├── global_mass_change_timeseries.py
│   ├── plot_mass_change_anomaly.py
│   └── export_mass_change_tables.py
├── RGI_Regions.py                  # RGI region bounding boxes (simplified)
├── data/                           # created at runtime (dataset + extracted files)
└── figures/                        # created at runtime (PNGs + tables/)
```
- The main ZIP and **all nested ZIPs** are extracted automatically.
- Original ZIPs are moved to `data/zipped_data/` (safe to delete).



---

## 2. To start

```bash
# 1) Create/activate an environment (conda recommended for Cartopy)
conda create -n venv python=3.11 -y
conda activate venv

# 2) Install dependencies
conda install -c conda-forge numpy xarray netcdf4 pandas matplotlib cartopy requests -y

# 3) Run everything
python run_all.py
```



---

## 3. What You Get

- **Figures**
  - `figures/global_plus_regions_timeseries.png`  
    Global cumulative (line + uncertainty band) and annual area-weighted mass balance bars per RGI region.

  - `figures/anomaly_<Region>_<Year>_vs_<Baseline>.png`  
    Spatial anomaly (m w.e.): `year` minus `baseline mean`.
    - Defaults: `Region=Greenland`, `Year=2023`, `Baseline=1991–2020`.
    - To change defaults, edit the constants at the top of `plots/plot_mass_change_anomaly.py`.

- **Tables (CSV)**
  - `figures/tables/global_annual_mass_change.csv`
  - `figures/tables/hemisphere_annual_mass_change.csv`
  - `figures/tables/top10_losses.csv`

---

## 4. Script Notes

### `plots/global_mass_change_timeseries.py`
- Loads `global-gridded-annual-glacier-mass-change.nc4`.
- Global series computed via area-weighted means (m w.e.) and spatial sums (Gt).
- Regional panels use simplified RGI bounding boxes from `RGI_Regions.py`.
- Uncertainties in Gt are added together using root-sum-of-squares. 

### `plots/plot_mass_change_anomaly.py`
- Subsets to a named RGI region (bounding box).
- Anomaly = mean of selected `Year` minus time-mean over `BaselineStart–BaselineEnd`.
- Uses **Cartopy** for basemap.

### `plots/export_mass_change_tables.py`
- Produces year-indexed tables for global + hemispheres.
- Also outputs **Top-10 most negative** mass-change years globally.

---

## 5. Requirements

- Python 3.9+
- Packages: `numpy`, `xarray`, `netcdf4`, `pandas`, `matplotlib`, `cartopy`, `requests`
---

