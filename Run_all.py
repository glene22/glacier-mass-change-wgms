"""
Run all - WGMS-AMCE downloader, nested-zip extractor, and plot/table runner
------------------------------------------
Author: Emily Glen
Date: 2025-08-16
------------------------------------------

This script:
- Downloads the WGMS dataset ZIP (default: WGMS-AMCE-2025-02b).
- Extracts the main ZIP and any nested ZIPs it contains.
- Moves all ZIP files into data/zipped_data/ to keep the workspace clean.
- Runs the three plotting exports:
    1) global_mass_change_timeseries
    2) plot_mass_change_anomaly
    3) export_mass_change_tables

Notes:
- If the dataset is already downloaded, it will not be re-downloaded.
- If you need to refresh the data, delete `data/wgms_dataset.zip` or the extracted folder(s).
"""

# ========= Imports =========
import argparse
import shutil
from pathlib import Path
from zipfile import ZipFile
import requests

# ========= User settings =========
URL_DEFAULT = "https://wgms.ch/downloads/wgms-amce-2025-02b.zip"
MAIN_ZIP_NAME = "wgms_dataset.zip"

# ========= Download, unzip, and plot =========
def run(url: str, data_dir: Path, fig_dir: Path):
    """End-to-end workflow: download -> extract (incl. nested zips) -> run plots."""
    # Ensure output folders exist
    data_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    # Stash for ZIPs after extraction (safe to delete anytime)
    to_delete_dir = data_dir / "zipped_data"
    to_delete_dir.mkdir(exist_ok=True)

    # 1) Download main ZIP if missing
    main_zip = data_dir / MAIN_ZIP_NAME
    if not main_zip.exists():
        print(f"Downloading {url} - {main_zip}")
        r = requests.get(url)
        r.raise_for_status()
        with open(main_zip, "wb") as fd:
            fd.write(r.content)

    # 2) Extract the main ZIP into data_dir
    with ZipFile(main_zip, "r") as zip_ref:
        zip_ref.extractall(data_dir)
    print(f"Extracted: {main_zip}")
    # Move the main ZIP into the stash folder after extracting
    shutil.move(str(main_zip), to_delete_dir / main_zip.name)

    # 3) Recursively extract *nested* ZIPs that came with the dataset
    processed = set()
    while True:
        # Find any ZIPs anywhere under data_dir that:
        #   - Haven't been processed yet, and
        #   - Are not located inside the zipped_data stash folder
        zips = [
            z for z in data_dir.rglob("*.zip")
            if z.resolve() not in processed and to_delete_dir not in z.parents
        ]
        if not zips:
            break

        for z in zips:
            with ZipFile(z, "r") as zip_ref:
                zip_ref.extractall(z.parent)
            print(f"Extracted: {z}")
            processed.add(z.resolve())
            shutil.move(str(z), to_delete_dir / z.name)

    # 4) Run plots/exports
    from plots.global_mass_change_timeseries import main as plot_one_main
    from plots.plot_mass_change_anomaly import main as plot_two_main
    from plots.export_mass_change_tables import main as plot_three_main

    print("\nRunning plots/exports …")
    plot_one_main(data_dir, fig_dir)    # saves global_plus_regions_timeseries.png
    plot_two_main(data_dir, fig_dir)    # saves anomaly_<Region>_<Year>_vs_<Baseline>.png
    plot_three_main(data_dir, fig_dir)  # saves CSVs under figures/tables/
    print(f"\nAll done. Figures/CSVs saved in: {fig_dir}")

# Anchor paths to the script’s folder
if __name__ == "__main__":
    repo_dir = Path(__file__).resolve().parent
    run(URL_DEFAULT, repo_dir / "data", repo_dir / "figures")

