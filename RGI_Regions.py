"""
RGI region bounding boxes
------------------------------------------
Author: Emily Glen
Date: 2025-08-16
------------------------------------------
-  (coarse lat/lon extents) for quick subsetting/plots.
- Region names follow the first-order Randolph Glacier Inventory (RGI) scheme.
- Boxes are **approximate**, not official polygons—use RGI region shapefiles for precise masks.
- Coordinates in degrees; longitude in [-180, 180].

Format: region_name: [lat_min, lat_max, lon_min, lon_max]
"""

# ========= Regions =========

RGI_REGIONS = {
    "Alaska": [50, 72, -170, -130],
    "Western Canada & US": [42, 70, -130, -110],
    "Arctic Canada North": [70, 85, -120, -60],
    "Arctic Canada South": [60, 75, -90, -50],
    "Greenland": [58, 83, -75, -10],
    "Iceland": [63, 67.5, -25, -12],
    "Svalbard": [74, 82, 5, 35],
    "Scandinavia": [58, 71, 5, 25],
    "Russian Arctic": [70, 82, 35, 180],
    "Central Europe": [44, 48, 5, 15],
    "Caucasus & Middle East": [39, 44, 40, 50],
    "Central Asia": [35, 50, 65, 95],
    "South Asia West": [30, 40, 70, 80],
    "South Asia East": [25, 35, 80, 100],
    "Low Latitudes": [-5, 5, -80, -60],
    "Southern Andes": [-56, -23, -80, -65],
    "New Zealand": [-47, -43, 166, 171],
    "Antarctic & Subantarctic": [-90, -60, -180, 180]
}
