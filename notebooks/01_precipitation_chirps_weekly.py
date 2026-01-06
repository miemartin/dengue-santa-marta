"""
Precipitation data extraction from CHIRPS (daily) aggregated by epidemiological week.

Study area:
- Santa Marta (AOI hosted in Google Earth Engine)

Description:
This script aggregates daily CHIRPS precipitation data into weekly sums
based on epidemiological weeks defined in an external table. Zonal statistics
are computed for each AOI polygon.

Inputs:
- Epidemiological weeks table (Excel): year, week (se), start date
- AOI FeatureCollection hosted in Google Earth Engine
- CHIRPS daily precipitation (UCSB-CHG/CHIRPS/DAILY)

Outputs:
- Weekly precipitation statistics per AOI (COUNT, MIN, MEAN, MAX, MEDIAN, MODE,
  STD, SUM, VARIANCE) exported as an Excel table

Notes:
- Requires a configured Google Earth Engine environment
- Authentication and project initialization are intentionally omitted
  for repository reproducibility
"""

import ee
import pandas as pd
import geopandas as gpd
import numpy as np
from datetime import datetime

# ---------------------------------------------------------------------
# USER-DEFINED PARAMETERS
# ---------------------------------------------------------------------

YEAR_LABEL = "2024"

EPIDEMIO_WEEKS_FILE = "data/SE_2024.xlsx"
AOI_ASSET = "projects/ee-insaurraldeiibyt/assets/Proyecto_IAI/AOI_split_SantaMarta"

OUTPUT_FILE = f"PRECIPITATION_{YEAR_LABEL}.xlsx"

SCALE = 1000

STATISTICS = [
    "COUNT", "MINIMUM", "MEAN", "MAXIMUM",
    "MEDIAN", "MODE", "STD", "SUM", "VARIANCE"
]

# ---------------------------------------------------------------------
# FUNCTIONS
# ---------------------------------------------------------------------

def format_timestamps(timestamps):
    return [ts.strftime("%Y-%m-%d") for ts in timestamps]


def zonal_stats(image, features, statistic, scale):
    reducers = {
        "COUNT": ee.Reducer.count(),
        "MEAN": ee.Reducer.mean(),
        "MAXIMUM": ee.Reducer.max(),
        "MEDIAN": ee.Reducer.median(),
        "MINIMUM": ee.Reducer.min(),
        "MODE": ee.Reducer.mode(),
        "STD": ee.Reducer.stdDev(),
        "SUM": ee.Reducer.sum(),
        "VARIANCE": ee.Reducer.variance()
    }

    stats = image.reduceRegions(
        collection=features,
        reducer=reducers[statistic],
        scale=scale
    ).getInfo()

    return gpd.GeoDataFrame.from_features(stats, crs="EPSG:4326")


# ---------------------------------------------------------------------
# DATA LOADING
# ---------------------------------------------------------------------

weeks_df = pd.read_excel(EPIDEMIO_WEEKS_FILE)

years = weeks_df["year"].tolist()
weeks = weeks_df["se"].tolist()
start_dates = weeks_df["fecha_ini"].tolist()

end_dates = start_dates[1:]
end_dates.append(start_dates[-1] + pd.DateOffset(days=7))

start_dates = format_timestamps(start_dates)
end_dates = format_timestamps(end_dates)

# ---------------------------------------------------------------------
# EARTH ENGINE DATA
# ---------------------------------------------------------------------

AOI = ee.FeatureCollection(AOI_ASSET)

chirps_daily = (
    ee.ImageCollection("UCSB-CHG/CHIRPS/DAILY")
    .select("precipitation")
)

results = []

# ---------------------------------------------------------------------
# PROCESSING LOOP
# ---------------------------------------------------------------------

for i in range(len(weeks)):

    collection = (
        chirps_daily
        .filterBounds(AOI)
        .filterDate(start_dates[i], end_dates[i])
    )

    if collection.size().getInfo() == 0:
        image = ee.Image.constant(0).clip(AOI)
    else:
        image = collection.sum().clip(AOI)

    for stat in STATISTICS:
        zs = zonal_stats(image, AOI, stat, SCALE)
        zs = zs.drop(columns=["geometry", "gid"])
        results.append(zs)

# ---------------------------------------------------------------------
# POST-PROCESSING
# ---------------------------------------------------------------------

values = [str(x).split()[2] for x in map(str, results)]
matrix = np.array(values).reshape(len(weeks), len(STATISTICS))

df_out = pd.DataFrame(
    matrix,
    index=weeks,
    columns=STATISTICS
)

df_out.to_excel(OUTPUT_FILE)
