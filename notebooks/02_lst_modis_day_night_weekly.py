"""
MODIS Land Surface Temperature (LST) extraction (Day and Night)
aggregated by epidemiological week.

Study area:
- Santa Marta (AOI hosted in Google Earth Engine)

Description:
This script processes MODIS MOD11A1 daily LST data (day and night),
applies quality control masks, converts temperatures from Kelvin
to Celsius, aggregates by epidemiological week (median), and computes
zonal statistics over the AOI.

Inputs:
- Epidemiological weeks table (Excel): year, week (se), start date
- AOI FeatureCollection hosted in Google Earth Engine
- MODIS MOD11A1 daily LST product

Outputs:
- Weekly LST statistics (day and night) per AOI polygon

Notes:
- Google Earth Engine authentication and initialization are intentionally
  omitted for repository reproducibility
"""

import ee
import pandas as pd
import geopandas as gpd
import numpy as np

# ---------------------------------------------------------------------
# USER PARAMETERS
# ---------------------------------------------------------------------

YEAR_LABEL = "2024"

EPIDEMIO_WEEKS_FILE = "data/SE_2024.xlsx"
AOI_ASSET = "projects/ee-insaurraldeiibyt/assets/Proyecto_IAI/AOI_split_SantaMarta"

OUTPUT_FILE = f"LST_DAY_NIGHT_{YEAR_LABEL}.xlsx"

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


def mask_day(image):
    qa = image.select("QC_Day")
    mask = qa.bitwiseAnd(1 << 2).eq(0).And(
           qa.bitwiseAnd(1 << 3).eq(0))
    return image.updateMask(mask)


def mask_night(image):
    qa = image.select("QC_Night")
    mask = qa.bitwiseAnd(1 << 2).eq(0).And(
           qa.bitwiseAnd(1 << 3).eq(0))
    return image.updateMask(mask)


def kelvin_to_celsius_day(img):
    return (
        img.select("LST_Day_1km")
        .multiply(0.02)
        .subtract(273.15)
        .copyProperties(img, ["system:time_start"])
    )


def kelvin_to_celsius_night(img):
    return (
        img.select("LST_Night_1km")
        .multiply(0.02)
        .subtract(273.15)
        .copyProperties(img, ["system:time_start"])
    )

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

AOI = ee.FeatureCollection(AOI_ASSET)

# ---------------------------------------------------------------------
# IMAGE COLLECTIONS
# ---------------------------------------------------------------------

MODIS_DAY = ee.ImageCollection("MODIS/061/MOD11A1").select(
    ["LST_Day_1km", "QC_Day"]
)

MODIS_NIGHT = ee.ImageCollection("MODIS/061/MOD11A1").select(
    ["LST_Night_1km", "QC_Night"]
)

results_day = []
results_night = []

# ---------------------------------------------------------------------
# PROCESSING LOOP
# ---------------------------------------------------------------------

for i in range(len(weeks)):

    day_col = (
        MODIS_DAY
        .filterBounds(AOI)
        .filterDate(start_dates[i], end_dates[i])
        .map(mask_day)
        .map(kelvin_to_celsius_day)
    )

    night_col = (
        MODIS_NIGHT
        .filterBounds(AOI)
        .filterDate(start_dates[i], end_dates[i])
        .map(mask_night)
        .map(kelvin_to_celsius_night)
    )

    day_img = (
        ee.Image.constant(0).clip(AOI)
        if day_col.size().getInfo() == 0
        else day_col.median().clip(AOI)
    )

    night_img = (
        ee.Image.constant(0).clip(AOI)
        if night_col.size().getInfo() == 0
        else night_col.median().clip(AOI)
    )

    for stat in STATISTICS:
        d = zonal_stats(day_img, AOI, stat, SCALE)
        n = zonal_stats(night_img, AOI, stat, SCALE)

        results_day.append(d.drop(columns=["geometry", "gid"]))
        results_night.append(n.drop(columns=["geometry", "gid"]))

# ---------------------------------------------------------------------
# POST-PROCESSING
# ---------------------------------------------------------------------

def reshape_results(results):
    values = [str(x).split()[2] for x in map(str, results)]
    matrix = np.array(values).reshape(len(weeks), len(STATISTICS))
    return pd.DataFrame(matrix, index=weeks, columns=STATISTICS)

df_day = reshape_results(results_day)
df_night = reshape_results(results_night)

df_out = pd.concat(
    [df_day.add_prefix("DAY_"), df_night.add_prefix("NIGHT_")],
    axis=1
)

df_out.to_excel(OUTPUT_FILE)
