"""
_____________________________________________
name: snowMapper 
doi: 10.5281/zenodo.17663731
licence: MIT
contact: konstantis.alexopoulos@gmail.com
_____________________________________________

__________________________________________________________________________________________
Description:
Builds the training dataset from individual station files, applies snow-phase logic, 
stratified sampling, attaches climate variables, generates histograms, and saves as CSV.

Input parameters:
- stations_ee (ee.Feature.Collection): Earth Engine Feature Collection of 'stations.csv'.
- output_dir (str): Path of folder where processed station files are stored.
- sample_size (int): Number of datapoints to be sampled in the dataset.
- output_root (str): Path where the final dataset will be stored.
- months_list (list): List of months of in the configuration timeframe.
- SC_THRES (int): Theshold for conferting snow depth (cm) to bianry snow cover.
- CRS (str): Desired CRS of grid.
- SCALE (int): Desired scale of grid.

Internal functions:
- gap_filter(): 

Output:
- Final training dataset (.csv).
__________________________________________________________________________________________

"""

#===============================================================================
# Load libraries
#===============================================================================
import os
import pandas as pd
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt
import geemap
import ee

#===============================================================================
# Get Earth Engine started
#===============================================================================
ee.Authenticate()
ee.Initialize()

#===============================================================================
# Train dataset
#===============================================================================
def train_dataset(output_dir, stations_ee, sample_size, CRS, SCALE,
                  months_list, SC_THRES, output_root, sc_probab_col,
                  slope, aspect, chili, mtpi, landcover):

    #------------------------------------------
    # Sample image at station points
    #------------------------------------------
    def img_sampling(img, stations_ee, SCALE, CRS):
        return img.sampleRegions(
            collection=stations_ee,
            scale=SCALE,
            projection=CRS,
            properties=["STAID", "lat", "lon", 'elev'],
            tileScale=4,
            geometries=False
        )
    
    #------------------------------------------
    # Sample monthly snow probability image collection
    #------------------------------------------
    stations_sc_probab_ee = sc_probab_col.map(
        lambda img: img_sampling(img, stations_ee, SCALE, CRS)
    ).flatten()
    
    stations_sc_probab_df = geemap.ee_to_df(stations_sc_probab_ee)
    stations_sc_probab_df = stations_sc_probab_df[["month", "STAID", "sc_probab", "lat", "lon", 'elev']]
    
    sampled_df = stations_sc_probab_df.copy()
    
    #------------------------------------------
    # Sample terrain images
    #------------------------------------------
    terrain_img = slope.addBands([aspect, chili, mtpi, landcover])
    
    stations_terrain_ee = img_sampling(terrain_img, stations_ee, SCALE, CRS)
    stations_terrain_df = geemap.ee_to_df(stations_terrain_ee)
    
    stations_meta = stations_sc_probab_df.merge(
        stations_terrain_df[["STAID", "slope", "aspect", "mtpi", "chili", "landcover"]],
        on=["STAID"], 
        how="left"
    )
    
    #------------------------------------------
    # Load station timeseries files
    #------------------------------------------
    meta_cols = ["lat", "lon", "elev", "slope", "aspect", "mtpi", "chili", "landcover", "sc_probab", "month"]
    
    # keep unique station–month metadata rows
    stations_meta = (
        stations_meta
        .drop_duplicates(subset=["STAID", "month"])[["STAID"] + meta_cols]
        .set_index(["STAID", "month"])
    )
    stations_meta
    
    df_list = []
    
    for staid in stations_meta.index.get_level_values("STAID").unique():
        staid_str = f"{int(staid):06d}"
        df = pd.read_csv(os.path.join(output_dir, f"STAID{staid_str}.csv"))
    
        df["STAID"] = int(staid)
        df["date"] = pd.to_datetime(df["date"])
        df["month"] = df["date"].dt.month
        df["year"] = df["date"].dt.year
        df["DOY"]   = df["date"].dt.dayofyear
    
        # merge on (STAID, month)
        df = df.merge(
            stations_meta.reset_index(),
            on=["STAID", "month"],
            how="left"
        )
    
        df_list.append(df)
    
    merged_df = pd.concat(df_list, ignore_index=True)
    
    #------------------------------------------
    # Snow-phase logic
    #------------------------------------------
    df = merged_df.copy()
    df["sc_phase"] = "unknown"
    df.loc[(df.sc == 0) & (df.sc_sum_prev == 0), "sc_phase"] = "no_snow"
    df.loc[(df.sc == 1) & (df.sc_sum_prev == 0), "sc_phase"] = "new_snow"
    df.loc[(df.sc == 1) & (df.sc_sum_prev >= 1), "sc_phase"] = "snow"
    df.loc[(df.sc == 0) & (df.sc_sum_prev >= 1), "sc_phase"] = "melted_snow"
    
    df = df.dropna()
    df = df[df["month"].isin(months_list)]
    
    df = df[~(
        ((df.sc_phase == 'new_snow')    & (df.cd_precip < 10.0)) |
        ((df.sc_phase == 'new_snow')    & (df.hd_precip > 0.0)) |
        ((df.sc_phase == 'new_snow')    & (df.hd > 0.0)) |
        ((df.sc_phase == 'new_snow')    & (df.cd == 0.0)) |
        ((df.sc_phase == 'new_snow')    & (df.t2mHres > 0.0)) |
        ((df.sc_phase == 'melted_snow') & (df.cd_precip > 10.0)) |
        ((df.sc_phase == 'melted_snow') & (df.cd > 0.0)) |
        ((df.sc_phase == 'melted_snow') & (df.hd == 0.0)) |
        ((df.sc_phase == 'melted_snow') & (df.t2mHres < 0.0))
    )]
    
    df = df.drop_duplicates()
    
    #------------------------------------------
    # Stratified sampling
    #------------------------------------------
    cols_interest = ["sc_phase", "month", "elev", "aspect", "slope", "mtpi", "chili", "lat", "lon", "sc_probab"]
    df["stratum"] = df[cols_interest].astype(str).agg("-".join, axis=1)
    
    max_per_phase = sample_size // 4
    
    strat_sizes = df.groupby("stratum").size()
    proportional_sizes = (strat_sizes / len(df) * sample_size).astype(int)
    
    sampled_list = []
    phase_counts = defaultdict(int)
    
    for stratum, size in proportional_sizes.items():
        group = df[df["stratum"] == stratum]
        phase = group["sc_phase"].iloc[0]
    
        remaining = max_per_phase - phase_counts[phase]
        if remaining <= 0:
            continue
    
        n = min(size, len(group), remaining)
        if n <= 0:
            continue
    
        samp = group.sample(n=n, random_state=42)
        sampled_list.append(samp)
        phase_counts[phase] += n
    
    sampled_df = pd.concat(sampled_list)
    
    # Fill remaining
    if len(sampled_df) < sample_size:
        remaining = sample_size - len(sampled_df)
        eligible = df.copy()
        eligible["current"] = eligible["sc_phase"].map(phase_counts)
        eligible["max_allowed"] = eligible["sc_phase"].apply(lambda p: max_per_phase - phase_counts[p])
        eligible = eligible[eligible["max_allowed"] > 0]
    
        extras = []
        for phase, group in eligible.groupby("sc_phase"):
            take = min(max_per_phase - phase_counts[phase], len(group), remaining)
            if take > 0:
                extra = group.sample(n=take, random_state=42)
                extras.append(extra)
                phase_counts[phase] += take
                remaining -= take
            if remaining <= 0:
                break
    
        sampled_df = pd.concat([sampled_df] + extras)
    
    #------------------------------------------
    # Clean output
    #------------------------------------------
    drop_cols = [
        "stratum","max_allowed","current","date","hd","cd","hd_precip","cd_precip",
        "sc_sum","cd_sum","hd_sum","hd_precip_sum","cd_precip_sum"
    ]
    
    sampled_df = sampled_df.drop(columns=drop_cols, errors="ignore")
    sampled_df = sampled_df.dropna(subset=["sc_probab"]).drop_duplicates()
    
    #------------------------------------------
    # Histogram plots
    #------------------------------------------
    numeric_cols = sampled_df.select_dtypes(include=[np.number]).columns.tolist()
    categorical_cols = sampled_df.select_dtypes(include=['object','string']).columns.tolist()
    all_cols = numeric_cols + categorical_cols
    
    n = len(all_cols)
    cols = 3
    rows = int(np.ceil(n / cols))
    
    plt.figure(figsize=(15, max(rows * 3, 4)))
    
    for i, col in enumerate(all_cols, start=1):
        plt.subplot(rows, cols, i)
        if col in numeric_cols:
            plt.hist(sampled_df[col].dropna(), bins=50, alpha=0.7)
        else:
            counts = sampled_df[col].astype(str).value_counts()
            plt.bar(counts.index, counts.values, alpha=0.7)
            plt.xticks(rotation=90)
        plt.yscale("log")
        plt.xlabel(col)
        plt.ylabel("Count (log scale)")
        plt.title(f"Distribution of {col}")
    
    plt.tight_layout()
    fig_path = os.path.join(output_root, "training_data_histograms.png")
    plt.savefig(fig_path, dpi=300)
    plt.close()
    
    #------------------------------------------
    # Final output
    #------------------------------------------
    training_df = sampled_df.drop(columns=["sc_phase"], errors="ignore").reset_index(drop=True)
    
    training_data_path = os.path.join(output_root, f"training_{SCALE}m_SD{SC_THRES}_{int(sample_size/1000)}k.csv")
    training_df.to_csv(training_data_path, index=False)
    
    training_ee = geemap.df_to_ee(training_df, longitude='lon', latitude='lat')

    return training_ee