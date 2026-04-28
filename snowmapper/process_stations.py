"""
_____________________________________________
name: snowMapper 
doi: 10.5281/zenodo.17663731
licence: MIT
contact: konstantis.alexopoulos@gmail.com
_____________________________________________

__________________________________________________________________________________________
Description:
Read station files, derive daily snow-cover flags at multiple thresholds, apply gap to 
keep only valid snow periods, compute cumulative features, and export tidy 
CSVs per threshold.

Input parameters:
- stations (str): Path for 'stations.csv' - a file containing the station metadata 
- input_dir (str): Path of folder where raw station files are stored.
- SC_THRES (int): Snow depth (cm) thresholds for binarization to snow cover.
- output_root (str): Path where the station output folder and processed data will be 
                     stored.

Internal functions:
- gap_filter(): One-pass gap filtering on the snow-cover flag 'sc'. 
                If a missing/1 appears in the middle of an sc=1 run that later resets 
                to 0, convert those interior NAs/1s to 0 to avoid artificial chaining 
                across gaps.
- remove_partial_snow_periods(): Remove snow periods (sc=1) that start or end next to 
                                 gaps in the DATE sequence. Marks invalid periods via a 
                                 boolean mask; current version returns the original df 
                                 (since final sc reset is commented out). Leave as-is if 
                                 you only want to *identify* partials; uncomment the 
                                 final assignment to zero them out.
                                 Logic:
                                 - Label contiguous runs of equal 'sc' with a period id.
                                 - For each snow run, require that the day before the 
                                   start and the day after the end exist and are 
                                   consecutive (±1 day). Otherwise mark it invalid.
- remove_all_snow_gaps(): If we have sc=1 ... [gap of >1 day] ... sc=1, drop the sc=1 
                          rows after the gap until the first sc=0 (or to the end if no 
                          zero occurs). This prevents a single snow period from spanning 
                          across date gaps.
- compute_features(): Compute degree-day and precipitation features during snow cover.
                      Uses simple forward iteration; assumes df sorted by DATE. 
                      Fixed: avoid calling .astype on scalars; use float() where needed.
                      Definitions:
                      - hd: heating degree-like metric for T>0  (positive temp)
                      - cd: cooling degree-like metric for T<=0 (absolute of 
                            negative temp)
                      - *_sum: running sums that *only* accumulate while sc==1
                      - *_precip: partition precipitation into hd_precip (T>0) and 
                                  cd_precip (T<=0)

Output:
- Path of folder with suffix '*_SD{SC_THRES}', containing post-processed station 
  data (.txt).
__________________________________________________________________________________________

"""
#===============================================================================
# Load libraries
#===============================================================================
import pandas as pd
import numpy as np
import glob
import os

#===============================================================================
# Utilities & loop
#===============================================================================
def process_stations(SC_THRES, stations, input_dir, output_root):
    SNOW_VAR = 'SD'
    TEMP_VAR = 'TG'
    DATE_VAR = 'DATE'
    LAT_VAR = 'lat'
    LON_VAR = 'lon'
    
    #------------------------------------------
    # Filter gaps
    #------------------------------------------
    def gap_filter(data_cube: pd.DataFrame) -> pd.DataFrame:
        df = data_cube.copy()
        reset_bool = False
        for i in range(1, len(df)):
            if df.loc[i-1, 'sc'] == 0:
                reset_bool = False
            elif pd.isna(df.loc[i, 'sc']) or df.loc[i, 'sc'] == 1:
                if reset_bool:
                    df.loc[i, 'sc'] = 0
            if df.loc[i-1, 'sc'] == 1 and pd.isna(df.loc[i, 'sc']):
                reset_bool = True
            if df.loc[i, 'sc'] == 0:
                reset_bool = False
        return df
    
    #------------------------------------------
    # Remove partial snow periods
    #------------------------------------------    
    def remove_partial_snow_periods(data_cube: pd.DataFrame) -> pd.DataFrame:
        df = data_cube.copy()
        df[DATE_VAR] = pd.to_datetime(df[DATE_VAR])
        df['sc_period'] = (df['sc'] != df['sc'].shift()).cumsum()
    
        valid = np.ones(len(df), dtype=bool)
    
        for period_id in df['sc_period'].unique():
            period_df = df[df['sc_period'] == period_id]
            if period_df['sc'].iloc[0] != 1:
                continue  # Only check snow periods
    
            start_idx = period_df.index[0]
            end_idx = period_df.index[-1]
    
            # Check previous date continuity
            if start_idx > 0:
                delta_prev = (df.loc[start_idx, DATE_VAR] - df.loc[start_idx - 1, DATE_VAR]).days
                if delta_prev > 1:
                    valid[start_idx:end_idx + 1] = False
                    continue
            else:
                valid[start_idx:end_idx + 1] = False
                continue
    
            # Check next date continuity
            if end_idx < len(df) - 1:
                delta_next = (df.loc[end_idx + 1, DATE_VAR] - df.loc[end_idx, DATE_VAR]).days
                if delta_next > 1:
                    valid[start_idx:end_idx + 1] = False
            else:
                valid[start_idx:end_idx + 1] = False
    
        # To actually remove partial snow periods from 'sc', uncomment:
        df.loc[~valid, 'sc'] = 0
        return df
    
    #------------------------------------------
    # Remove all snow gaps
    #------------------------------------------  
    def remove_all_snow_gaps(data_cube: pd.DataFrame, date_col: str = 'DATE') -> pd.DataFrame:
        df = data_cube.copy()
        df[date_col] = pd.to_datetime(df[date_col])
        df = df.sort_values(date_col).reset_index(drop=True)
    
        to_remove = set()
        i = 1
        while i < len(df):
            gap_days = (df.loc[i, date_col] - df.loc[i - 1, date_col]).days
            if df.loc[i - 1, 'sc'] == 1 and df.loc[i, 'sc'] == 1 and gap_days > 1:
                gap_start = i
    
                # Find first sc==0 after the gap start
                zero_found = False
                for j in range(gap_start, len(df)):
                    if df.loc[j, 'sc'] == 0:
                        zero_found = True
                        break
    
                if zero_found:
                    # Remove the sc=1 rows bridging the gap up to (but not including) that zero
                    to_remove.update(range(gap_start, j))
                    i = j  # Continue scanning from the zero
                else:
                    # No zero after the gap → drop everything from gap_start to end
                    to_remove.update(range(gap_start, len(df)))
                    break
            else:
                i += 1
    
        df_filtered = df.drop(index=to_remove).reset_index(drop=True)
        return df_filtered
    
    #------------------------------------------
    # Compute features
    #------------------------------------------    
    def compute_features(data_cube: pd.DataFrame) -> pd.DataFrame:
        df = data_cube.copy()
    
        # Base daily metrics
        df['sc_sum_prev'] = 0  # Initialized; will be replaced later by a shift
        df['hd'] = np.round(np.maximum(0, df[TEMP_VAR]).astype(float), 1)
        df['cd'] = np.round(np.maximum(0, 0 - df[TEMP_VAR]).astype(float), 1)
    
        # Running sums while sc==1
        df['hd_sum'] = 0.0
        df['cd_sum'] = 0.0
        df['sc_sum'] = 0.0
    
        # Precip partitions and their sums while sc==1
        df['precip_sum'] = 0.0
        df['precip_sum_prev'] = 0.0
        df['hd_precip'] = 0.0
        df['cd_precip'] = 0.0
        df['hd_precip_sum'] = 0.0
        df['cd_precip_sum'] = 0.0
    
        # Ensure DATE sorted for cumulative logic
        df[DATE_VAR] = pd.to_datetime(df[DATE_VAR])
        df = df.sort_values(DATE_VAR).reset_index(drop=True)
    
        for i in range(1, len(df)):
            # Split precipitation by temperature sign
            if df.loc[i, TEMP_VAR] <= 0:
                df.loc[i, 'cd_precip'] = df.loc[i, 'RR']
            else:
                df.loc[i, 'hd_precip'] = df.loc[i, 'RR']
    
            # Accumulate only while sc==1 (multiply by 0/1)
            sc_flag = float(df.loc[i, 'sc'])
    
            df.loc[i, 'sc_sum'] = (df.loc[i-1, 'sc_sum'] + df.loc[i, 'sc']) * sc_flag
    
            df.loc[i, 'hd_sum'] = np.round((df.loc[i-1, 'hd_sum'] + df.loc[i, 'hd']) * sc_flag, 1)
            df.loc[i, 'hd_precip_sum'] = np.round((df.loc[i-1, 'hd_precip_sum'] + df.loc[i, 'hd_precip']) * sc_flag, 1)
    
            df.loc[i, 'cd_sum'] = np.round((df.loc[i-1, 'cd_sum'] + df.loc[i, 'cd']) * sc_flag, 1)
            df.loc[i, 'cd_precip_sum'] = np.round((df.loc[i-1, 'cd_precip_sum'] + df.loc[i, 'cd_precip']) * sc_flag, 1)
    
        return df
    
    #------------------------------------------
    # Main processing loop
    #------------------------------------------
    # SC_THRES output folder (ensure it exists)
    output_dir = os.path.join(output_root, f"stations_daily_SD{SC_THRES}")
    os.makedirs(output_dir, exist_ok=True)

    for staid in stations["STAID"]:
        # Input file follows naming STAIDXXXXXX.txt in MERGED_DATA
        filename = f"STAID{int(staid):06d}.txt"
        file_path = os.path.join(input_dir, filename)
    
        # Read station data
        data_cube = pd.read_csv(file_path)
    
        # For saving names later
        base_name_no_ext = os.path.splitext(filename)[0]  # e.g., "STAID000123"

        # Work on a copy per SC_THRES
        df = data_cube.copy()

        # Binary snow-cover flag at the chosen SC_THRES
        df['sc'] = (df[SNOW_VAR] >= SC_THRES).astype(int)

        # Apply gap logic
        df = gap_filter(df)
        df = remove_partial_snow_periods(df)
        df = remove_all_snow_gaps(df)

        # Compute cumulative features
        df = compute_features(df)

        # Time handling & previous-day features
        df[DATE_VAR] = pd.to_datetime(df[DATE_VAR])
        df = df.set_index(DATE_VAR).sort_index()

        # Previous-day temperature (only if previous row is exactly previous calendar day)
        df['TG_prev'] = df[TEMP_VAR].shift(1)
        exact_prev_day = (df.index.to_series().diff() == pd.Timedelta(days=1))
        df['TG_prev'] = df['TG_prev'].where(exact_prev_day)

        # Lags of cumulative sums (prev day)
        for col in ['cd_sum', 'cd_precip_sum', 'sc_sum', 'hd_sum', 'hd_precip_sum']:
            df[f'{col}_prev'] = df[col].shift(1)

        # Drop rows where we don't have a strict previous day
        df = df.dropna(subset=['TG_prev']).reset_index()

        # Round numeric columns to 2 decimals for tidy output
        df = df.round(2)

        # Final column selection & renaming
        training_vars = [
            'DATE', 'TG', 'RR', 'sc',
            'sc_sum', 'sc_sum_prev',
            'hd', 'hd_sum', 'hd_sum_prev',
            'cd', 'cd_sum', 'cd_sum_prev',
            'hd_precip', 'hd_precip_sum', 'hd_precip_sum_prev',
            'cd_precip', 'cd_precip_sum', 'cd_precip_sum_prev'
        ]

        # Ensure all requested vars exist even if some are missing due to earlier filtering
        missing = [c for c in training_vars if c not in df.columns]
        for c in missing:
            df[c] = np.nan

        # Rename for consistency with typical downstream naming
        df = df[training_vars].rename(columns={"TG": "t2mHres", "RR": "precipHres", "DATE": "date"})

        # Compose output filename per station & SC_THRES
        out_file = os.path.join(output_dir, f"{base_name_no_ext}.csv")

        # Write CSV
        df.to_csv(out_file, index=False)

    return output_dir