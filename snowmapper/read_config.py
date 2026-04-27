"""
_____________________________________________
name: snowMapper 
doi: 10.5281/zenodo.17663731
_____________________________________________

__________________________________________________________________________________________
Description:
Reads a YAML configuration file (config.yml) and returns its values as a dictionary for 
use within a notebook.

Input parameters:
- config_file (yml): Configuration file.

Output:
- Dictionary containing all configuration settings.
__________________________________________________________________________________________

"""
#===============================================================================
# Load libraries
#===============================================================================
import yaml
import os
import geopandas as gpd
import pandas as pd
from datetime import datetime, timedelta
import ee
import geemap
import re

#===============================================================================
# Get Earth Engine started
#===============================================================================
ee.Authenticate()
ee.Initialize()

#===============================================================================
# Apply configuration settings
#===============================================================================
def read_config(config_file):
    
    # Ensure the config file exists
    if not os.path.exists(config_file):
        raise FileNotFoundError(f"The configuration file {config_file} was not found.")
    
    # Load YAML configuration file into a dictionary
    with open(config_file, "r") as file:
        cfg = yaml.safe_load(file)

    # Assign config values
    config_data = {}

    #------------------------------------------------------------------------------------------------------------------------------
    # EE project name
    #------------------------------------------------------------------------------------------------------------------------------
    config_data["ee_project"] = cfg.get("ee_project")
    
    #------------------------------------------------------------------------------------------------------------------------------
    # Geographic parameters
    #------------------------------------------------------------------------------------------------------------------------------
    config_data["SCALE"] = cfg.get("SCALE")
    config_data["CRS"] = cfg.get("CRS")

    #------------------------------------------------------------------------------------------------------------------------------
    # Domain & subdomain shapefiles
    #------------------------------------------------------------------------------------------------------------------------------
    config_data["roi_name"] = cfg.get("roi_name")
    config_data["roi_path"] = re.sub(r"\s+", "_", config_data["roi_name"].strip().lower())
    
    # Check if EE folders exist and create if not
    def asset_exists(asset_id):
        try:
            ee.data.getAsset(asset_id)
            return True
        except Exception:
            return False
    
    def ensure_folder(path):
        if not asset_exists(path):
            print("Creating folder:", path)
            ee.data.createFolder(path)
    
    root = f"projects/{config_data['ee_project']}/assets/"
    snowMapper_root = root + "snowMapper"
    roi_folder = snowMapper_root + f"/{config_data['roi_path']}"
    
    ensure_folder(snowMapper_root)
    ensure_folder(roi_folder)
    
    # config_data["stations"] = pd.read_csv(cfg["stations"], delimiter=",")
    # config_data["stations_ee"] = geemap.df_to_ee(config_data["stations"], latitude="lat", longitude="lon")
    
    stations_path = cfg.get("stations", "")
    # If CSV file provided, load and convert to EE asset
    if stations_path and stations_path.lower().endswith(".csv"):
        stations_df = pd.read_csv(stations_path)
        stations_df = stations_df[stations_df.is_valid]
        config_data["stations_ee"] = geemap.df_to_ee(stations_df, latitude="lat", longitude="lon")
    # If non-CSV ending and non-empty path, treat as EE asset
    elif stations_path:
        config_data["stations_ee"] = (ee.FeatureCollection(stations_path))
    # If path empty, then pass
    else:
        pass

    
    # Load domain local SHP path & convert to EE asset or load from EE asset path
    domain_path = cfg["domain_path"]
    if domain_path.lower().endswith(".shp"):    
        domain_gdf = gpd.read_file(domain_path).to_crs(epsg=4326)
        domain_gdf = domain_gdf[domain_gdf.is_valid]
        domain_ee = geemap.gdf_to_ee(domain_gdf)
        domain_ee_path = f"projects/{config_data['ee_project']}/assets/snowMapper/domain"
        task = ee.batch.Export.table.toAsset(
            collection=domain_ee,
            description="Export_Domain",
            assetId=domain_ee_path,
        )
        task.start()
        config_data["domain_ee"] = domain_ee.filter(ee.Filter.eq('MapName', config_data["roi_name"]))
        # config_data["domain_ee"] = ee.FeatureCollection(domain_ee_path).filter(ee.Filter.eq('MapName', config_data["roi_name"]))
    else:
        config_data["domain_ee"] = ee.FeatureCollection(domain_path).filter(ee.Filter.eq('MapName', config_data["roi_name"]))

    # Load subdomain local SHP path & convert to EE asset or load from EE asset path
    subdomain_path = cfg["subdomain_path"]
    if subdomain_path.lower().endswith(".shp"):    
        subdomain_gdf = gpd.read_file(subdomain_path).to_crs(epsg=4326)
        subdomain_gdf = subdomain_gdf[subdomain_gdf.is_valid]
        subdomain_ee = geemap.gdf_to_ee(subdomain_gdf)
        subdomain_ee_path = f"projects/{config_data['ee_project']}/assets/snowMapper/subdomain"
        task = ee.batch.Export.table.toAsset(
            collection=subdomain_ee,
            description="Export_Subdomain",
            assetId=subdomain_ee_path,
        )
        task.start()
        config_data["subdomain_ee"] = subdomain_ee.filter(ee.Filter.eq('MapName', config_data["roi_name"]))
        # config_data["subdomain_ee"] = ee.FeatureCollection(subdomain_ee_path).filter(ee.Filter.eq('MapName', config_data["roi_name"]))
    else:
        config_data["subdomain_ee"] = ee.FeatureCollection(subdomain_path).filter(ee.Filter.eq('MapName', config_data["roi_name"]))

    # Create domain bounds
    config_data["domain_bounds_ee"] = config_data["domain_ee"].geometry().buffer(config_data["SCALE"]).bounds(proj=config_data["CRS"], maxError=1)
    bounds = config_data["domain_ee"].geometry().bounds(proj=config_data["CRS"], maxError=1)
    coordList = ee.Array.cat(bounds.coordinates(), 1)
    xCoords = coordList.slice(1, 0, 1)
    yCoords = coordList.slice(1, 1, 2)
    xMin = xCoords.reduce('min', [0]).get([0,0])
    yMax = yCoords.reduce('max', [0]).get([0,0])
    config_data["CRS_TRANSFORM"] = ee.List([
        config_data["SCALE"],   # xScale
        0,                      # xShearing
        xMin,                   # xTranslation
        0,                      # yShearing
        -config_data["SCALE"],  # yScale
        yMax                    # yTranslation
    ]).getInfo()
    
    #------------------------------------------------------------------------------------------------------------------------------
    # Training data
    #------------------------------------------------------------------------------------------------------------------------------
    config_data['SC_THRES'] =  cfg["SC_THRES"]

    config_data['input_dir'] =  cfg["input_dir"]
    config_data['output_root'] =  cfg["output_root"]
    os.makedirs(config_data['output_root'], exist_ok=True)
    
    config_data['sample_size'] =  cfg["sample_size"]


    training_data_path = cfg.get("training_dataset", "")
    # If CSV file provided, load and convert to EE asset
    if training_data_path and training_data_path.lower().endswith(".csv"):
        training_df = pd.read_csv(training_data_path)
        training_df = training_df[training_df.is_valid]
        config_data["training_ee"] = geemap.df_to_ee(training_df, latitude="lat", longitude="lon")
    # If non-CSV ending and non-empty path, treat as EE asset
    elif training_data_path:
        config_data["training_ee"] = (
            ee.FeatureCollection(training_data_path)
              .filter(ee.Filter.eq('MapName', config_data["roi_name"]))
        )
    # If path empty, then pass
    else:
        pass

    
    classifier_path = cfg.get("classifier_path", "")
    # If classifier provided, load from EE asset path, otherwise pass.
    if classifier_path:
        config_data["classifier"] = ee.Classifier.load(classifier_path)
    else:
        pass
    
    #------------------------------------------------------------------------------------------------------------------------------
    # Timeframe parameters
    #------------------------------------------------------------------------------------------------------------------------------
    config_data["START_MONTH"] = cfg.get("START_MONTH")
    config_data["START_YEAR"] = cfg.get("START_YEAR")
    config_data["END_MONTH"] = cfg.get("END_MONTH")
    config_data["END_YEAR"] = cfg.get("END_YEAR")
    config_data["START_HOUR"] = cfg.get("START_HOUR")
    config_data["SPIN_UP_PERIOD"] = cfg.get("SPIN_UP_PERIOD")
    end_today = cfg.get("END_TODAY", False)
    config_data["SC_PROBAB_START_YEAR"] = cfg.get("SC_PROBAB_START_YEAR")
    config_data["SC_PROBAB_END_YEAR"] = cfg.get("SC_PROBAB_END_YEAR")

    # BUILD RECONSTRUCTION TIMEFRAME
    # Build start date
    start_year = config_data["START_YEAR"]
    start_month = config_data["START_MONTH"]
    start_date = pd.Timestamp(f"{start_year}-{start_month:02d}-01")
    
    # Determine end date
    if not end_today:
        end_year = config_data["END_YEAR"]
        end_month = config_data["END_MONTH"]
        if end_year is None or end_month is None:
            raise ValueError("END_YEAR and END_MONTH must be provided when TODAY=False.")
        base_end = pd.Timestamp(f"{end_year}-{end_month:02d}-01") + pd.offsets.MonthEnd(0)
    
        # Enforce max 12-month rule
        month_diff = (end_year - start_year) * 12 + (end_month - start_month)
        if month_diff > 11:
            raise ValueError(
                "Configuration cannot exceed one hydrological year / 12-month period (maximum)."
            )
        elif month_diff < 0:
            raise ValueError("END_YEAR/END_MONTH must be after START_YEAR/START_MONTH.")
    
    else:
        # TODAY=True then use today
        today = pd.Timestamp.today().normalize()
        base_end = today
    
        # Enforce max 12 months from start_date
        month_diff = (today.year - start_year) * 12 + (today.month - start_month)
        if month_diff > 11:
            raise ValueError(
                "TODAY=True would exceed one hydrological year / 12-month period (maximum)."
            )
        elif month_diff < 0:
            raise ValueError(
                "TODAY=True results in a date before START_YEAR/START_MONTH. Please re-configure."
            )
    
    # Make end_date exclusive (+1 day)
    end_date = base_end + pd.Timedelta(days=1)
    
    # Build list of selected months
    season_start = start_date
    season_end = base_end
    month_range = pd.date_range(start=season_start, end=season_end, freq="MS")
    months_list = [d.month for d in month_range]
    
    # Store results
    config_data["start_date"] = start_date.strftime("%Y-%m-%d")
    config_data["end_date"] = end_date.strftime("%Y-%m-%d")
    config_data["months_list"] = months_list

    
    # BUILD SNOW COVER PROBABILITY CALCULATION TIMEFRAME
    # Build start date
    probab_start_year = config_data["SC_PROBAB_START_YEAR"]
    probab_start_date = pd.Timestamp(f"{probab_start_year}-{start_month:02d}-01")
    
    # Determine end date
    probab_end_year = config_data["SC_PROBAB_END_YEAR"]
    probab_base_end = pd.Timestamp(f"{probab_end_year}-{end_month:02d}-01") + pd.offsets.MonthEnd(0)

    # Enforce max 12-month rule
    probab_year_diff = (probab_end_year - probab_start_year)
    if probab_year_diff < 10:
        raise ValueError(
            "Timeframe configuration for the calculation of multi-year monthly snow cover probabilities must be at least 10 years long."
        )
    
    # Make end_date exclusive (+1 day)
    probab_end_date = probab_base_end + pd.Timedelta(days=1)
    
    # Store results
    config_data["scprobab_start_date"] = probab_start_date.strftime("%Y-%m-%d")
    config_data["scprobab_end_date"] = probab_end_date.strftime("%Y-%m-%d")
    
    #------------------------------------------------------------------------------------------------------------------------------
    # Satellite missions
    #------------------------------------------------------------------------------------------------------------------------------
    config_data["missions"] = cfg["missions"]
    if not (
        config_data["missions"]["Sentinel_2"]
        or config_data["missions"]["Landsat_9"]
        or config_data["missions"]["Landsat_8"]
        or config_data["missions"]["Landsat_7"]
        or config_data["missions"]["Landsat_5"]
        or config_data["missions"]["Landsat_4"]
    ):
        raise ValueError("At least one satellite mission must be activated.")  
        
    #------------------------------------------------------------------------------------------------------------------------------
    # Snow classification parameters
    #------------------------------------------------------------------------------------------------------------------------------
    sc_method = cfg.get("sc_method", {})
    active_method = None
    thresholds = None

    for method_name, method_data in sc_method.items():
        if method_data.get("active", False):
            if active_method is not None:
                raise ValueError(f"Multiple active methods found: '{active_method}' and '{method_name}'. Only one method can be active.")
            active_method = method_name
            thresholds = method_data.get("thresholds", {})

    if active_method is None:
        raise ValueError("No active method found in 'sc_method'")

    # Store in config_data
    config_data["active_method"] = active_method
    config_data["thresholds"] = thresholds
    
    #------------------------------------------------------------------------------------------------------------------------------
    # Decision tree booleans & thresholds
    #------------------------------------------------------------------------------------------------------------------------------
    config_data["decision_tree_settings"] = cfg["decision_tree_settings"]
    if not (
        config_data["decision_tree_settings"]["temperature"]
        or config_data["decision_tree_settings"]["modis"]
        or config_data["decision_tree_settings"]["sc_probab"]
    ):
        raise ValueError("At least one decision-tree gap-filling method needs to be actived.")
        
    #------------------------------------------------------------------------------------------------------------------------------
    # Auxiliary data
    #------------------------------------------------------------------------------------------------------------------------------
    dem = ee.Image(cfg["dem"]).resample('bicubic')
    
    config_data["elev"] = dem.select(["elevation"], ["elev"])
    config_data["slope"] = ee.Terrain.slope(config_data["elev"])
    config_data["aspect"] = ee.Terrain.aspect(config_data["elev"])
    
    config_data["chili"] = (
        ee.Image(cfg["chili"])
          .select(["constant"], ["chili"])
          .resample("bicubic")
    )

    config_data["mtpi"] = (
        ee.Image(cfg["mtpi"])
          .select(["elevation"], ["mtpi"])
          .resample("bicubic")
    )
    
    config_data["landcover"] = (
        ee.Image(ee.ImageCollection(cfg["landcover"]).first())
          .remap([10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 100], ee.List.sequence(0, 10)).select(["remapped"], ["landcover"]).toByte()
    )

    config_data["clim"] = (
        ee.ImageCollection(cfg["clim"])
          .select(
              ["temperature_2m", "total_precipitation_hourly"],  # , "snow_cover"
              ["t2m", "precip"]                                  # , "fsc"
          )
    )
    
    config_data["clim_elev"] = ee.Image(cfg["clim_elev"])  

    config_data["constants"] = ee.Dictionary(cfg["constants"])

    config_data["sc_probab_col"] = ee.ImageCollection(cfg["sc_probab_col"])
    
    #------------------------------------------------------------------------------------------------------------------------------
    # Masks
    #------------------------------------------------------------------------------------------------------------------------------
    config_data["masks"] = cfg["masks"]
    if not (
        config_data["masks"]["glaciers"]
        or config_data["masks"]["water"]
        or config_data["masks"]["forest"]
        or config_data["masks"]["urban"]
    ):
        print("IMPORTANT: No masks, other than the domain (default), have been selected.")
    
    #------------------------------------------------------------------------------------------------------------------------------
    # Machine learning settings
    #------------------------------------------------------------------------------------------------------------------------------
    classifiers = cfg.get("classifiers", {})
    active_classifier = None
    settings = None

    for classifier_name, classifier_data in classifiers.items():
        if classifier_data.get("active", False):
            if active_classifier is not None:
                raise ValueError(f"Multiple active classifiers found: '{active_classifier}' and '{classifier_name}'. Only one classifier can be active.")
            active_classifier = classifier_name
            settings = classifier_data.get("settings", {})

    if active_classifier is None:
        raise ValueError("No active classifier found in 'classifiers'")

    # Store in config_data
    config_data["active_classifier"] = active_classifier
    config_data["settings"] = settings
    config_data["classifier_path"] = f"projects/{config_data['ee_project']}/assets/snowMapper/{config_data['active_classifier']}_classifier"

    
    config_data["input_vars"] = [k for k, v in cfg["input_vars"].items() if v]
    config_data["output_vars"] = [k for k, v in cfg["output_vars"].items() if v]
    
    return config_data