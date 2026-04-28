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
import xarray as xr
import geopandas as gpd
import pandas as pd
from datetime import datetime, timedelta
import ee
import geemap
import re
import rasterio
import time

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
        raise FileNotFoundError(f"Configuration file {config_file} not found")
    
    # Load YAML configuration file into a dictionary
    print(f" --- START -- Configuring snowMappar settings...")
    with open(config_file, 'r') as file:
        cfg = yaml.safe_load(file)

    # Assign config values
    config_data = {}

    #------------------------------------------------------------------------------------------------------------------------------
    # EE project name
    #------------------------------------------------------------------------------------------------------------------------------
    config_data['ee_project'] = cfg.get('ee_project')
    
    #------------------------------------------------------------------------------------------------------------------------------
    # Geographic parameters
    #------------------------------------------------------------------------------------------------------------------------------
    config_data['SCALE'] = cfg.get('SCALE')
    print(f"              Spatial resolution: {config_data['SCALE']}")
    
    config_data['CRS'] = cfg.get('CRS')
    print(f"              CRS: {config_data['CRS']}")

    #------------------------------------------------------------------------------------------------------------------------------
    # Domain & subdomain shapefiles
    #------------------------------------------------------------------------------------------------------------------------------
    config_data['roi_name'] = cfg.get('roi_name')
    if cfg.get('custom_path_name'):
        config_data['roi_path'] = cfg.get('custom_path_name')
    else:
        config_data['roi_path'] = re.sub(r'\s+', '_', config_data['roi_name'].strip().lower())
    print(f"              ROI name: {config_data['roi_name']}")
    print(f"              ROI path: {config_data['roi_path']}")

    # Check if EE folders exist and create if not
    def asset_exists(asset_id):
        try:
            ee.data.getAsset(asset_id)
            return True
        except Exception:
            return False
    
    def ensure_folder(path):
        if not asset_exists(path):
            print('              Creating folder:', path)
            ee.data.createFolder(path)
    
    root = f"projects/{config_data['ee_project']}/assets/"
    snowMapper_root = root + 'snowMapper'
    roi_folder = snowMapper_root + f"/{config_data['roi_path']}"
    
    ensure_folder(snowMapper_root)
    ensure_folder(roi_folder)
    
    # config_data['stations'] = pd.read_csv(cfg['stations'], delimiter=',')
    # config_data['stations_ee'] = geemap.df_to_ee(config_data['stations'], latitude='lat', longitude='lon')
    
    stations_path = cfg.get('stations', '')
    # If CSV file provided, load and convert to EE asset
    if stations_path and stations_path.lower().endswith('.csv'):
        stations_df = pd.read_csv(stations_path)
        stations_df = stations_df[stations_df.is_valid]
        config_data['stations_ee'] = geemap.df_to_ee(stations_df, latitude='lat', longitude='lon')
    # If non-CSV ending and non-empty path, treat as EE asset
    elif stations_path:
        config_data['stations_ee'] = (ee.FeatureCollection(stations_path))
    # If path empty, then pass
    else:
        pass

    
    # Load domain local SHP path & convert to EE asset or load from EE asset path
    domain_path = cfg['domain_path']
    if domain_path.lower().endswith('.shp'):    
        domain_gdf = gpd.read_file(domain_path).to_crs(epsg=4326)
        domain_gdf = domain_gdf[domain_gdf.is_valid]
        domain_ee = geemap.gdf_to_ee(domain_gdf)
        domain_ee_path = f"projects/{config_data['ee_project']}/assets/snowMapper/domain"
        task = ee.batch.Export.table.toAsset(
            collection=domain_ee,
            description='Export_Domain',
            assetId=domain_ee_path,
        )
        task.start()
        config_data['domain_ee'] = domain_ee.filter(ee.Filter.eq('MapName', config_data['roi_name']))
    else:
        config_data['domain_ee'] = ee.FeatureCollection(domain_path).filter(ee.Filter.eq('MapName', config_data['roi_name']))

    # Load subdomain local SHP path & convert to EE asset or load from EE asset path
    subdomain_path = cfg['subdomain_path']
    if subdomain_path.lower().endswith(".shp"):    
        subdomain_gdf = gpd.read_file(subdomain_path).to_crs(epsg=4326)
        subdomain_gdf = subdomain_gdf[subdomain_gdf.is_valid]
        subdomain_ee = geemap.gdf_to_ee(subdomain_gdf)
        subdomain_ee_path = f"projects/{config_data['ee_project']}/assets/snowMapper/subdomain"
        task = ee.batch.Export.table.toAsset(
            collection=subdomain_ee,
            description='Export_Subdomain',
            assetId=subdomain_ee_path,
        )
        task.start()
        config_data['subdomain_ee'] = subdomain_ee.filter(ee.Filter.eq('MapName', config_data['roi_name']))
    else:
        config_data['subdomain_ee'] = ee.FeatureCollection(subdomain_path).filter(ee.Filter.eq('MapName', config_data['roi_name']))

    # Create domain bounds
    config_data['domain_bounds_ee'] = config_data['domain_ee'].geometry().buffer(config_data['SCALE']).bounds(proj=config_data['CRS'], maxError=1)
    bounds = config_data['domain_ee'].geometry().bounds(proj=config_data['CRS'], maxError=1)
    coordList = ee.Array.cat(bounds.coordinates(), 1)
    xCoords = coordList.slice(1, 0, 1)
    yCoords = coordList.slice(1, 1, 2)
    xMin = xCoords.reduce('min', [0]).get([0,0])
    yMax = yCoords.reduce('max', [0]).get([0,0])
    config_data['CRS_TRANSFORM'] = ee.List([
        config_data['SCALE'],   # xScale
        0,                      # xShearing
        xMin,                   # xTranslation
        0,                      # yShearing
        -config_data['SCALE'],  # yScale
        yMax                    # yTranslation
    ]).getInfo()
    
    #------------------------------------------------------------------------------------------------------------------------------
    # Training data
    #------------------------------------------------------------------------------------------------------------------------------
    config_data['SC_THRES'] =  cfg['SC_THRES']

    config_data['input_dir'] =  cfg['input_dir']
    config_data['output_root'] =  cfg['output_root']
    os.makedirs(config_data['output_root'], exist_ok=True)
    
    config_data['sample_size'] =  cfg['sample_size']


    training_data_path = cfg.get('training_dataset', '')
    # If CSV file provided, load and convert to EE asset
    if training_data_path and training_data_path.lower().endswith('.csv'):
        training_df = pd.read_csv(training_data_path)
        training_df = training_df[training_df.is_valid]
        config_data['training_ee'] = geemap.df_to_ee(training_df, latitude='lat', longitude='lon')
    # If non-CSV ending and non-empty path, treat as EE asset
    elif training_data_path:
        config_data['training_ee'] = (
            ee.FeatureCollection(training_data_path)
              .filter(ee.Filter.eq('MapName', config_data['roi_name']))
        )
    # If path empty, then pass
    else:
        pass

    
    classifier_path = cfg.get('classifier_path', '')
    # If classifier provided, load from EE asset path, otherwise pass.
    if classifier_path:
        config_data['classifier'] = ee.Classifier.load(classifier_path)
    else:
        pass
    
    #------------------------------------------------------------------------------------------------------------------------------
    # Timeframe parameters
    #------------------------------------------------------------------------------------------------------------------------------
    start_date_cfg = cfg['start_date']
    end_date_cfg = cfg['end_date']
    
    config_data['START_HOUR'] = cfg.get('start_hour', 10)
    
    # Parse dates
    start_date = pd.Timestamp(start_date_cfg)
    base_end = pd.Timestamp(end_date_cfg)
    
    # Derive year/month for compatibility outputs
    start_year, start_month = start_date.year, start_date.month
    end_year, end_month = base_end.year, base_end.month
    
    # Enforce max 12-month window
    month_diff = (end_year - start_year) * 12 + (end_month - start_month)
    if month_diff > 11:
        raise ValueError(
            "Configuration cannot exceed one hydrological year / 12-month period (maximum)."
        )
    elif month_diff < 0:
        raise ValueError("end_date must be after start_date.")
    
    # Add a 1 day buffer in start & end
    initialisation_date = start_date + pd.Timedelta(days=-1)
    end_date = base_end + pd.Timedelta(days=1)
    
    # Build list of selected months
    month_range = pd.date_range(start=start_date, end=base_end, freq='MS')
    months_list = [d.month for d in month_range]
    
    # Store results
    config_data['initialisation_date'] = initialisation_date.strftime('%Y-%m-%d')
    config_data['start_date'] = start_date.strftime('%Y-%m-%d')
    config_data['end_date'] = end_date.strftime('%Y-%m-%d')
    config_data['months_list'] = months_list

    # STILL EXPOSE FOR DOWNSTREAM COMPATIBILITY
    config_data['START_YEAR'] = start_year
    config_data['END_YEAR'] = end_year

    print(f"              Start date: {config_data['start_date']}")
    print(f"              End date: {config_data['end_date']}")
    print(f"              Daily meteo aggregation time: {config_data['START_HOUR']}:00")
    
    # BUILD SNOW COVER PROBABILITY CALCULATION TIMEFRAME
    # Build start date
    config_data['SC_PROBAB_START_YEAR'] = cfg.get('SC_PROBAB_START_YEAR')
    config_data['SC_PROBAB_END_YEAR'] = cfg.get('SC_PROBAB_END_YEAR')
    probab_start_year = config_data['SC_PROBAB_START_YEAR']
    probab_start_date = pd.Timestamp(f"{probab_start_year}-{start_month:02d}-01")
    
    # Determine end date
    probab_end_year = config_data['SC_PROBAB_END_YEAR']
    probab_base_end = pd.Timestamp(f"{probab_end_year}-{end_month:02d}-01") + pd.offsets.MonthEnd(0)

    # Enforce max 12-month rule
    probab_year_diff = (probab_end_year - probab_start_year)
    if probab_year_diff < 10:
        raise ValueError(
            'Timeframe configuration for the calculation of multi-year monthly snow cover probabilities must be at least 10 years long.'
        )
    
    # Make end_date exclusive (+1 day)
    probab_end_date = probab_base_end + pd.Timedelta(days=1)
    
    # Store results
    config_data['scprobab_start_date'] = probab_start_date.strftime('%Y-%m-%d')
    config_data['scprobab_end_date'] = probab_end_date.strftime('%Y-%m-%d')
    
    #------------------------------------------------------------------------------------------------------------------------------
    # Satellite missions
    #------------------------------------------------------------------------------------------------------------------------------
    config_data['missions'] = cfg['missions']
    if not (
        config_data['missions']['Sentinel_2']
        or config_data['missions']['Landsat_9']
        or config_data['missions']['Landsat_8']
        or config_data['missions']['Landsat_7']
        or config_data['missions']['Landsat_5']
        or config_data['missions']['Landsat_4']
    ):
        raise ValueError('At least one satellite mission must be activated.')
        
    #------------------------------------------------------------------------------------------------------------------------------
    # Snow classification parameters
    #------------------------------------------------------------------------------------------------------------------------------
    sc_method = cfg.get('sc_method', {})
    active_method = None
    thresholds = None

    for method_name, method_data in sc_method.items():
        if method_data.get('active', False):
            if active_method is not None:
                raise ValueError(f"Multiple active methods found: '{active_method}' and '{method_name}'. Only one method can be active.")
            active_method = method_name
            thresholds = method_data.get('thresholds', {})

    if active_method is None:
        raise ValueError("No active method found in 'sc_method'")

    # Store in config_data
    config_data['active_method'] = active_method
    config_data['thresholds'] = thresholds

    print(f"              Snow mapping method: {config_data['active_method']}")
    if active_method == 'Otsu_ndsi' or active_method == 'Clustering_ndsi':
        print(f"-- WARNING -- {config_data['active_method']} is not yet fully operational")

    #------------------------------------------------------------------------------------------------------------------------------
    # Decision tree booleans & thresholds
    #------------------------------------------------------------------------------------------------------------------------------
    config_data['decision_tree_settings'] = cfg['decision_tree_settings']
    if not (
        config_data['decision_tree_settings']['temperature']
        or config_data['decision_tree_settings']['modis']
        or config_data['decision_tree_settings']['sc_probab']
    ):
        raise ValueError('At least one decision-tree gap-filling method needs to be actived.')
        
    #------------------------------------------------------------------------------------------------------------------------------
    # Auxiliary data
    #------------------------------------------------------------------------------------------------------------------------------
    metadata = cfg.get('metadata', {})
    
    # DEM (required)
    def load_dem(entry):
        if not entry or not entry.get('dir'):
            raise ValueError('DEM must be provided in config')
    
        asset_type = entry.get('type', 'image')  # default to image
    
        if asset_type == 'collection':
            return ee.Image(ee.ImageCollection(entry['dir']).mosaic()).resample('bicubic')
        elif asset_type == 'image':
            return ee.Image(entry['dir']).resample('bicubic')
        else:
            raise ValueError(f"Unknown asset type: {asset_type}")
            
    dem_entry = metadata.get('dem')
    dem = load_dem(dem_entry)
    config_data['elev'] = dem.select([dem_entry['band_name']], ['elev'])
    config_data['slope'] = ee.Terrain.slope(config_data['elev'])
    config_data['aspect'] = ee.Terrain.aspect(config_data['elev'])
    
    # Optional layers
    def load_image(entry, rename_to, resample=True):
        if not entry or not entry.get('dir'):
            return None
        
        img = ee.Image(entry['dir'])
        
        if entry.get('band_name'):
            img = img.select([entry['band_name']], [rename_to])
        
        if resample:
            img = img.resample('bicubic')
        
        return img
    
    config_data['chili'] = load_image(metadata.get('chili'), 'chili')
    config_data['mtpi']  = load_image(metadata.get('mtpi'), 'mtpi')  


    # Landcover
    def load_landcover(entry):
        if not entry or not entry.get('dir'):
            return None
    
        asset_type = entry.get('type', 'image')  # default to image
    
        if asset_type == 'collection':
            return ee.Image(ee.ImageCollection(entry['dir']).first())
        elif asset_type == 'image':
            return ee.Image(entry['dir'])
        else:
            raise ValueError(f'Unknown asset type: {asset_type}')
        
    lc_entry = metadata.get('landcover')
    
    if lc_entry and lc_entry.get('dir'):
        img = load_landcover(lc_entry)
    
        classes = lc_entry.get('classes', [])
        new_values = list(range(len(classes)))
    
        config_data['landcover'] = (
            img.select([lc_entry.get('band_name')])
               .remap(classes, new_values)
               .rename('landcover')
               .toByte()
        )
    else:
        config_data['landcover'] = None
    
    # Snow cover probability collection
    sc_entry = metadata.get('sc_probab_col')
    if sc_entry and sc_entry.get('dir'):
        config_data['sc_probab_col'] = ee.ImageCollection(sc_entry['dir'])
    else:
        config_data['sc_probab_col'] = None
    
    #------------------------------------------
    # Meteorological forcing data
    #------------------------------------------
    meteo = cfg.get('meteo', {})

    config_data['meteo_format'] = meteo.get('meteo_format')
    config_data['meteo_crs'] = meteo.get('meteo_crs')
    config_data['meteo_scale'] = meteo.get('meteo_scale')
    config_data['meteo_elev'] = ee.Image(f"{meteo.get('meteo_elev_dir')}").rename('meteo_elev')
    config_data['constants'] = ee.Dictionary(cfg['constants'])
   
    temperature = meteo.get('temperature', {})
    temperature_variable_name = temperature.get('var_name')
    temperature_units = temperature.get('units')
    precipitation = meteo.get('precipitation', {})
    precipitation_variable_name = precipitation.get('var_name')
    precipitation_units = precipitation.get('units')

    if config_data['meteo_format'] == 'ee':
        print('              Retrieving meteo from Earth Engine Data/Community Catalog...')
        meteo = (
            ee.ImageCollection(meteo.get('meteo_dir'))
            .filterDate(
                ee.Date(config_data['initialisation_date']).advance(config_data['START_HOUR'], 'hour'), 
                ee.Date(config_data['end_date']).advance(config_data['START_HOUR'], 'hour')
            )
            .select(
                [temperature_variable_name, precipitation_variable_name],
                ['t2m', 'precip']
            )

        )
        print('              Corresponding Earth Engine meteo assets have been found')

        # Convert temperature units
        if temperature_units == 'C':
            meteo_t2m = meteo.select('t2m')
            # pass
        elif temperature_units == 'K':
            print("              Converting temperature units to '℃'...")
            meteo_t2m = meteo.select('t2m').map(
                lambda img: img.select('t2m')
                               .subtract(ee.Image.constant(273.15))
                               .copyProperties(img, ["system:time_start"])
            )
        else:
            raise ValueError('Cannot process temperature: please provide in Celcius or Kelvin')

        # Convert precipitation units
        if precipitation_units == 'mm':
            meteo_precip = meteo.select('precip')
            # pass
        elif precipitation_units == 'm':
            print("              Converting precipitation units to 'mm'...")
            meteo_precip = meteo.select('precip').map(
                lambda img: img.select('precip')
                               .multiply(ee.Image.constant(1000))
                               .copyProperties(img, ["system:time_start"])
            )
        else:
            raise ValueError('Cannot process precipitation: please provide in millimetres or metres')

        config_data['meteo'] = ee.ImageCollection(meteo_t2m.combine(meteo_precip))
        print('              Meteo preprocessing successful')
    
    elif config_data['meteo_format'] == 'nc':
        print('              Retrieving meteo from local drive...')
        print('              Checking if already uploaded as Earth Engine assets...')
        try:
            # Earth Engine paths were data is meant to be stored
            meteo_t2m_ee_path = root + f"snowMapper/{config_data['roi_path']}/meteo_t2m_{config_data['initialisation_date']}_{config_data['end_date']}"
            meteo_precip_ee_path = root + f"snowMapper/{config_data['roi_path']}/meteo_precip_{config_data['initialisation_date']}_{config_data['end_date']}"
            
            ee.data.getAsset(meteo_t2m_ee_path)
            ee.data.getAsset(meteo_precip_ee_path)
            
            config_data['meteo'] = ''
            print('              Corresponding Earth Engine meteo assets have been found')
        
        except ee.EEException:
            print('              No corresponding Earth Engine meteo assets were found')
            print('              Preprocessing netCDF meteo files...')
                         
            # Open and prepare dataset
            ds = xr.open_dataset(meteo.get('meteo_dir'), engine='netcdf4')
    
            for var in ds.data_vars:
                da = ds[var]
                extra_dims = [d for d in da.dims if d not in {'lat', 'lon', 'time'}]
            
                if extra_dims:
                    ds[var] = da.mean(dim=extra_dims)
            ds = ds.squeeze(drop=True)       
    
            ds = ds[[temperature_variable_name, precipitation_variable_name]] \
                    .rename({temperature_variable_name: 't2m',
                             precipitation_variable_name: 'precip'})
            # Output directory
            out_dir_nc = './input_data/meteo/nc'
            os.makedirs(out_dir_nc, exist_ok=True)
            
            # Read config
            start_hour = int(config_data['START_HOUR'])
    
            ds_sel = ds.sel(time=slice((initialisation_date+pd.Timedelta(days=1)).strftime('%Y-%m-%d'), config_data['end_date']))
            
            # Shift time so your "day" starts at start_hour
            ds_shifted = ds_sel.assign_coords(time=ds_sel.time - pd.Timedelta(hours=start_hour))
            
            # Resample to daily
            ds_daily = xr.Dataset({
                't2m': ds_shifted['t2m'].resample(time='1D').mean(),
                'precip': ds_shifted['precip'].resample(time='1D').sum()
            })
    
            # Convert temperature units
            if temperature_units == 'C':
                pass
            elif temperature_units == 'K':
                print("              Converting temperature units to '℃'...")
                ds_daily['t2m'] = ds_daily['t2m'] - 273.15
            else:
                raise ValueError('Cannot process temperature: please provide in Celcius or Kelvin')
    
            # Convert precipitation units
            if precipitation_units == 'mm':
                pass
            elif precipitation_units == 'm':
                print("              Converting precipitation units to 'mm'...")
                ds_daily['precip'] = ds_daily['precip'] * 1000
            else:
                raise ValueError('Cannot process precipitation: please provide in millimetres or metres')
            
            # Shift time back
            ds_daily = ds_daily.assign_coords(time=ds_daily.time + pd.Timedelta(hours=start_hour))
            print("              Meteo preprocessing successful")
                
            # Save
            ds_daily.to_netcdf(os.path.join(out_dir_nc, f"meteo_{config_data['initialisation_date']}_{config_data['end_date']}.nc"))
            print(f"              Saved: meteo_{config_data['initialisation_date']}_{config_data['end_date']}.nc under '{out_dir_nc}'")
    
            # Output directory
            out_dir_tif = './input_data/meteo/geotif'
            os.makedirs(out_dir_tif, exist_ok=True)
            
            # Make sure spatial dimensions are set correctly
            ds_daily = ds_daily.rio.set_spatial_dims(x_dim='lon', y_dim='lat')
            ds_daily = ds_daily.rio.write_crs(config_data['meteo_crs'])
            
            def save_multiband(ds, var_name, out_path):
                da = ds[var_name]
            
                # Ensure ordering: (time, lat, lon)
                da = da.transpose("time", 'lat', 'lon')
            
                # Convert time to milliseconds since epoch
                time_ms = (da["time"].values.astype('datetime64[ms]')
                           .astype("int64"))
            
                # Write GeoTIFF
                da.rio.to_raster(
                    out_path,
                    driver='GTiff',
                    dtype=str(da.dtype),
                    count=da.sizes['time']
                )
            
                # Add band descriptions (timestamps)
                with rasterio.open(out_path, 'r+') as dst:
                    for i, t in enumerate(time_ms, start=1):
                        dst.set_band_description(i, str(t))
    
            # Save both variables
            t2m_assetId = f"meteo_t2m_{config_data['initialisation_date']}_{config_data['end_date']}"
            precip_assetId = f"meteo_precip_{config_data['initialisation_date']}_{config_data['end_date']}"
            save_multiband(ds_daily, 't2m', os.path.join(out_dir_tif, t2m_assetId + '.tif'))
            save_multiband(ds_daily, 'precip', os.path.join(out_dir_tif, precip_assetId + '.tif'))
            print(f"              Saved: meteo_{config_data['initialisation_date']}_{config_data['end_date']}.tif under '{out_dir_tif}'")
        
            # # Convert xarray to ee.Image
            # def xarray_to_ee_image(ds, var_name):
            #     da = ds[var_name].transpose("time", "lat", "lon")
                
            #     arr = da.values  # safe after size check
                
            #     n_bands = da.sizes["time"]
            #     band_names = [f"b{i}" for i in range(n_bands)]
                
            #     img = geemap.numpy_to_ee(
            #         arr,
            #         crs=config_data['meteo_crs'],
            #         # transform=config_data['CRS_TRANSFORM']
            #     )
                
            #     return img.rename(band_names)
            
            # # Export function
            # def export_to_ee(image, description, asset_id, config):
            #     task = ee.batch.Export.image.toAsset(
            #         image=image,
            #         description=description,
            #         assetId=asset_id,
            #         # region=config_data['domain_ee'],
            #         crs=config_data['meteo_crs'],
            #         # crsTransform=config_data['CRS_TRANSFORM'],
            #         scale=config_data['meteo_scale'],
            #         maxPixels=1e13
            #     )
            #     task.start()
            #     print(f"              {asset_id.split('/')[-1]} upload started")
            #     return task
            
            
            # # Check if asset exists
            # def asset_exists(asset_path):
            #     try:
            #         ee.data.getAsset(asset_path)
            #         return True
            #     except ee.EEException:
            #         return False
            
            # # Asset paths
            # base_asset_path = f"projects/{config_data['ee_project']}/assets/snowMapper/{config_data['roi_path']}"
            # roi_folder = base_asset_path
            
            # t2m_asset_full = f"{base_asset_path}/{t2m_assetId}"
            # precip_asset_full = f"{base_asset_path}/{precip_assetId}"
            
            
            # # Check size before loading into memory
            # n_elements = ds_daily['t2m'].size + ds_daily['precip'].size
            # MAX_ELEMENTS = 5e7  # adjust if needed (~50M elements)
            
            # if n_elements < MAX_ELEMENTS:
            #     print(f"              Meteo dataset size OK for automatic upload ({n_elements:.2e} elements)")
            
            #     # Convert xarray to ee.Image
            #     t2m_ee = xarray_to_ee_image(ds_daily, 't2m')
            #     precip_ee = xarray_to_ee_image(ds_daily, 'precip')
            
            #     # Launch exports
            #     exports = [
            #         (t2m_ee, 'export_meteo_t2m', t2m_asset_full),
            #         (precip_ee, 'export_meteo_precip', precip_asset_full),
            #     ]
            
            #     tasks = [export_to_ee(img, desc, aid, config_data) for img, desc, aid in exports]
            
            #     # Wait loop
            #     print("              Meteo files: uploading to Earth Engine...")
            
            #     while True:
            #         ready = [asset_exists(aid) for _, _, aid in exports]
                    
            #         if all(ready):
            #             print(f"              Meteo files successfully uploaded at {roi_folder}")
            #             break
                    
            #         print("              Meteo files: uploading to Earth Engine...")
            #         time.sleep(20)
            
            # else:
            #     print(
            #         f"-- WARNING -- Meteo GeoTIFF files too large for automatic upload. "
            #         f"              Before moving on to next steps you must first manually upload the meteo*.tif's "
            #         f".             on Earth Engine under: '{roi_folder}'"
            #     )

            config_data['meteo'] = ''
    
            print(f"-- WARNING -- Before proceeding to following steps, you must first manually upload the meteo*.tif's on Earth Engine under:")
            print(f"              '{roi_folder}'")
            
    else:
        raise ValueError(f"Meteo format must be 'nc' for netCDF or 'ee' for Earth Engine Image Collection.")

    #------------------------------------------------------------------------------------------------------------------------------
    # Masks
    #------------------------------------------------------------------------------------------------------------------------------
    config_data['masks'] = cfg['masks']
    if not (
        config_data['masks']['glaciers']
        or config_data['masks']['water']
        or config_data['masks']['forest']
        or config_data['masks']['urban']
    ):
        print('              IMPORTANT: No masks, other than the domain (default), have been selected.')
    
    #------------------------------------------------------------------------------------------------------------------------------
    # Machine learning settings
    #------------------------------------------------------------------------------------------------------------------------------
    classifiers = cfg.get('classifiers', {})
    active_classifier = None
    settings = None

    for classifier_name, classifier_data in classifiers.items():
        if classifier_data.get('active', False):
            if active_classifier is not None:
                raise ValueError(f"Multiple active classifiers found: '{active_classifier}' and '{classifier_name}'. Only one classifier can be active.")
            active_classifier = classifier_name
            settings = classifier_data.get('settings', {})

    if active_classifier is None:
        raise ValueError("No active classifier found in 'classifiers'")

    # Store in config_data
    config_data['active_classifier'] = active_classifier
    config_data['settings'] = settings
    config_data['classifier_path'] = f"projects/{config_data['ee_project']}/assets/snowMapper/{config_data['active_classifier']}_classifier"

    
    config_data['input_vars'] = [k for k, v in cfg['input_vars'].items() if v]
    config_data['output_vars'] = [k for k, v in cfg['output_vars'].items() if v]
    
    print(f" --- DONE --- Configuration successfully passed!")
    return config_data