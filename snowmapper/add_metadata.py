"""
_____________________________________________
name: snowMapper 
doi: 10.5281/zenodo.17663731
_____________________________________________

__________________________________________________________________________________________
Description: 
Adds metadata to the image collection as image bands. The metadata includes variables
that are later used for the reconstruction of snow cover.

Input parameters:
- collection (ee.ImageCollection): Preprocessed collection & synthetic images.
- domain_ee (ee.FeatureCollection): Region of interest (e.g. mountain range).
- constants (float): Constants for MicroMet meteo downscaling (Liston & Elder, 2006).
- input_vars (list): List of variables to be used as predictors for the PIML classifier.
- decision_tree_settings (dict): Contains thresholds for temperature-based 
                                 differentiation between rain and snow, and for
                                 binarising downsampled fractional snow cover
                                 from MODIS Terra NDSI.
- meteo_format (str): netCDF or ee.ImageCollection.
- meteo (ee.ImageCollection): Meteo data.
- roi_path (str): ROI name for path use.
- ee_project (str): The name of the Earth Engine project.
- meteo_elev (ee.Image): Coarse elevation used for downscaling.
- elev (ee.Image): Elevation.
- sc_probab (ee.ImageCollection): Monthly multiyear snow cover probability.
- slope (ee.Image): Slope.
- aspect (ee.Image): Aspect.
- chili (ee.Image): CHILI.
- mtpi (ee.Image): Multi-topographic position index.
- landcover (ee.Image): Landcover.
- mask (ee.ImageCollection): Mask of water/glaciers/forest/outside roi.
- initialisation_date (str): Initialisation ('%Y-%m-%d') 1 day before the start-date.
- end_date (str): End-date of the season ('%Y-%m-%d').
- START_HOUR (int): It sets the time around which a 24h period ir wrapped for the 
                    daily aggregation of meteo data. It is usually set ~10:00 since
                    that is when more satellite overpasses occur.

Internal functions:
- date_bands(): Adds bands of year, month, day of year (DOY), and date.
- terrain_bands(): Adds bands of elevation, slope, aspect, and landforms.
- meteo_bands(): Downscales precipitation and 2m air temperature via MicroMet,
                   (Liston & Elder, 2006), and adds them as bands.
- degreeday_bands(): Calculates heating & cooling degree days (HD & CD) and adds bands.
- sum_bands(): Adds placeholder cumulative snow cover, hd, and cd bands.
- flag_band(): Adds flag band (clear sky pixel vs masked pixel).

Output:
- ee.ImageCollection with metadata.
__________________________________________________________________________________________

"""
#===============================================================================
# Load libraries
#===============================================================================
import ee
import geemap
import pandas as pd

#===============================================================================
# Get Earth Engine started
#===============================================================================
ee.Authenticate()
ee.Initialize()

#===============================================================================
# Add metadata
#===============================================================================
def add_metadata(collection, domain_ee, constants, input_vars, decision_tree_settings, meteo_format, meteo, roi_path, ee_project,
                 elev, meteo_elev, sc_probab_col, slope, aspect, chili, mtpi, landcover, mask, initialisation_date, end_date, START_HOUR):

    if meteo_format == 'nc':
        # Setup Time
        t0 = ee.Date(initialisation_date).advance(START_HOUR, 'hour')
        n_days = ee.Date(end_date).advance(START_HOUR, 'hour').difference(t0, 'day')
        date_list = ee.List.sequence(0, n_days).map(lambda d: t0.advance(d, 'day').millis())
        
        # Turn a multi-band Image into a timed ImageCollection
        def to_collection(path, name):
            img = ee.Image(path)
            # Map over band names to create a list of single-band images with timestamps
            imgs = img.bandNames().zip(date_list).map(lambda pair: 
                img.select([ee.List(pair).get(0)]).rename(name).set('system:time_start', ee.List(pair).get(1))
            )
            return ee.ImageCollection.fromImages(imgs)
        
        # Apply the above & combine to a single ee.ImageCollection
        meteo_t2m_col = to_collection(
            f'projects/{ee_project}/assets/snowMapper/{roi_path}/meteo_t2m_{initialisation_date}_{end_date}',
            't2m'
        )
        meteo_precip_col = to_collection(
            f'projects/{ee_project}/assets/snowMapper/{roi_path}/meteo_precip_{initialisation_date}_{end_date}',
            'precip'
        )
        meteo = meteo_t2m_col.combine(meteo_precip_col).sort('system:time_start')
    else:
        pass

    def img_metadata(img, domain_ee, constants, input_vars, decision_tree_settings, meteo_format, meteo, roi_path,
                     elev, meteo_elev, sc_probab_col, slope, aspect, chili, mtpi, landcover, mask, START_HOUR):
        # Extract the acquisition date from the image
        acquisition_date = ee.Date(img.get('system:time_start'))
        
        # Get the acquisition date in milliseconds and add 'START_HOUR' number of hours in milliseconds (ms)
        date_millis = acquisition_date.advance(START_HOUR, 'hour').millis()

        #------------------------------------------
        # Add acquisition date info as bands
        #------------------------------------------
        def date_bands(img, acquisition_date, START_HOUR):
            # Extract year, month, and day of year (DOY)
            year = acquisition_date.get('year')
            month = acquisition_date.get('month')
            # doy = acquisition_date.getRelative('day', 'year')
            
            img = img.set('month', month)  # Required for linking collection to the monthly 'sc_probab' img.
            # img = img.set('year', year)
        
            # Convert year, month, DOY, and date in millis to images
            year_img = ee.Image.constant(year).toUint16().rename('year')           # Cast to Int16
            month_img = ee.Image.constant(month).toUint8().rename('month')         # Cast to Int8
            # doy_img = ee.Image.constant(doy).toUint16().rename('doy')            # Cast to Int16
            # date_img = ee.Image.constant(date_millis).toUint64().rename('date')  # Cast to Int64 for large range
        
            return ee.Image(
                img.addBands([
                    year_img, month_img  # , doy_img,  date_img
                ])
            )

        #------------------------------------------
        # Add terrain variables
        #------------------------------------------
        def terrain_bands(img, input_vars, elev, slope, aspect, chili, mtpi, landcover):
            bands_to_add = []
        
            for var in input_vars:
                if var in ['lat', 'lon']:
                    coords = ee.Image.pixelLonLat()
                    if var == 'lat':
                        band = coords.select(['latitude'], ['lat']).toFloat()
                    else:
                        band = coords.select(['longitude'], ['lon']).toFloat()
                else:
                    band = {'elev':elev, 
                            'slope':slope,
                            'aspect':aspect,
                            'chili':chili, 
                            'mtpi':mtpi,
                            'landcover':landcover
                           }.get(var)
        
                if band is not None:
                    bands_to_add.append(band)
        
            return img.addBands(bands_to_add)
                    
        #------------------------------------------
        # Add meteo variables
        #------------------------------------------
        def meteo_bands(img, meteo_format, meteo, elev, meteo_elev, constants, acquisition_date, START_HOUR, roi_path):         
            month = acquisition_date.get('month')

            # Get the indices of the month & define constants
            month = ee.Number(acquisition_date.get('month'))
            index = month.subtract(1)  # because lists are 0-based
            lr = ee.Number(
                ee.List(
                    constants.get('lr')
                ).get(index)
            )
            x = ee.Number(
                ee.List(
                    constants.get('x')
                ).get(index)
            )
    
            end_time = ee.Date(acquisition_date.advance(START_HOUR, 'hour').advance(1, 'hour'))  # e.g. START_HOUR + 1 the current day
            start_time = ee.Date(end_time.advance(-25, 'hour'))                                  # e.g. 24 + 1 h since last known snow cover state

            # Filter to daily & aggregate temperature using mean and precipitation using sum
            if meteo_format == 'ee':
                meteo_filtered = meteo.select(['t2m', 'precip']).filterDate(start_time, end_time)
                daily_mean_temp = ee.Image(meteo_filtered.select('t2m').mean())     # Calculate daily average
                daily_sum_precip = ee.Image(meteo_filtered.select('precip').sum())  # Calculate daily cumulative
                meteo_img = daily_mean_temp.addBands([daily_sum_precip])
            elif meteo_format == 'nc':
                meteo_img = ee.Image(meteo.select(['t2m', 'precip']).filterDate(start_time, end_time).first())
            
            # Calculate elevation difference
            elev_difference = elev.subtract(meteo_elev)

            #------------------------------------------
            # Downscale temperature via MicroMet (Liston & Elder, 2006)
            #------------------------------------------
            temp_micromet = (
                meteo_img.select(['t2m'], ['t2mHres'])
                        .subtract(ee.Image.constant(lr).multiply(elev_difference))
                        .toFloat()
            )

            #------------------------------------------
            # Downscale precipitation via MicroMet (Liston & Elder, 2006)
            #------------------------------------------
            # Calculate the correction factor for MicroMet precipitation
            correction_factor = (
                ee.Image.constant(1).add(elev_difference.multiply(x))
                  .divide(ee.Image.constant(1).subtract(elev_difference.multiply(x)))
            )
            
            # Apply the correction factor to precipitation
            precip_micromet = (
                meteo_img.select(['precip'], ['precipHres'])
                        .multiply(correction_factor)
                        .toFloat()
            )
            
            #------------------------------------------
            # Return downscaled/downsampled vars
            #------------------------------------------
            return ee.Image(
                img.addBands([
                    temp_micromet, precip_micromet
                ])
            )
            
        #------------------------------------------
        # Degree-day vars
        #------------------------------------------
        def degreeday_bands(img, decision_tree_settings):
            precip_micromet = img.select('precipHres')
            temp_micromet = img.select('t2mHres')
        
            # HD = max(0, T - 0) = max(0, T)
            hd = temp_micromet.max(0).rename('hd').toFloat()
        
            # CD = max(0, 0 - T)
            cd = ee.Image(0).subtract(temp_micromet).max(0).rename('cd').toFloat()
        
            # Precipitation when T > 0
            hd_precip = precip_micromet.where(temp_micromet.lte(decision_tree_settings["FREEZING_TEMP"]), 0).rename('hd_precip').toFloat()
            
            # Precipitation when T <= 0
            cd_precip = precip_micromet.where(temp_micromet.gt(decision_tree_settings["FREEZING_TEMP"]), 0).rename('cd_precip').toFloat()
        
            return ee.Image(
                img.addBands([
                    hd, cd, hd_precip, cd_precip
                ])
            )
        
        #------------------------------------------
        # SC, cumulative, previous day cumulative, and evaluation bands
        #------------------------------------------
        def sum_bands(img):
            return ee.Image(
                img.addBands([
                    ee.Image(2).toUint8().rename('sc'),  # Value of 2 corresponds to empty
                    ee.Image(0).toUint8().rename('tp'),
                    ee.Image(0).toUint8().rename('tn'),
                    ee.Image(0).toUint8().rename('fp'),
                    ee.Image(0).toUint8().rename('fn'),
                    ee.Image(0).toUint16().rename('sc_sum'), 
                    ee.Image(0).toFloat().rename('hd_sum'), 
                    ee.Image(0).toFloat().rename('cd_sum'), 
                    ee.Image(0).toFloat().rename('hd_precip_sum'),
                    ee.Image(0).toFloat().rename('cd_precip_sum'),
                    ee.Image(0).toFloat().rename('sc_sum_prev'), 
                    ee.Image(0).toFloat().rename('hd_sum_prev'),
                    ee.Image(0).toFloat().rename('cd_sum_prev'), 
                    ee.Image(0).toFloat().rename('hd_precip_sum_prev'),
                    ee.Image(0).toFloat().rename('cd_precip_sum_prev')
                ])
            )
            
        #------------------------------------------
        # Data flag
        #------------------------------------------
        def flag_band(img):
            # 1 where sc_obs is 2, 0 otherwise
            flag_all = img.select('sc_obs').eq(2).rename('flag_all').toUint8()
            flag_dt = img.select('sc_obs').eq(2).rename('flag_dt').toUint8()
            
            return ee.Image(
                img.addBands([
                    flag_all, flag_dt
                ])
            )

        #------------------------------------------
        # Add multi-year snow cover probability
        #------------------------------------------
        def sc_probab_band(img, sc_probab_col):        
            monthly_sc_probab = img.linkCollection(
                imageCollection = sc_probab_col,
                linkedBands = ['sc_probab'],
                linkedProperties = ['month'],
                matchPropertyName = 'month'
            )
            return ee.Image(
                monthly_sc_probab
            )

        #------------------------------------------
        # Add MODIS Terra snow cover
        #------------------------------------------
        def sc_modis_band(img, acquisition_date, decision_tree_settings, domain_ee):
            start_time = ee.Date(acquisition_date)
            end_time = ee.Date(start_time.advance(1, 'day'))
            
            modis_sc = (
                ee.ImageCollection('MODIS/061/MOD10A1')
                  .filterDate(start_time, end_time)
                  .map(lambda img: img.updateMask(img.select('NDSI_Snow_Cover_Basic_QA').lte(1)))
                  .first()
                  .resample('bicubic')
                  .select(['NDSI_Snow_Cover'], ['sc_modis']).gt(decision_tree_settings['MODIS_FSC_THRES'])
                  .toUint8()
            )

            return ee.Image(
                img.addBands([modis_sc])
            )

        #------------------------------------------
        # Apply the above
        #------------------------------------------
        img = date_bands(img, acquisition_date, START_HOUR)
        img = terrain_bands(img, input_vars, elev, slope, aspect, chili, mtpi, landcover)
        img = meteo_bands(img, meteo_format, meteo, elev, meteo_elev, constants, acquisition_date, START_HOUR, roi_path)
        img = degreeday_bands(img, decision_tree_settings)
        img = sum_bands(img)
        img = flag_band(img)
        if 'sc_probab' in input_vars:
            img = sc_probab_band(img, sc_probab_col)
        if decision_tree_settings["modis"]:
            img = sc_modis_band(img, acquisition_date, decision_tree_settings, domain_ee)
        img = img.updateMask(mask.Not())
        return ee.Image(
            img.set('system:time_start', date_millis)
        )

    #------------------------------------------
    # Apply the above
    #------------------------------------------
    return ee.ImageCollection(
        collection.map(
            lambda img: img_metadata(img, domain_ee, constants, input_vars, decision_tree_settings, meteo_format, meteo, roi_path,
                                     elev, meteo_elev, sc_probab_col, slope, aspect, chili, mtpi, landcover, mask, START_HOUR)
        )
    )