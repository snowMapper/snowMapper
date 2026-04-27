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
- constants (float): Constants for MicroMet climate downscaling (Liston & Elder, 2006).
- clim (ee.ImageCollection): Hourly climate data.
- clim_elev (ee.Image): Coarse elevation used for downscaling.
- elev (ee.Image): Elevation.
- sc_probab (ee.ImageCollection): Monthly multiyear snow cover probability.
- slope (ee.Image): Slope.
- aspect (ee.Image): Aspect.
- chili (ee.Image): CHILI.
- mtpi (ee.Image): Multi-topographic position index.
- landcover (ee.Image): Landcover.
- mask (ee.ImageCollection): Mask of water/glaciers/forest/outside roi.
- decision_tree_settings (dict): Contains thresholds for temperature-based 
                                 differentiation between rain and snow, and for
                                 binarising downsampled fractional snow cover
                                 from MODIS Terra NDSI.
- START_HOUR (int): It sets the time around which a 24h period ir wrapped for the 
                    daily aggregation of climate data. It is usually set ~10:00 since
                    that is when more satellite overpasses occur.

Internal functions:
- date_bands(): Adds bands of year, month, day of year (DOY), and date.
- terrain_bands(): Adds bands of elevation, slope, aspect, and landforms.
- climate_bands(): Downscales precipitation and 2m air temperature via MicroMet,
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

#===============================================================================
# Get Earth Engine started
#===============================================================================
ee.Authenticate()
ee.Initialize()

#===============================================================================
# Add metadata
#===============================================================================
def add_metadata(collection, domain_ee, constants, decision_tree_settings, clim, elev, clim_elev, 
                 sc_probab_col, slope, aspect, chili, mtpi, landcover, mask, START_HOUR):

    def img_metadata(img, domain_ee, constants, decision_tree_settings, clim, elev, clim_elev, 
                     sc_probab_col, slope, aspect, chili, mtpi, landcover, mask, START_HOUR):
        # Extract the acquisition date from the image
        acquisition_date = ee.Date(img.get('system:time_start'))
        
        # Get the acquisition date in milliseconds and add 'START_HOUR' number of hours in milliseconds (ms)
        date_millis = acquisition_date.millis().add(START_HOUR * 60 * 60 * 1000)

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
        def terrain_bands(img, elev, slope, aspect, chili, mtpi, landcover):
            coords = ee.Image.pixelLonLat()
            # lat = coords.select(['latitude'], ['lat']).toFloat()
            # lon = coords.select(['longitude'], ['lon']).toFloat()
            return ee.Image(
                img.addBands([
                    chili, mtpi  #, elev, slope, aspect, landcover, lat, lon
                ])
            )

        #------------------------------------------
        # Add climate variables
        #------------------------------------------
        def climate_bands(img, clim, elev, clim_elev, constants, acquisition_date, START_HOUR):         
            month = acquisition_date.get('month')
    
            # Get the indices of the month
            month_list = ee.List(constants.get('month'))
            lr_list = ee.List(constants.get('lr'))
            x_list = ee.List(constants.get('x'))
                                
            # Define constants
            index = month_list.indexOf(month)
            lr = ee.Number(lr_list.get(index))
            x = ee.Number(x_list.get(index))
    
            # Aggregate daily & filter
            end_time = ee.Date(acquisition_date.advance(START_HOUR, 'hour'))  # e.g. 10:00 of the current day
            start_time = ee.Date(end_time.advance(-24, 'hour'))               # e.g. 24h since last known snow cover state
            
            clim_filtered = clim.select([
                't2m', 'precip'  # , 'fsc'
            ]).filterDate(start_time, end_time)
    
            # Aggregate temperature using mean and precipitation using sum
            daily_mean_temp = ee.Image(clim_filtered.select('t2m').mean()).subtract(273.15)   # Calculate daily average & convert from K to degC
            daily_sum_precip = ee.Image(clim_filtered.select('precip').sum()).multiply(1000)  # Calculate daily cumulative & convert from m to mm
            # daily_mean_fsc = ee.Image(clim_filtered.select('fsc').mean()).divide(100)       # Calculate daily mean & convert percentage to fraction
            clim_img = daily_mean_temp.addBands([daily_sum_precip])  # , daily_mean_fsc
            
            # Calculate elevation difference
            elev_difference = elev.subtract(clim_elev)

            #------------------------------------------
            # Downscale temperature via MicroMet (Liston & Elder, 2006)
            #------------------------------------------
            temp_micromet = (
                clim_img.select(['t2m'], ['t2mHres'])
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
                clim_img.select(['precip'], ['precipHres'])
                        .multiply(correction_factor)
                        .toFloat()
            )
            
            #------------------------------------------
            # Downsample fractional snow cover
            #------------------------------------------
            # fsc_downsampled = (
            #     clim_img.select(['fsc'], ['fscLres'])
            #             .toFloat()
            # )
            
            #------------------------------------------
            # Return downscaled/downsampled vars
            #------------------------------------------
            return ee.Image(
                img.addBands([
                    temp_micromet, precip_micromet  # , fsc_downsampled
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
        img = terrain_bands(img, elev, slope, aspect, chili, mtpi, landcover)
        img = climate_bands(img, clim, elev, clim_elev, constants, acquisition_date, START_HOUR)
        img = degreeday_bands(img, decision_tree_settings)
        img = sum_bands(img)
        img = flag_band(img)
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
            lambda img: img_metadata(img, domain_ee, constants, decision_tree_settings, clim, elev, clim_elev, 
                                     sc_probab_col, slope, aspect, chili, mtpi, landcover, mask, START_HOUR)
        )
    )