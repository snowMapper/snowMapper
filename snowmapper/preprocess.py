"""
_____________________________________________
name: snowMapper 
doi: 10.5281/zenodo.17663731
licence: MIT
contact: konstantis.alexopoulos@gmail.com
_____________________________________________

__________________________________________________________________________________________
Description:
Applies a number of preprocessing algorithms necessary for the following steps of the
model. 'If' statements at the end use the booleans from the 'missions' dictionary to 
select images from one specific sensor at a time (from the ones inputed in config.yml,
to then apply the relavant preprocessing function.

Input parameters:
- domain_ee (ee.FeatureCollection): Region of interest.
- start_date (str): Start-date of the season ('%Y-%m-%d').
- end_date (str): End-date of the season ('%Y-%m-%d').
- months_list (list): List of months included in the configured reconstruction timeframe.
- missions (dict): Dictionary of booleans, defining the selected missions to be used.

Internal functions:
- preprocess_*(): Different preprocessing functions for satellite missions with similar
                  characteristics. Preprocessing functions may perform cloud-masking
                  using different algorithms, band reflectance scaling, band renaming,
                  properties calculations (e.g. solar zenith) and renaming.

Output:
- Preprocessed ee.ImageCollection.
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
# Preprocessing
#===============================================================================
def preprocess(domain_ee, start_date, end_date, months_list, missions):
    
    #------------------------------------------
    # Sentinel-2 MSI preprocessing
    #------------------------------------------
    def preprocess_S2(img):
            
        #------------------------------------------
        # Mask cloudy pixels based on QA60 band
        #------------------------------------------
        # def cloudmask_QA60(img):
        #     qa = img.select('QA60')
            
        #     # Bits 10 and 11 are clouds and cirrus, respectively.
        #     cloud_bit_mask = 1 << 10
            
        #     # Both flags should be set to zero, indicating clear conditions
        #     mask = qa.bitwiseAnd(cloud_bit_mask).eq(0)
        #     cloudmasked_img = img.updateMask(mask).select(['B2', 'B3', 'B4', 'B8', 'B11', 'B12']).copyProperties(img, ["system:time_start"])
        #     return ee.Image(cloudmasked_img)

        #------------------------------------------
        # Mask cloudy pixels based on Cloudscore+ collection
        #------------------------------------------
        def cloudmask_CS_CDF(img):
            # The threshold for masking; values between 0.50 and 0.65 generally work well.
            # Higher values will remove thin clouds, haze & cirrus shadows.
            CLEAR_THRESHOLD = 0.60
            
            # Use 'cs' or 'cs_cdf', depending on your use case; see docs for guidance.            
            cloud_free_img = img.updateMask(img.select('cs_cdf').gte(CLEAR_THRESHOLD)) \
                                .copyProperties(img, ['system:time_start', 'MEAN_SOLAR_AZIMUTH_ANGLE', 'MEAN_SOLAR_ZENITH_ANGLE'])

            return ee.Image(cloud_free_img)
        
        #------------------------------------------
        # Apply scaling factors to Sentinel-2 bands
        #------------------------------------------
        def scale_factors_S2(img):
            opticalBands = img.select(['B2', 'B3', 'B4', 'B8', 'B11', 'B12']).divide(10000).toDouble()
            scaled_img = img.addBands(opticalBands, None, True).copyProperties(img, ['system:time_start'])
            return ee.Image(scaled_img)

        #------------------------------------------
        # Rename Sentinel-2 bands
        #------------------------------------------
        def rename_bands_S2(img):
            bands = ['B2', 'B3', 'B4', 'B8', 'B11', 'B12']
            new_bands = ['blue', 'green', 'red', 'nir', 'swir1', 'swir2']
            renamed_img = img.select(bands).rename(new_bands).copyProperties(img, ['system:time_start'])
            return ee.Image(renamed_img) 
            
        #------------------------------------------
        # Add a null temperature band
        #------------------------------------------
        def add_null_temp_band(img):
            null_temp = ee.Image(0).rename('temp').updateMask(ee.Image.constant(0)).toDouble()
            return img.addBands(null_temp)
            
        #------------------------------------------
        # Rename properties to match Lansdat properties
        #------------------------------------------
        def rename_properties_S2(img):
            return (
                img.set('SPACECRAFT_ID', img.get('SPACECRAFT_NAME'))
                   .set('SUN_AZIMUTH', img.get('MEAN_SOLAR_AZIMUTH_ANGLE'))
                   .set('SUN_ZENITH', img.get('MEAN_SOLAR_ZENITH_ANGLE'))
            )

        #------------------------------------------
        # Apply the above
        #------------------------------------------
        img = cloudmask_CS_CDF(img)  # Options: cloudmask_CS_CDF(img), cloudmask_QA60(img)
        img = scale_factors_S2(img)
        img = rename_bands_S2(img)
        img = add_null_temp_band(img)
        img = rename_properties_S2(img)
        return ee.Image(img)

    #------------------------------------------
    # Landsat 8-9 OLI/TIRS preprocessing
    #------------------------------------------
    def preprocess_L89(img):       

        #------------------------------------------
        # Mask cloudy pixels based on QA60 band
        #------------------------------------------
        def cloudmask_QA_PIXEL(img):
            qa = img.select('QA_PIXEL') 
            cloudBitMask = 1 << 3
            cloudShadowBitMask = 1 << 4
            snowBitMask = 1 << 5
            
            cloudMask = qa.bitwiseAnd(cloudBitMask).eq(0)              # Cloud-free pixels
            cloudShadowMask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)  # Cloud shadow-free pixels
            snowMask = qa.bitwiseAnd(snowBitMask).neq(0)               # Snow-covered pixels

            # Mask
            mask = cloudMask  # .And(cloudShadowMask) .Or(snowMask)          
            
            cloudmasked_img = (
                img.updateMask(mask)
                   .select(['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7', 'ST_B10'])
                   .copyProperties(img, ['system:time_start'])
            )
            return ee.Image(cloudmasked_img)
            
        #------------------------------------------
        # Apply scaling factors to Landsat 8-9 image bands
        #------------------------------------------
        def scale_factors_L89(img):
            optical_bands = img.select(['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7']).multiply(0.0000275).add(-0.2).toDouble()
            thermal_band = img.select('ST_B10').multiply(0.00341802).add(149.0).toDouble()
            scaled_img = (
                img.addBands(optical_bands, None, True)
                   .addBands(thermal_band, None, True)
                   .copyProperties(img, ['system:time_start'])
            )
            return ee.Image(scaled_img)

        #------------------------------------------
        # Rename Landsat 8-9 image bands
        #------------------------------------------
        def rename_bands_L89(img):
            bands = ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7', 'ST_B10']
            new_bands = ['blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'temp']
            renamed_img = img.select(bands).rename(new_bands).copyProperties(img, ['system:time_start'])
            return ee.Image(renamed_img)

        #------------------------------------------
        # Calculate solar zenith
        #------------------------------------------
        def set_zenith_L89(img):
            sun_elev = ee.Number(img.get('SUN_ELEVATION'))
            sun_zenith = ee.Number(90).subtract(sun_elev)
            return img.set('SUN_ZENITH', sun_zenith)
            
        #------------------------------------------
        # Apply the above
        #------------------------------------------
        img = cloudmask_QA_PIXEL(img)
        img = scale_factors_L89(img)
        img = rename_bands_L89(img)
        img = set_zenith_L89(img)
        return ee.Image(img)

    #------------------------------------------
    # Landsat 7 ETM+ & 4-5 TM preprocessing
    #------------------------------------------
    def preprocess_L457(img):          
        
        #------------------------------------------
        # Mask cloudy pixels based on SR_CLOUD_QA band
        #------------------------------------------
        # def cloudmask_SR_CLOUD_QA(img):
        #     qa = img.select('SR_CLOUD_QA')
        #     cloudBitMask = 1 << 1  # Bit position for cloud flag
        
        #     # Both flags should be set to zero, indicating clear conditions
        #     mask = qa.bitwiseAnd(cloudBitMask).eq(0)
        
        #     cloudmasked_img = (
        #         img.updateMask(mask)
        #            .select(['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7', 'ST_B6'])
        #            .copyProperties(img, ['system:time_start'])
        #     )
        #     return ee.Image(cloudmasked_img)

        #------------------------------------------
        # Mask cloudy pixels based on QA60 band
        #------------------------------------------
        def cloudmask_QA_PIXEL(img):
            qa = img.select('QA_PIXEL') 
            cloudBitMask = 1 << 3
            cloudShadowBitMask = 1 << 4
            snowBitMask = 1 << 5
            
            cloudMask = qa.bitwiseAnd(cloudBitMask).eq(0)              # Cloud-free pixels
            cloudShadowMask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)  # Cloud shadow-free pixels
            snowMask = qa.bitwiseAnd(snowBitMask).neq(0)               # Snow-covered pixels

            # Mask
            mask = cloudMask  # .And(cloudShadowMask) .Or(snowMask)          
            
            cloudmasked_img = (
                img.updateMask(mask)
                   .select(['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7', 'ST_B6'])
                   .copyProperties(img, ['system:time_start'])
            )
            return ee.Image(cloudmasked_img)
        
        #------------------------------------------
        # Apply scaling factors to Landsat 4-5-7 image bands
        #------------------------------------------
        def scale_factors_L457(img):
            optical_bands = img.select(['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7']).multiply(0.0000275).add(-0.2).toDouble()
            thermal_band = img.select('ST_B6').multiply(0.00341802).add(149.0).toDouble()
            scaled_img = (
                img.addBands(optical_bands, None, True)
                   .addBands(thermal_band, None, True)
                   .copyProperties(img, ['system:time_start'])
            )
            return ee.Image(scaled_img)

        #------------------------------------------
        # Rename Landsat 4-5-7 image bands
        #------------------------------------------
        def rename_bands_L457(img):
            bands = ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7', 'ST_B6']
            new_bands = ['blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'temp']
            renamed_img = img.select(bands).rename(new_bands).copyProperties(img, ['system:time_start'])
            return ee.Image(renamed_img)

        #------------------------------------------
        # Calculate solar zenith
        #------------------------------------------
        def set_zenith_L457(img):
            sun_elev = ee.Number(img.get('SUN_ELEVATION'))
            sun_zenith = ee.Number(90).subtract(sun_elev)
            return img.set('SUN_ZENITH', sun_zenith)
            
        #------------------------------------------
        # Apply the above
        #------------------------------------------
        img = cloudmask_QA_PIXEL(img)  # Options: cloudmask_QA_PIXEL(img), cloudmask_SR_CLOUD_QA(img)
        img = scale_factors_L457(img)
        img = rename_bands_L457(img)
        img = set_zenith_L457(img)
        return ee.Image(img)

    #------------------------------------------
    # Run based on Spacecraft ID
    #------------------------------------------
    raw_col = ee.ImageCollection([])
    
    start_month = months_list[0]
    end_month = months_list[-1]
    
    if missions["Sentinel_2"]:
        csPlus = (
            ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED')
              .filterDate(start_date, end_date)
              .filter(ee.Filter.calendarRange(start_month, end_month, 'month'))
        )
        
        S2 = (
            ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
              .filterBounds(domain_ee)
              .filterDate(start_date, end_date)
              .filter(ee.Filter.calendarRange(start_month, end_month, 'month'))
              .linkCollection(csPlus, ['cs_cdf'])
              .map(preprocess_S2)
        )
        raw_col = raw_col.merge(S2)
    
    if missions["Landsat_9"]:
        L9 = (
            ee.ImageCollection("LANDSAT/LC09/C02/T1_L2")
              .filterBounds(domain_ee)
              .filterDate(start_date, end_date)
              .filter(ee.Filter.calendarRange(start_month, end_month, 'month'))
              .map(preprocess_L89)
        )
        raw_col = raw_col.merge(L9)               
    
    if missions["Landsat_8"]:
        L8 = (
            ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
              .filterBounds(domain_ee)
              .filterDate(start_date, end_date)
              .filter(ee.Filter.calendarRange(start_month, end_month, 'month'))
              .map(preprocess_L89)
        )
        raw_col = raw_col.merge(L8)
    
    if missions["Landsat_7"]:
        L7 = (
            ee.ImageCollection("LANDSAT/LE07/C02/T1_L2")
              .filterBounds(domain_ee)
              .filterDate(start_date, end_date)
              .filter(ee.Filter.calendarRange(start_month, end_month, 'month'))
              .map(preprocess_L457)
        )
        raw_col = raw_col.merge(L7)
    
    if missions["Landsat_5"]:
        L5 = (
            ee.ImageCollection("LANDSAT/LT05/C02/T1_L2")
              .filterBounds(domain_ee)
              .filterDate(start_date, end_date)
              .filter(ee.Filter.calendarRange(start_month, end_month, 'month'))
              .map(preprocess_L457)
        )
        raw_col = raw_col.merge(L5)
    
    if missions["Landsat_4"]:
        L4 = (
            ee.ImageCollection("LANDSAT/LT04/C02/T1_L2")
              .filterBounds(domain_ee)
              .filterDate(start_date, end_date)
              .filter(ee.Filter.calendarRange(start_month, end_month, 'month'))
              .map(preprocess_L457)
        )
        raw_col = raw_col.merge(L4)
                
    return ee.ImageCollection(raw_col)