"""
_____________________________________________
name: snowMapper 
doi: 10.5281/zenodo.17663731
_____________________________________________

__________________________________________________________________________________________
Description: 
Gap-fills data on the first 'initial state' image. Currently, the random forest option is 
deactivated, and all gaps are filled with no-snow values.

Input parameters:
- collection (ee.ImageCollection): Image collection with metadata.
- domain_ee (ee.Feature): Region of interest (e.g. mountain range).
- START_MONTH (int): The first month of the preprocessed image collection, which will be
                     composited and used as an initial state for subsequent snow cover
                     gap-filling and reconstructing routines.

Internal functions:
- obs_img(): Replace the empty 'sc' with the 'sc_obs'.
- gapfill_img(): Passive gapfill (i.e. fills all gaps with 'no-snow' values).
- starting_vars(): Sets the starting point for cumulative variables.

Output:
- ee.ImageCollection with fully gap-filled initial state (first ee.Image).
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
# Gap-fill first binary snow image to use for spin-up
#===============================================================================
def initial_state(collection, domain_ee, START_MONTH):

    img = collection.first()  # This is the initialisation image
    
    #------------------------------------------
    # Replace 'sc' band with available 'sc_obs' band
    #------------------------------------------
    def obs_img(img):
        sc_band = img.select('sc_obs').rename('sc')
        updated_img = img.addBands(sc_band, overwrite=True)

        return ee.Image(
            updated_img.copyProperties(img, ['system:time_start'])
        )
    
    #------------------------------------------
    # Spatial Gapfill
    #------------------------------------------
    def gapfill_img(img):
        sc = img.select('sc')
        sc_filled = sc.where(sc.eq(2), 0)
        updated_img = img.addBands(sc_filled.rename('sc'), overwrite=True)
        
        return ee.Image(
            updated_img.copyProperties(img, ['system:time_start'])
        )
    
    #------------------------------------------
    # Add starting cumulative snow cover value (=sc)
    #------------------------------------------
    def starting_vars(img):
        sc = img.select('sc')
        hd = img.select('hd')
        cd = img.select('cd')
        hd_precip = img.select('hd_precip')
        cd_precip = img.select('cd_precip')

        # Compute new bands
        sc_sum = sc.rename('sc_sum')
        hd_sum = hd.multiply(sc).rename('hd_sum')
        cd_sum = cd.multiply(sc).rename('cd_sum')
        hd_precip_sum = hd_precip.multiply(sc).rename('hd_precip_sum')
        cd_precip_sum = cd_precip.multiply(sc).rename('cd_precip_sum')        
        
        return ee.Image(
            img.addBands([sc_sum, hd_sum, cd_sum, hd_precip_sum, cd_precip_sum], overwrite=True)
        ).copyProperties(img, ['system:time_start'])
    
    #------------------------------------------
    # Apply the above
    #------------------------------------------
    filtered_col = collection.filter(ee.Filter.calendarRange(START_MONTH, START_MONTH, 'month').Not())

    updated_img = obs_img(img)
    gapfilled_img = gapfill_img(updated_img)
    gapfilled_img = starting_vars(gapfilled_img)

    updated_col = filtered_col.merge(ee.ImageCollection([gapfilled_img]))
    
    return ee.ImageCollection(updated_col).sort('system:time_start')