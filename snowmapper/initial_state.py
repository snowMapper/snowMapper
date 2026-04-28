"""
_____________________________________________
name: snowMapper 
doi: 10.5281/zenodo.17663731
_____________________________________________

__________________________________________________________________________________________
Description: 
Gap-fills data on the 'initial state' image. All gaps are filled with no-snow values.

Input parameters:
- collection (ee.ImageCollection): Image collection with metadata.
- domain_ee (ee.Feature): Region of interest (e.g. mountain range).
- initialisation_date (str): Initialisation ('%Y-%m-%d') 1 day before the start-date.
- start_date (str): Start-date of the season ('%Y-%m-%d').

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
def initial_state(collection, domain_ee, initialisation_date, start_date):
    
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
    img = collection.first()  # Initialisation image
    post_initialisation_col = ee.ImageCollection(collection.filter(ee.Filter.date(initialisation_date, start_date).Not()))
    updated_img = obs_img(img)
    updated_img = gapfill_img(updated_img)
    updated_img = starting_vars(updated_img)
    updated_col = post_initialisation_col.merge(ee.ImageCollection([updated_img]))
    
    return ee.ImageCollection(updated_col).sort('system:time_start')