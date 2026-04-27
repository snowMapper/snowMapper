"""
_____________________________________________
name: snowMapper 
doi: 10.5281/zenodo.17663731
_____________________________________________

__________________________________________________________________________________________
Description:
Unmasks pixels masked during preprocessing (does not include roi preprocessing - these
remain permanently masked), and also creates synthetic images for days where no images
exist at all. This is to create a base upon which daily metadata can be added, for the
gapfilling and reconstruction algorithms to be then performed.

Input parameters:
- collection (ee.ImageCollection): Preprocessed image collection.
- domain_ee (ee.FeatureCollection): Region of interest (e.g. mountain range).
- start_date (ee.Date): Start-date of the season.
- end_date (ee.Date): End-date of the season.

Internal functions:
- synthesise(): Creates synthetic images for days within the timeframe when no images 
                were captured from any of the existing satellite sensors.
- unmask_null_pixels(): Unmasks the pixels of daily images, using a value of '2' to
                        represent pixels with no snow cover value ('0' or '1').

Output:
- ee.ImageCollection of daily pictures and no masked pixels.
__________________________________________________________________________________________

"""
#===============================================================================
# Load libraries
#===============================================================================
import ee
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta

#===============================================================================
# Get Earth Engine started
#===============================================================================
ee.Authenticate()
ee.Initialize()

#===============================================================================
# Image synthesis
#===============================================================================
def img_synth(collection, domain_ee, start_date, end_date):
    
    #------------------------------------------
    # Prerequisites for synthesise()
    #------------------------------------------
    # Recalculate start_date based on the output of set_first (e.g. 31 August)
    start_date = (datetime.strptime(start_date, "%Y-%m-%d") + relativedelta(months=1) - timedelta(days=1)).strftime("%Y-%m-%d")
    
    # Ensure start_date is earlier than end_date by comparing ee.Date objects
    ordered_dates = sorted([start_date, end_date])
    start_date_ee = ee.Date(ordered_dates[0])
    end_date_ee = ee.Date(ordered_dates[1])
    
    # Get the number of days in the date range
    num_days = end_date_ee.difference(start_date_ee, 'day')
    
    # Generate a list of all dates in the range
    date_list = ee.List.sequence(0, num_days.subtract(1)).map(
        lambda day: start_date_ee.advance(day, 'day')
    )
    
    # Get the dates available in the processed collection
    existing_dates = collection.aggregate_array('system:time_start').map(
        lambda time: ee.Date(time)
    )
    
    # Find missing dates
    missing_dates = date_list.filter(ee.Filter.inList('item', existing_dates).Not())
    
    #------------------------------------------
    # Create synthetic images for missing dates
    #------------------------------------------
    def synthesise(date):
        return ee.Image(
            ee.Image(0) \
            .rename('sc_obs') \
            .updateMask(ee.Image.constant(0)) \
            .set('system:time_start', ee.Date(date).millis())
        )
    
    #------------------------------------------
    # Unmask and fill missing `sc` pixels
    #------------------------------------------
    def unmask_null_pixels(img):
        return ee.Image(
            img.select('sc_obs') \
            .unmask(2) \
            .set('system:time_start', img.get('system:time_start'))
        )
        
    #------------------------------------------
    # Apply the above
    #------------------------------------------
    synthetic_imgs = ee.ImageCollection(missing_dates.map(synthesise))
    synthetic_col = collection.merge(synthetic_imgs)
    unmasked_col = synthetic_col.map(unmask_null_pixels)
    sorted_col = unmasked_col.sort('system:time_start')
    return ee.ImageCollection(sorted_col)