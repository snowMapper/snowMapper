"""
_____________________________________________
name: snowMapper 
doi: 10.5281/zenodo.17663731
_____________________________________________

__________________________________________________________________________________________
Description:
Using a predefined start month, here we gather all the images from that month (usually
August due to limited cloud cover & minimum/no snow cover) and create a median composite
of the preprocessed spectral bands. This is used as an 'initial state' for later processes
of the model. The date of the composite image is then also modified to represent the last 
date of the start month, after which we will eventually begin reconstruction.

Input parameters:
- collection (ee.ImageCollection): Preprocessed image collection.
- START_MONTH (int): The first year.
- START_MONTH (int): The first month of the preprocessed image collection, which will be
                     composited and used as an initial state for subsequent snow cover
                     gap-filling and reconstructing routines.

Output:
- ee.ImageCollection where the first image is a median composite of all images within the 
  starting month (initial state).
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
# Filter collection for the month of August and create composite
#===============================================================================
def set_first(collection, START_MONTH, START_YEAR):
    month_col = collection.filter(ee.Filter.calendarRange(START_MONTH, START_MONTH, 'month'))
    
    start_date = ee.Date.fromYMD(START_YEAR, START_MONTH, 1)
    end_date = start_date.advance(1, 'month').advance(-1, 'day')  # Last day of the month
    
    empty_img = ee.Image(0).toDouble().updateMask(ee.Image(0))
    month_img = ee.Image.cat([
        empty_img.rename('blue'),
        empty_img.rename('green'),
        empty_img.rename('red'),
        empty_img.rename('nir'),
        empty_img.rename('swir1'),
        empty_img.rename('swir2'),
        empty_img.rename('temp')
    ]).set('system:time_start', end_date.millis())
    
    non_month_col = collection.filter(ee.Filter.calendarRange(START_MONTH, START_MONTH, 'month').Not())
    updated_col = non_month_col.merge(ee.ImageCollection([month_img]))

    return ee.ImageCollection(updated_col).sort('system:time_start')