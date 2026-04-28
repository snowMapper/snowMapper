"""
_____________________________________________
name: snowMapper 
doi: 10.5281/zenodo.17663731
licence: MIT
contact: konstantis.alexopoulos@gmail.com
_____________________________________________

__________________________________________________________________________________________
Description:
This function creates a composite of all available images that fall within the region of
interest within a day, using the median of each band across the available images.

Input parameters:
- collection (ee.ImageCollection): Earth Engine image collection of preprocessed 
                                   satellite images.

Output:
- ee.ImageCollection of daily composites from clear-sky pixels, where each day is the max
  composite of all available images for that day.
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
# Create daily composites
#===============================================================================
def daily_composites(collection):
    # Convert system:time_start to dates with no time component
    dates = (
        collection.aggregate_array('system:time_start')
                  .map(lambda t: ee.Date(t).format('YYYY-MM-dd'))
                  .distinct()
    )
    
    # Function to create a composite image for a specific date
    def daily_composite(date_str):
        date = ee.Date.parse('YYYY-MM-dd', date_str)
        daily_imgs = collection.filterDate(date, date.advance(1, 'day'))
        
        # Use median across images of the same day
        composite_img = daily_imgs.median().set({
            'system:time_start': date.millis(),
            'SUN_ZENITH': daily_imgs.aggregate_mean('SUN_ZENITH'),
            'SUN_AZIMUTH': daily_imgs.aggregate_mean('SUN_AZIMUTH')
        })
        
        return composite_img
    
    # Map over the dates to compute daily averages
    composites_col = ee.ImageCollection(dates.map(daily_composite, dropNulls = True))
    return composites_col.sort('system:time_start')