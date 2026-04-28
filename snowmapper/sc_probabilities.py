"""
_____________________________________________
name: snowMapper 
doi: 10.5281/zenodo.17663731
licence: MIT
contact: konstantis.alexopoulos@gmail.com
_____________________________________________

__________________________________________________________________________________________
Description:
Performs temporal aggregation of snow cover, in order to calculate the multiyear 
probability of snow cover in each given pixel, across the entire period of available
satellite observations from Landsat missions and Sentinel-2. Apply without running the
sm.daily_composites() which will produce a memory error due to the number of aggregations.
To run this, re-configure model timeframe to include multiple years (usually >30). Check
that all available satellite missions set as active.

IMPORTANT: If the function returns an error, it is likely due to the size of the dataset.
           Try with fewer months at a time by adapting line 43 and re-runing the script
           after re-loading the snowMapper package.
 
Input parameters:
- collection (ee.ImageCollection): Binary snow collection.
- mask (ee.ImageCollection): Mask of water/glaciers/forest/outside roi.
- months_list (list): List of months included in the configured reconstruction timeframe.

Output:
- ee.ImageCollection of monthly snow cover probabilities ('sc_probab').
__________________________________________________________________________________________

"""
#===============================================================================
# Load libraries
#===============================================================================
import ee
from datetime import datetime

#===============================================================================
# Get Earth Engine started
#===============================================================================
ee.Authenticate()
ee.Initialize()

#===============================================================================
# Aggregate into monthly collection
#===============================================================================
def sc_probabilities(collection, months_list, mask):  
        
    # ------------------------------------------
    # Extract the month from the date of each image
    # ------------------------------------------
    def add_month(image):
        date = image.date()
        return image.set({'month': date.get('month')})
    
    collection = collection.map(add_month)
    
    # ------------------------------------------
    # Group by month and reduce
    # ------------------------------------------
    def monthly_composite(month, mask):
        month = ee.Number(month)
        monthly_imgs = collection.filter(ee.Filter.eq('month', month))
        probability_img = (
            monthly_imgs.mean()
                        .select(['sc_obs'], ['sc_probab'])
                        .set({'month': month})
        )

        probability_img = probability_img.addBands(
            [ee.Image.constant(month).toInt8().rename('month')]
        )

        probability_img = probability_img.updateMask(mask.Not())
        return ee.Image(probability_img)
        
    # Create the monthly composites
    return ee.ImageCollection(
        ee.List(months_list).map(lambda month: monthly_composite(month, mask))
    ).sort('system:time_start')