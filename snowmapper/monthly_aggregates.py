"""
_____________________________________________
name: snowMapper 
doi: 10.5281/zenodo.17663731
licence: MIT
contact: konstantis.alexopoulos@gmail.com
_____________________________________________

__________________________________________________________________________________________
Description:
Performs temporal aggregation of results, using the mean value of a pixel across a of its
instances (daily) during a given month. This is being performed separately for every month 
of a list, defined using the start and end date of the season, and skipping the months
that formed the spin-up period.
 
Input parameters:
- reconstructed_col (ee.ImageCollection): Final reconstructed collection.
- meta_col (ee.ImageCollection): Metadata collection from which to source downscaled
                                 climate variables.
- months_list (list): List of months of in the configuration timeframe.
- SPIN_UP_PERIOD (int): Number of months to be discarted from the start of the season.

Output:
- ee.ImageCollection of monthly means of snow cover and other metadata bands.
__________________________________________________________________________________________

"""
#===============================================================================
# Load libraries
#===============================================================================
import ee
import pandas as pd

#===============================================================================
# Get Earth Engine started
#===============================================================================
ee.Authenticate()
ee.Initialize()

#===============================================================================
# Aggregate into monthly collection
#===============================================================================
def monthly_aggregates(reconstructed_col, meta_col, months_list, SPIN_UP_PERIOD):
    
    months_list = months_list[SPIN_UP_PERIOD:]
    months_list = ee.List(months_list)

    #------------------------------------------
    # Link reconstructed & meta cols
    #------------------------------------------
    collection = reconstructed_col.linkCollection(
        imageCollection = meta_col,
        linkedBands = ['t2mHres', 'precipHres'],
        linkedProperties = ['system:time_start'],
        matchPropertyName = 'system:time_start'
    )
    
    #------------------------------------------
    # Extract the year and month from the date of each image
    #------------------------------------------
    def add_month(image):
        date = image.date()
        return image.set({
            'year': date.get('year'),
            'month': date.get('month')
        })
    
    collection = collection.map(add_month)
    
    #------------------------------------------
    # Group by year and month, then reduce each group
    #------------------------------------------
    def monthly_composite(month):
        month = ee.Number(month)
        monthly_imgs = collection.filter(ee.Filter.eq('month', month))
        year = ee.Number(monthly_imgs.first().get('year'))
        output_date = ee.Date.fromYMD(year, month, 1)

        # Combine reducers with prefixes to distinguish output bands
        combined_reducer = ee.Reducer.mean().combine(
            reducer2=ee.Reducer.sum(),
            sharedInputs=True
        )

        # Reduce image collection
        reduced_img = monthly_imgs.reduce(combined_reducer) \
                                      .set({
                                          'month': month, 
                                          'year': year, 
                                          'system:time_start': output_date.millis()
                                      })
        reduced_img = reduced_img.select(
            ['sc_mean', 'flag_all_mean', 'flag_dt_mean', 't2mHres_mean', 'precipHres_sum', 'year_mean', 'month_mean', 'tp_sum', 'tn_sum', 'fp_sum', 'fn_sum'],
            ['sc_mean', 'flag_all_mean', 'flag_dt_mean', 't2mHres_mean', 'precipHres_sum', 'year', 'month', 'tp_sum', 'tn_sum', 'fp_sum', 'fn_sum']
        )  
        
        return ee.Image(reduced_img)

    # Create monthly composites by applying the monthly_composite function
    monthly_col = ee.ImageCollection(months_list.map(monthly_composite))

    # Mask pixles that were not reconstructed. This will happen if the gridded 
    # climate dataset used is masked to land, and the domain is close to sea.
    monthly_col = monthly_col.map(lambda img: img.updateMask(img.select('sc_mean').lte(1)))
    
    return monthly_col