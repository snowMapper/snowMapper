"""
_____________________________________________
name: snowMapper 
doi: 10.5281/zenodo.17663731
licence: MIT
contact: konstantis.alexopoulos@gmail.com
_____________________________________________

__________________________________________________________________________________________
Description:
Spatially aggregates results, using the mean value of pixels falling within a sub-region
of interest, to derive one representative value of snow cover (as well as other metabands)
for that area, on a given reconstructed image.

Input parameters:
- collection (ee.ImageCollection): Final reconstructed and temporally aggregated collection.
- subdomain_ee (ee.FeatureCollection): Sub-region of interest (e.g. massif instead of 
                                       mountain range) in which spatial aggregation is
                                       intended.
- SCALE (int): Desired scale of grid.
- CRS (str): Desired CRS of grid.

Output:
- ee.FeatureCollection of polygon geometries with a mean value of snow cover and other
  metabands, for a given image (date).
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
# Aggregate (reduce) data over the subdomain_ee (e.g. massif)
#===============================================================================
def subdomain_reduction(collection, subdomain_ee, SCALE, CRS):
    # Combine reducers with prefixes
    combined_reducer = ee.Reducer.mean().combine(
        reducer2=ee.Reducer.sum(),
        sharedInputs=True
    )
    
    # Reduce each image over the subdomain
    def img_reduction(img):
        return img.reduceRegions(
            collection=subdomain_ee,
            reducer=combined_reducer,
            tileScale=16,
            scale=SCALE,
            crs=CRS
        )
    
    # Map over collection and flatten the result
    reduced = collection.map(img_reduction)
    reduced_flat = reduced.flatten()
    
    # Optionally, select only relevant properties
    return ee.FeatureCollection(reduced_flat)