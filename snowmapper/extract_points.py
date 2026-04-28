"""
_____________________________________________
name: snowMapper 
doi: 10.5281/zenodo.17663731
licence: MIT
contact: konstantis.alexopoulos@gmail.com
_____________________________________________

__________________________________________________________________________________________
Description:
Extracts pixel values from the bands of a reconstructed image. These could be used for
point-scale analysis or validation purposes.

Input parameters:
- collection (ee.ImageCollection): Final reconstructed collection.
- points_ee (ee.FeatureCollection): Earth Engine feature collection of points where
                                        in-situ measurements exist.
- SCALE (int): Desired scale of sampling for validation.
- CRS (str): Desired CRS of grid.

Output:
- ee.FeatureCollection of points with the bands and values of the image.
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
# Extract data points for validation
#===============================================================================
def extract_points(collection, points_ee, SCALE, CRS):

    # Define the function to sample points and map it over the collection
    def img_sampling(img, points_ee, SCALE, CRS):
        return img.sampleRegions(
            collection=points_ee,
            scale=SCALE,
            projection=CRS,
            # properties=[],
            tileScale=4,
            geometries=False
        )
    
    # Map the function and flatten the resulting collection of collections
    sampled_ee = collection.map(lambda img: img_sampling(img, points_ee, SCALE, CRS)).flatten()

    return ee.FeatureCollection(sampled_ee)