"""
_____________________________________________
name: snowMapper 
doi: 10.5281/zenodo.17663731
_____________________________________________

__________________________________________________________________________________________
Description:
Function to save an ee.ImageCollection as an asset.

Input parameters:
- collection (ee.ImageCollection): Image collection at any stage.
- ee_project (str): The name of the Earth Engine project.
- name (str): 'meta', 'reconstructed', 'monthly'.
- roi_path (str): ROI name for path use.
- domain_ee (ee.FeatureCollection): Region of interest (e.g. mountain range).
- CRS (str): Desired CRS of grid.
- SCALE (int): Desired scale of grid.
- CRS_TRANSFORM (str): CRS transform for that domain.

Internal functions:
-

Output:
- ee.ImageCollection asset saved.
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
# Reconstruct all binary snow images
#===============================================================================
def export_ee(collection, ee_project, name, roi_path, START_YEAR, END_YEAR, domain_ee, CRS, SCALE, CRS_TRANSFORM):
    path = f"projects/{ee_project}/assets/snowMapper/{roi_path}/{name}_{roi_path}_{START_YEAR}_{END_YEAR}"
    image_ids = collection.aggregate_array('system:index').getInfo()
    
    # Create ImageCollection asset
    ee.data.createAsset({'type': 'ImageCollection'}, path)

    tasks = []
    for image_id in image_ids:
        image = ee.Image(collection.filter(ee.Filter.eq("system:index", image_id)).first())
        
        task = ee.batch.Export.image.toAsset(
            image=image,
            description=f"Export_{image_id}",
            assetId=f"{path}/{image_id}",
            region=domain_ee,
            scale=SCALE,
            crs=CRS,
            crsTransform=CRS_TRANSFORM,
            maxPixels=1e10
        )
        task.start()
        tasks.append(task)    
    return tasks