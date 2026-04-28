'''
_____________________________________________
name: snowMapper 
doi: 10.5281/zenodo.17663731
_____________________________________________

__________________________________________________________________________________________
Description: 
Creates a binary image of water, glacier, and forest-covered areas, which is later used 
to mask the snow maps.

Input parameters:
- masks (dict): Dictionary of booleans, defining the selected masks to be used.
- domain_ee (ee.FeatureCollection): Region of interest (e.g. mountain range).
- CRS (str): Desired CRS of grid.
- SCALE (int): Desired scale of grid.

Output:
- ee.Image (binary) with '1' set for mask pixels.
__________________________________________________________________________________________

'''
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
# Create mask image
#===============================================================================
def create_mask(masks, domain_ee, CRS, SCALE):
    mask_imgs = []

    domain_img = ee.Image().byte().paint(featureCollection=domain_ee, color=0).unmask(1).rename('mask')
    mask_imgs.append(domain_img)
    
    if masks['elevation']:
        elevation_img = ee.Image('USGS/SRTMGL1_003').select('elevation')
        elevation_mask = elevation_img.lt(masks['elevation']).rename('mask')
        mask_imgs.append(elevation_mask)

    if masks['glaciers'] == 'glims':
        glaciers = ee.FeatureCollection('GLIMS/20230607').filterBounds(domain_ee)
        glacier_img = ee.Image().byte().paint(featureCollection=glaciers, color=1).unmask(0).rename('mask')
        mask_imgs.append(glacier_img)
             
    if masks['water'] == 'jrc_global_surface_water':
        water_img = ee.Image('JRC/GSW1_4/GlobalSurfaceWater').select('max_extent')
        water_mask = water_img.eq(1).byte().rename('mask')
        mask_imgs.append(water_mask)
    
    if masks['water'] == 'esa_worldcover':
        landcover_img = ee.ImageCollection('ESA/WorldCover/v100').first().select('Map')
        water_mask = landcover_img.eq(80).byte().rename('mask')
        mask_imgs.append(water_mask)
    
    if masks['forest'] == 'copernicus_tree_cover_density':
        tree_cover_density_img = ee.Image('projects/snowmapper/assets/tree_cover_density_2015_3035_100m').select('b1')
        forest_mask = tree_cover_density_img.updateMask(tree_cover_density_img.neq(255)).gt(50).byte().rename('mask')
        mask_imgs.append(forest_mask)

    if masks['forest'] == 'esa_worldcover':
        landcover_img = ee.ImageCollection('ESA/WorldCover/v100').first().select('Map')
        forest_mask = landcover_img.eq(10).byte().rename('mask')
        mask_imgs.append(forest_mask)

    if masks['forest'] == 'jrc_global_forest_cover':
        forest_img = ee.ImageCollection('JRC/GFC2020/V2').filterBounds(domain_ee).mosaic().select('Map')
        forest_mask = forest_img.eq(1).byte().rename('mask')
        mask_imgs.append(forest_mask)

    if masks['urban'] == 'esa_worldcover':
        landcover_img = ee.ImageCollection('ESA/WorldCover/v100').first().select('Map')
        urban_mask = landcover_img.eq(50).byte().rename('mask')
        mask_imgs.append(urban_mask)
    
    combined_mask = (
        ee.ImageCollection(mask_imgs)
          .sum()
          .gt(0)
          .unmask(0)
          .clip(
              domain_ee.bounds(maxError=1, proj=CRS)
                       .buffer(SCALE)
                       .bounds(maxError=1, proj=CRS)
               )
          .rename('mask')
    )

    return combined_mask