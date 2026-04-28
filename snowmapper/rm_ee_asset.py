"""
_____________________________________________
name: snowMapper 
doi: 10.5281/zenodo.17663731
licence: MIT
contact: konstantis.alexopoulos@gmail.com
_____________________________________________
__________________________________________________________________________________________
Description: Automated deletion of image collections.

Input parameters:
- collection_path (str): path for the collection to be deleted.

Output:
- Will delete the given collection.
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
# Method registry
#===============================================================================
def rm_ee_asset(collection_path):
    # Get children (images) under the collection
    children = ee.data.listAssets({'parent': collection_path}).get('assets', [])
    print(f"Found {len(children)} images")

    for child in children:
        asset_id = child['id']
        try:
            print(f"Deleting: {asset_id}")
            ee.data.deleteAsset(asset_id)
        except Exception as e:
            print(f"Failed to delete {asset_id}: {e}")

    # Delete the now-empty collection
    try:
        ee.data.deleteAsset(collection_path)
        print(f"Deleted collection: {collection_path}")
    except Exception as e:
        print(f"Failed to delete collection: {e}")