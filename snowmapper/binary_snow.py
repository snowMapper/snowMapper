"""
_____________________________________________
name: snowMapper 
doi: 10.5281/zenodo.17663731
_____________________________________________
__________________________________________________________________________________________
Description: Converts reflectance from clear-sky pixels into binary snow cover
(i.e. '0': no snow, '1': snow).

Input parameters:
- img (ee.Image): Preprocessed & daily-composited image with spectral bands.
- thresholds (dict): Thresholds for the snow cover mapping approaches.
- active_method (str): Method that has been selected in the configuration file for
                       binarising snow.
- SCALE (int): Desired scale of grid.
- domain_ee (ee.FeatureCollection): Region of interest (e.g. mountain range).

Internal functions:
- register_method(): Finds the method that needs to be applied.
- raw_binary(): Applies the method selected for retrieving binary snow cover.

Output:
- ee.ImageCollection of binary snow cover.
__________________________________________________________________________________________

"""
#===============================================================================
# Load libraries
#===============================================================================
import ee
import math

#===============================================================================
# Get Earth Engine started
#===============================================================================
ee.Authenticate()
ee.Initialize()

#===============================================================================
# Method registry
#===============================================================================
METHOD_REGISTRY = {}

def register_method(name):
    def decorator(func):
        METHOD_REGISTRY[name] = func
        return func
    return decorator

#===============================================================================
# Method: NDSI
#===============================================================================
@register_method("ndsi")
def method_ndsi(img, thresholds, domain_ee=None, SCALE=None):
    ndsi = img.normalizedDifference(["green", "swir1"])  # Normalised Difference Snow Index
    ndsi_condition = ndsi.gt(thresholds["NDSI_THRES"])
    
    sc_obs = ndsi_condition
    return img.addBands(sc_obs.rename("sc_obs")).select("sc_obs").toUint8().copyProperties(img, ["system:time_start"])

#===============================================================================
# Method: Otsu NDSI, based on Otsu (1979)
#===============================================================================
@register_method("Otsu_ndsi")
def method_otsu_ndsi(img, thresholds, domain_ee, SCALE):
    ndsi = img.normalizedDifference(["green", "swir1"]).rename('NDSI')  # Normalised Difference Snow Index    
    
    # -----------------------------
    # Otsu threshold function
    # -----------------------------
    def otsu(hist_dict):
        counts = ee.List(hist_dict.get('histogram'))
        means = ee.List(hist_dict.get('bucketMeans'))
    
        size = counts.length()
        total = ee.Number(counts.reduce(ee.Reducer.sum()))
    
        # total mean
        sum_ = ee.Number(
            ee.List(means.zip(counts)).map(
                lambda x: ee.Number(ee.List(x).get(0))
                .multiply(ee.Number(ee.List(x).get(1)))
            ).reduce(ee.Reducer.sum())
        )
    
        indices = ee.List.sequence(1, size.subtract(1))
    
        def bss_func(i):
            i = ee.Number(i)
    
            counts1 = counts.slice(0, i)
            counts2 = counts.slice(i)
    
            means1 = means.slice(0, i)
            means2 = means.slice(i)
    
            w1 = ee.Number(counts1.reduce(ee.Reducer.sum())).divide(total)
            w2 = ee.Number(1).subtract(w1)
    
            mu1 = ee.Number(
                ee.List(means1.zip(counts1)).map(
                    lambda x: ee.Number(ee.List(x).get(0))
                    .multiply(ee.Number(ee.List(x).get(1)))
                ).reduce(ee.Reducer.sum())
            ).divide(ee.Number(counts1.reduce(ee.Reducer.sum())))
    
            mu2 = ee.Number(
                ee.List(means2.zip(counts2)).map(
                    lambda x: ee.Number(ee.List(x).get(0))
                    .multiply(ee.Number(ee.List(x).get(1)))
                ).reduce(ee.Reducer.sum())
            ).divide(ee.Number(counts2.reduce(ee.Reducer.sum())))
    
            return w1.multiply(w2).multiply(mu1.subtract(mu2).pow(2))
    
        bss = indices.map(bss_func)
    
        maxval = ee.Number(bss.reduce(ee.Reducer.max()))
        idx = bss.indexOf(maxval)
    
        return ee.Number(means.get(idx))

    # -----------------------------
    # Compute histogram
    # -----------------------------
    hist = ndsi.reduceRegion(
        reducer=ee.Reducer.histogram(maxBuckets=256),
        geometry=domain_ee,
        scale=SCALE,
        bestEffort=True
    )
    hist_dict = ee.Dictionary(hist.get('NDSI'))
    
    # -----------------------------
    # Calculate & apply optimised NDSI threshold
    # -----------------------------    
    otsu_ndsi_thres = otsu(hist_dict)
    ndsi_condition = ndsi.gt(ee.Image.constant(otsu_ndsi_thres))
    
    sc_obs = ndsi_condition
    return img.addBands(sc_obs.rename("sc_obs")).select("sc_obs").toUint8().copyProperties(img, ["system:time_start"])

#===============================================================================
# Method: Gascoin et al. (2015)
#===============================================================================
@register_method("Gascoin_2015")
def method_gascoin_2015(img, thresholds, domain_ee=None, SCALE=None):
    ndsi = img.normalizedDifference(["green", "swir1"])  # Normalised Difference Snow Index
    ndsi_condition = ndsi.gt(thresholds["NDSI_THRES"])
    red_condition = img.select("red").gt(thresholds["RED_THRES"])
    swir_condition = img.select("swir1").lt(thresholds["SWIR_THRES"])

    sc_obs = ndsi_condition.And(red_condition).And(swir_condition)
    return img.addBands(sc_obs.rename("sc_obs")).select("sc_obs").toUint8().copyProperties(img, ["system:time_start"])

#===============================================================================
# Method: Gascoin et al. (2019)
#===============================================================================
@register_method("Gascoin_2019")
def method_gascoin_2019(img, thresholds, domain_ee=None, SCALE=None):
    ndsi = img.normalizedDifference(["green", "swir1"])  # Normalized Difference Snow Index
    ndsi_condition = ndsi.gt(thresholds["NDSI_THRES"])
    red_condition = img.select("red").gt(thresholds["RED_THRES"])

    sc_obs = ndsi_condition.And(red_condition)
    return img.addBands(sc_obs.rename("sc_obs")).select("sc_obs").toUint8().copyProperties(img, ["system:time_start"])

#===============================================================================
# Method: Koehler et al. (2022)
#===============================================================================
@register_method("Koehler_2022")
def method_koehler_2022(img, thresholds, domain_ee=None, SCALE=None):
    ndsi = img.normalizedDifference(["green", "swir1"])  # Normalised Difference Snow Index
    ndvi = img.normalizedDifference(["nir", "red"])      # Normalised Difference Vegetation Index
    ndwi = img.normalizedDifference(["green", "nir"])    # Normalised Difference Water Index
    ndli = img.normalizedDifference(["blue", "nir"])     # Normalised Difference Lignin Index

    # Snow binarization conditions
    ndsi_condition = ndsi.gt(thresholds["NDSI_THRES"])
    ndvi_condition = ndvi.lt(thresholds["NDVI_THRES"])
        
    water_condition = (
        ndwi.gt(thresholds["WATER_NDWI_THRES"])
            .And(img.select("green").lt(thresholds["WATER_GREEN_THRES"]))
    )
    shadow_condition = (
        img.select("swir1").lt(thresholds["SHADOW_SWIR_THRES"])
           .And(img.select("green").lt(thresholds["SHADOW_GREEN_THRES"]))
           .And(ndvi.lt(thresholds["SHADOW_NDVI_THRES"]))
    )
    ndli_condition = ndli.gt(thresholds["NDLI_THRES"])

    sc_obs = (
        ndsi_condition.And(ndvi_condition)
                      .And(water_condition.Not())
                      .And(shadow_condition.Not())
                      .And(ndli_condition.Not())
    )
    
    # Optional temperature condition
    temp_thres = thresholds.get("TEMP_THRES")
    if temp_thres is not None:
        temp_condition = img.select("temp").lt(temp_thres)
        sc_obs = sc_obs.And(temp_condition)
    
    return img.addBands(sc_obs.rename("sc_obs")).select("sc_obs").toUint8().copyProperties(img, ["system:time_start"])

#===============================================================================
# Method: Wang et al. (2025)
#===============================================================================
@register_method("Wang_2025")
def method_wang_2025(img, thresholds, domain_ee=None, SCALE=None):
    ndsi = img.normalizedDifference(["green", "swir1"])  # Normalised Difference Snow Index
    ndfsi = img.normalizedDifference(["nir", "swir1"])   # Normalised Difference Forest Snow Index
    ndwi = img.normalizedDifference(["green", "nir"])    # Normalised Difference Water Index

    # Forest cover masks
    landcover = ee.ImageCollection("ESA/WorldCover/v100").first().select(["Map"], ["forest_cover"]).byte()
    forest = landcover.eq(10)
    non_forest = landcover.neq(10)

    # Terrain data
    dem = ee.Image("USGS/SRTMGL1_003")
    terrain = ee.Algorithms.Terrain(dem)
    slope = terrain.select("slope").multiply(math.pi / 180)    # radians
    aspect = terrain.select("aspect").multiply(math.pi / 180)  # radians

    # Solar angles from metadata
    zenith_deg = ee.Number(img.get("SUN_ZENITH"))
    azimuth_deg = ee.Number(img.get("SUN_AZIMUTH"))
    zenith_rad = zenith_deg.multiply(math.pi / 180)
    azimuth_rad = azimuth_deg.multiply(math.pi / 180)

    # Convert to constant images for pixelwise operations
    zenith_img = ee.Image.constant(zenith_rad)
    azimuth_img = ee.Image.constant(azimuth_rad)

    # Solar incidence angle (in degrees)
    cos_i = (
        zenith_img.cos().multiply(slope.cos())
                  .add(zenith_img.sin().multiply(slope.sin()).multiply((aspect.subtract(azimuth_img)).cos()))
    )
    incidence_angle_rad = cos_i.acos()
    incidence_angle_deg = incidence_angle_rad.multiply(180 / math.pi)

    # Shadow mask
    shadow = ee.Terrain.hillShadow(
        image=dem.select("elevation"),
        azimuth=azimuth_deg,  # degrees
        zenith=zenith_deg,    # degrees
        neighborhoodSize=200,
        hysteresis=True
    )

    # Illumination condition masks
    shaded = shadow.eq(1)
    poorly_lit = incidence_angle_deg.gt(45).And(incidence_angle_deg.lt(90)).And(shaded.Not())
    non_shaded = incidence_angle_deg.lte(45).And(shaded.Not())

    # Snow binarization paths
    p1 = non_forest.And(ndsi.gt(thresholds["NDSI_THRES_P1"])).And(img.select("blue").gt(thresholds["BLUE_THRES_P1"]))
    p2 = (
        non_forest.And(poorly_lit).And(ndsi.gt(thresholds["NDSI_THRES_P2"]))
                  .And(img.select("blue").gt(thresholds["BLUE_THRES_P2"]))
                  .And(ndwi.lt(thresholds["NDWI_THRES_P2"]))
    )
    p3 = non_forest.And(shaded).And(ndsi.gt(thresholds["NDSI_THRES_P3"])).And(img.select("blue").gt(thresholds["BLUE_THRES_P3"]))
    p4 = (
        forest.And(non_shaded).And(ndsi.gt(thresholds["NDSI_THRES_P4"]))
              .And(ndfsi.gt(thresholds["NDFSI_THRES_P4"]))
              .And(img.select("blue").gt(thresholds["BLUE_THRES_P4"]))
    )
    p5 = (
        forest.And(shaded).And(ndsi.gt(thresholds["NDSI_THRES_P5"]))
              .And(ndfsi.lt(thresholds["NDFSI_THRES_P5"]))
    )

    # Combine all classification paths
    sc_obs = p1.Or(p2).Or(p3).Or(p4).Or(p5)

    # Optional temperature condition
    temp_thres = thresholds.get("TEMP_THRES")
    if temp_thres is not None:
        temp_condition = img.select("temp").lt(temp_thres)
        sc_obs = sc_obs.And(temp_condition)

    return (
        img.addBands(sc_obs.rename("sc_obs"))
           .select("sc_obs")
           .toUint8()
           .copyProperties(img, ["system:time_start"])
    )

#===============================================================================
# Main binary snow classification function
#===============================================================================
def binary_snow(collection, thresholds, active_method, domain_ee, SCALE):
    def raw_binary(img):
        return METHOD_REGISTRY[active_method](img, thresholds, domain_ee, SCALE)
    return collection.map(raw_binary).sort("system:time_start")