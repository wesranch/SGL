### Calculating depth and volume from 2021 landsat scenes
### Wesley Rancher, Nathan Rowley, Lily Bechina
### 10 July 2024
### Goals for the script:

# 1. read in masked lakes from NDWI at different dates alonside all scenes
        # composite scenes with B2, B3, B4
        # B3 scenes are read in seperately
# 2. Use Pope method to estimate lake volume for these 10 lakes for each landsat date

# Subsequent steps:

# 3. Sneed/Hamilton do the same (shall we compare two depth-reflectance methods?)
# 4. Using the delineated lake perimeter (from NDWI), extract lake elevations from the DEM (for 10 lakes for each Landsat date)
# 5. Watersheds have been delineated → identify a surface melt model (PDD and ?? energy balance model??) to melt water within the watershed 
#################################################################################################################################################

################################################################
# landsat data acquisition

# start session
import ee
import geemap
ee.Authenticate()
ee.Initialize()

# add lakes for filter landsat collection
prospective_lakes = ee.FeatureCollection('users/wesranch/WV01_2021_0721_selected_watersheds_shp')
Map = geemap.Map(center=[70, -50], zoom=5, basemap='Esri.WorldGrayCanvas')
Map.addLayer(prospective_lakes, {}, 'wtrshds from NR')
Map


#########################################################################################################################################
# attempting cloud mask
begin = '2022-05-15'
end = '2022-09-30'
# filter the landsat 8 collection where overlapping with watersheds of interest 
l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2') \
    .filterDate(begin, end) \
    .filterBounds(prospective_lakes) \
    .filter(ee.Filter.lte('CLOUD_COVER', 10)) \
    .sort('CLOUD_COVER')
l8_as_list = l8.toList(l8.size())
tcc = {'bands': ['SR_B4', 'SR_B3', 'SR_B2']}
img0 = ee.Image(l8_as_list.get(0))
img1 = ee.Image(l8_as_list.get(1))
img2 = ee.Image(l8_as_list.get(2))
img3 = ee.Image(l8_as_list.get(3))
img4 = ee.Image(l8_as_list.get(4))
img5 = ee.Image(l8_as_list.get(5))
img6 = ee.Image(l8_as_list.get(6))
img7 = ee.Image(l8_as_list.get(7))
img8 = ee.Image(l8_as_list.get(8))
img9 = ee.Image(l8_as_list.get(9))
Map.addLayer(img0, tcc, 'tcc 0')
Map.addLayer(img1, tcc, 'tcc 1')
Map.addLayer(img2, tcc, 'tcc 2')
Map.addLayer(img3, tcc, 'tcc 3')
Map.addLayer(img4, tcc, 'tcc 4')
Map.addLayer(img5, tcc, 'tcc 5')
Map.addLayer(img6, tcc, 'tcc 6')
Map.addLayer(img7, tcc, 'tcc 7')
Map.addLayer(img8, tcc, 'tcc 8')
Map.addLayer(img9, tcc, 'tcc 9')

def cloud_mask_landsat8(img):
    img = ee.Image(img)
    quality_band = img.select('QA_PIXEL')
    shadow = quality_band.bitwiseAnd(8).neq(0)  # gets the areas containing shadows (bitwiseAnd checks whether the quality band contains the bit 3 (2^3 = 8))
    cloud = quality_band.bitwiseAnd(32).neq(0)  # gets the areas containing clouds
    # cloud confidence is comprised of bits 6-7.
    # add the two bits and interpolate them to a range from 0-3.
    # 0 = None, 1 = Low, 2 = Medium, 3 = High.
    cloud_confidence = quality_band.bitwiseAnd(64).add(quality_band.bitwiseAnd(128)).interpolate([0, 64, 128, 192], [0, 1, 2, 3], 'clamp').int()
    cloud_confidence_medium_high = cloud_confidence.gte(2)
    cloudM = shadow.Or(cloud).Or(cloud_confidence_medium_high).select([0], ['cloudM'])
    # add cirrus confidence to cloud mask (cloudM) for Landsat 8
    cirrus_confidence = quality_band.bitwiseAnd(256).add(quality_band.bitwiseAnd(512)).interpolate([0, 256, 512, 768], [0, 1, 2, 3], 'clamp').int()
    cirrus_confidence_medium_high = cirrus_confidence.gte(2)
    cloudM = cloudM.Or(cirrus_confidence_medium_high)
    cloudM = cloudM.Not()  # Not required to swap 0 and 1 (so clouds have number 0)
    # mask image with cloud mask
    image_cloud_masked = img.updateMask(cloudM)  # second mask removes places where cloudM has zeroes (needed to avoid strips between sentinel scenes)
    # add cloud mask as band
    image_cloud_masked = image_cloud_masked.addBands(cloudM)
    return image_cloud_masked

l8_cmasked = l8.map(cloud_mask_landsat8) # apply cloud mask

# apply ndwi
def compute_ndwi_ice(image):
    ndwi_ice = image.normalizedDifference(['SR_B2', 'SR_B4']).rename('NDWI_Ice')
    return image.addBands(ndwi_ice)

l8_cmasked_w_ndwi = l8_cmasked.map(compute_ndwi_ice)
l8_cmasked_ndwi_list = l8_cmasked_w_ndwi.toList(l8_cmasked_w_ndwi.size())

# visualize
ndwi_params = {'bands': ['NDWI_Ice'],'min': -1,'max': 1,'palette': ['gray', 'blue']}
image_0 = ee.Image(l8_cmasked_ndwi_list.get(0))
image_1 = ee.Image(l8_cmasked_list.get(1))
image_2 = ee.Image(l8_cmasked_list.get(2))
image_3 = ee.Image(l8_cmasked_list.get(3))
image_4 = ee.Image(l8_cmasked_list.get(4))

Map.addLayer(image_0, ndwi_params, 'ndwi image 0')

# export cmasked ndwi images
def export_ndwi(image_collection, label):
    image_list = image_collection.toList(image_collection.size())
    image_count = image_collection.size().getInfo()

    for i in range(image_count):
        image = ee.Image(image_list.get(i))
        image_id = image.get('system:index').getInfo()  

        export_task = ee.batch.Export.image.toDrive(
            image=image.select('NDWI_Ice'),
            description=f'{image_id}_NDWI',
            folder='gee_SGL_exports',
            scale=30,
            region=image.geometry().bounds(),
            fileFormat='GeoTIFF'
        )
        export_task.start()
        print(f"exporting {label}, image {i+1} of {image_count} with ID: {image_id}")

export_ndwi(l8_cmasked_w_ndwi, 'ndwi ice')
###########################################################################################################################################
# individual images
# dates of interest
    # 05/28/21
    # 6/22/21
    # 7/08/21
    # 7/31/21
    # 8/23/21 009011 and 009012
    # 09/01/21

may28 = ee.Image('LANDSAT/LC08/C02/T1_TOA/LC08_008012_20210528').select('B2', 'B3', 'B4').multiply(10000.0).toUint16()
june22 = ee.Image('LANDSAT/LC08/C02/T1_TOA/LC08_007012_20210622').select('B2', 'B3', 'B4')
july8 = ee.Image('LANDSAT/LC08/C02/T1_TOA/LC08_007012_20210708').select('B2', 'B3', 'B4')
july31 = ee.Image('LANDSAT/LC08/C02/T1_TOA/LC08_008012_20210731').select('B2', 'B3', 'B4')
aug23_1 = ee.Image('LANDSAT/LC08/C02/T1_TOA/LC08_009012_20210823').select('B2', 'B3', 'B4')
aug23_2 = ee.Image('LANDSAT/LC08/C02/T1_TOA/LC08_009011_20210823').select('B2', 'B3', 'B4')
sep1 = ee.Image('LANDSAT/LC08/C02/T1_TOA/LC08_008012_20210901').select('B2', 'B3', 'B4')

# specific images for TOA
may28 = ee.Image('LANDSAT/LC08/C02/T1_TOA/LC08_008012_20210528').select('B8')
june22 = ee.Image('LANDSAT/LC08/C02/T1_TOA/LC08_007012_20210622').select('B8')
july8 = ee.Image('LANDSAT/LC08/C02/T1_TOA/LC08_007012_20210708').select('B8')
july31 = ee.Image('LANDSAT/LC08/C02/T1_TOA/LC08_008012_20210731').select('B8')
aug23_1 = ee.Image('LANDSAT/LC08/C02/T1_TOA/LC08_009012_20210823').select('B8')
aug23_2 = ee.Image('LANDSAT/LC08/C02/T1_TOA/LC08_009011_20210823').select('B8')
sep1 = ee.Image('LANDSAT/LC08/C02/T1_TOA/LC08_008012_20210901').select('B8')
toa_images = [may28, june22, july8, july31, aug23_1, aug23_2, sep1]

# convert to surface reflectance
mult = 0.00002
additive = -0.100000

# obtained by inspecting image properties (sun elevation)
se_may28 = 43.1666093
se_june22 = 44.98071892
se_july8 = 43.91102106
se_july31 = 39.6100356
se_aug23_1 = 32.79109399
se_aug23_2 = 31.53688709
se_sep1 = 29.64967597

# math
import math
#may28_b8_sr = (may28_b8*0.00002-0.100000) / sin(43.1666093 * pi / 180)

may28_b8_sr = (may28.multiply(mult).subtract(additive)).divide(math.sin(math.radians(se_may28)))
june22_b8_sr = (june22.multiply(mult).subtract(additive)).divide(math.sin(math.radians(se_june22)))
july8_b8_sr = (july8.multiply(mult).subtract(additive)).divide(math.sin(math.radians(se_july8)))
july31_b8_sr = (july31.multiply(mult).subtract(additive)).divide(math.sin(math.radians(se_july31))).toUint16()
aug23_1_b8_sr = (aug23_1.multiply(mult).subtract(additive)).divide(math.sin(math.radians(se_aug23_1)))
aug23_2_b8_sr = (aug23_2.multiply(mult).subtract(additive)).divide(math.sin(math.radians(se_aug23_2)))
sep1_b8_sr = (sep1.multiply(mult).subtract(additive)).divide(math.sin(math.radians(se_sep1)))

# export images to google drive
may28_task = ee.batch.Export.image.toDrive(
    image=may28,
    description='TOA_B8_008012_20210528',
    folder='gee_SGL_exports',
    scale=15,            
    region=may28.geometry().bounds(),
    fileFormat='GeoTIFF',
    maxPixels = 1e13)

june22_task = ee.batch.Export.image.toDrive(
    image=june22,
    description='TOA_B8_007012_20210622',
    folder='gee_SGL_exports',
    scale=15,            
    region=june22.geometry().bounds(),
    fileFormat='GeoTIFF',
    maxPixels = 1e13)

july8_task = ee.batch.Export.image.toDrive(
    image=july8,
    description='TOA_B8_007012_20210708',
    folder='gee_SGL_exports',
    scale=15,            
    region=july8.geometry().bounds(),
    fileFormat='GeoTIFF',
    maxPixels = 1e13)

july31_task = ee.batch.Export.image.toDrive(
    image=july31,
    description='TOA_B8_008012_20210731',
    folder='gee_SGL_exports',
    scale=15,            
    region=july31.geometry().bounds(),
    fileFormat='GeoTIFF',
    maxPixels = 1e13)

aug23_1_task = ee.batch.Export.image.toDrive(
    image=aug23_1,
    description='TOA_B8_009012_20210823',
    folder='gee_SGL_exports',
    scale=15,            
    region=aug23_1.geometry().bounds(),
    fileFormat='GeoTIFF',
    maxPixels = 1e13)

aug23_2_task = ee.batch.Export.image.toDrive(
    image=aug23_2,
    description='TOA_B8_009011_20210823',
    folder='gee_SGL_exports',
    scale=15,            
    region=aug23_2.geometry().bounds(),
    fileFormat='GeoTIFF',
    maxPixels = 1e13)

sep1_task = ee.batch.Export.image.toDrive(
    image=sep1,
    description='TOA_B8_008012_20210901',
    folder='gee_SGL_exports',
    scale=15,            
    region=sep1.geometry().bounds(),
    fileFormat='GeoTIFF',
    maxPixels = 1e13)

may28_task.start()
june22_task.start()
july8_task.start()
july31_task.start()
aug23_1_task.start()
aug23_2_task.start()
sep1_task.start()

####################################################################

# read in normalized surface reflectance images 
import rasterio
import numpy as np
import os
import glob
import matplotlib.pyplot as plt
from rasterio.plot import show

# directory where landsat surface reflectance images are stored
dir = 'E:/rsl_owu/SGL/data/output/gee'
os.chdir(dir)

# read in scenes
landsat_images = glob.glob("LC08*.tif")
landsat_images = [raster for raster in landsat_images if not raster.endswith('_NDWI.tif')]

#############################################################################
# print indices for images of interest
for i, image in enumerate(landsat_images):
    print(f"{i}: {image}")

indices = [16, 6, 7, 12, 8, 4, 15]
# dates of interest
    # 05/28/21
    # 6/22/21
    # 7/08/21
    # 7/31/21
    # 8/23/21 009011 and 009012
    # 09/01/21

# attempt to plot
landsat_rasters = [rasterio.open(landsat_images[i]) for i in indices]
for raster in landsat_rasters:
    show(raster)

landsat_rasters = {
    "may28_008012": landsat_images[16],
    "june22_007012": landsat_images[6],
    "july8_007012": landsat_images[7],
    "july31_008012": landsat_images[12],
    "sep1_008012": landsat_images[8],
    "aug23_009011": landsat_images[4],
    "aug23_009012": landsat_images[15]
}

