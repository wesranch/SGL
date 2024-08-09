# Inspecting potential ways to delineate lakes using NDWI
# Wesley Rancher
# 4 July 2024

# libs
library(terra)
library(sf)

# directory
setwd("E:/rsl_owu/SGL/data/output/gee/")

# read in a good image date of NDWI// starting to think we want mndwi
july_31_ndwi_p8_r12 <- rast("LC08_008012_20210731_NDWI.tif")
print(crs(july_31_ndwi_p8_r12))

# inspect statistics -- how does this change when clipped
summary(july_31_ndwi_p8_r12)

# read in watersheds for clipping
wtrshds <- st_read("E:/rsl_owu/SGL/data/output/shp/WV01_2021_0721_selected_watersheds.shp")
plot(wtrshds)

# reproject watersheds so crs match landsat scene
wtrshd_vec <- vect(wtrshds)
wtrshd_reproj <- project(wtrshd_vec, crs(july_31_ndwi_p8_r12))

ext(wtrshd_reproj)
plot(wtrshd_reproj)

# mask the raster by the watersheds vector
clipped_raster <- crop(july_31_ndwi_p8_r12, wtrshd_reproj) #clips to the extent not the geom
summary(clipped_raster)
plot(clipped_raster) #alright this works

# what should the threshold be???
lakes <- july_31_ndwi_p8_r12 >= 0.483 #value from symbology in arc
plot(lakes)
#writeRaster(lakes, "test_lake_delineation.tif")


# quantile method of 90th pctl for determining water body threshold
# alternatively if there is a literature value // 0.25 on ndwi (Williamson et al 2018)
# we could download one mndwi to compare to this test
setwd("E:/rsl_owu/SGL/data/output/raster/")

lakes_25_thres <- july_31_ndwi_p8_r12 >= 0.25
plot(lakes_25_thres)

# calculate ndwi ice and compare
july31_mndwi <- rast("LC08_008012_20210731_MNDWI.tif")
july_31_p8_r12 <- rast("LC08_008012_20210731.tif")
july31_ndwi_ice <- (july_31_p8_r12$SR_B2 - july_31_p8_r12$SR_B4) / (july_31_p8_r12$SR_B2 + july_31_p8_r12$SR_B4)

# plotting comparisons
plot(july31_ndwi_ice) #it looks the best
plot(july_31_ndwi_p8_r12) #ndwi using b3 and b5
plot(july31_mndwi) #mndwi

names(july31_ndwi_ice) <- "NDWI_Ice"
summary(july31_ndwi_ice)
summary(july31_mndwi)

lakes_using_ice_ndwi <- july31_ndwi_ice >= 0.25
lakes_mndwi <- july31_mndwi >= 0.25

plot(lakes_using_ice_ndwi) #best for replication
plot(lakes_mndwi)
plot(lakes2) #no sir -- using threshold on regular ndwi
plot(lakes) # assigned the value based on symbology

summary(lakes)
summary(lakes_using_ice_ndwi) #compares pretty well to visual thresholding
summary(lakes_mndwi)

# clip to wtrshd
clipped_lakes <- crop(lakes_using_ice_ndwi, wtrshd_reproj)
clipped_mndwi_lakes <- crop(lakes_mndwi, wtrshd_reproj)

plot(clipped_lakes) #lets use ndwi_ice 
plot(clipped_mndwi_lakes)

##############################
# set NA to anything less than zero in lake mask for sieving in Q
july31 <- rast("thres_mask_008012_20210731.tif")
plot(july31)

july31_binary <- july31$NDWI_Ice > 0 
plot(july31_binary)
summary(july31_binary)

# should be good
writeRaster(july31_binary, "july31_thres_binary.tif")
  
  
  
  
