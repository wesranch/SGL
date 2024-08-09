# Creating NDWI maps and Delineating Lakes in watersheds
# Wesley Rancher, Nathan Rowley, Lily Bechina, Chris Karmosky
# 4 July 2024

## DROPBOX file designations for this script:
## landsat_2021_images are in DATA/GEE_Exports
## wtrshds is in DATA/DEM
## files exported in the end can be found in DATA/Raster

# directory and list files
dir <- "E:/rsl_owu/SGL/data/output/gee/"
setwd(dir)
library(terra)
library(sf)

# this will fetch ALL the tiffs in the dir and might require cleaning
list_of_files <- basename(list.files(pattern = "^LC08.*\\.tif$", full.names = TRUE))

# rm NDWI tiffs and the one MNDWI
list_of_files <- list_of_files[!grepl("_NDWI|_MNDWI", list_of_files)]
print(list_of_files)

# read each file in the list as a raster
landsat_2021_images <- lapply(list_of_files, rast)
plot(landsat_2021_images[[1]])

# function for calculating ndwi
# blue and red
ndwi_ice_calc <- function(x) {
  ndwi_ice <- (x$SR_B2 - x$SR_B4) / (x$SR_B2 + x$SR_B4)
  names(ndwi_ice) <- "NDWI_Ice"
  return(ndwi_ice)
}

# apply the ndwi_ice function to each raster
ndwi_ice_images <- lapply(landsat_2021_images, ndwi_ice_calc)
summary(ndwi_ice_images[[1]]) #make sure band name is ndwi_ice
plot(ndwi_ice_images[[17]])
# assign threshold from Williamson et al., 2018
threshold <- 0.25
gte_thres <- function(x) {
  x[x >= threshold] <- 1
  x[x < threshold] <- 0
  return(x)
}

# apply it to each image
ndwi_ice_images_thres <- lapply(ndwi_ice_images, gte_thres)
for (i in seq_along(ndwi_ice_images_thres)) {
  plot(ndwi_ice_images_thres[[i]])
}

################################################################################
# read in aoi (watersheds from Rowley)
wtrshds <- vect("E:/rsl_owu/SGL/data/output/shp/WV01_2021_0721_selected_watersheds.shp")
ref_rast <- ndwi_ice_images[[1]]
wtrshd_reproj <- project(wtrshds, crs(ref_rast)) #reproject
wtrshd_sf <- st_as_sf(wtrshd_reproj)

# read in ndwi files
setwd("E:/rsl_owu/SGL/data/output/raster/")
ndwi_files <- basename(list.files(pattern = "NDWI_.*\\.tif$", full.names = TRUE))
ndwi_raster_list <- lapply(ndwi_files, rast)


# use this if using lapply method in lieu of for loop above
#crop_to_wtrshd_ext <- function(x) {
#  crop(x, wtrshd_reproj)
#}
#ndwi_ice_thres_cropped <- lapply(ndwi_ice_images_thres, crop_to_wtrshd_ext)

# inspect that the clipping worked with some plots
for (i in seq_along(ndwi_ice_thres_cropped)) {
  plot(ndwi_ice_thres_cropped[[i]])
}

# masking function within the wtrshd geometries
mask_to_wtrshd <- function(x) {
  mask(x, wtrshd_sf)
}
ndwi_ice_thres_masked <- lapply(ndwi_ice_images_thres, mask_to_wtrshd)

# see how it looks // if it takes up a whole wtrshd then don't use or change thres
for (i in seq_along(ndwi_ice_thres_masked)) {
  plot(ndwi_ice_thres_masked[[i]])
}

################################################################################
# function for raster to polygon
raster_to_polygons <- function(raster_chunk) {
  pols <- as.polygons(raster_chunk, dissolve = TRUE)
  st_as_sf(pols)
}
ndwi_images_as_polygons <- lapply(ndwi_ice_thres_masked, raster_to_polygons)

for (i in seq_along(ndwi_images_as_polygons)) {
  plot(ndwi_images_as_polygons[[i]])
}

one <- ndwi_images_as_polygons[[1]]
plot(ndwi_images_as_polygons[[5]])

out_dir <- "E:/rsl_owu/SGL/data/output/shp/"
write_raster_polygons <- function(ndwi_images_as_polygons, ndwi_files, out_dir) {
  for (i in seq_along(ndwi_images_as_polygons)) {
    shp <- ndwi_images_as_polygons[[i]]
    filename <- paste0(out_dir, "lake_polygons_", extract_date_n_prow(ndwi_files[i]), ".shp")
    st_write(shp, filename)
  }
}
write_raster_polygons(ndwi_images_as_polygons, ndwi_files, out_dir) 

###### save the processed rasters
out_dir <- "E:/rsl_owu/SGL/data/output/raster/"

# call the individual rasters from list and assign naming convention
extract_date_n_prow <- function(filename) {
  substr(filename, 10, 24) #position of date in file nomenclature
}
test_file_names <- extract_date_n_prow(ndwi_files)

# function to save each raster

# this one for clipped and masked
write_clipped_masked <- function(ndwi_ice_thres_masked, list_of_files, out_dir) {
  for (i in seq_along(ndwi_ice_thres_masked)) {
    raster <- ndwi_ice_thres_masked[[i]]
    filename <- paste0(out_dir, "delin_lakes_", extract_date_n_prow(list_of_files[i]), ".tif")
    writeRaster(raster, filename)
  }
}
write_clipped_masked(ndwi_ice_thres_masked, list_of_files, out_dir) 

# this one for pre threshold ndwi ice imagery
write_ndwi_ice <- function(ndwi_ice_images, list_of_files, out_dir) {
  for (i in seq_along(ndwi_ice_images)) {
    raster <- ndwi_ice_images[[i]]
    filename <- paste0(out_dir, "NDWI_ice_", extract_date_n_prow(list_of_files[i]), ".tif")
    writeRaster(raster, filename, overwrite = TRUE)
  }
}
write_ndwi_ice(ndwi_ice_images, list_of_files, out_dir) 





