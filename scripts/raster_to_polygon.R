# shpfile from delineated lake rasters
# fed in from delineate_lakes_from_ndwi.R
# Wesley Rancher, Nathan Rowley, Lily Bechina
# 4 July 2024

# libraries
library(terra)
library(sf)
library(RColorBrewer)

# dir
dir <- "E:/rsl_owu/SGL/data/output/raster/"
setwd(dir)

# read in rasters
raster_paths <- basename(list.files(pattern = "^delin_lakes_.*\\.tif$", full.names = TRUE))
ndwi_paths <- basename(list.files(pattern = "^NDWI_ice_.*\\.tif$", full.names = TRUE))
raster_list <- lapply(raster_paths, rast)
ndwi_list <- lapply(ndwi_paths, rast)

# plot
for (i in seq_along(raster_list)) {
  plot(raster_list[[i]])
}

for (i in seq_along(ndwi_list)) {
  plot(ndwi_list[[i]])
}

# again for ndwi tiffs but at a given threshold
threshold <- 0.25 
gte_thres <- function(x) {
  x[x >= threshold] <- 1
  x[x < threshold] <- 0
  return(x)
}
# apply it to each image
ndwi_thres_list <- lapply(ndwi_list, gte_thres)
color_palette <- c("gray", "blue") 
for (i in seq_along(ndwi_thres_list)) {
  plot(ndwi_thres_list[[i]], col = color_palette)
}

### conversion
# loop over rasters and polygon them
poly_list <- list()
for (i in seq_along(raster_list)){
  poly <- as.polygons(raster_list[[i]])
  poly_sf <- st_as_sf(poly)
  poly_list[[i]] <- poly_sf
}

# plot
plot(poly_list[[1]])
plot(poly_list[[2]])
plot(poly_list[[3]])
plot(poly_list[[4]])
plot(poly_list[[5]])
plot(poly_list[[6]])
plot(poly_list[[7]])

# save // I keep changing this so careful of output
out_dir <- "E:/rsl_owu/SGL/data/output/shp/"
out_dir <- "E:/rsl_owu/SGL/data/output/raster/"

# retain path, row, and date
extract_date_n_prow <- function(filename) {
  substr(filename, 10, 24) #position of date in ndwi_paths
}
#test_file_names <- extract_date_n_prow(raster_paths)
test_file_names <- extract_date_n_prow(ndwi_paths)

# function to save delineated shapes 
write_shp_files <- function(poly_list, raster_paths, out_dir) {
  for (i in seq_along(poly_list)) {
    shp <- poly_list[[i]]
    filename <- paste0(out_dir, "delin_lakes_", extract_date_n_prow(raster_paths[i]), ".shp")
    st_write(shp, filename)
  }
}
write_shp_files(poly_list, raster_paths, out_dir) 

# function to save thresholded ndwi images
write_ndwi_thres <- function(ndwi_thres_list, ndwi_paths, out_dir) {
  for (i in seq_along(ndwi_thres_list)) {
    raster <- ndwi_thres_list[[i]]
    filename <- paste0(out_dir, "ndwi_gte25_", extract_date_n_prow(ndwi_paths[i]), ".tif")
    writeRaster(raster, filename)
  }
}
write_ndwi_thres(ndwi_thres_list, ndwi_paths, out_dir) 
