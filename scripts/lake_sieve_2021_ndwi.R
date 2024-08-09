# lake sieving ndwi delineated scenes 
# Wesley Rancher, Nathan Rowley, Chris Karmosky, Lily Bechina
# source: https://github.com/allenpope/Landsat8_SupraglacialLakes/blob/SermeqKujalleq2015/Depth/Jakobshavn_lake_sieve.m
# adapted to R

# libs
library(terra)

# dir and outdir
dir <- "E:/rsl_owu/SGL/data/output/raster/"
setwd(dir)

# list of all landsat scenes
ndwi_paths_thres <- basename(list.files(pattern = "^thres_.*\\.tif$", full.names = TRUE))
ndwi_paths <- basename(list.files(pattern = "NDWI_.*\\.tif$", full.names = TRUE))
threshold_rast_list <- lapply(ndwi_paths_thres, rast)
ndwi_rast_list <- lapply(ndwi_paths, rast)


################################################################################
# function to remove outlier pixels from ndwi threshold images
lake_sieve <- function(ndwi_thres_raster) {
  # get connected components
  connected_comp <- patches(ndwi_thres_raster, directions = 4, zeroAsNA = TRUE)
  components <- unique(values(connected_comp), na.rm = TRUE)
  components <- components[components != 0 & !is.na(components)]
  cell_indices <- unique(values(connected_comp, na.rm = TRUE))
  # empty mask for valid lakes
  valid_lake_mask <- rast(ndwi_thres_raster)
  values(valid_lake_mask) <- 0
  min_pixel <- 10
  min_width <- 1
  #loop over each ndwi scene
  for (comp_id in components) {
    # get cell indices for the current cc
    cell_indices <- which(values(connected_comp) == comp_id)
    if (length(cell_indices) < min_pixel) next
    # convert cell indices to coordinates
    coords <- xyFromCell(ndwi_thres_raster, cell_indices)
    width <- length(unique(coords[, "x"]))
    height <- length(unique(coords[, "y"]))
    # check if the component meets width/height requirements
    if (width <= min_width || height <= min_width) next
    values(valid_lake_mask)[cell_indices] <- 1 #update valid lake
  }
  return(valid_lake_mask)
}

sieved_lake_list <- lapply(threshold_rast_list, lake_sieve)
plot(sieved_lake_list[[3]])


################################################################################
# save
extract_date_n_prow <- function(filename) {
  substr(filename, 12, 26) #position of date in file nomenclature
}
test_file_names <- extract_date_n_prow(ndwi_paths_thres)

write_clipped_masked <- function(sieved_lake_list, ndwi_paths_thres, dir) {
  for (i in seq_along(sieved_lake_list)) {
    raster <- sieved_lake_list[[i]]
    filename <- paste0(dir, "SIEVED_", extract_date_n_prow(ndwi_paths_thres[i]), ".tif")
    writeRaster(raster, filename)
  }
}
write_clipped_masked(sieved_lake_list, ndwi_paths_thres, dir) 

