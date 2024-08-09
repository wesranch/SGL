### Calculating depth and volume from 2021 landsat scenes
### Wesley Rancher, Nathan Rowley, Lily Bechina, Chris Karmosky
### 25 July 2024

### Citation: Pope, A. (2016). Reproducibly estimating and evaluating 
### supraglacial lake depth with Landsat 8 and other multispectral sensors. 
### Earth and Space Science, 3(4), 176–188. https://doi.org/10.1002/2015EA000125


### Goals for the script:
# 1. read in masked lakes from NDWI at different dates
# 2. Use Pope method to estimate lake depth for each landsat date


################################################################################
# libs
library(terra)
library(pbapply)

# dirs
raster_dir <- "E:/rsl_owu/SGL/data/output/raster"
output_dir <- "E:/rsl_owu/SGL/data/output/"
gee_dir <- "E:/rsl_owu/SGL/data/output/gee"

# read in sieved images
setwd(raster_dir)
sieve_paths <- basename(list.files(pattern = "^SIEVED_.*\\.tif$", full.names = TRUE))
sieved_images <- lapply(sieve_paths, rast)
names(sieved_images) <- sieve_paths


################################################################################
# read in landsat images
setwd(gee_dir)
# file paths
land_paths <- basename(list.files(pattern = "^TOA_.*\\.tif$", full.names = TRUE))

# remove B8 scenes and read in rasters for individual bands
land_paths <- land_paths[!grepl("_B8_.", land_paths)]
red_scenes <- lapply(land_paths, rast, lyr = 3)
names(red_scenes) <- land_paths

# select for B8 scenes from file path and read in
land_paths <- basename(list.files(pattern = "^TOA_.*\\.tif$", full.names = TRUE))
land_paths <- land_paths[grepl("_B8_.", land_paths)]
pan_scenes <- lapply(land_paths, rast)
names(pan_scenes) <- land_paths


################################################################################
# calculate depth // Rinf values acquired from manual inspection of scenes in Arc
Rinf.csv <- read.csv("E:/rsl_owu/SGL/data/output/csv/Rinf.csv")
RinfB8.csv <- read.csv("E:/rsl_owu/SGL/data/output/csv/Rinf_B8.csv")

#g values calculated according to Pope et al., 2015 -- need to add this paper to LITERATURE
## G = [0.0178, 0.0341, 0.1413, 0.7507, 0.3817] # OLI 1, 2, 3, 4, 8
#g.blue <- 0.0341
#g.green <- 0.1413
g.red <- 0.7507
g.pan <- 0.3817

################################################################################
# function to pull rinf values from applicable scenes
get_B8_rinf <- function(scene_id, RinfB8.csv) {
  return(RinfB8.csv[RinfB8.csv$scene == scene_id, "b8"])  
}

# band 8 depth
depth_calc_b8 <- function(pan_scenes, i, sieved_images, RinfB8.csv, g.pan) {
  # extract scene id and rinf
  scene_id <- names(pan_scenes[i])
  rinf <- get_B8_rinf(scene_id, RinfB8.csv) #call the function
  g <- g.pan #panchromatic g-constant
  
  # read in B8 img as double from float
  doubled <- as.double(pan_scenes[[i]])
  masked <- sieved_images[[i]]
  resampled <- resample(doubled, masked) #15m to 30m 
  #lakes <- mask(resampled, masked) #mask by sieved scene
  lakes <- ifel(masked, resampled, NA) #since zeros are not NA
  
  # obtain lake specific bottom albedo
  Ad <- rast(lakes)
  values(Ad) <- 0 #init Ad
  # extract unique connected components from the lakes raster
  temp_file <- tempfile(pattern = "cc_", fileext = ".tif")
  connected_comp <- patches(lakes, directions = 4, zeroAsNA = TRUE, filename = temp_file)
  components <- unique(values(connected_comp))
  components <- components[components != 0 & !is.na(components)]
  for (comp_id in components) {
    cell_indices <- which(values(connected_comp) == comp_id) #pixel index for cc
    Ad_lake <- mean(values(lakes)[cell_indices], na.rm = TRUE) #Ad_lake
    values(Ad)[cell_indices] <- Ad_lake #assign Ad to app. cell indices
  }
   unlink(temp_file)
  rm(temp_file, connected_comp, Ad_lake, components, cell_indices)
  
  # Pope equation
  #Ad2 <- mean(values(lakes), na.rm = TRUE) #not a static value
  temp <- as.numeric(values(lakes) - rinf) #no neg.log
  valid_index <- which(temp > 0 & (values(Ad) - rinf) > 0)
  z <- rep(NA_real_, ncell(masked)) #init Z as numeric
  z[valid_index] <- (log(values(Ad)[valid_index] - rinf) - log(temp[valid_index])) / g
  z_raster <- masked #create a z raster for storing the values
  values(z_raster) <- as.numeric(z)
  z_raster[z_raster < 0] <- 0 #no negatives
  z_raster[1] <- 0 #adjusting the first value
  names(z_raster) <- "Depth"

  return(z_raster)  
}

b8_depth_list <- pblapply(seq_along(pan_scenes), function(i) {
  depth_calc_b8(pan_scenes, i, sieved_images, RinfB8.csv, g.pan)
})


################################################################################
# function to pull rinf values from applicable scenes
get_B4_rinf <- function(scene_id, Rinf.csv) {
  return(Rinf.csv[Rinf.csv$scene == scene_id, "b4"])
}

# band 4 depth
depth_calc_b4 <- function(red_scenes, i, sieved_images, Rinf.csv, g.red) {
  # extract scene id and rinf
  scene_id <- names(red_scenes[i])
  rinf <- get_B4_rinf(scene_id, Rinf.csv) #call the function
  g <- g.red #red g-constant
  
  # read in B8 img as double from float
  doubled <- as.double(red_scenes[[i]])
  masked <- sieved_images[[i]]
  resampled <- resample(doubled, masked) #15m to 30m 
  #lakes <- mask(resampled, masked) #mask by sieved scene
  lakes <- ifel(masked, resampled, NA) #since zeros are not NA

  # obtain lake specific bottom albedo
  Ad <- rast(lakes)
  values(Ad) <- 0 #init Ad
  # extract unique connected components from the lakes raster
  temp_file <- tempfile(pattern = "cc_", fileext = ".tif")
  connected_comp <- patches(lakes, directions = 4, zeroAsNA = TRUE, filename = temp_file)
  components <- unique(values(connected_comp))
  components <- components[components != 0 & !is.na(components)]
  for (comp_id in components) {
    cell_indices <- which(values(connected_comp) == comp_id) #pixel index for cc
    Ad_lake <- mean(values(lakes)[cell_indices], na.rm = TRUE) #Ad_lake
    values(Ad)[cell_indices] <- Ad_lake #assign Ad to app. cell indices
  }
  unlink(temp_file)
  rm(temp_file, connected_comp, Ad_lake, components, cell_indices)
  
  # Pope equation
  temp <- as.numeric(values(lakes) - rinf) #no neg.log
  valid_index <- which(temp > 0 & (values(Ad) - rinf) > 0)
  z <- rep(NA_real_, ncell(masked)) #init Z as numeric
  z[valid_index] <- (log(values(Ad)[valid_index] - rinf) - log(temp[valid_index])) / g
  z[valid_index] <- (log(Ad - rinf) - log(temp[valid_index])) / g
  z_raster <- masked #create a z raster for storing the values
  values(z_raster) <- as.numeric(z)
  z_raster[z_raster < 0] <- 0 #no negatives
  z_raster[1] <- 0 #adjusting the first value
  names(z_raster) <- "Depth"
  
  return(z_raster)  
}

b4_depth_list <- pblapply(seq_along(red_scenes), function(i) {
  depth_calc_b4(red_scenes, i, sieved_images, Rinf.csv, g.red)
})


################################################################################
# average band 4 and 8 depth and saving
average_depth <- function (b8_depth_raster, b4_depth_raster) {
  avg_depth <- (b8_depth_raster + b4_depth_raster) / 2
  names(avg_depth) <- "Depth"
  return(avg_depth)
}
average_depth_list <- mapply(average_depth, b8_depth_list, b4_depth_list, SIMPLIFY = FALSE)


################################################################################
# get path row from scene id
extract_date_n_prow <- function(filename) {
  substr(filename, 8, 22) #position of date in file nomenclature
}
test_file_names <- extract_date_n_prow(sieve_paths)

# save b8 depth
write_b8_depth <- function(b8_depth_list) {
  for (i in seq_along(b8_depth_list)) {
    raster <- b8_depth_list[[i]]
    filename <- paste0(raster_dir, "/V2_B8_DEPTH_", extract_date_n_prow(sieve_paths[i]), ".tif")
    writeRaster(raster, filename, overwrite = TRUE)
  }
}
write_b8_depth(b8_depth_list) 

# save b4 depth
write_b4_depth <- function(b4_depth_list) {
  for (i in seq_along(b4_depth_list)) {
    raster <- b4_depth_list[[i]]
    filename <- paste0(raster_dir, "/V2_B4_DEPTH_", extract_date_n_prow(sieve_paths[i]), ".tif")
    writeRaster(raster, filename, overwrite = TRUE)
  }
}
write_b4_depth(b4_depth_list) 

# save avg depth
write_avg_depth <- function(average_depth_list) {
  for (i in seq_along(average_depth_list)) {
    raster <- average_depth_list[[i]]
    filename <- paste0(raster_dir, "/V2_AVG_DEPTH_", extract_date_n_prow(sieve_paths[i]), ".tif")
    writeRaster(raster, filename, overwrite = TRUE)
  }
}
write_avg_depth(average_depth_list) 




