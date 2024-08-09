# calculate volume
# Wesley Rancher, Nathan Rowley, Lily Bechina
# 21 July 2024

# goals for the script:
# 3. clip to watershed ext and build attribute table


################################################################################
# libs
library(terra)
library(sf)
library(dplyr)


################################################################################
# read in lake perimeters from Rowley // raster to polygon of sieved scenes
# would allow for volume calculations of the whole scene

# read in depth scenes
# as.data.frame depth scenes
# then you would get x y cell depths
# clip this to each lake each date
# multiply depth value x y by 900m (this would be vol of ea. pixel)
# sum this across ncell in the clip

# dirs
shpdir <- "E:/rsl_owu/SGL/data/output/shp/lake_perims_v2/"
depth_dir <- "E:/rsl_owu/SGL/data/output/raster/"

# files 
lake_perim_paths <- list.files(shpdir, pattern = "\\.shp$", full.names = TRUE)
depth_files <- list.files(depth_dir, pattern = "AVG_DEPTH_.*\\.tif$", full.names = TRUE)
lake_perim_paths <- lake_perim_paths[!grepl(".xml", lake_perim_paths)]
depth_files <- depth_files[!grepl("FULLSCENE|V2", depth_files)]

# read rasters
depth_rasters <- lapply(depth_files, rast)
names(depth_rasters) <- basename(depth_files)
# read shp files
lake_perims <- lapply(lake_perim_paths, st_read)
names(lake_perims) <- basename(lake_perim_paths)

# selecting for lake perimeters on the same day
june22_perims <- lake_perims[grep("622", lake_perim_paths)]
july8_perims <- lake_perims[grep("708", lake_perim_paths)]
july31_perims <- lake_perims[grep("731", lake_perim_paths)]
aug23_perims <- lake_perims[grep("823", lake_perim_paths)]
sep1_perims <- lake_perims[grep("901", lake_perim_paths)]


################################################################################
# cropping rasters by all lake perims of respective day
crop_raster_by_perimeter <- function(depth_raster, lake_perimeter) {
  cropped_lake <- crop(depth_raster, vect(lake_perimeter))
  cropped_lake_df <- as.data.frame(cropped_lake, xy = TRUE)
  return(cropped_lake_df)
}


cropped_lakes_june22 <- vector("list", length(june22_perims)) 
cropped_lakes_july8 <- vector("list", length(july8_perims)) 
cropped_lakes_july31 <- vector("list", length(july31_perims)) 
cropped_lakes_aug23 <- vector("list", length(aug23_perims)) 
cropped_lakes_sep1 <- vector("list", length(sep1_perims)) 
for (i in seq_along(june22_perims)) {
  cropped_lakes_june22[[i]] <- crop_raster_by_perimeter(depth_rasters[[1]], june22_perims[[i]])
}
for (i in seq_along(july8_perims)) {
  cropped_lakes_july8[[i]] <- crop_raster_by_perimeter(depth_rasters[[2]], july8_perims[[i]])
}
for (i in seq_along(july31_perims)) {
  cropped_lakes_july31[[i]] <- crop_raster_by_perimeter(depth_rasters[[4]], july31_perims[[i]])
}
for (i in seq_along(aug23_perims)) {
  cropped_lakes_aug23[[i]] <- crop_raster_by_perimeter(depth_rasters[[7]], aug23_perims[[i]])
}
for (i in seq_along(sep1_perims)) {
  cropped_lakes_sep1[[i]] <- crop_raster_by_perimeter(depth_rasters[[5]], sep1_perims[[i]])
}

names(cropped_lakes_june22) <- names(june22_perims)
names(cropped_lakes_july8) <- names(july8_perims)
names(cropped_lakes_july31) <- names(july31_perims)
names(cropped_lakes_aug23) <- names(aug23_perims)
names(cropped_lakes_sep1) <- names(sep1_perims)


################################################################################
# function to calculate volume
lake_volume_each_lake <- function(df) {
  df$Pix.Vol <- df$Depth * 900 #900sqm and depth (m) = m3
  Lake.Vol <- sum(df$Pix.Vol, na.rm = TRUE)
}
lake_volumes_june22 <- lapply(cropped_lakes_june22, lake_volume_each_lake)
lake_volumes_july8 <- lapply(cropped_lakes_july8, lake_volume_each_lake) #lake A in Arc says 729900 // seems to track
lake_volumes_july31 <- lapply(cropped_lakes_july31, lake_volume_each_lake)
lake_volumes_aug23 <- lapply(cropped_lakes_aug23, lake_volume_each_lake)
lake_volumes_sep1 <- lapply(cropped_lakes_sep1, lake_volume_each_lake)

# rename
names(lake_volumes_june22) <- names(cropped_lakes_june22)
names(lake_volumes_july8) <- names(cropped_lakes_july8)
names(lake_volumes_july31) <- names(cropped_lakes_july31)
names(lake_volumes_aug23) <- names(cropped_lakes_aug23)
names(lake_volumes_sep1) <- names(cropped_lakes_sep1)

# printing
lake_volumes_june22 <- as.data.frame(unlist(lake_volumes_june22))
lake_volumes_july8 <- as.data.frame(unlist(lake_volumes_july8))
lake_volumes_july31 <- as.data.frame(unlist(lake_volumes_july31))
lake_volumes_aug23 <- as.data.frame(unlist(lake_volumes_aug23))
lake_volumes_sep1 <- as.data.frame(unlist(lake_volumes_sep1))
colnames(lake_volumes_june22)[1] <- "Volume"
colnames(lake_volumes_july8)[1] <- "Volume"
colnames(lake_volumes_july31)[1] <- "Volume"
colnames(lake_volumes_aug23)[1] <- "Volume"
colnames(lake_volumes_sep1)[1] <- "Volume"
lake_volumes_june22$Date <- "2021-06-22"
lake_volumes_july8$Date <- "2021-07-08"
lake_volumes_july31$Date <- "2021-07-31"
lake_volumes_aug23$Date <- "2021-08-23"
lake_volumes_sep1$Date <- "2021-09-01"
all_volumes_df <- bind_rows(lake_volumes_june22, lake_volumes_july8, lake_volumes_july31, lake_volumes_aug23, lake_volumes_sep1)
all_volumes_df$Volume <- format(all_volumes_df$Volume, scientific = FALSE)
write.csv(all_volumes_df, "E:/rsl_owu/SGL/data/output/csv/pope_daily_lake_vol.csv")




