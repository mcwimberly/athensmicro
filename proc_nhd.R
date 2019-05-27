library(raster)
library(sf)

ga_line1 <- st_read(dsn = "./nhd/Shape/NHDFlowline.shp")
ga_line2 <- st_read(dsn = "./nhd/Shape/NHDFlowline2.shp")
ga_line3 <- st_read(dsn = "./nhd/Shape/NHDFlowline3.shp")

ga_wbody <- st_read(dsn = "./nhd/Shape/NHDWaterbody.shp")

ga_line <- rbind(ga_line1, ga_line2, ga_line3)

athens_rast <- raster("treecov30.tif")

ga_line_reproj <- st_transform(ga_line, st_crs(athens_rast))
ga_line_reproj <- st_zm(ga_line_reproj)
ga_line_crop <- st_crop(ga_line_reproj, athens_rast)
ga_line_rast <- rasterize(as_Spatial(ga_line_crop), athens_rast, field = "InNetwork")
writeRaster(ga_line_rast, "ga_line_rast.tif", overwrite =T)

ga_wbody_reproj <- st_transform(ga_wbody, st_crs(athens_rast))
ga_wbody_reproj <- st_zm(ga_wbody_reproj)
ga_wbody_crop <- st_crop(ga_wbody_reproj, athens_rast)
ga_wbody_rast <- rasterize(as_Spatial(ga_wbody_crop), athens_rast, field = "Resolution")
writeRaster(ga_wbody_rast, "ga_wbody_rast.tif", overwrite=T)

reclassify(r, cbind(NA, NA, 0), right=FALSE)

ga_water_rast <- reclassify(ga_line_rast, cbind(NA, NA, 0), right=FALSE) + reclassify(ga_wbody_rast, cbind(NA, NA, 0), right=FALSE)
writeRaster(ga_water_rast, "ga_water_rast.tif", overwrite=T)
ga_water_rast <- raster("ga_water_rast.tif")


ga_water_utm <- projectRaster(ga_water_rast, crs = CRS("+proj=utm +zone=17 +datum=WGS84"), method="ngb")

ga_water_utm1 <- ga_water_utm
ga_water_utm1[ga_water_utm1 == 0] <- NA
ga_water_utmdist <- distance(ga_water_utm1)
ga_water_dist <- projectRaster(ga_water_utmdist, ga_water_rast)
writeRaster(ga_water_dist, "ga_water_dist.tif", overwrite=T)

ga_water_utm2 <- ga_water_utm
ga_water_utm2[ga_water_utm2 > 0] <- 1
fwts500 <- focalWeight(ga_water_utm2, d=500, "circle")
fwts500 <- ifelse(fwts500 > 0, 1, 0)
fwts1000 <- focalWeight(ga_water_utm2, d=1000, "circle")
fwts1000 <- ifelse(fwts1000 > 0, 1, 0)
fwts1500 <- focalWeight(ga_water_utm2, d=1500, "circle")
fwts1500 <- ifelse(fwts1500 > 0, 1, 0)
fwts2000 <- focalWeight(ga_water_utm2, d=2000, "circle")
fwts2000 <- ifelse(fwts2000 > 0, 1, 0)
ga_water_utm500m <- focal(ga_water_utm2, w=fwts500, fun=mean, na.rm = TRUE, pad = TRUE)
ga_water_utm1000m <- focal(ga_water_utm2, w=fwts1000, fun=mean, na.rm = TRUE, pad = TRUE)
ga_water_utm1500m <- focal(ga_water_utm2, w=fwts1500, fun=mean, na.rm = TRUE, pad = TRUE)
ga_water_utm2000m <- focal(ga_water_utm2, w=fwts2000, fun=mean, na.rm = TRUE, pad = TRUE)

ga_water_500m <- projectRaster(ga_water_utm500m, ga_water_rast)
ga_water_1000m <- projectRaster(ga_water_utm1000m, ga_water_rast)
ga_water_1500m <- projectRaster(ga_water_utm1500m, ga_water_rast)
ga_water_2000m <- projectRaster(ga_water_utm2000m, ga_water_rast)


writeRaster(ga_water_500m, "ga_water_500m.tif", overwrite=T)
writeRaster(ga_water_1000m, "ga_water_1000m.tif", overwrite=T)
writeRaster(ga_water_1500m, "ga_water_1500m.tif", overwrite=T)
writeRaster(ga_water_2000m, "ga_water_2000m.tif", overwrite=T)




