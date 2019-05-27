library(tidyverse)
library(lubridate)
library(raster)
library(rgdal)
library(sp)
library(snow)

# Read in data files
# Daily microclimate summaries 
micro_daily <- read_csv("micro_2018_daily.csv")
# Microclimate logger locations
loggers <- read_csv("Logger_latlon.csv")
coordinates(loggers) <- c("lon", "lat")

# Set raster options - increase memory and chunk sizes
rasterOptions(maxmemory = 1e+10)
rasterOptions(chunksize = 1e+9)

# Read in land cover data - treecover, imprevious surface, and ndwi
treecov <- raster("Athens_treecover_v2.tif")
imperv <- raster("Athens_impervious_v2.tif")
ndwi <- raster("Athens_NDWI.tif")
waterdist <- raster("ga_water_dist.tif")
water500 <- raster("ga_water_500m.tif")
water1000 <- raster("ga_water_1000m.tif")
water1500 <- raster("ga_water_1500m.tif")
water2000 <- raster("ga_water_2000m.tif")

# Set tree cover and impervious to zero for water pixels
treecov[ndwi > 0] <- 0
imperv[ndwi > 0] <- 0

# Scale up original 10 m resolution raster to 30 m
treecov30 <- aggregate(treecov, fact=3)
imperv30 <- aggregate(imperv, fact=3)

writeRaster(treecov30, "treecov30.tif", overwrite=T)
writeRaster(imperv30, "imperv30.tif", overwrite=T)

# Generate focal means
# 3 x 3 and 5 x 5 windows for tree cover
forest3 <- focal(treecov30, w=matrix(1/9, nc=3, nr=3))
forest5 <- focal(treecov30, w=matrix(1/25, nc=5, nr=5))

# Larger 1 and 2 km windows for tree cover
forest500 <- focal(treecov30, focalWeight(treecov30, 0.005, "circle"))
forest1000 <- focal(treecov30, focalWeight(treecov30, 0.01, "circle"))
forest1500 <- focal(treecov30, focalWeight(treecov30, 0.015, "circle"))
forest2000 <- focal(treecov30, focalWeight(treecov30, 0.02, "circle"))

# 3 x 3 and 5 x 5 for impervious surface
imperv3 <- focal(imperv30, w=matrix(1/9, nc=3, nr=3))
imperv5 <- focal(imperv30, w=matrix(1/25, nc=5, nr=5))

# Larger 1 and 2 km windows for impervious
imperv500 <- focal(imperv30, focalWeight(imperv30, 0.005, "circle"))
imperv1000 <- focal(imperv30, focalWeight(imperv30, 0.01, "circle"))
imperv1500 <- focal(imperv30, focalWeight(imperv30, 0.015, "circle"))
imperv2000 <- focal(imperv30, focalWeight(imperv30, 0.02, "circle"))


# Extract land cover data from rasters at logger locations
rsbrick <- brick(treecov30, forest3, forest5,
                 imperv30, imperv3, imperv5,
                 imperv500, imperv1000, imperv1500, imperv2000,
                 forest500, forest1000, forest1500, forest2000,
                 waterdist, water500, water1000, water1500, water2000)
names(rsbrick) <- c("treecov30", "forest3", "forest5",
                    "imperv30", "imperv3", "imperv5",
                    "imperv500", "imperv1000", "imperv1500", "imperv2000",
                    "forest500", "forest1000", "forest1500", "forest2000",
                    "waterdist", "water500", "water1000", "water1500", "water2000")
rsextract <- raster::extract(rsbrick, loggers, df = TRUE)
rstable <- rsextract %>%
  mutate(Site = loggers$name)

# Join the extracted data to the microclimate data
micro_tab <- micro_daily %>%
  left_join(rstable, by = "Site") %>%
  filter(num_temp == 144) #%>%

# Read in daily meterological grids
tminstack <- stack("./gridmet/tmmn_2018.nc")
tmaxstack <- stack("./gridmet/tmmx_2018.nc")
rminstack <- stack("./gridmet/rmin_2018.nc")
rmaxstack <- stack("./gridmet/rmax_2018.nc")
sradstack <- stack("./gridmet/srad_2018.nc")
prstack <- stack("./gridmet/pr_2018.nc")
vsstack <- stack("./gridmet/vs_2018.nc")

# Generate a larger grid extend for the meteorological data
# Do this to avoid edge effects when resampling the met data to 30m
big30mgrid <- extend(treecov30, c(150, 150))
athense <- extent(big30mgrid)

# Generate a stack for each meterological variable using the following steps:
# 1. Crop to expanded Athens area
# 2. Extract the desired dates (June 15 - October 10 )
# 3. Resample to the 30 m grid
# 4. Crop to the original extent 

beginCluster(12)

tmincrop <- crop(tminstack, athense, filename = "tmincrop.tif", overwrite = T)
tmincrop <- tmincrop[[166: 283]]
tmincrop2 <- resample(tmincrop, big30mgrid, method = "bilinear", filename = "tmincrop2.tif", overwrite = T)
tmincrop3 <- crop(tmincrop2, extent(treecov30), filename = "tmincrop3.tif", overwrite = T)
tmincrop3 <- tmincrop3 - 273.15

tmaxcrop <- crop(tmaxstack, athense, filename = "tmaxcrop.tif", overwrite = T)
tmaxcrop <- tmaxcrop[[166: 283]]
tmaxcrop2 <- resample(tmaxcrop, big30mgrid, method = "bilinear", filename = "tmaxcrop2.tif", overwrite = T)
tmaxcrop3 <- crop(tmaxcrop2, extent(treecov30), filename = "tmaxcrop3.tif", overwrite = T)
tmaxcrop3 <- tmaxcrop3 - 273.15

rmincrop <- crop(rminstack, athense, filename = "rmincrop.tif", overwrite = T)
rmincrop <- rmincrop[[166: 283]]
rmincrop2 <- resample(rmincrop, big30mgrid, method = "bilinear", filename = "rmincrop2.tif", overwrite = T)
rmincrop3 <- crop(rmincrop2, extent(treecov30), filename = "rmincrop3.tif", overwrite = T)

rmaxcrop <- crop(rmaxstack, athense, filename = "rmaxcrop.tif", overwrite = T)
rmaxcrop <- rmaxcrop[[166: 283]]
rmaxcrop2 <- resample(rmaxcrop, big30mgrid, method = "bilinear", filename = "rmaxcrop2.tif", overwrite = T)
rmaxcrop3 <- crop(rmaxcrop2, extent(treecov30), filename = "rmaxcrop3.tif", overwrite = T)

sradcrop <- crop(sradstack, athense, filename = "sradcrop.tif", overwrite = T)
sradcrop <- sradcrop[[166: 283]]
sradcrop2 <- resample(sradcrop, big30mgrid, method = "bilinear", filename = "sradcrop2.tif", overwrite = T)
sradcrop3 <- crop(sradcrop2, extent(treecov30), filename = "sradcrop3.tif", overwrite = T)

prcrop <- crop(prstack, athense, filename = "prcrop.tif", overwrite = T)
prcrop <- prcrop[[166: 283]]
prcrop2 <- resample(prcrop, big30mgrid, method = "bilinear", filename = "prcrop2.tif", overwrite = T)
prcrop3 <- crop(prcrop2, extent(treecov30), filename = "prcrop3.tif", overwrite = T)

vscrop <- crop(vsstack, athense, filename = "vscrop.tif", overwrite = T)
vscrop <- vscrop[[166: 283]]
vscrop2 <- resample(vscrop, big30mgrid, method = "bilinear", filename = "vscrop2.tif", overwrite = T)
vscrop3 <- crop(vscrop2, extent(treecov30), filename = "vscrop3.tif", overwrite = T)

endCluster()


###############################################################################
tmincrop3 <- brick("tmincrop3.tif")
tmaxcrop3 <- brick("tmaxcrop3.tif")
rmincrop3 <- brick("rmincrop3.tif")
rmaxcrop3 <- brick("rmaxcrop3.tif")
sradcrop3 <- brick("sradcrop3.tif")
prcrop3 <- brick('prcrop3.tif')
vscrop3 <- brick("vscrop3.tif")
###############################################################################


# Assign a date name to each layer in the stacks
microdates <- seq(ymd('2018-06-15'), ymd('2018-10-10'), by = 'days')
names(tmincrop3) <- microdates
names(tmaxcrop3) <- microdates
names(rmincrop3) <- microdates
names(rmaxcrop3) <- microdates
names(sradcrop3) <- microdates
names(prcrop3) <- microdates
names(vscrop3) <- microdates

tmincrop3 <- tmincrop3 - 273.15
tmaxcrop3 <- tmaxcrop3 - 273.15

# Extract the gridded meterological data at the logger locations
tmin_ex <- raster::extract(tmincrop3, loggers)
tmax_ex <- raster::extract(tmaxcrop3, loggers)
rmin_ex <- raster::extract(rmincrop3, loggers)
rmax_ex <- raster::extract(rmaxcrop3, loggers)
srad_ex <- raster::extract(sradcrop3, loggers)
pr_ex <- raster::extract(prcrop3, loggers)
vs_ex <- raster::extract(vscrop3, loggers)

#colnames(tmin_ex) <- as.character(microdates)
#colnames(tmax_ex) <- as.character(microdates)
#colnames(rmin_ex) <- as.character(microdates)
#colnames(rmax_ex) <- as.character(microdates)
#colnames(srad_ex) <- as.character(microdates)
#colnames(pr_ex) <- as.character(microdates)
#colnames(vs_ex) <- as.character(microdates)

tmin_ex <- data.frame(Site = loggers$name, tmin_ex)
tmax_ex <- data.frame(Site = loggers$name, tmax_ex)
rmin_ex <- data.frame(Site = loggers$name, rmin_ex)
rmax_ex <- data.frame(Site = loggers$name, rmax_ex)
srad_ex <- data.frame(Site = loggers$name, srad_ex)
pr_ex <- data.frame(Site = loggers$name, pr_ex)
vs_ex <- data.frame(Site = loggers$name, vs_ex)

#names(tmin_ex)[1] <- "Site"
#names(tmax_ex)[1] <- "Site"
#names(rmin_ex)[1] <- "Site"
#names(rmax_ex)[1] <- "Site"
#names(srad_ex)[1] <- "Site"
#names(pr_ex)[1] <- "Site"
#names(vs_ex)[1] <- "Site"

tmin_ex2 <- tmin_ex %>%
  gather(Date, tmin, -Site) %>%
  mutate(Date = ymd(substr(Date, 2, 11))) 

tmax_ex2 <- tmax_ex %>%
  gather(Date, tmax, -Site) %>%
  mutate(Date = ymd(substr(Date, 2, 11))) 

rmin_ex2 <- rmin_ex %>%
  gather(Date, rmin, -Site) %>%
  mutate(Date = ymd(substr(Date, 2, 11))) 

rmax_ex2 <- rmax_ex %>%
  gather(Date, rmax, -Site) %>%
  mutate(Date = ymd(substr(Date, 2, 11))) 

srad_ex2 <- srad_ex %>%
  gather(Date, srad, -Site) %>%
  mutate(Date = ymd(substr(Date, 2, 11))) 

pr_ex2 <- pr_ex %>%
  gather(Date, pr, -Site) %>%
  mutate(Date = ymd(substr(Date, 2, 11))) 

vs_ex2 <- vs_ex %>%
  gather(Date, vs, -Site) %>%
  mutate(Date = ymd(substr(Date, 2, 11))) 

tall <- tmin_ex2 %>%
  left_join(micro_tab, by = (c("Site", "Date"))) %>%
  left_join(tmax_ex2, by = (c("Site", "Date"))) %>%
  left_join(rmin_ex2, by = (c("Site", "Date"))) %>%
  left_join(rmax_ex2, by = (c("Site", "Date"))) %>%
  left_join(srad_ex2, by = (c("Site", "Date"))) %>%
  left_join(pr_ex2, by = (c("Site", "Date"))) %>%
  left_join(vs_ex2, by = (c("Site", "Date")))

# Only include data through October 10 in the model fit
fitdata <- filter(tall, !is.na(mean_temp) & !is.na(min_temp) & !is.na(max_temp))

write_csv(fitdata, "fitdata.csv")
fitdata <- read_csv("fitdata.csv")



library(leaps)
mintemp_sub <- regsubsets(min_temp ~ tmin + treecov30 + forest3 + forest5 + imperv30 + imperv3 + imperv5 + 
                            imperv500 + imperv1000 + imperv1500 + imperv2000 +
                            rmin + rmax + srad + log(pr + 1) + vs + 
                            forest500 + forest1000 + forest1500 + forest2000 +
                            sqrt(waterdist) + water500 + water1000 + water1500 + water2000, data = fitdata)

maxtemp_sub <- regsubsets(max_temp ~ tmax + treecov30 + forest3 + forest5 + imperv30 + imperv3 + imperv5 + 
                            imperv500 + imperv1000 + imperv1500 + imperv2000 +
                            rmin + rmax + srad + log(pr + 1) + vs + 
                            forest500 + forest1000 + forest1500 + forest2000 +
                            sqrt(waterdist) + water500 + water1000 + water1500 + water2000, data = fitdata)

minrh_sub <- regsubsets(min_rh ~  tmax + treecov30 + forest3 + forest5 + imperv30 + imperv3 + imperv5 + 
                          imperv500 + imperv1000 + imperv1500 + imperv2000 +
                          rmin + srad + log(pr + 1) + vs + 
                          forest500 + forest1000 + forest1500 + forest2000 +
                          sqrt(waterdist) + water500 + water1000 + water1500 + water2000, data = fitdata)

meanrh_sub <- regsubsets(mean_rh ~ tmin + tmax + treecov30 + forest3 + forest5 + imperv30 + imperv3 + imperv5 + 
                           imperv500 + imperv1000 + imperv1500 + imperv2000 +
                           rmin + rmax + srad + log(pr + 1) + vs + 
                           forest500 + forest1000 + forest1500 + forest2000 +
                           sqrt(waterdist) + water500 + water1000 + water1500 + water2000, data = fitdata)



dscale_minrh <- lm(min_rh ~ rmin + imperv1000 + treecov30 + tmin + log(pr + 1), data = fitdata)

dscale_meanrh <- lm(mean_rh ~ rmin + imperv1000 + log(pr + 1) + sqrt(waterdist), data = fitdata)

dscale_min <- lm(min_temp ~ tmin + forest5 + imperv1000 + rmax, data = fitdata)

dscale_max <- lm(max_temp ~ tmax + treecov30 + srad + rmin, data = fitdata)


mean(abs(fitdata$min_temp - predict(dscale_min)))
mean(abs(fitdata$min_temp - fitdata$tmin))
mean(fitdata$min_temp - predict(dscale_min))

cor(fitdata$min_temp, predict(dscale_min))
cor(fitdata$min_temp, fitdata$tmin)

mean(abs(fitdata$max_temp - predict(dscale_max)))
mean(abs(fitdata$max_temp - fitdata$tmax))
mean(fitdata$max_temp - predict(dscale_max))

cor(fitdata$max_temp, predict(dscale_max))
cor(fitdata$max_temp, fitdata$tmax)

mean(abs(fitdata$min_rh - predict(dscale_minrh)))
mean(abs(fitdata$min_rh - fitdata$rmin))
mean(fitdata$min_rh - predict(dscale_minrh))

cor(fitdata$min_rh, predict(dscale_minrh))
cor(fitdata$min_rh, fitdata$rmin)


water1000_v2 <- water1000
water1000_v2[water1000_v2 > max(fitdata$water1000)] <- max(fitdata$water1000)

water2000_v2 <- water2000
water2000_v2[water2000_v2 > max(fitdata$water2000)] <- max(fitdata$water2000)

beginCluster(8)
tempminstack <- brick()
for(x in 1:nlayers(tmincrop3)) {
  print(x)
  print(Sys.time())
  tempras <- stack(tmincrop3[[x]], forest5, imperv1000, rmaxcrop3[[x]])
  names(tempras) <- c("tmin", "forest5", "imperv1000", "rmax")
  #tempout <- predict(tempras, dscale_lm)
  tempout <- clusterR(tempras, predict, args=list(dscale_min))
  tempminstack <-addLayer(tempminstack, tempout)
}
endCluster()

beginCluster(8)
tempmaxstack <- brick()
for(x in 1:nlayers(tmaxcrop3)) {
  print(x)
  print(Sys.time())
  tempras <- stack(tmaxcrop3[[x]], treecov30, sradcrop3[[x]], rmincrop3[[x]])
  names(tempras) <- c("tmax", "treecov30", "srad", "rmin")
  #tempout <- predict(tempras, dscale_lm)
  tempout <- clusterR(tempras, predict, args=list(dscale_max))
  tempmaxstack <-addLayer(tempmaxstack, tempout)
}
endCluster()

beginCluster(8)
rhminstack <- brick()
for(x in 1:nlayers(rmincrop3)) {
  print(x)
  print(Sys.time())
  tempras <- stack(rmincrop3[[x]], imperv1000, treecov30, tmincrop3[[x]], prcrop3[[x]])
  names(tempras) <- c("rmin", "imperv1000", "treecov30", "tmin", "pr")
  #tempout <- predict(tempras, dscale_lm)
  tempout <- clusterR(tempras, predict, args=list(dscale_minrh))
  rhminstack <-addLayer(rhminstack, tempout)
}
endCluster()



writeRaster(tempminstack, "tempminstack.tif", overwrite=T)
writeRaster(tempmaxstack, "tempmaxstack.tif", overwrite=T)
writeRaster(rhminstack, "rhminstack.tif", overwrite=T)


rm(list = ls())

