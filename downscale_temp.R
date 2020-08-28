#################################################################################
# downscale_temp.r
# Generate downscaled microclimate maps using GridMET meterological grids, 
# microclimate measurements, and land cover data
# written by Mike Wimberly, University of Oklahoma
#################################################################################

library(tidyverse)
library(lubridate)
library(raster)
library(rgdal)
library(sp)
library(snow)
library(mgcv)
library(splines)
library(leaps)
library(caret)

# Read in data files
# Daily microclimate summaries 
micro_daily <- read_csv("./outdata/micro_2018_daily.csv")
# Microclimate logger locations
loggers <- read_csv("c:/Users/wimb0002/Google Drive/AthensData/Logger_latlon.csv")
coordinates(loggers) <- c("lon", "lat")

# Set raster options - increase memory and chunk sizes
rasterOptions(maxmemory = 1e+10)
rasterOptions(chunksize = 1e+9)

# Read in land cover data - treecover, impervious surface, and ndwi
treecov <- raster("c:/Users/wimb0002/Google Drive/AthensData/Athens_treecover.tif")
imperv <- raster("c:/Users/wimb0002/Google Drive/AthensData/Athens_impervious.tif")
herb <- raster("c:/Users/wimb0002/Google Drive/AthensData/Athens_herbaceous.tif")
ndwi <- raster("c:/Users/wimb0002/Google Drive/AthensData/Athens_NDWI.tif")

# Set tree cover and impervious to zero for water pixels
treecov[ndwi > 0] <- 0
imperv[ndwi > 0] <- 0
herb[ndwi > 0] <- 0

# Scale up original 10 m resolution raster to 30 m
treecov30 <- aggregate(treecov, fact=3)
imperv30 <- aggregate(imperv, fact=3)
herb30 <- aggregate(herb, fact=3)

writeRaster(treecov30, "./outdata/treecov30.tif", overwrite=T)
writeRaster(imperv30, "./outdata/imperv30.tif", overwrite=T)
writeRaster(herb30, "./outdata/herb30.tif", overwrite=T)

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

# 3 x 3 and 5 x 5 for herbaceous
herb3 <- focal(herb30, w=matrix(1/9, nc=3, nr=3))
herb5 <- focal(herb30, w=matrix(1/25, nc=5, nr=5))

# Larger 1 and 2 km windows for herbaceous
herb500 <- focal(herb30, focalWeight(herb30, 0.005, "circle"))
herb1000 <- focal(herb30, focalWeight(herb30, 0.01, "circle"))
herb1500 <- focal(herb30, focalWeight(herb30, 0.015, "circle"))
herb2000 <- focal(herb30, focalWeight(herb30, 0.02, "circle"))

# Extract land cover data from rasters at logger locations
rsbrick <- brick(treecov30, forest3, forest5,
                 imperv30, imperv3, imperv5,
                 herb30, herb3, herb5,
                 imperv500, imperv1000, imperv1500, imperv2000,
                 forest500, forest1000, forest1500, forest2000,
                 herb500, herb1000, herb1500, herb2000)
names(rsbrick) <- c("treecov30", "forest3", "forest5",
                    "imperv30", "imperv3", "imperv5",
                    "herb30", "herb3", "herb5",
                    "imperv500", "imperv1000", "imperv1500", "imperv2000",
                    "forest500", "forest1000", "forest1500", "forest2000",
                    "herb500", "herb1000", "herb1500", "herb2000")
rsextract <- raster::extract(rsbrick, loggers, df = TRUE)
rstable <- rsextract %>%
  mutate(Site = loggers$name)

# Join the extracted data to the microclimate data
micro_tab <- micro_daily %>%
  left_join(rstable, by = "Site") %>%
  filter(num_temp == 144) #%>%

# Read in daily meterological grids
tminstack <- stack("./gridmet_old/tmmn_2018.nc")
tmaxstack <- stack("./gridmet_old/tmmx_2018.nc")
rminstack <- stack("./gridmet_old/rmin_2018.nc")
rmaxstack <- stack("./gridmet_old/rmax_2018.nc")
sradstack <- stack("./gridmet_old/srad_2018.nc")
prstack <- stack("./gridmet_old/pr_2018.nc")
vsstack <- stack("./gridmet_old/vs_2018.nc")

# Generate a larger grid extent for the meteorological data
# Do this to avoid edge effects when resampling the met data to 30m
big30mgrid <- extend(treecov30, c(150, 150))
athense <- extent(big30mgrid)

# Generate a stack for each meterological variable using the following steps:
# 1. Crop to expanded Athens area
# 2. Extract the desired dates (June 15 - October 10 )
# 3. Resample to the 30 m grid
# 4. Crop to the original extent 

beginCluster(8)

tmincrop <- crop(tminstack, athense, filename = "./outdata/tmincrop.tif", overwrite = T)
tmincrop <- tmincrop[[166: 283]]
tmincrop2 <- resample(tmincrop, big30mgrid, method = "bilinear", filename = "./outdata/tmincrop2.tif", overwrite = T)
tmincrop3 <- crop(tmincrop2, extent(treecov30), filename = "./outdata/tmincrop3.tif", overwrite = T)
tmincrop3 <- tmincrop3 - 273.15

tmaxcrop <- crop(tmaxstack, athense, filename = "./outdata/tmaxcrop.tif", overwrite = T)
tmaxcrop <- tmaxcrop[[166: 283]]
tmaxcrop2 <- resample(tmaxcrop, big30mgrid, method = "bilinear", filename = "./outdata/tmaxcrop2.tif", overwrite = T)
tmaxcrop3 <- crop(tmaxcrop2, extent(treecov30), filename = "./outdata/tmaxcrop3.tif", overwrite = T)
tmaxcrop3 <- tmaxcrop3 - 273.15

rmincrop <- crop(rminstack, athense, filename = "./outdata/rmincrop.tif", overwrite = T)
rmincrop <- rmincrop[[166: 283]]
rmincrop2 <- resample(rmincrop, big30mgrid, method = "bilinear", filename = "./outdata/rmincrop2.tif", overwrite = T)
rmincrop3 <- crop(rmincrop2, extent(treecov30), filename = "./outdata/rmincrop3.tif", overwrite = T)

rmaxcrop <- crop(rmaxstack, athense, filename = "./outdata/rmaxcrop.tif", overwrite = T)
rmaxcrop <- rmaxcrop[[166: 283]]
rmaxcrop2 <- resample(rmaxcrop, big30mgrid, method = "bilinear", filename = "./outdata/rmaxcrop2.tif", overwrite = T)
rmaxcrop3 <- crop(rmaxcrop2, extent(treecov30), filename = "./outdata/rmaxcrop3.tif", overwrite = T)

sradcrop <- crop(sradstack, athense, filename = "./outdata/sradcrop.tif", overwrite = T)
sradcrop <- sradcrop[[166: 283]]
sradcrop2 <- resample(sradcrop, big30mgrid, method = "bilinear", filename = "./outdata/sradcrop2.tif", overwrite = T)
sradcrop3 <- crop(sradcrop2, extent(treecov30), filename = "./outdata/sradcrop3.tif", overwrite = T)

prcrop <- crop(prstack, athense, filename = "./outdata/prcrop.tif", overwrite = T)
prcrop <- prcrop[[166: 283]]
prcrop2 <- resample(prcrop, big30mgrid, method = "bilinear", filename = "./outdata/prcrop2.tif", overwrite = T)
prcrop3 <- crop(prcrop2, extent(treecov30), filename = "./outdata/prcrop3.tif", overwrite = T)

vscrop <- crop(vsstack, athense, filename = "./outdata/vscrop.tif", overwrite = T)
vscrop <- vscrop[[166: 283]]
vscrop2 <- resample(vscrop, big30mgrid, method = "bilinear", filename = "./outdata/vscrop2.tif", overwrite = T)
vscrop3 <- crop(vscrop2, extent(treecov30), filename = "./outdata/vscrop3.tif", overwrite = T)

endCluster()

# Read in saved TIF files to create raster brick objects
###############################################################################
tmincrop3 <- brick("./outdata/tmincrop3.tif")
tmaxcrop3 <- brick("./outdata/tmaxcrop3.tif")
rmincrop3 <- brick("./outdata/rmincrop3.tif")
rmaxcrop3 <- brick("./outdata/rmaxcrop3.tif")
sradcrop3 <- brick("./outdata/sradcrop3.tif")
prcrop3 <- brick('./outdata/prcrop3.tif')
vscrop3 <- brick("./outdata/vscrop3.tif")
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

# Convert temperatures from Kevin to Celcius
tmincrop3 <- tmincrop3 - 273.15
tmaxcrop3 <- tmaxcrop3 - 273.15

# Extract the gridded meteorological data at the logger locations
tmin_ex <- raster::extract(tmincrop3, loggers)
tmax_ex <- raster::extract(tmaxcrop3, loggers)
rmin_ex <- raster::extract(rmincrop3, loggers)
rmax_ex <- raster::extract(rmaxcrop3, loggers)
srad_ex <- raster::extract(sradcrop3, loggers)
pr_ex <- raster::extract(prcrop3, loggers)
vs_ex <- raster::extract(vscrop3, loggers)

# Add column with logger IDs
tmin_ex <- data.frame(Site = loggers$name, tmin_ex)
tmax_ex <- data.frame(Site = loggers$name, tmax_ex)
rmin_ex <- data.frame(Site = loggers$name, rmin_ex)
rmax_ex <- data.frame(Site = loggers$name, rmax_ex)
srad_ex <- data.frame(Site = loggers$name, srad_ex)
pr_ex <- data.frame(Site = loggers$name, pr_ex)
vs_ex <- data.frame(Site = loggers$name, vs_ex)

# Convert GridMET data to long format
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

# Join evertyhing into a single table
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

write_csv(fitdata, "./outdata/fitdata.csv")
fitdata <- read_csv("./outdata/fitdata.csv")

# Use best subsets regression to guide initial variable selection
mintemp_sub <- regsubsets(min_temp ~ tmin + treecov30 + forest3 + forest5 + imperv30 + imperv3 + imperv5 + 
                            imperv500 + imperv1000 + imperv1500 + imperv2000 +
                            rmin + rmax + srad + log(pr + 1) + vs + 
                            forest500 + forest1000 + forest1500 + forest2000, data = fitdata)

maxtemp_sub <- regsubsets(max_temp ~ tmax + treecov30 + forest3 + forest5 + imperv30 + imperv3 + imperv5 + 
                            imperv500 + imperv1000 + imperv1500 + imperv2000 +
                            rmin + rmax + srad + log(pr + 1) + vs + 
                            forest500 + forest1000 + forest1500 + forest2000, data = fitdata)

minrh_sub <- regsubsets(min_rh ~  tmax + treecov30 + forest3 + forest5 + imperv30 + imperv3 + imperv5 + 
                          imperv500 + imperv1000 + imperv1500 + imperv2000 +
                          rmin + srad + log(pr + 1) + vs + 
                          forest500 + forest1000 + forest1500 + forest2000, data = fitdata)

meanrh_sub <- regsubsets(mean_rh ~ tmin + tmax + treecov30 + forest3 + forest5 + imperv30 + imperv3 + imperv5 + 
                           imperv500 + imperv1000 + imperv1500 + imperv2000 +
                           rmin + rmax + srad + log(pr + 1) + vs + 
                           forest500 + forest1000 + forest1500 + forest2000, data = fitdata)


# Fit final models
dscale_minrh <- lm(min_rh ~ rmin + imperv1000 + treecov30 + tmin + log(pr + 1), data = fitdata)

dscale_min <- lm(min_temp ~ tmin + forest5 + imperv1000 + rmax, data = fitdata)

dscale_max <- lm(max_temp ~ tmax + treecov30 + srad + rmin, data = fitdata)

# Examine distribution of residuals by site
boxplot(resid(dscale_min) ~ fitdata$Site)
boxplot(resid(dscale_max) ~ fitdata$Site)

# Cross-validate the final models
logger_sites <- groupKFold(fitdata$Site)
train_control <- trainControl(method="cv", number = 51, index= logger_sites)

tmin_cv <- train(min_temp ~ tmin + forest5 + imperv1000 + rmax, data=fitdata, trControl=train_control, method="lm")
tmax_cv <- train(max_temp ~ tmax + treecov30 + srad + rmin, data=fitdata, trControl=train_control, method="lm")

# Apply model to downscale minimum temperature
beginCluster(8)
tempminstack <- brick()
for(x in 1:nlayers(tmincrop3)) {
  print(x)
  print(Sys.time())
  tempras <- stack(tmincrop3[[x]], forest5, imperv1000, rmaxcrop3[[x]])
  names(tempras) <- c("tmin", "forest5", "imperv1000", "rmax")
  tempout <- clusterR(tempras, predict, args=list(dscale_min))
  tempminstack <-addLayer(tempminstack, tempout)
}
endCluster()

# Apply model to downscale maximum temperature
beginCluster(8)
tempmaxstack <- brick()
for(x in 1:nlayers(tmaxcrop3)) {
  print(x)
  print(Sys.time())
  tempras <- stack(tmaxcrop3[[x]], treecov30, sradcrop3[[x]], rmincrop3[[x]])
  names(tempras) <- c("tmax", "treecov30", "srad", "rmin")
  tempout <- clusterR(tempras, predict, args=list(dscale_max))
  tempmaxstack <-addLayer(tempmaxstack, tempout)
}
endCluster()

# Apply model to downscale minimum relative humidity
beginCluster(8)
rhminstack <- brick()
for(x in 1:nlayers(rmincrop3)) {
  print(x)
  print(Sys.time())
  tempras <- stack(rmincrop3[[x]], imperv1000, treecov30, tmincrop3[[x]], prcrop3[[x]])
  names(tempras) <- c("rmin", "imperv1000", "treecov30", "tmin", "pr")
  tempout <- clusterR(tempras, predict, args=list(dscale_minrh))
  rhminstack <-addLayer(rhminstack, tempout)
}
endCluster()

# Write the output grids
writeRaster(tempminstack, "./outdata/tempminstack.tif", overwrite=T)
writeRaster(tempmaxstack, "./outdata/tempmaxstack.tif", overwrite=T)
writeRaster(rhminstack, "./outdata/rhminstack.tif", overwrite=T)

#####################################################################
# Generate 7-day rolling averages of minimum and maximum temperature
#####################################################################

# First, we have to process the GridMET data for the six days at the beginning 
# of the time series
tmincrop <- crop(tminstack, athense)
tmincrop <- tmincrop[[160: 165]]
tmincrop2 <- resample(tmincrop, big30mgrid, method = "bilinear")
tmincrop4 <- crop(tmincrop2, extent(treecov30), filename = "./outdata/tmincrop4.tif", overwrite = T)
tmincrop4 <- tmincrop4 - 273.15

tmaxcrop <- crop(tmaxstack, athense)
tmaxcrop <- tmaxcrop[[160: 165]]
tmaxcrop2 <- resample(tmaxcrop, big30mgrid, method = "bilinear")
tmaxcrop4 <- crop(tmaxcrop2, extent(treecov30), filename = "./outdata/tmaxcrop4.tif", overwrite = T)
tmaxcrop4 <- tmaxcrop4 - 273.15

rmincrop <- crop(rminstack, athense)
rmincrop <- rmincrop[[160: 165]]
rmincrop2 <- resample(rmincrop, big30mgrid, method = "bilinear")
rmincrop4 <- crop(rmincrop2, extent(treecov30), filename = "./outdata/rmincrop4.tif", overwrite = T)

rmaxcrop <- crop(rmaxstack, athense)
rmaxcrop <- rmaxcrop[[160: 165]]
rmaxcrop2 <- resample(rmaxcrop, big30mgrid, method = "bilinear")
rmaxcrop4 <- crop(rmaxcrop2, extent(treecov30), filename = "./outdata/rmaxcrop4.tif", overwrite = T)

sradcrop <- crop(sradstack, athense)
sradcrop <- sradcrop[[160: 165]]
sradcrop2 <- resample(sradcrop, big30mgrid, method = "bilinear")
sradcrop4 <- crop(sradcrop2, extent(treecov30), filename = "./outdata/sradcrop4.tif", overwrite = T)

# Generate a new minimum temperature stack
beginCluster(8)
tempminstack2 <- brick()
for(x in 1:nlayers(tmincrop4)) {
  print(x)
  print(Sys.time())
  tempras <- stack(tmincrop4[[x]], forest5, imperv1000, rmaxcrop4[[x]])
  names(tempras) <- c("tmin", "forest5", "imperv1000", "rmax")
  tempout <- clusterR(tempras, predict, args=list(dscale_min))
  tempminstack2 <-addLayer(tempminstack2, tempout)
}
endCluster()

# Generate a new maximum temeprature stack
beginCluster(8)
tempmaxstack2 <- brick()
for(x in 1:nlayers(tmaxcrop4)) {
  print(x)
  print(Sys.time())
  tempras <- stack(tmaxcrop4[[x]], treecov30, sradcrop4[[x]], rmincrop4[[x]])
  names(tempras) <- c("tmax", "treecov30", "srad", "rmin")
  tempout <- clusterR(tempras, predict, args=list(dscale_max))
  tempmaxstack2 <-addLayer(tempmaxstack2, tempout)
}
endCluster()

# Add new data to the front of the min and max temperature stacks
tempminstack2 <- stack(tempminstack2, tempminstack)
tempmaxstack2 <- stack(tempmaxstack2, tempmaxstack)

# Generate a new stack of 7-day rolling averages
tempminstack_7d <- calc(tempminstack2, function(x) movingFun(x, n=7, fun=mean, type='to'))
tempmaxstack_7d <- calc(tempmaxstack2, function(x) movingFun(x, n=7, fun=mean, type='to'))

# Strip off the first six dates in the rolling average stack, which should all be NODATA
tempminstack_7d <- tempminstack_7d[[-(1:6)]]
tempmaxstack_7d <- tempmaxstack_7d[[-(1:6)]]

# Write seven-day moving average stacks to disk
writeRaster(tempminstack_7d, "./outdata/tempminstack_7d.tif", overwrite=T)
writeRaster(tempmaxstack_7d, "./outdata/tempmaxstack_7d.tif", overwrite=T)

rm(list = ls())

