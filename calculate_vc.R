#################################################################################
# calculate_vc.r
# Predict vectorial capacity based on downscaled temperatures
# written by Mike Wimberly, University of Oklahoma
#################################################################################

library(tidyverse)
library(lubridate)
library(raster)
library(rgdal)
library(sp)
library(snow)
library(RcppRoll)

# Set raster options - increase memory and chunk sizes
rasterOptions(maxmemory = 1e+10)
rasterOptions(chunksize = 1e+9)

# Function for temperature-trait relationships using a quadratic curve
quadratic <- function(temp) {
  resp <- c1 * (temp - t0) * (temp - tm)
  resp[temp >= tm] <- 0
  resp[temp <= t0] <- 0
  resp
}

# function for temperature-trait relationships using a Briere durve
briere <- function(temp) {
  resp <-c1 * temp * (temp - t0) * (tm - temp)^(1/2)
  resp[temp >= tm] <- 0
  resp[temp <= t0] <- 0
  resp
}

# Simple function to compute hourly temperature based on a sigmoidal surve
hourly_temp <- function(minmax) {
  htemp <- (minmax[[2]] - minmax[[1]])/2 * sin(2 * pi * hour/24) + (minmax[[1]] + minmax[[2]])/2
  htemp
}

# Function to compute sunrise, sunset, and day length for a given day of year
calcsun <- function(doy) {
  Gamma <- 2 * pi/365 * ((doy) - 1)
  Delta <- 180/pi * (0.006918 - 0.399912 * cos(Gamma) + 0.070257 * 
                       sin(Gamma) - 0.006758 * cos(Gamma) + 0.000907 * sin(Gamma) - 
                       0.002697 * cos(3 * (Gamma)) + 0.00148 * sin(3 * (Gamma)))
  CosWo <- (sin(-0.8333/360 * 2 * pi) - sin(latitude/360 * 2 * pi) * 
              sin(Delta/360 * 2 * pi))/(cos(latitude/360 * 
                                              2 * pi) * cos(Delta/360 * 2 * pi))
  Sunrise <- 12 - acos(CosWo)/(15/360 * 2 * pi)
  Sunset <- 12 + acos(CosWo)/(15/360 * 2 * pi)
  Daylength <- Sunset - Sunrise
  results <- c(Sunrise, Sunset, Daylength)
}

# Function to compute hourly temperatures
hourly_temp2 <- function(tstack) {
  Tmin <- tstack[1]
  Tmax <- tstack[2]
  prev_Tmin <- tstack[3]
  prev_Tmax <- tstack[4]
  next_Tmin <- tstack[5]
  
  cur_Tsunset <- Tmin + (Tmax - Tmin) * sin(pi * (cur_sun[2] - cur_sun[1]) / (cur_sun[3] + 4 ))
  prev_Tsunset <- prev_Tmin + (prev_Tmax - prev_Tmin) * sin(pi * (prev_sun[2] - prev_sun[1]) / (prev_sun[3] + 4 ))
  
  if(hour <= cur_sun[1]) {
    Temp <- prev_Tsunset - (prev_Tsunset - Tmin) / log(24 - (prev_sun[2] - cur_sun[1]) + 1) * log(hour + 24 - prev_sun[2] + 1)
  } else if (hour > cur_sun[1] & hour <= cur_sun[2]) {
    Temp <- Tmin + (Tmax - Tmin) * sin(pi * (hour - cur_sun[1]) / (cur_sun[3] + 4 ))
  } else if (hour > cur_sun[2]) {
    Temp <- cur_Tsunset - (cur_Tsunset - next_Tmin) / log(24 - (cur_sun[2] - next_sun[1]) + 1) * log(hour - cur_sun[2] + 1)
  }
  Temp
}

###############################
# Generate four monthly VC maps
###############################

# Read temperature data stacks
tempminstack <- brick("./outdata/tempminstack.tif")
tempmaxstack <- brick("./outdata/tempmaxstack.tif")
tempmaxstack_7d <- brick("./outdata/tempmaxstack_7d.tif")
tempminstack_7d <- brick("./outdata/tempminstack_7d.tif")

# Read empirical mosquito abundance model
minmaxmod <- readRDS("./outdata/minmaxmod.rds")

# Beginning and end data (position in stack for each of four months)
begdays <- c(1, 31, 63, 94)
enddays <- c(30, 62, 93, 118)

# Latitude for diurnal temperature model
latitude <- 38

# 4 Loops - One for Each Month
for(mnum in 1:4) {
  
  begday <- begdays[mnum]
  endday <- enddays[mnum]
  
  ###############
  beginCluster(12)
  ###############
  astack <- stack()
  TFDstack <- stack()
  pEAstack <- stack()
  MDRstack <- stack()
  lfstack <- stack()
  bstack <- stack()
  cstack <- stack()
  PDRstack <- stack()
  eMstack <- stack()
  # Daily Loop
  for(x in begday:endday) {
    print(x)
    # Set up variables for the diurnal temeprature equation
    tmin <- tempminstack[[x]]
    tmax <- tempmaxstack[[x]]
    if(x > 1) {
      prev_tmin <- tempminstack[[x - 1]]
      prev_tmax <- tempmaxstack[[x - 1]]
      
    } else {
      prev_tmin <- tempminstack[[x]]
      prev_tmax <- tempmaxstack[[x]]
    }
    if(x < 118) {
      next_tmin <- tempminstack[[x + 1]]
    } else {
      next_tmin <- tempminstack[[x]]
    }
    tempras <- stack(tmin, tmax, prev_tmin, prev_tmax, next_tmin)
    cur_sun <- calcsun(165 + x)
    prev_sun <- calcsun(165 + x - 1)
    next_sun <- calcsun(165 + x + 1)
    daystack <- stack()
    # Day Loop
    for(hour in 1:24) {
      print(hour)
      print(Sys.time())
      # Generate hourly temperature estimate
      hourlytemp <-clusterR(tempras, calc, args=list(fun=hourly_temp2), export=c('hour', 'cur_sun', 'prev_sun', 'next_sun'))
      # Calcualte trait values based on temperature
      c1 <- 0.000193
      t0 <- 10.25
      tm <- 38.32
      a = clusterR(hourlytemp, calc, args=list(fun=briere), export=c('c1', 't0', 'tm'))
      c1 <- 0.0488
      t0 <- 8.02
      tm <- 35.65    
      TFD = clusterR(hourlytemp, calc, args=list(fun=briere), export=c('c1', 't0', 'tm'))
      c1 <- -0.00361
      t0 <- 9.04
      tm <- 33.93     
      pEA = clusterR(hourlytemp, calc, args=list(fun=quadratic), export=c('c1', 't0', 'tm'))
      c1 <- 0.0000638
      t0 <- 8.6
      tm <- 39.66      
      MDR = clusterR(hourlytemp, calc, args=list(fun=briere), export=c('c1', 't0', 'tm'))
      c1 <- -1.43
      t0 <- 13.41
      tm <- 31.51          
      lf = clusterR(hourlytemp, calc, args=list(fun=quadratic), export=c('c1', 't0', 'tm'))
      c1 <- 0.000735
      t0 <- 15.84
      tm <- 36.4
      lf[lf <= 1] <- 1
      b = clusterR(hourlytemp, calc, args=list(fun=briere), export=c('c1', 't0', 'tm'))
      c1 <- 0.000439
      t0 <- 3.62
      tm <- 36.82        
      c2 = clusterR(hourlytemp, calc, args=list(fun=briere), export=c('c1', 't0', 'tm'))
      c1 <- 0.000109
      t0 <- 10.39
      tm <- 43.05  
      PDR = clusterR(hourlytemp, calc, args=list(fun=briere), export=c('c1', 't0', 'tm'))
      astack <- addLayer(astack, a)
      TFDstack <- addLayer(TFDstack, TFD)
      pEAstack <- addLayer(pEAstack, pEA)
      MDRstack <- addLayer(MDRstack, MDR)
      lfstack <- addLayer(lfstack, lf)
      bstack <- addLayer(bstack, b)
      cstack <- addLayer(cstack, c2)
      PDRstack <- addLayer(PDRstack, PDR)
    } # End the hourly loop

    # Each day, generate empirical estimate of M based on rolling average min and max temps
    tempras <- stack(tempminstack_7d[[x]], tempmaxstack_7d[[x]])
    names(tempras) <- c("MA_tmin", "MA_tmax")
    eM <- clusterR(tempras, predict, args=list(minmaxmod, type = "response"))
    eMstack <- addLayer(eMstack, eM)

  } # End the daily loop
  
  # Compute the monthly mean for each trait (based on the hourly estimates)
  amean <- mean(astack)
  TFDmean <- mean (TFDstack)
  pEAmean <- mean(pEAstack)
  MDRmean <- mean(MDRstack)
  lfmean <- mean(lfstack)
  bmean <- mean(bstack)
  cmean <- mean(cstack)
  PDRmean <- mean(PDRstack)
  # Convert lifespan to mortality rate
  mumean <- 1/lfmean
  # Theoretical estimate of M
  M = TFDmean * pEAmean * MDRmean / (mumean)^2
  # Empirical estimate of M
  M2 = 0.5 * mean(eMstack)
  # VC based on theoretical estimate
  VC = M * amean^2 * bmean * cmean * exp(-1 * mumean/PDRmean)/mumean
  # VC based on empirical estimate
  VC2 = M2 * amean^2 * bmean * cmean * exp(-1 * mumean/PDRmean)/mumean
  # Stack up and output the monthly results
  outstack <- stack(amean, TFDmean, pEAmean, MDRmean, lfmean, bmean, cmean, PDRmean, mumean, M, M2, VC, VC2)
  outname = paste0("./outdata/vcstack_", mnum, ".tif")
  writeRaster(x = outstack, filename = outname, overwrite = TRUE) 
  
  ############
  endCluster()
  ############
} # End the monthly loop

junestack <- brick("./outdata/vcstack_1.tif")
julystack <- brick("./outdata/vcstack_2.tif")
augstack <- brick("./outdata/vcstack_3.tif")
septstack <- brick("./outdata/vcstack_4.tif")


##########################################
# Generate four monthly VC point estimates
##########################################

# Read data file from UGA weather station
uga_met <- read_csv("x:/Users/wimb0002/Google Drive/AthensData/uga_ggy_weather.csv")

# Process data to produce daily min and max temps and rolling mean temps
uga_2018 <- uga_met %>%
  mutate(Date = mdy(Date),
         Year = year(Date),
         Month = month(Date),
         Doy = yday(Date),
         temp = (temp - 32) * 5/9,
         temp_lo = (temp_lo - 32) * 5/9,
         temp_hi = (temp_hi - 32) * 5/9,
         temp_lo_7d = roll_meanr(temp_lo, 7),
         temp_hi_7d = roll_meanr(temp_hi, 7)) %>%
  filter(Year == 2018 & Doy >= 166 & Doy <= (166 + 117))

# Process using the same code as for the gridded temperatures
latitude <- 38
outdata <- data.frame()

for(mnum in 1:4) {
  
  begday <- begdays[mnum]
  endday <- enddays[mnum]
  
  astack <- numeric()
  TFDstack <- numeric()
  pEAstack <- numeric()
  MDRstack <- numeric()
  lfstack <- numeric()
  bstack <- numeric()
  cstack <- numeric()
  PDRstack <- numeric()
  eMstack <- numeric()
  #begday <- 1
  #endday <- 30
  for(x in begday:endday) {
    print(x)
    # Set up variables for the diurnal temperature equation
    tmin <- as.numeric(uga_2018[x, "temp_lo"])
    tmax <- as.numeric(uga_2018[x, "temp_hi"])
    tmin_7d <- as.numeric(uga_2018[x, "temp_lo_7d"])
    tmax_7d <- as.numeric(uga_2018[x, "temp_hi_7d"])
    if(x > 1) {
      prev_tmin <- as.numeric(uga_2018[x - 1, "temp_lo"])
      prev_tmax <- as.numeric(uga_2018[x - 1, "temp_hi"])
      
    } else {
      prev_tmin <- as.numeric(uga_2018[x, "temp_lo"])
      prev_tmax <- as.numeric(uga_2018[x, "temp_hi"])
    }
    if(x < 118) {
      next_tmin <- as.numeric(uga_2018[x + 1, "temp_lo"])
    } else {
      next_tmin <- as.numeric(uga_2018[x, "temp_lo"])
    }
    tempvec <- c(tmin, tmax, prev_tmin, prev_tmax, next_tmin)
    cur_sun <- calcsun(165 + x)
    prev_sun <- calcsun(165 + x - 1)
    next_sun <- calcsun(165 + x + 1)
    daystack <- stack()
    for(hour in 1:24) {
      print(hour)
      print(Sys.time())
      # Generate hourly temperature estimate
      hourlytemp <- hourly_temp2(tempvec)
      # Calculate trait values based on temperature
      c1 <- 0.000193
      t0 <- 10.25
      tm <- 38.32
      a = briere(hourlytemp)
      c1 <- 0.0488
      t0 <- 8.02
      tm <- 35.65    
      TFD = briere(hourlytemp)
      c1 <- -0.00361
      t0 <- 9.04
      tm <- 33.93     
      pEA = quadratic(hourlytemp)
      c1 <- 0.0000638
      t0 <- 8.6
      tm <- 39.66      
      MDR = briere(hourlytemp)
      c1 <- -1.43
      t0 <- 13.41
      tm <- 31.51          
      lf = quadratic(hourlytemp)
      c1 <- 0.000735
      t0 <- 15.84
      tm <- 36.4
      lf[lf <= 1] <- 1
      b = briere(hourlytemp)
      c1 <- 0.000439
      t0 <- 3.62
      tm <- 36.82        
      c2 = briere(hourlytemp)
      c1 <- 0.000109
      t0 <- 10.39
      tm <- 43.05  
      PDR = briere(hourlytemp)
      astack <- c(astack, a)
      TFDstack <- c(TFDstack, TFD)
      pEAstack <- c(pEAstack, pEA)
      MDRstack <- c(MDRstack, MDR)
      lfstack <- c(lfstack, lf)
      bstack <- c(bstack, b)
      cstack <- c(cstack, c2)
      PDRstack <- c(PDRstack, PDR)
    }

    # Generate empirical estimates of M for each day based on rolling min and max temperature
    templist <- list("MA_tmin" = tmin_7d, "MA_tmax" = tmax_7d)
    eM <- predict(minmaxmod, newdata = templist, type = "response")
    eMstack <- c(eMstack, eM)   
 
  }
  amean <- mean(astack)
  TFDmean <- mean (TFDstack)
  pEAmean <- mean(pEAstack)
  MDRmean <- mean(MDRstack)
  lfmean <- mean(lfstack)
  bmean <- mean(bstack)
  cmean <- mean(cstack)
  PDRmean <- mean(PDRstack)
  # Compute mortality rate from lifespan
  mumean <- 1/lfmean
  # Theoretical estimate of M
  M = TFDmean * pEAmean * MDRmean / (mumean)^2
  # Empirical estimate of M
  M2 = 0.5 * mean(eMstack)
  # Vectorial capacity based on the theoretical estimate
  VC = M * amean^2 * bmean * cmean * exp(-1 * mumean/PDRmean)/mumean
  # Vectorial capacity based on the empirical estimate
  VC2 = M2 * amean^2 * bmean * cmean * exp(-1 * mumean/PDRmean)/mumean
  outvect <- c(amean, TFDmean, pEAmean, MDRmean, lfmean, bmean, cmean, PDRmean, mumean, M, M2, VC, VC2)
  outdata <- rbind(outdata, outvect)
}

# Write results to disk
uga_met_vc <- outdata
names(uga_met_vc) <- c("amean", "TFDmean", "pEAmean", "MDRmean", "lfmean", "bmean", "cmean", "PDRmean", "mumean",
                       "M", "M2", "VC", "VC2")
write_csv(uga_met_vc, "./outdata/uga_met_vc.csv")


