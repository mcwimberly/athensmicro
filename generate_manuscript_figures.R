#################################################################################
# generate_manuscript_figures.r
# Summarize results of micro-scale vectorial capacity modeling
# written by Mike Wimberly, University of Oklahoma
#################################################################################

library(raster)
library(rasterVis)
library(ggplot2)
library(rgdal)
library(RColorBrewer)
library(colorRamps)
library(RStoolbox)
library(tidyverse)
library(gridExtra)
library(scales)
library(lubridate)

# Read in saved outputs
junestack <- stack("./outdata/vcstack_1.tif")
julystack <- stack("./outdata/vcstack_2.tif")
augstack <- stack("./outdata/vcstack_3.tif")
septstack <- stack("./outdata/vcstack_4.tif")
gacnty_sp <- readOGR(layer = "GA_SHP", dsn="x:/Users/wimb0002/Google Drive/AthensData")
gacnty_f <- fortify(gacnty_sp, regions="GEOID10")
tempminstack <- brick("./outdata/tempminstack.tif")
tempmaxstack <- brick("./outdata/tempmaxstack.tif")
rhminstack <- brick("./outdata/rhminstack.tif")
tmingmet <- stack("./outdata/tmincrop.tif")
tmaxgmet <- stack("./outdata/tmaxcrop.tif")
vcplot <- stack(junestack[[13]], julystack[[13]], augstack[[13]], septstack[[13]])
names(vcplot) <- c("June-July", "July-Aug", "Aug-Sept", "Sept-Oct")
vcjune_fm <- focal(junestack[[13]], focalWeight(junestack[[13]], 0.001, "circle"))
vcjuly_fm <- focal(julystack[[13]], focalWeight(julystack[[13]], 0.001, "circle"))
vcaug_fm <- focal(augstack[[13]], focalWeight(augstack[[13]], 0.001, "circle"))
vcsept_fm <- focal(septstack[[13]], focalWeight(septstack[[13]], 0.001, "circle"))
vcplot_fm <- stack(vcjune_fm, vcjuly_fm, vcaug_fm, vcsept_fm)
names(vcplot_fm) <- c("June-July", "July-Aug", "Aug-Sept", "Sept-Oct")

treecov <- raster("x:/Users/wimb0002/Google Drive/AthensData/Athens_treecover.tif")
imperv <- raster("x:/Users/wimb0002/Google Drive/AthensData/Athens_impervious.tif")
ndwi <- raster("x:/Users/wimb0002/Google Drive/AthensData/Athens_NDWI.tif")
sentinel <- stack("x:/Users/wimb0002/Google Drive/AthensData/SentinelComposite.tif")
popdens <- raster("x:/Users/wimb0002/Google Drive/AthensData/athens_popdens.tif")

ndwi2 <- resample(ndwi, vcplot)

###############################################################################
# Temperature time series figure
###############################################################################
uga_met <- read_csv("x:/Users/wimb0002/Google Drive/AthensData/uga_ggy_weather.csv")

uga_2018 <- uga_met %>%
  mutate(Date = mdy(Date),
         Year = year(Date),
         Month = month(Date),
         Doy = yday(Date),
         temp = (temp - 32) * 5/9,
         temp_lo = (temp_lo - 32) * 5/9,
         temp_hi = (temp_hi - 32) * 5/9) %>%
  filter(Year == 2018 & Doy >= 166 & Doy <= (166 + 117))


micro <- read_csv("./outdata/micro_2018_daily.csv")

microsum <- micro %>%
  filter(Date >= ymd("2018-06-15"),
         Date <= ymd("2018-10-10"),
         num_temp > 95) %>%
  group_by(Date) %>%
  summarise(dailymin_max = quantile(min_temp, 0.9),
            dailymin_min = quantile(min_temp, 0.1),
            dailymax_max = quantile(max_temp, 0.9),
            dailymax_min = quantile(max_temp, 0.1),
            dailymin_med = median(min_temp, na.rm=T),
            dailymax_med = median(max_temp, na.rm=T),
            dailycount = n())

microsum1 <- microsum %>%
  filter(Date < ymd("2018-06-22"))

microsum2 <- microsum %>%
  filter(Date > ymd("2018-06-22") & Date < ymd("2018-09-15"))

microsum3 <- microsum %>%
  filter(Date > ymd("2018-09-20"))

#png("./outfigs/tempts.png", width = 6, height = 3.5, units = "in", res = 300)
tiff("./outfigs/Fig_2.tif", width = 6, height = 3.5, units = "in", res = 300, compression="lzw")

ggplot() +
  geom_ribbon(data = filter(microsum1, dailycount >= 0), aes(x = Date, ymin = dailymin_min, ymax = dailymin_max), fill = "orange") +
  geom_ribbon(data = filter(microsum2, dailycount >= 0), aes(x = Date, ymin = dailymin_min, ymax = dailymin_max), fill = "orange") +
  geom_ribbon(data = filter(microsum3, dailycount >= 0), aes(x = Date, ymin = dailymin_min, ymax = dailymin_max), fill = "orange") +
  geom_ribbon(data = filter(microsum1, dailycount >= 0), aes(x = Date, ymin = dailymax_min, ymax = dailymax_max), fill = "indianred1") +
  geom_ribbon(data = filter(microsum2, dailycount >= 0), aes(x = Date, ymin = dailymax_min, ymax = dailymax_max), fill = "indianred1") +
  geom_ribbon(data = filter(microsum3, dailycount >= 0), aes(x = Date, ymin = dailymax_min, ymax = dailymax_max), fill = "indianred1") +
  geom_line(data = filter(microsum1, dailycount >= 0), aes(x = Date, y = dailymax_med), color = "red3", size = 0.5) +
  geom_line(data = filter(microsum2, dailycount >= 0), aes(x = Date, y = dailymax_med), color = "red3", size = 0.5) +
  geom_line(data = filter(microsum3, dailycount >= 0), aes(x = Date, y = dailymax_med), color = "red3", size = 0.5) +
  geom_line(data = filter(microsum1, dailycount >= 0), aes(x = Date, y = dailymin_med), color = "orange3", size = 0.5) +
  geom_line(data = filter(microsum2, dailycount >= 0), aes(x = Date, y = dailymin_med), color = "orange3", size = 0.5) +
  geom_line(data = filter(microsum3, dailycount >= 0), aes(x = Date, y = dailymin_med), color = "orange3", size = 0.5) +
  geom_line(data = uga_2018, aes(x = Date, y = temp_lo), color = "black", size = 0.5, linetype = 2) +
  geom_line(data = uga_2018, aes(x = Date, y = temp_hi), color = "black", size = 0.5, linetype = 2)  +
  labs(y = "Temperature (\u00B0C)", x = "") +
  theme_bw()
dev.off()

###############################################################################
# Temperature Map Figure
###############################################################################

#png("./outfigs/tempmaps.png", width = 8, height = 4.5, units = "in", res = 300)
tiff("./outfigs/Fig_3.tif", width = 8, height = 4.5, units = "in", res = 300, compression="lzw")

tempminstack[ndwi2 > 0] <- NA
tempmaxstack[ndwi2 > 0] <- NA

fig2a <- gplot(tempminstack[[1 + 29]], maxpixels = 5000000) +
  geom_raster(aes(fill = value)) +
  geom_polygon(data = gacnty_f, aes(x = long, y = lat, group = group), color="black", fill=NA, size=0.5) +
  scale_fill_gradientn(name = "Temp (\u00B0C)", colors=brewer.pal(9, "Oranges")) +
  coord_fixed(xlim = c(-83.55, -83.24), ylim = c(33.847, 34.044), expand = F, ratio = 1) +
  labs(title = "A") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12, hjust = 0, face = "bold"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        panel.border = element_rect(fill = NA, color = "black"),
        panel.background = element_rect(fill = "white", color = NA))

fig2b <- gplot(tempmaxstack[[1 + 29]], maxpixels = 5000000) +
  geom_raster(aes(fill = value)) +
  geom_polygon(data = gacnty_f, aes(x = long, y = lat, group = group), color="black", fill=NA, size=0.5) +
  scale_fill_gradientn(name = "Temp (\u00B0C)", colors=brewer.pal(9, "Reds")) +
  coord_fixed(xlim = c(-83.55, -83.24), ylim = c(33.847, 34.044), expand = F, ratio = 1) +
  labs(title = "B") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12, hjust = 0, face = "bold"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        panel.border = element_rect(fill = NA, color = "black"),
        panel.background = element_rect(fill = "white", color = NA))

fig2c <- gplot(tmingmet[[166 + 29]] - 273.15, maxpixels = 5000000) +
  geom_raster(aes(fill = value)) +
  geom_polygon(data = gacnty_f, aes(x = long, y = lat, group = group), color="black", fill=NA, size=0.5) +
  scale_fill_gradientn(name = "Temp (\u00B0C)", colors=brewer.pal(9, "Oranges")) +
  coord_fixed(xlim = c(-83.55, -83.24), ylim = c(33.847, 34.044), expand = F, ratio = 1) +
  labs(title = "C") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12, hjust = 0, face = "bold"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        panel.border = element_rect(fill = NA, color = "black"),
        panel.background = element_rect(fill = "white", color = NA))

fig2d <- gplot(tmaxgmet[[166 + 29]] - 273.15, maxpixels = 5000000) +
  geom_raster(aes(fill = value)) +
  geom_polygon(data = gacnty_f, aes(x = long, y = lat, group = group), color="black", fill=NA, size=0.5) +
  scale_fill_gradientn(name = "Temp (\u00B0C)", colors=brewer.pal(9, "Reds")) +
  coord_fixed(xlim = c(-83.55, -83.24), ylim = c(33.847, 34.044), expand = F, ratio = 1) +
  labs(title = "D") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12, hjust = 0, face = "bold"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        panel.border = element_rect(fill = NA, color = "black"),
        panel.background = element_rect(fill = "white", color = NA))

grid.arrange(fig2a, fig2b, fig2c, fig2d, ncol = 2)

dev.off()

#################################
# Seasonal Temperature Map Figure
#################################

tempminstack[ndwi2 > 0] <- NA
tempmaxstack[ndwi2 > 0] <- NA

tmin_june <- mean(tempminstack[[1:31]])
tmin_july <- mean(tempminstack[[32:62]])
tmin_aug <- mean(tempminstack[[63:93]])
tmin_sept <- mean(tempminstack[[94:118]])

tmin_months <- stack(tmin_june, tmin_july, tmin_aug, tmin_sept)
names(tmin_months) <- c("A", "B", "C", "D")

tmax_june <- mean(tempmaxstack[[1:31]])
tmax_july <- mean(tempmaxstack[[32:62]])
tmax_aug <- mean(tempmaxstack[[63:93]])
tmax_sept <- mean(tempmaxstack[[94:118]])

tmax_months <- stack(tmax_june, tmax_july, tmax_aug, tmax_sept)
names(tmax_months) <- c("A", "B", "C", "D")

tiff("./outfigs/tmin_month.tif", width = 8, height = 5, units = "in", res = 300, compression = "lzw")

gplot(tmin_months, maxpixels=5000000) +
  geom_raster(aes(fill = value)) +
  geom_polygon(data = gacnty_f, aes(x = long, y = lat, group = group), color="black", fill=NA, size=0.5) +
  scale_fill_distiller(name = "Temp (\u00B0C)", palette = "Oranges", direction = 1) +
  coord_fixed(xlim = c(-83.55, -83.24), ylim = c(33.847, 34.044), expand = F, ratio = 1) +
  facet_wrap(~ variable, ncol = 2) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 14, hjust = 0),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        strip.background = element_rect(fill = NA, color = NA),
        strip.text = element_text(hjust = 0, face = "bold"),
        panel.border = element_rect(fill = NA, color = "black"),
        panel.background = element_rect(fill = "white", color = NA))

dev.off()

tiff("./outfigs/tmax_month.tif", width = 8, height = 5, units = "in", res = 300, compression = "lzw")

gplot(tmax_months, maxpixels=5000000) +
  geom_raster(aes(fill = value)) +
  geom_polygon(data = gacnty_f, aes(x = long, y = lat, group = group), color="black", fill=NA, size=0.5) +
  scale_fill_distiller(name = "Temp (\u00B0C)", palette = "Reds", direction = 1) +
  coord_fixed(xlim = c(-83.55, -83.24), ylim = c(33.847, 34.044), expand = F, ratio = 1) +
  facet_wrap(~ variable, ncol = 2) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 14, hjust = 0),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        strip.background = element_rect(fill = NA, color = NA),
        strip.text = element_text(hjust = 0, face = "bold"),
        panel.border = element_rect(fill = NA, color = "black"),
        panel.background = element_rect(fill = "white", color = NA))

dev.off()

#############################################
# Monthly Temperature Vs. Mosquito Abundance
#############################################

clarke <- gacnty_sp[gacnty_sp@data$COUNTYFP10 == "059",]
clarke_r <- rasterize(clarke, vcplot, field = "COUNTYFP10")
clarke_r2 <- projectRaster(clarke_r, vcplot)

tempminstack[ndwi2 > 0] <- NA
tempmaxstack[ndwi2 > 0] <- NA

tmin_june <- mean(tempminstack[[1:31]])
tmin_july <- mean(tempminstack[[32:62]])
tmin_aug <- mean(tempminstack[[63:93]])
tmin_sept <- mean(tempminstack[[94:118]])

tmin_months <- stack(tmin_june, tmin_july, tmin_aug, tmin_sept)
names(tmin_months) <- c("A", "B", "C", "D")
tmin_months <- mask(tmin_months, clarke_r2)


tmax_june <- mean(tempmaxstack[[1:31]])
tmax_july <- mean(tempmaxstack[[32:62]])
tmax_aug <- mean(tempmaxstack[[63:93]])
tmax_sept <- mean(tempmaxstack[[94:118]])

tmax_months <- stack(tmax_june, tmax_july, tmax_aug, tmax_sept)
names(tmax_months) <- c("A", "B", "C", "D")
tmax_months <- mask(tmax_months, clarke_r2)

mplot <- stack(junestack[[10]], julystack[[10]], augstack[[10]], septstack[[10]])
mplot2 <- stack(junestack[[11]], julystack[[11]], augstack[[11]], septstack[[11]])
vcplot <- stack(junestack[[12]], julystack[[12]], augstack[[12]], septstack[[12]])
vcplot2 <- stack(junestack[[13]], julystack[[13]], augstack[[13]], septstack[[13]])
mplot[ndwi2 > 0] <- NA
mplot2[ndwi2 > 0] <- NA
vcplot[ndwi2 > 0] <- NA
vcplot2[ndwi2 > 0] <- NA
mplot <- mask(mplot, clarke_r2)
mplot2 <- mask(mplot2, clarke_r2)
vcplot <- mask(vcplot, clarke_r2)
vcplot2 <- mask(vcplot2, clarke_r2)

tmin <- as.vector(tmin_months)
tmax <- as.vector(tmax_months)
mtheor <- as.vector(mplot)
memp <- as.vector(mplot2)
vctheor <- as.vector(vcplot)
vcemp <- as.vector(vcplot2)

tmabund <- data.frame(tmin, tmax, mtheor, memp, vctheor, vcemp)
sampdata <- sample_n(tmabund, 200000) 

loess1 <- loess(mtheor ~ tmin * tmax, data = sampdata)
loess2 <- loess(memp ~ tmin * tmax, data = sampdata)
loess3 <- loess(vctheor ~ tmin * tmax, data = sampdata)
loess4 <- loess(vcemp ~ tmin * tmax, data = sampdata)

tmaxgrid <-  seq(29, 33.5, 0.05)
tmingrid <-  seq(19, 22.5, 0.05)
# Generate a dataframe with every possible combination of wt and hp
data.fit <-  expand.grid(tmin = tmingrid, tmax = tmaxgrid)

# Feed the dataframe into the loess model and receive a matrix output with estimates of 
# acceleration for each combination of wt and hp
loessp1 <-  predict(loess1, newdata = data.fit)
loessp2 <-  predict(loess2, newdata = data.fit)
loessp3 <-  predict(loess3, newdata = data.fit)
loessp4 <-  predict(loess4, newdata = data.fit)

plotdat <- data.frame(data.fit, 
                      mtheor = as.numeric(loessp1), 
                      memp = as.numeric(loessp2),
                      vctheor = as.numeric(loessp3),
                      vcemp = as.numeric(loessp4))

maxdiurn <- max(sampdata$tmax - sampdata$tmin, na.rm=TRUE)
mindiurn <- min(sampdata$tmax - sampdata$tmin, na.rm=TRUE)

plotdat$mtheor[plotdat$tmin > (plotdat$tmax - mindiurn) | plotdat$tmax > (plotdat$tmin + maxdiurn)] <- NA
plotdat$memp[plotdat$tmin > (plotdat$tmax - mindiurn) | plotdat$tmax > (plotdat$tmin + maxdiurn)] <- NA
plotdat$vctheor[plotdat$tmin > (plotdat$tmax - mindiurn) | plotdat$tmax > (plotdat$tmin + maxdiurn)] <- NA
plotdat$vcemp[plotdat$tmin > (plotdat$tmax - mindiurn) | plotdat$tmax > (plotdat$tmin + maxdiurn)] <- NA

#plotdat <- bind_rows(plotdat1, plotdat2, plotdat3, plotdat4)

#png("./outfigs/lcmaps.png", width = 7, height = 6, units = "in", res = 300)
tiff("./outfigs/temp_m.tif", width = 7, height = 6, units = "in", res = 300, compression="lzw")

ggplot(plotdat, aes(x = tmin, y = tmax, z = mtheor)) +
  stat_contour(geom = "polygon", aes(fill = ..level..)) +
  geom_tile(aes(fill = mtheor)) +
  stat_contour(binwidth = 1000) +
  scale_fill_distiller(name = expression(italic("M(T)")), palette = "RdYlBu",
                       na.value = NA) +
  scale_x_continuous(limits = c(19, 22.2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(29.2, 33.5), expand = c(0, 0)) +
  labs(x="Minimum Temperature (\u00B0C)", y = "Maximum Temperature (\u00B0C)") +
  coord_fixed(ratio = 1) + 
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 12))


dev.off()


#png("./outfigs/lcmaps.png", width = 7, height = 6, units = "in", res = 300)
tiff("./outfigs/temp_m2.tif", width = 7, height = 6, units = "in", res = 300, compression="lzw")

ggplot(plotdat, aes(x = tmin, y = tmax, z = memp)) +
  stat_contour(geom = "polygon", aes(fill = ..level..)) +
  geom_tile(aes(fill = memp)) +
  stat_contour(binwidth = 0.5) +
  scale_fill_distiller(name = expression(italic(M(T[min], T[max]))), palette = "RdYlBu",
                       na.value = NA) +
  scale_x_continuous(limits = c(19, 22.2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(29.2, 33.5), expand = c(0, 0)) +
  labs(x="Minimum Temperature (\u00B0C)", y = "Maximum Temperature (\u00B0C)") +
  coord_fixed(ratio = 1) + 
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 12))


dev.off()

tiff("./outfigs/temp_vc.tif", width = 7, height = 6, units = "in", res = 300, compression="lzw")

ggplot(plotdat, aes(x = tmin, y = tmax, z = vctheor)) +
  stat_contour(geom = "polygon", aes(fill = ..level..)) +
  geom_tile(aes(fill = vctheor)) +
  stat_contour(binwidth = 2000) +
  scale_fill_distiller(name = expression(italic("VC(T)")), palette = "RdYlBu",
                       na.value = NA) +
  scale_x_continuous(limits = c(19, 22.2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(29.2, 33.5), expand = c(0, 0)) +
  labs(x="Minimum Temperature (\u00B0C)", y = "Maximum Temperature (\u00B0C)") +
  coord_fixed(ratio = 1) + 
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 12))


dev.off()

tiff("./outfigs/temp_vc2.tif", width = 7, height = 6, units = "in", res = 300, compression="lzw")

ggplot(plotdat, aes(x = tmin, y = tmax, z = vcemp)) +
  stat_contour(geom = "polygon", aes(fill = ..level..)) +
  geom_tile(aes(fill = vcemp)) +
  stat_contour(binwidth = 1) +
  scale_fill_distiller(name = expression(italic("VC(T)")), palette = "RdYlBu",
                       na.value = NA) +
  scale_x_continuous(limits = c(19, 22.2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(29.2, 33.5), expand = c(0, 0)) +
  labs(x="Minimum Temperature (\u00B0C)", y = "Maximum Temperature (\u00B0C)") +
  coord_fixed(ratio = 1) + 
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 12))


dev.off()


###############################################################################
# Temperature-Humidity Map Figure
###############################################################################

png("tempmaps2.png", width = 8, height = 3, units = "in", res = 300)

fig2a <- gplot(tempminstack[[60]], maxpixels = 5000000) +
  geom_raster(aes(fill = value)) +
  geom_polygon(data = gacnty_f, aes(x = long, y = lat, group = group), color="black", fill=NA, size=0.5) +
  scale_fill_gradientn(name = "Temp (C)", colors=brewer.pal(9, "Oranges")) +
  coord_fixed(xlim = c(-83.55, -83.24), ylim = c(33.847, 34.044), expand = F, ratio = 1) +
  labs(title = "a) Minimum Temperature") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 10, hjust = 0),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.key.size =  unit(0.15, "in"),
        legend.position = "bottom",
        panel.border = element_rect(fill = NA, color = "black"),
        panel.background = element_rect(fill = "white", color = NA))

fig2b <- gplot(tempmaxstack[[60]], maxpixels = 5000000) +
  geom_raster(aes(fill = value)) +
  geom_polygon(data = gacnty_f, aes(x = long, y = lat, group = group), color="black", fill=NA, size=0.5) +
  scale_fill_gradientn(name = "Temp (C)", colors=brewer.pal(9, "Reds")) +
  coord_fixed(xlim = c(-83.55, -83.24), ylim = c(33.847, 34.044), expand = F, ratio = 1) +
  labs(title = "b) Maximum Temperature") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 10, hjust = 0),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.key.size =  unit(0.15, "in"),
        legend.position = "bottom",
        panel.border = element_rect(fill = NA, color = "black"),
        panel.background = element_rect(fill = "white", color = NA))

fig2c <- gplot(rhminstack[[60]], maxpixels = 5000000) +
  geom_raster(aes(fill = value)) +
  geom_polygon(data = gacnty_f, aes(x = long, y = lat, group = group), color="black", fill=NA, size=0.5) +
  scale_fill_gradientn(name = "RH (%)  ", colors=brewer.pal(9, "Blues")) +
  coord_fixed(xlim = c(-83.55, -83.24), ylim = c(33.847, 34.044), expand = F, ratio = 1) +
  labs(title = "c) Minimum Relative Humidity") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 10, hjust = 0),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.key.size =  unit(0.15, "in"),
        legend.position = "bottom",
        panel.border = element_rect(fill = NA, color = "black"),
        panel.background = element_rect(fill = "white", color = NA))


grid.arrange(fig2a, fig2b, fig2c, ncol = 3)

dev.off()

###############################################################################
# Temperature-Humidity Map Figure #2
###############################################################################

png("tempmaps2.png", width = 8, height = 4.5, units = "in", res = 300)

july_m2 <- julystack[[11]]  
july_m2r <- july_m2 / max(as.matrix(july_m2), na.rm=T)
july_vc2 <- julystack[[13]]  
july_vc2r <- july_vc2 / max(as.matrix(july_vc2), na.rm=T)

writeRaster(july_m2, "july_m2.tif")
writeRaster(july_vc2, "july_vc2.tif")
writeRaster(tempmaxstack[[31]], "July15tmax.tif")
writeRaster(rhminstack[[31]], "July15rhmin.tif")


fig2a <- gplot(tempmaxstack[[31]], maxpixels = 5000000) +
  geom_raster(aes(fill = value)) +
  geom_polygon(data = gacnty_f, aes(x = long, y = lat, group = group), color="black", fill=NA, size=0.5) +
  scale_fill_gradientn(name = "Temp (C)", colors=brewer.pal(9, "Reds")) +
  coord_fixed(xlim = c(-83.55, -83.24), ylim = c(33.847, 34.044), expand = F, ratio = 1) +
  labs(title = "A) Maximum Temperature - July 15, 2018") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 10, hjust = 0),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        #legend.key.size =  unit(0.15, "in"),
        #legend.position = "bottom",
        panel.border = element_rect(fill = NA, color = "black"),
        panel.background = element_rect(fill = "white", color = NA))

fig2b <- gplot(rhminstack[[31]], maxpixels = 5000000) +
  geom_raster(aes(fill = value)) +
  geom_polygon(data = gacnty_f, aes(x = long, y = lat, group = group), color="black", fill=NA, size=0.5) +
  scale_fill_gradientn(name = "RH (%)  ", colors=brewer.pal(9, "Blues")) +
  coord_fixed(xlim = c(-83.55, -83.24), ylim = c(33.847, 34.044), expand = F, ratio = 1) +
  labs(title = "B) Minimum Relative Humidity - July 15, 2018") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 10, hjust = 0),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        #legend.key.size =  unit(0.15, "in"),
        #legend.position = "bottom",
        panel.border = element_rect(fill = NA, color = "black"),
        panel.background = element_rect(fill = "white", color = NA))

sa3sfb <- gplot(july_m2r, maxpixels=5000000) +
  geom_raster(aes(fill = value)) +
  geom_polygon(data = gacnty_f, aes(x = long, y = lat, group = group), color="black", fill=NA, size=0.5) +
  #scale_fill_gradient2(name = "M", low = "#313695", mid = "#ffffbf", high = "#a50026", midpoint=mmed2) +
  scale_fill_distiller(name = "M", palette = "YlOrRd", direction = 1) +
  coord_fixed(xlim = c(-83.55, -83.24), ylim = c(33.847, 34.044), expand = F, ratio = 1) +
  #facet_wrap(~ variable, ncol = 2) +
  labs(title = "C) Mosquito Abundance - July 2018") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 10, hjust = 0),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        strip.background = element_rect(fill = NA, color = NA),
        strip.text = element_text(hjust = 0, face = "bold"),
        panel.border = element_rect(fill = NA, color = "black"),
        panel.background = element_rect(fill = "white", color = NA))

sa3sfc <- gplot(july_vc2r, maxpixels=5000000) +
  geom_raster(aes(fill = value)) +
  #geom_polygon(data = blocks_f, aes(x = long, y = lat, group = group), color="blue", fill=NA, size=1) +
  geom_polygon(data = gacnty_f, aes(x = long, y = lat, group = group), color="black", fill=NA, size=0.5) +
  #scale_fill_gradient2(name = "VCr", low = "#313695", mid = "#ffffbf", high = "#a50026", midpoint=vcmed) +
  scale_fill_distiller(name = "VC", palette = "YlOrRd", direction = 1) +
  coord_fixed(xlim = c(-83.55, -83.24), ylim = c(33.847, 34.044), expand = F, ratio = 1) +
  #facet_wrap(~ variable, ncol = 2) +
  labs(title = "D) Vectorial Capacity - July 2018") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 10, hjust = 0),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        strip.background = element_rect(fill = NA, color = NA),
        strip.text = element_text(hjust = 0, face = "bold"),
        panel.border = element_rect(fill = NA, color = "black"),
        panel.background = element_rect(fill = "white", color = NA))  


grid.arrange(fig2a, fig2b, sa3sfb, sa3sfc, ncol = 2)

dev.off()

###############################################################################
# VC Histogram Plot
###############################################################################
uga_met_vc <- read_csv("./outdata/uga_met_vc.csv")

vcplot <- stack(junestack[[13]], julystack[[13]], augstack[[13]], septstack[[13]])
names(vcplot) <- c("June-July", "July-Aug", "Aug-Sept", "Sept-Oct")

clarke <- gacnty_sp[gacnty_sp@data$COUNTYFP10 == "059",]
clarke_r <- rasterize(clarke, vcplot, field = "COUNTYFP10")
clarke_r2 <- projectRaster(clarke_r, vcplot)
vcplot <- mask(vcplot, clarke_r2)

vcmat <- as.matrix(vcplot)
vcmax <- max(vcmat, na.rm=T)
vcmin <- min(vcmat, na.rm=T)

vcplot <- (vcplot - vcmin) / (vcmax - vcmin)

# Note - taking all data for the histograms, not just a sample
vcsample <- as.data.frame(as.matrix(vcplot))

vcsamplot <- vcsample %>%
  gather(key = "month", value = "vc", June.July, July.Aug, Aug.Sept, Sept.Oct) %>%
  mutate(month = factor(month, levels = c("June.July", "July.Aug", "Aug.Sept", "Sept.Oct"), 
                        labels = c("June-July", "July-Aug", "Aug-Sept", "Sept-Oct"))) %>%
  drop_na()

vline.data <- data.frame(z = uga_met_vc$VC2, month=c("June-July","July-Aug","Aug-Sept","Sept-Oct"))
vline.data <- vline.data %>%
  mutate(z = z/vcmax)

#png("./outfigs/vchist.png", width = 4, height = 6, units = "in", res = 300)
tiff("./outfigs/Fig_4.tif", width = 4, height = 6, units = "in", res = 300, compression="lzw")

ggplot(data = vcsamplot) +
  geom_histogram(aes(x = vc), fill = "skyblue2", bins = 20) +
  #geom_histogram(aes(x = vc), vcsamplot2, fill = "indianred2", bins = 20) +
  geom_vline(aes(xintercept = z), vline.data, color = "black") +
  scale_y_continuous(labels =  function(x) scientific(x, digits = 2)) +
  facet_wrap(~ month, ncol = 1, scales = "free_y") +
  theme_bw() + 
  labs(x=expression(italic("VC(T)")), y = "Count") 

dev.off()

#######################################
# Vectorial Capacity Maps - Empirical
#######################################

vcplot <- stack(junestack[[13]], julystack[[13]], augstack[[13]], septstack[[13]])
vcplot[ndwi2 > 0] <- NA
names(vcplot) <- c("A", "B", "C", "D")
vcmat <- as.matrix(vcplot)
vcmax <- max(vcmat, na.rm=T)
#vcplot <- vcplot / vcmax
vcmed <- median(as.matrix(vcplot), na.rm=T)

tiff("./outfigs/vcmap_empm.tif", width = 8, height = 5, units = "in", res = 300, compression = "lzw")

gplot(vcplot, maxpixels=5000000) +
  geom_raster(aes(fill = value)) +
  #geom_polygon(data = blocks_f, aes(x = long, y = lat, group = group), color="blue", fill=NA, size=1) +
  geom_polygon(data = gacnty_f, aes(x = long, y = lat, group = group), color="black", fill=NA, size=0.5) +
  #scale_fill_gradient2(name = "VCr", low = "#313695", mid = "#ffffbf", high = "#a50026", midpoint=vcmed) +
  scale_fill_distiller(name = expression(italic("VC(T)")), palette = "YlOrRd", direction = 1) +
  coord_fixed(xlim = c(-83.55, -83.24), ylim = c(33.847, 34.044), expand = F, ratio = 1) +
  facet_wrap(~ variable, ncol = 2) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 14, hjust = 0),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        strip.background = element_rect(fill = NA, color = NA),
        strip.text = element_text(hjust = 0, face = "bold"),
        panel.border = element_rect(fill = NA, color = "black"),
        panel.background = element_rect(fill = "white", color = NA))

dev.off()



#######################################
# Vectorial Capacity Maps - Theoretical
#######################################

vcplot <- stack(junestack[[12]], julystack[[12]], augstack[[12]], septstack[[12]])
vcplot[ndwi2 > 0] <- NA
names(vcplot) <- c("A", "B", "C", "D")
vcmat <- as.matrix(vcplot_fm)
vcmax <- max(vcmat, na.rm=T)
#vcplot <- vcplot / vcmax
vcmed <- median(as.matrix(vcplot), na.rm=T)

tiff("./outfigs/vcmap_theom.tif", width = 8, height = 5, units = "in", res = 300, compression = "lzw")

gplot(vcplot, maxpixels=5000000) +
  geom_raster(aes(fill = value)) +
  #geom_polygon(data = blocks_f, aes(x = long, y = lat, group = group), color="blue", fill=NA, size=1) +
  geom_polygon(data = gacnty_f, aes(x = long, y = lat, group = group), color="black", fill=NA, size=0.5) +
  #scale_fill_gradient2(name = "VCr", low = "#313695", mid = "#ffffbf", high = "#a50026", midpoint=vcmed) +
  scale_fill_distiller(name = expression(italic("VC(T)")), palette = "YlOrRd", direction = 1) +
  coord_fixed(xlim = c(-83.55, -83.24), ylim = c(33.847, 34.044), expand = F, ratio = 1) +
  facet_wrap(~ variable, ncol = 2) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 14, hjust = 0),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        strip.background = element_rect(fill = NA, color = NA),
        strip.text = element_text(hjust = 0, face = "bold"),
        panel.border = element_rect(fill = NA, color = "black"),
        panel.background = element_rect(fill = "white", color = NA))

dev.off()

########################################################
# Vectorial Capacity Maps - separate legend for each map
########################################################

names(vcplot) <- c("A", "B", "C", "D")
vcmat <- as.matrix(vcplot_fm)
vcmax <- max(vcmat, na.rm=T)
vcmin <- min(vcmat, na.rm=T)
vcplot <- (vcplot_fm - vcmin) / (vcmax - vcmin)
vcmed <- median(as.matrix(vcplot), na.rm=T)

#png("./outfigs/vcmap_v2.png", width = 8, height = 4.5, units = "in", res = 300)
tiff("./outfigs/Fig_5.tif", width = 8, height = 4.5, units = "in", res = 300, compression="lzw")

vcplot[ndwi2 > 0] <- NA

june_vc <- gplot(vcplot[[1]], maxpixels=5000000) +
  geom_raster(aes(fill = value)) +
  #geom_polygon(data = blocks_f, aes(x = long, y = lat, group = group), color="blue", fill=NA, size=1) +
  geom_polygon(data = gacnty_f, aes(x = long, y = lat, group = group), color="black", fill=NA, size=0.5) +
  #scale_fill_gradient2(name = "VCr", low = "#313695", mid = "#ffffbf", high = "#a50026", midpoint=vcmed) +
  scale_fill_distiller(name = expression(italic("VC(T)")), palette = "RdYlBu", direction = -1) +
  coord_fixed(xlim = c(-83.55, -83.24), ylim = c(33.847, 34.044), expand = F, ratio = 1) +
  labs(title = "A") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12, hjust = 0, face = "bold"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        panel.border = element_rect(fill = NA, color = "black"),
        panel.background = element_rect(fill = "white", color = NA))

july_vc <- gplot(vcplot[[2]], maxpixels=5000000) +
  geom_raster(aes(fill = value)) +
  #geom_polygon(data = blocks_f, aes(x = long, y = lat, group = group), color="blue", fill=NA, size=1) +
  geom_polygon(data = gacnty_f, aes(x = long, y = lat, group = group), color="black", fill=NA, size=0.5) +
  #scale_fill_gradient2(name = "VCr", low = "#313695", mid = "#ffffbf", high = "#a50026", midpoint=vcmed) +
  scale_fill_distiller(name = expression(italic("VC(T)")), palette = "RdYlBu", direction = -1) +
  coord_fixed(xlim = c(-83.55, -83.24), ylim = c(33.847, 34.044), expand = F, ratio = 1) +
  labs(title = "B") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12, hjust = 0, face = "bold"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        panel.border = element_rect(fill = NA, color = "black"),
        panel.background = element_rect(fill = "white", color = NA))

aug_vc <- gplot(vcplot[[3]], maxpixels=5000000) +
  geom_raster(aes(fill = value)) +
  #geom_polygon(data = blocks_f, aes(x = long, y = lat, group = group), color="blue", fill=NA, size=1) +
  geom_polygon(data = gacnty_f, aes(x = long, y = lat, group = group), color="black", fill=NA, size=0.5) +
  #scale_fill_gradient2(name = "VCr", low = "#313695", mid = "#ffffbf", high = "#a50026", midpoint=vcmed) +
  scale_fill_distiller(name = expression(italic("VC(T)")), palette = "RdYlBu", direction = -1) +
  coord_fixed(xlim = c(-83.55, -83.24), ylim = c(33.847, 34.044), expand = F, ratio = 1) +
  labs(title = "C") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12, hjust = 0, face = "bold"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        panel.border = element_rect(fill = NA, color = "black"),
        panel.background = element_rect(fill = "white", color = NA))

sept_vc <- gplot(vcplot[[4]], maxpixels=5000000) +
  geom_raster(aes(fill = value)) +
  #geom_polygon(data = blocks_f, aes(x = long, y = lat, group = group), color="blue", fill=NA, size=1) +
  geom_polygon(data = gacnty_f, aes(x = long, y = lat, group = group), color="black", fill=NA, size=0.5) +
  #scale_fill_gradient2(name = "VCr", low = "#313695", mid = "#ffffbf", high = "#a50026", midpoint=vcmed) +
  scale_fill_distiller(name = expression(italic("VC(T)")), palette = "RdYlBu", direction = -1) +
  coord_fixed(xlim = c(-83.55, -83.24), ylim = c(33.847, 34.044), expand = F, ratio = 1) +
  labs(title = "D") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12, hjust = 0, face = "bold"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        panel.border = element_rect(fill = NA, color = "black"),
        panel.background = element_rect(fill = "white", color = NA))

grid.arrange(june_vc, july_vc, aug_vc, sept_vc, ncol = 2)

dev.off()

#######################################
# Mosquito Abundance Maps - Theoretical
#######################################

  mplot <- stack(junestack[[10]], julystack[[10]], augstack[[10]], septstack[[10]])
  mplot[ndwi2 > 0] <- NA
  names(mplot) <- c("A", "B", "C", "D")
  mmed <- median(as.matrix(mplot), na.rm=T)
  
  tiff("./outfigs/mmaps.tif", width = 8, height = 5, units = "in", res = 300, compression="lzw")
  
  gplot(mplot, maxpixels=5000000) +
    geom_raster(aes(fill = value)) +
    geom_polygon(data = gacnty_f, aes(x = long, y = lat, group = group), color="black", fill=NA, size=0.5) +
    #scale_fill_gradient2(name = "M", low = "#313695", mid = "#ffffbf", high = "#a50026", midpoint=mmed) +
    scale_fill_distiller(name = expression(italic("M(T)")), palette = "YlOrRd", direction = 1) +
    coord_fixed(xlim = c(-83.55, -83.24), ylim = c(33.847, 34.044), expand = F, ratio = 1) +
    facet_wrap(~ variable, ncol = 2) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 14, hjust = 0),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 8),
          strip.background = element_rect(fill = NA, color = NA),
          strip.text = element_text(hjust = 0, face = "bold"),
          panel.border = element_rect(fill = NA, color = "black"),
          panel.background = element_rect(fill = "white", color = NA))
  
  dev.off()


  #######################################
  # Mosquito Abundance Maps - Empirical
  #######################################
 
  mplot2 <- stack(junestack[[11]], julystack[[11]], augstack[[11]], septstack[[11]])
  mplot2[ndwi2 > 0] <- NA
  names(mplot2) <- c("A", "B", "C", "D")
  mmed2 <- median(as.matrix(mplot2), na.rm=T)
  
  tiff("./outfigs/mmaps2.tif", width = 8, height = 5, units = "in", res = 300, compression = "lzw")
  
  gplot(mplot2, maxpixels=5000000) +
    geom_raster(aes(fill = value)) +
    geom_polygon(data = gacnty_f, aes(x = long, y = lat, group = group), color="black", fill=NA, size=0.5) +
    #scale_fill_gradient2(name = "M", low = "#313695", mid = "#ffffbf", high = "#a50026", midpoint=mmed2) +
    scale_fill_distiller(name = expression(italic(M(T[min], T[max]))), palette = "YlOrRd", direction = 1) +
    coord_fixed(xlim = c(-83.55, -83.24), ylim = c(33.847, 34.044), expand = F, ratio = 1) +
    facet_wrap(~ variable, ncol = 2) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 14, hjust = 0),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 8),
          strip.background = element_rect(fill = NA, color = NA),
          strip.text = element_text(hjust = 0, face = "bold"),
          panel.border = element_rect(fill = NA, color = "black"),
          panel.background = element_rect(fill = "white", color = NA))
  
  dev.off()
  
  
  
####################################################
# Composite map with M, M2, and VC2 for the proposal  
####################################################

july_m <- julystack[[10]]  
july_mr <- july_m / max(as.matrix(july_m), na.rm=T)
july_m2 <- julystack[[11]]  
july_m2r <- july_m2 / max(as.matrix(july_m2), na.rm=T)
july_vc2 <- julystack[[13]]  
july_vc2r <- july_vc2 / max(as.matrix(july_vc2), na.rm=T)

png("./outfigs/sa3fig.png", width = 4, height = 6, units = "in", res = 300)

sa3sfa <- gplot(july_mr, maxpixels=5000000) +
  geom_raster(aes(fill = value)) +
  geom_polygon(data = gacnty_f, aes(x = long, y = lat, group = group), color="black", fill=NA, size=0.5) +
  #scale_fill_gradient2(name = "M", low = "#313695", mid = "#ffffbf", high = "#a50026", midpoint=mmed) +
  scale_fill_distiller(name = "M", palette = "YlOrRd", direction = 1) +
  coord_fixed(xlim = c(-83.55, -83.24), ylim = c(33.847, 34.044), expand = F, ratio = 1) +
  #facet_wrap(~ variable, ncol = 2) +
  labs(title = "A) Mosquito Abundance from VC(T) Model") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 10, hjust = 0),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        strip.background = element_rect(fill = NA, color = NA),
        strip.text = element_text(hjust = 0, face = "bold"),
        panel.border = element_rect(fill = NA, color = "black"),
        panel.background = element_rect(fill = "white", color = NA))
  
sa3sfb <- gplot(july_m2r, maxpixels=5000000) +
  geom_raster(aes(fill = value)) +
  geom_polygon(data = gacnty_f, aes(x = long, y = lat, group = group), color="black", fill=NA, size=0.5) +
  #scale_fill_gradient2(name = "M", low = "#313695", mid = "#ffffbf", high = "#a50026", midpoint=mmed2) +
  scale_fill_distiller(name = "M", palette = "YlOrRd", direction = 1) +
  coord_fixed(xlim = c(-83.55, -83.24), ylim = c(33.847, 34.044), expand = F, ratio = 1) +
  #facet_wrap(~ variable, ncol = 2) +
  labs(title = "B) Mosquito Abundance from Empirical Model") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 10, hjust = 0),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        strip.background = element_rect(fill = NA, color = NA),
        strip.text = element_text(hjust = 0, face = "bold"),
        panel.border = element_rect(fill = NA, color = "black"),
        panel.background = element_rect(fill = "white", color = NA))

sa3sfc <- gplot(july_vc2r, maxpixels=5000000) +
  geom_raster(aes(fill = value)) +
  #geom_polygon(data = blocks_f, aes(x = long, y = lat, group = group), color="blue", fill=NA, size=1) +
  geom_polygon(data = gacnty_f, aes(x = long, y = lat, group = group), color="black", fill=NA, size=0.5) +
  #scale_fill_gradient2(name = "VCr", low = "#313695", mid = "#ffffbf", high = "#a50026", midpoint=vcmed) +
  scale_fill_distiller(name = "VC", palette = "YlOrRd", direction = 1) +
  coord_fixed(xlim = c(-83.55, -83.24), ylim = c(33.847, 34.044), expand = F, ratio = 1) +
  #facet_wrap(~ variable, ncol = 2) +
  labs(title = "C) Vectorial Capacity from Combined Model") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 10, hjust = 0),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        strip.background = element_rect(fill = NA, color = NA),
        strip.text = element_text(hjust = 0, face = "bold"),
        panel.border = element_rect(fill = NA, color = "black"),
        panel.background = element_rect(fill = "white", color = NA))  
    
grid.arrange(sa3sfa, sa3sfb, sa3sfc, ncol = 1)

dev.off()  





mstack <- stack(junestack[[10]], junestack[[11]], julystack[[10]], julystack[[11]],
                augstack[[10]], augstack[[11]], septstack[[10]], septstack[[11]])
names(mstack) <- c("junem", "junem2", "julym", "julym2", "augm", "augm2", "septm", "septm2")
mframe <- data.frame(na.omit(values(mstack)))

junelm <- lm(junem2 ~ junem, data = mframe)
julylm <- lm(julym2 ~ julym, data = mframe)
auglm <- lm(augm2 ~ augm, data = mframe)
septlm <- lm(septm2 ~ septm, data = mframe)

junep <- predict(mstack, junelm)
julyp <- predict(mstack, julylm)
augp <- predict(mstack, auglm)
septp <- predict(mstack, septlm)

juner <- mstack$junem2 - junep
julyr <- mstack$julym2 - julyp
augr <- mstack$augm2 - augp
septr <- mstack$septm2 - septp


f <- stack(system.file("external/rlogo.grd", package="raster")) 
s <- stack(f)
names(s)
v <- data.frame(na.omit(values(s)))
m <- lm(red ~ green, data=v)
m

p <- predict(s, m)
residuals <- s$red - p



##############
# Overlay Maps
##############
vcjune_fm <- focal(junestack[[13]], focalWeight(junestack[[13]], 0.001, "circle"))
vcjuly_fm <- focal(julystack[[13]], focalWeight(julystack[[13]], 0.001, "circle"))
vcaug_fm <- focal(augstack[[13]], focalWeight(augstack[[13]], 0.001, "circle"))
vcsept_fm <- focal(septstack[[13]], focalWeight(septstack[[13]], 0.001, "circle"))
vcplot_fm <- stack(vcjune_fm, vcjuly_fm, vcaug_fm, vcsept_fm)
vcmax <- quantile(as.matrix(vcplot_fm), probs=0.99, na.rm=T)
vcmin <- quantile(as.matrix(vcplot_fm), probs = 0.01, na.rm=T)
vctrunc_fm <- vcplot_fm
vctrunc_fm[vctrunc_fm < vcmin] <- vcmin
vctrunc_fm[vctrunc_fm > vcmax] <- vcmax
vcindex_fm <- (vctrunc_fm - vcmin)/(vcmax - vcmin)

popdens_fm <- focal(popdens, focalWeight(popdens, 0.001, "circle"))
logpop_fm <- log(popdens_fm + 1)
popmax <- quantile(as.matrix(logpop_fm), probs=0.99, na.rm=T)
poptrunc_fm <- logpop_fm
poptrunc_fm[poptrunc_fm > popmax] <- popmax
pindex_fm <- poptrunc_fm / popmax

vcrisk <- vcindex_fm * pindex_fm

png("./outfigs/vcpopover.png", width = 8, height = 4.5, units = "in", res = 300)

gplot(vcrisk, maxpixels=5000000) +
  geom_raster(aes(fill = value)) +
  #geom_polygon(data = blocks_f, aes(x = long, y = lat, group = group), color="blue", fill=NA, size=1) +
  geom_polygon(data = gacnty_f, aes(x = long, y = lat, group = group), color="black", fill=NA, size=0.5) +
  #scale_fill_gradient2(name = "VC", low = "#ffffcc", mid = "#fd8d3c", high = "#800026", midpoint=vcmed) +
  scale_fill_gradientn(name = "VCr", colors=brewer.pal(9, "YlOrBr")) +
  coord_fixed(xlim = c(-83.55, -83.24), ylim = c(33.847, 34.044), expand = F, ratio = 1) +
  facet_wrap(~ variable, ncol = 2) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 14, hjust = 0),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        strip.background = element_rect(fill = NA, color = NA),
        strip.text = element_text(hjust = 0, face = "bold"),
        panel.border = element_rect(fill = NA, color = "black"),
        panel.background = element_rect(fill = "white", color = NA))



dev.off()




############################
# Vectorial Capacity Maps 2
############################

vczoom_fm <- vcplot_fm
vczoom_fm[popdens < 10] <- NA
#vczoom_fm <- vczoom_fm / vcmax

png("vcmaps2.png", width = 8, height = 4.5, units = "in", res = 300)

gplot(vczoom_fm, maxpixels=5000000) +
  geom_raster(aes(fill = value)) +
  #geom_polygon(data = blocks_f, aes(x = long, y = lat, group = group), color="blue", fill=NA, size=1) +
  geom_polygon(data = gacnty_f, aes(x = long, y = lat, group = group), color="black", fill=NA, size=0.5) +
  #scale_fill_gradient2(name = "VC", low = "#ffffcc", mid = "#fd8d3c", high = "#800026", midpoint=vcmed) +
  scale_fill_gradientn(name = "VCr", colors=brewer.pal(9, "YlOrBr")) +
  coord_fixed(xlim = c(-83.48, -83.31), ylim = c(33.89, 34.01), expand = F, ratio = 1) +
  facet_wrap(~ variable, ncol = 2) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 14, hjust = 0),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        strip.background = element_rect(fill = NA, color = NA),
        strip.text = element_text(hjust = 0, face = "bold"),
        panel.border = element_rect(fill = NA, color = "black"),
        panel.background = element_rect(fill = "white", color = NA))

dev.off()

####################################
# Compute and Map Force of Infection
####################################



junefoi <- junestack[[9]] * junestack[[1]] * (1 - junestack[[11]] ^ (1 / junestack[[8]]))
julyfoi <- julystack[[9]] * julystack[[1]] * (1 - julystack[[11]] ^ (1 / julystack[[8]]))
augfoi <- augstack[[9]] * augstack[[1]] * (1 - augstack[[11]] ^ (1 / augstack[[8]]))
septfoi <- septstack[[9]] * septstack[[1]] * (1 - septstack[[11]] ^ (1 / septstack[[8]]))
junefoi_fm <- focal(junefoi, focalWeight(junefoi, 0.001, "circle"))
julyfoi_fm <- focal(julyfoi, focalWeight(julyfoi, 0.001, "circle"))
augfoi_fm <- focal(augfoi, focalWeight(augfoi, 0.001, "circle"))
septfoi_fm <- focal(septfoi, focalWeight(septfoi, 0.001, "circle"))
foistack_fm <- stack(junefoi_fm, julyfoi_fm, augfoi_fm, septfoi_fm)
names(foistack_fm) <- c("A", "B", "C", "D")


png("foimaps.png", width = 8, height = 4.5, units = "in", res = 300)

gplot(foistack_fm, maxpixels=5000000) +
  geom_raster(aes(fill = value)) +
  #geom_polygon(data = blocks_f, aes(x = long, y = lat, group = group), color="blue", fill=NA, size=1) +
  geom_polygon(data = gacnty_f, aes(x = long, y = lat, group = group), color="black", fill=NA, size=0.5) +
  #scale_fill_gradient2(name = "VC", low = "#ffffcc", mid = "#fd8d3c", high = "#800026", midpoint=vcmed) +
  scale_fill_gradientn(name = "FOI", colors=brewer.pal(9, "YlOrBr")) +
  coord_fixed(xlim = c(-83.55, -83.24), ylim = c(33.847, 34.044), expand = F, ratio = 1) +
  facet_wrap(~ variable, ncol = 2) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 14, hjust = 0),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        strip.background = element_rect(fill = NA, color = NA),
        strip.text = element_text(hjust = 0, face = "bold"),
        panel.border = element_rect(fill = NA, color = "black"),
        panel.background = element_rect(fill = "white", color = NA))

dev.off()

############################
# Land cover versus VC plots
############################

# Set tree cover and impervious to zero for water pixels
treecov[ndwi > 0] <- 0
imperv[ndwi > 0] <- 0

# Scale up original 10 m resolution raster to 30 m
treecov30 <- aggregate(treecov, fact=3)
imperv30 <- aggregate(imperv, fact=3)

vcplot <- stack(junestack[[13]], julystack[[13]], augstack[[13]], septstack[[13]])
names(vcplot) <- c("June-July", "July-Aug", "Aug-Sept", "Sept-Oct")
vcmat <- as.matrix(vcplot)
vcmax <- max(vcmat, na.rm=T)
vcmin <- min(vcmat, na.rm=T)
vcplot <- (vcplot - vcmin) / (vcmax - vcmin)

imperv1k <- focal(imperv30, focalWeight(imperv30, 0.01, "circle"))

lcstack <- stack(treecov30, imperv1k, vcplot)


clarke <- gacnty_sp[gacnty_sp@data$COUNTYFP10 == "059",]
clarke_r <- rasterize(clarke, vcplot, field = "COUNTYFP10")
clarke_r2 <- projectRaster(clarke_r, vcplot)
lcstack <- mask(lcstack, clarke_r2)

sampdata <- sampleRandom(lcstack, 100000)
sampdata <- as.data.frame(sampdata)
names(sampdata)[1:2] <- c("treecov", "imperv")
summary(sampdata)

loess1 <- loess(June.July ~ treecov * imperv, data = sampdata)
loess2 <- loess(July.Aug ~ treecov * imperv, data = sampdata)
loess3 <- loess(Aug.Sept ~ treecov * imperv, data = sampdata)
loess4 <- loess(Sept.Oct ~ treecov * imperv, data = sampdata)

treegrid <-  seq(1, 99, 0.5)
impervgrid <-  seq(2, 48, 0.5)
# Generate a dataframe with every possible combination of wt and hp
data.fit <-  expand.grid(treecov = treegrid, imperv = impervgrid)
# Feed the dataframe into the loess model and receive a matrix output with estimates of 
# acceleration for each combination of wt and hp
loessp1 <-  predict(loess1, newdata = data.fit)
loessp2 <-  predict(loess2, newdata = data.fit)
loessp3 <-  predict(loess3, newdata = data.fit)
loessp4 <-  predict(loess4, newdata = data.fit)


plotdat1 <- data.frame(data.fit, vc = as.numeric(loessp1), month = rep("A", nrow(data.fit)))
plotdat2 <- data.frame(data.fit, vc = as.numeric(loessp2), month = rep("B", nrow(data.fit)))
plotdat3 <- data.frame(data.fit, vc = as.numeric(loessp3), month = rep("C", nrow(data.fit)))
plotdat4 <- data.frame(data.fit, vc = as.numeric(loessp4), month = rep("D", nrow(data.fit)))

plotdat <- bind_rows(plotdat1, plotdat2, plotdat3, plotdat4)

#png("./outfigs/lcmaps.png", width = 7, height = 6, units = "in", res = 300)
tiff("./outfigs/Fig_6.tif", width = 7, height = 6, units = "in", res = 300, compression="lzw")

ggplot(plotdat, aes(x = treecov, y = imperv, z = vc)) +
  stat_contour(geom = "polygon", aes(fill = ..level..)) +
  geom_tile(aes(fill = vc)) +
  stat_contour(binwidth = 0.05) +
  scale_fill_distiller(name = expression(italic("VC(T)")), palette = "RdYlBu") +
  scale_x_continuous(limits = c(0, 98), expand = c(0, 0)) +
  scale_y_continuous(limits = c(2, 48), expand = c(0, 0)) +
  labs(x="Tree Cover (%)", y = "Impervious Surface (%)") +
  facet_wrap(~ month) +
  theme(strip.background = element_rect(fill = NA, color = NA),
        strip.text = element_text(hjust = 0, face = "bold", size = 14),
        panel.border = element_rect(fill = NA, color = "black"),
        panel.background = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 12))

  
dev.off()

g <- ggplot(faithfuld, aes(waiting, eruptions))

g + geom_raster(aes(fill = density))


ggplot(faithfuld, aes(waiting, eruptions))

g + geom_raster(aes(fill = density))





vczmean_fm <- mean(vczoom_fm)

gplot(vczmean_fm, maxpixels=5000000) +
  geom_raster(aes(fill = value)) +
  #geom_polygon(data = blocks_f, aes(x = long, y = lat, group = group), color="blue", fill=NA, size=1) +
  geom_polygon(data = gacnty_f, aes(x = long, y = lat, group = group), color="black", fill=NA, size=0.5) +
  #scale_fill_gradient2(name = "VC", low = "#ffffcc", mid = "#fd8d3c", high = "#800026", midpoint=vcmed) +
  scale_fill_gradientn(name = "VC", colors=brewer.pal(9, "YlOrBr")) +
  coord_fixed(xlim = c(-83.48, -83.30), ylim = c(33.90, 34.01), expand = F, ratio = 1) 
  #facet_wrap(~ variable, ncol = 2)







vcmean_fm <- mean(vcplot)
gplot(vcmean_fm, maxpixels=5000000) +
  geom_raster(aes(fill = value)) +
  #geom_polygon(data = blocks_f, aes(x = long, y = lat, group = group), color="blue", fill=NA, size=1) +
  geom_polygon(data = gacnty_f, aes(x = long, y = lat, group = group), color="black", fill=NA, size=0.5) +
  #scale_fill_gradient2(name = "VC", low = "#ffffcc", mid = "#fd8d3c", high = "#800026", midpoint=vcmed) +
  scale_fill_gradientn(name = "VC", colors=brewer.pal(9, "YlOrBr")) +
  coord_fixed(xlim = c(-83.55, -83.24), ylim = c(33.847, 34.044), expand = F, ratio = 1) +
  facet_wrap(~ variable, ncol = 2)


#################
# Single Risk Map
#################
vcmean_fm <- mean(vcplot_fm)

vcmax <- max(as.matrix(vcmean_fm), na.rm=T)
vcmin <- min(as.matrix(vcmean_fm), na.rm=T)
vcindex_fm <- (vcmean_fm - vcmin)/(vcmax - vcmin)

popdens_fm <- focal(popdens, focalWeight(popdens, 0.001, "circle"))

pmax <- max(as.matrix(log(popdens_fm + 1)), na.rm=T)
pmin <- min(as.matrix(log(popdens_fm + 1)), na.rm=T)
pindex_fm <- (log(popdens_fm + 1) - pmin)/(pmax - pmin)

vcrisk <- sqrt(vcindex_fm * pindex_fm)




vc_matrix <- as.matrix(vcmean_fm)
vc_vector <- as.numeric(vc_matrix)[!is.na(as.numeric(vc_matrix))]
vc_quantile <- quantile(vc_vector, seq(0, 1, 0.1))

vc_reclass <- matrix(c(vc_quantile[1:10], vc_quantile[2:11], 1:10), ncol = 3)

vcindex_fm <- reclassify(vcmean_fm, vc_reclass, include.lowest = T)


pop_matrix <- as.matrix(popdens_fm)
pop_vector <- as.numeric(pop_matrix)[!is.na(as.numeric(pop_matrix))]
pop_quantile <- quantile(pop_vector, seq(0, 1, 0.1))

pop_reclass <- matrix(c(pop_quantile[1:10], pop_quantile[2:11], 1:10), ncol = 3)
pindex_fm <- reclassify(popdens_fm, pop_reclass, include.lowest = T)





#pindex_fm <- logpop_fm / max(as.matrix(logpop_fm), na.rm=T)
#pindex_fm <- scale(logpop_fm)
temp <- popdens
temp[temp > 20] <- 20
pindex_fm <- temp / 20

vcrisk <- sqrt(vcindex_fm * pindex_fm)
vcrisk[ndwi2 > 0] <- NA

#png("./outfigs/riskmap.png", width = 8, height = 4.5, units = "in", res = 300)
tiff("./outfigs/Fig_7.tif", width = 8, height = 4.5, units = "in", res = 300, compression="lzw")

gplot(vcrisk, maxpixels=5000000) +
  geom_raster(aes(fill = value)) +
  #geom_polygon(data = blocks_f, aes(x = long, y = lat, group = group), color="blue", fill=NA, size=1) +
  geom_polygon(data = gacnty_f, aes(x = long, y = lat, group = group), color="black", fill=NA, size=0.5) +
  scale_fill_gradientn(name = "Index", colors=(brewer.pal(9, "YlOrRd"))) +
  coord_fixed(xlim = c(-83.55, -83.24), ylim = c(33.847, 34.044), expand = F, ratio = 1) +
  #facet_wrap(~ variable, ncol = 2)
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 10, hjust = 0),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        strip.background = element_rect(fill = NA, color = NA),
        strip.text = element_text(hjust = 0, face = "bold"),
        panel.border = element_rect(fill = NA, color = "black"),
        panel.background = element_rect(fill = "white", color = NA))  
dev.off()




