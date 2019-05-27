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

junestack <- stack("vcstack_1.tif")
julystack <- stack("vcstack_2.tif")
augstack <- stack("vcstack_3.tif")
septstack <- stack("vcstack_4.tif")
gacnty_sp <- readOGR(layer = "GA_SHP", dsn=".")
gacnty_f <- fortify(gacnty_sp, regions="GEOID10")
#tminplot <- stack("tminplot.tif")
#tmaxplot <- stack("tmaxplot.tif")
tempminstack <- brick("tempminstack.tif")
tempmaxstack <- brick("tempmaxstack.tif")
rhminstack <- brick("rhminstack.tif")
tmingmet <- stack("tmincrop.tif")
tmaxgmet <- stack("tmaxcrop.tif")

treecov <- raster("Athens_treecover_v2.tif")
imperv <- raster("Athens_impervious_v2.tif")
ndwi <- raster("Athens_NDWI.tif")

sentinel <- stack("SentinelComposite.tif")
popdens <- raster("athens_popdens.tif")

treecov[ndwi > 0] <- NA
treefig <- gplot(treecov, maxpixels = 5000000) +
  geom_raster(aes(fill = value)) +
  geom_polygon(data = gacnty_f, aes(x = long, y = lat, group = group), color="black", fill=NA, size=1.0) +
  scale_fill_gradientn(name = "% Cover", colors=brewer.pal(9, "Greens")[1:8]) +
  coord_fixed(xlim = c(-83.55, -83.24), ylim = c(33.847, 34.044), expand = F, ratio = 1) +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "white", color = NA))
imperv[ndwi > 0] <- NA
impervfig <- gplot(imperv, maxpixels = 5000000) +
  geom_raster(aes(fill = value)) +
  geom_polygon(data = gacnty_f, aes(x = long, y = lat, group = group), color="black", fill=NA, size=1.0) +
  scale_fill_gradientn(name = "% Cover", colors=brewer.pal(9, "Reds")) +
  coord_fixed(xlim = c(-83.55, -83.24), ylim = c(33.847, 34.044), expand = F, ratio = 1) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "white", color = NA))

senfig <- ggRGB(sentinel, r=1, g=2, b=3, stretch = 'lin', maxpixels = 5000000) +
  geom_polygon(data = gacnty_f, aes(x = long, y = lat, group = group), color="black", fill=NA, size=1.0) +
  coord_fixed(xlim = c(-83.55, -83.24), ylim = c(33.847, 34.044), expand = F, ratio = 1) +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "white", color = NA))

popdens2 <- log10(popdens + 1)
popfig <- gplot(popdens2, maxpixels = 5000000) +
  geom_raster(aes(fill = value)) +
  geom_polygon(data = gacnty_f, aes(x = long, y = lat, group = group), color="black", fill=NA, size=1.0) +
  scale_fill_gradientn(name = "Density", colors=brewer.pal(9, "YlOrRd")) +
  coord_fixed(xlim = c(-83.55, -83.24), ylim = c(33.847, 34.044), expand = F, ratio = 1) +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "white", color = NA))

grid.arrange(senfig, treefig, impervfig, popfig, ncol = 2)

###############################################################################
# Temperature time series figure
###############################################################################
uga_met <- read_csv("uga_ggy_weather.csv")


uga_2018 <- uga_met %>%
  mutate(Date = mdy(Date),
         Year = year(Date),
         Month = month(Date),
         Doy = yday(Date),
         temp = (temp - 32) * 5/9,
         temp_lo = (temp_lo - 32) * 5/9,
         temp_hi = (temp_hi - 32) * 5/9) %>%
  filter(Year == 2018 & Doy >= 166 & Doy <= (166 + 117))


micro <- read_csv("micro_2018_daily.csv")

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

png("tempts.png", width = 6, height = 3.5, units = "in", res = 300)


ggplot() +
  geom_ribbon(data = filter(microsum, dailycount >= 0), aes(x = Date, ymin = dailymin_min, ymax = dailymin_max), fill = "orange") +
  geom_ribbon(data = filter(microsum, dailycount >= 0), aes(x = Date, ymin = dailymax_min, ymax = dailymax_max), fill = "indianred1") +
  geom_line(data = filter(microsum, dailycount >= 0), aes(x = Date, y = dailymax_med), color = "red3", size = 0.5) +
  geom_line(data = filter(microsum, dailycount >= 0), aes(x = Date, y = dailymin_med), color = "orange3", size = 0.5) +
  geom_line(data = uga_2018, aes(x = Date, y = temp_lo), color = "black", size = 0.5, linetype = 2) +
  geom_line(data = uga_2018, aes(x = Date, y = temp_hi), color = "black", size = 0.5, linetype = 2)  +
  labs(y = "Temperature (\u00B0C)", x = "") +
  theme_bw()

dev.off()


###############################################################################
# Temperature Map Figure
###############################################################################

png("tempmaps.png", width = 8, height = 4.5, units = "in", res = 300)

fig2a <- gplot(tempminstack[[1 + 29]], maxpixels = 5000000) +
  geom_raster(aes(fill = value)) +
  geom_polygon(data = gacnty_f, aes(x = long, y = lat, group = group), color="black", fill=NA, size=0.5) +
  scale_fill_gradientn(name = "Temp (C)", colors=brewer.pal(9, "Oranges")) +
  coord_fixed(xlim = c(-83.55, -83.24), ylim = c(33.847, 34.044), expand = F, ratio = 1) +
  labs(title = "a") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 14, hjust = 0),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        panel.border = element_rect(fill = NA, color = "black"),
        panel.background = element_rect(fill = "white", color = NA))

fig2b <- gplot(tempmaxstack[[1 + 29]], maxpixels = 5000000) +
  geom_raster(aes(fill = value)) +
  geom_polygon(data = gacnty_f, aes(x = long, y = lat, group = group), color="black", fill=NA, size=0.5) +
  scale_fill_gradientn(name = "Temp (C)", colors=brewer.pal(9, "Reds")) +
  coord_fixed(xlim = c(-83.55, -83.24), ylim = c(33.847, 34.044), expand = F, ratio = 1) +
  labs(title = "b") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 14, hjust = 0),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        panel.border = element_rect(fill = NA, color = "black"),
        panel.background = element_rect(fill = "white", color = NA))

fig2c <- gplot(tmingmet[[166 + 29]] - 273.15, maxpixels = 5000000) +
  geom_raster(aes(fill = value)) +
  geom_polygon(data = gacnty_f, aes(x = long, y = lat, group = group), color="black", fill=NA, size=0.5) +
  scale_fill_gradientn(name = "Temp (C)", colors=brewer.pal(9, "Oranges")) +
  coord_fixed(xlim = c(-83.55, -83.24), ylim = c(33.847, 34.044), expand = F, ratio = 1) +
  labs(title = "c") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 14, hjust = 0),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        panel.border = element_rect(fill = NA, color = "black"),
        panel.background = element_rect(fill = "white", color = NA))

fig2d <- gplot(tmaxgmet[[166 + 29]] - 273.15, maxpixels = 5000000) +
  geom_raster(aes(fill = value)) +
  geom_polygon(data = gacnty_f, aes(x = long, y = lat, group = group), color="black", fill=NA, size=0.5) +
  scale_fill_gradientn(name = "Temp (C)", colors=brewer.pal(9, "Reds")) +
  coord_fixed(xlim = c(-83.55, -83.24), ylim = c(33.847, 34.044), expand = F, ratio = 1) +
  labs(title = "d") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 14, hjust = 0),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        panel.border = element_rect(fill = NA, color = "black"),
        panel.background = element_rect(fill = "white", color = NA))

grid.arrange(fig2a, fig2b, fig2c, fig2d, ncol = 2)

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
# VC Histogram Plot
###############################################################################
uga_met_vc <- read_csv("uga_met_vc.csv")

vcplot <- stack(junestack[[10]], julystack[[10]], augstack[[10]], septstack[[10]])
names(vcplot) <- c("June-July", "July-Aug", "Aug-Sept", "Sept-Oct")

clarke <- gacnty_sp[gacnty_sp@data$COUNTYFP10 == "059",]
clarke_r <- rasterize(clarke, vcplot, field = "COUNTYFP10")
clarke_r2 <- projectRaster(clarke_r, vcplot)
vcplot <- mask(vcplot, clarke_r2)

vczoom <- vcplot
vczoom[popdens < 10] <- NA

vcmat <- as.matrix(vcplot)
vcmax <- max(vcmat, na.rm=T)

vcplot <- vcplot / vcmax
vczoom <- vczoom / vcmax

vcsample <- sampleRandom(vcplot, size = 1000000)
vcsample <- as.data.frame(vcsample)

vcsample2 <- sampleRandom(vczoom, size = 1000000)
vcsample2 <- as.data.frame(vcsample2)

vcsamplot <- vcsample %>%
  gather(key = "month", value = "vc", June.July, July.Aug, Aug.Sept, Sept.Oct) %>%
  mutate(month = factor(month, levels = c("June.July", "July.Aug", "Aug.Sept", "Sept.Oct"), 
                        labels = c("June-July", "July-Aug", "Aug-Sept", "Sept-Oct")))


vcsamplot2 <- vcsample2 %>%
  gather(key = "month", value = "vc", June.July, July.Aug, Aug.Sept, Sept.Oct) %>%
  mutate(month = factor(month, levels = c("June.July", "July.Aug", "Aug.Sept", "Sept.Oct"), 
                        labels = c("June-July", "July-Aug", "Aug-Sept", "Sept-Oct")))


vline.data <- data.frame(z = uga_met_vc$VC, month=c("June-July","July-Aug","Aug-Sept","Sept-Oct"))
vline.data[,1] <- vline.data[,1] / vcmax

png("vchist.png", width = 4, height = 6, units = "in", res = 300)

ggplot(data = vcsamplot) +
  geom_histogram(aes(x = vc), fill = "skyblue2", bins = 20) +
  geom_histogram(aes(x = vc), vcsamplot2, fill = "indianred2", bins = 20) +
  geom_vline(aes(xintercept = z), vline.data, color = "black") +
  #scale_y_continuous(labels = percent_format(accuracy = 2)) +
  facet_wrap(~ month, ncol = 1, scales = "free_y") +
  theme_bw() + 
  labs(x="VCr", y = "Count") 

dev.off()




  
#########################
# Vectorial Capacity Maps
#########################

vcjune_fm <- focal(junestack[[10]], focalWeight(junestack[[10]], 0.001, "circle"))
vcjuly_fm <- focal(julystack[[10]], focalWeight(julystack[[10]], 0.001, "circle"))
vcaug_fm <- focal(augstack[[10]], focalWeight(augstack[[10]], 0.001, "circle"))
vcsept_fm <- focal(septstack[[10]], focalWeight(septstack[[10]], 0.001, "circle"))
vcplot_fm <- stack(vcjune_fm, vcjuly_fm, vcaug_fm, vcsept_fm)
names(vcplot_fm) <- c("A", "B", "C", "D")
vcplot_fm <- vcplot_fm / vcmax
#vcmed <- median(as.matrix(vcplot), na.rm=T)

png("vcmaps.png", width = 8, height = 4.5, units = "in", res = 300)

gplot(vcplot_fm, maxpixels=5000000) +
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


# Land cover versus CV plots

# Set tree cover and impervious to zero for water pixels
treecov[ndwi > 0] <- 0
imperv[ndwi > 0] <- 0

# Scale up original 10 m resolution raster to 30 m
treecov30 <- aggregate(treecov, fact=3)
imperv30 <- aggregate(imperv, fact=3)

imperv1k <- focal(imperv30, focalWeight(imperv30, 0.01, "circle"))
lcstack <- stack(treecov30, imperv1k, vcplot)

lcstack <- mask(lcstack, clarke_r2)

sampdata <- sampleRandom(lcstack, 25000)
sampdata <- as.data.frame(sampdata)
names(sampdata)[1:2] <- c("treecov", "imperv")

loess1 <- loess(June.July ~ treecov * imperv, data = sampdata)
loess2 <- loess(July.Aug ~ treecov * imperv, data = sampdata)
loess3 <- loess(Aug.Sept ~ treecov * imperv, data = sampdata)
loess4 <- loess(Sept.Oct ~ treecov * imperv, data = sampdata)

treegrid <-  seq(0, 98, 0.5)
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

png("lcmaps.png", width = 8, height = 4.5, units = "in", res = 300)

ggplot(plotdat, aes(x = treecov, y = imperv, z = vc)) +
  stat_contour(geom = "polygon", aes(fill = ..level..)) +
  geom_tile(aes(fill = vc)) +
  stat_contour(bins = 15) +
  scale_fill_distiller(name = "VCr", palette = "RdYlBu") +
  scale_x_continuous(limits = c(0, 98), expand = c(0, 0)) +
  scale_y_continuous(limits = c(2, 48), expand = c(0, 0)) +
  labs(x="Tree Cover (%)", y = "Impervious Surface (%)") +
  facet_wrap(~ month) +
  theme(strip.background = element_rect(fill = NA, color = NA),
        strip.text = element_text(hjust = 0, face = "bold"),
        panel.border = element_rect(fill = NA, color = "black"),
        panel.background = element_rect(fill = "white", color = NA))

  
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







vcmean_fm <- mean(vcplot_fm)
gplot(vcmean_fm, maxpixels=5000000) +
  geom_raster(aes(fill = value)) +
  #geom_polygon(data = blocks_f, aes(x = long, y = lat, group = group), color="blue", fill=NA, size=1) +
  geom_polygon(data = gacnty_f, aes(x = long, y = lat, group = group), color="black", fill=NA, size=0.5) +
  #scale_fill_gradient2(name = "VC", low = "#ffffcc", mid = "#fd8d3c", high = "#800026", midpoint=vcmed) +
  scale_fill_gradientn(name = "VC", colors=brewer.pal(9, "YlOrBr")) +
  coord_fixed(xlim = c(-83.55, -83.24), ylim = c(33.847, 34.044), expand = F, ratio = 1) +
  facet_wrap(~ variable, ncol = 2)


vcmax <- max(as.matrix(vcmean_fm), na.rm=T)
vcmin <- min(as.matrix(vcmean_fm), na.rm=T)
vcindex_fm <- (vcmean_fm - vcmin)/(vcmax - vcmin)

popdens_fm <- focal(popdens, focalWeight(popdens, 0.001, "circle"))
logpop_fm <- log(popdens_fm + 1)


pindex_fm <- logpop_fm / max(as.matrix(logpop_fm), na.rm=T)

vcrisk <- vcindex_fm * pindex_fm

riskmed <- median(as.matrix(vcrisk), na.rm=T)
gplot(vcrisk, maxpixels=5000000) +
  geom_raster(aes(fill = value)) +
  #geom_polygon(data = blocks_f, aes(x = long, y = lat, group = group), color="blue", fill=NA, size=1) +
  geom_polygon(data = gacnty_f, aes(x = long, y = lat, group = group), color="black", fill=NA, size=0.5) +
  scale_fill_gradientn(name = "Risk", colors=rev(brewer.pal(11, "RdYlBu"))) +
  coord_fixed(xlim = c(-83.55, -83.24), ylim = c(33.847, 34.044), expand = F, ratio = 1) +
  facet_wrap(~ variable, ncol = 2)





