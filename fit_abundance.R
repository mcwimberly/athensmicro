rm(list=ls())

library(ggplot2)
library(mgcv)
library(tidyr)
library(dplyr)
library(reshape2)
library(MASS)
library(quantreg)
library(splines)
library(raster)
library(pracma)
library(caret)

library(dplyr)
library(zoo)

albo <- read.csv(".\\michelle_regression\\albo-abundance.csv")
tlog  <- read.csv(".\\michelle_regression\\loggers-site-daily.csv")


#test2<-arrange(test,ID,YEAR_VISIT) %>% group_by(subject)%>%
#  mutate(ma2=rollapply(BLOOD_PRESSURE,2,mean,align='right',fill=NA))


albo <- albo %>%
  mutate(date = as.Date(Date))

tsum <- tlog %>%
  mutate(date = as.Date(day)) %>%
  arrange(Site, day) %>%
  group_by(Site) %>%
  mutate(MA_tmin = rollapply(Temp_min, 7, mean, align='right', fill=NA),
         MA_tmax = rollapply(Temp_max, 7, mean, align='right', fill=NA))

albo <- albo %>%
  left_join(tsum,
            by=c("date"="date",
                 "Site"="Site")) %>%
  dplyr::select(one_of("Site", "date", "count", "MA_tmin", "MA_tmax"))
  
tegam    <- gam(count ~ s(MA_tmax, MA_tmin), data=albo, family=nb())
adgam    <- gam(count ~ s(MA_tmax, bs="tp") + s(MA_tmin, bs="tp"), data=albo, family=nb())

tegam
adgam

saveRDS(adgam, "./outdata/minmaxmod.rds")

#Cross-validation
sites <- unique(albo$Site)
xval <- data.frame()
for (site in sites) {
  fitframe <- filter(albo, Site != site)
  predframe <- filter(albo, Site == site)
  curfit <- gam(count ~ s(MA_tmax, bs="tp") + s(MA_tmin, bs="tp"), data=fitframe, family=nb())
  curpred <- predict(curfit, predframe, type="response")
  curout <- data.frame(predframe$count/2, curpred/2)
  names(curout) <- c("obs", "pred")
  xval <- rbind(xval, curout)
}
mean(xval$obs - xval$pred)
mean(abs(xval$obs - xval$pred))
median(abs(xval$obs - xval$pred))
cor(xval$obs, xval$pred, method="spearman")


logger_sites <- groupKFold(albo$Site)
train_control <- trainControl(method="cv", number = 9, index= logger_sites)

albo_cv <- train(count ~ MA_tmax + MA_tmin, data=albo, family=nb(), trControl=train_control, method="gam", bs = "tp")



predframe <- expand.grid(MA_tmax=seq(from=min(albo$MA_tmax, na.rm=TRUE),
                                     to  =max(albo$MA_tmax, na.rm=TRUE),
                                     length.out=500),
                         MA_tmin=seq(from=min(albo$MA_tmin, na.rm=TRUE),
                                     to  =max(albo$MA_tmin, na.rm=TRUE),
                                     length.out=500))
predframe$preds <- predict(adgam, newdata=predframe, type="response")

maxdiurn <- max(albo$MA_tmax - albo$MA_tmin, na.rm=TRUE)
mindiurn <- min(albo$MA_tmax - albo$MA_tmin, na.rm=TRUE)

predframe$preds[predframe$MA_tmin > (predframe$MA_tmax - mindiurn)] <- NA
predframe$preds[predframe$MA_tmax > (predframe$MA_tmin + maxdiurn)] <- NA

# predframe$cut <- cut(predframe$preds,
#                      #breaks=c(0,5,10,20,100),
#                      breaks=seq(from=0, to=100, by=2),
#                      include.lower=TRUE)

tiff("./outfigs/mosquito_abundance.tif", width = 6, height = 6, units = "in", res = 300, compression="lzw")

ggplot() + geom_tile(data=predframe, aes(x=MA_tmin, y=MA_tmax, fill=preds/2)) +
  geom_abline(slope=1, intercept=mindiurn, linetype=2, color="white", size=1) + 
  geom_abline(slope=1, intercept=maxdiurn, linetype=2, color="white", size=1) + 
  geom_point(data=albo, aes(x=MA_tmin, y=MA_tmax, size=count), color="black") + 
  #xlab(expression(paste(italic(T[min]), " ", ( degree*C)))) + 
  #ylab(expression(paste(italic(T[max]), " ", ( degree*C)))) + 
  xlab("Minimum Temperature (\u00B0C)") + 
  ylab("Maximum Temperature (\u00B0C)") + 
  scale_fill_distiller(name = "M", palette = "RdYlBu",
                       na.value = NA) +
  #scale_fill_gradient2(low="lightblue", high="red", mid="lightyellow", midpoint=10,
  #                     name="M",
  #                     na.value=NA) +
  scale_size_continuous(guide=FALSE) + 
  scale_x_continuous(expand = c(0.012, 0.012)) +
  scale_y_continuous(expand = c(0.012, 0.012)) +
  theme(legend.position="bottom",
        legend.key.width=unit(2, "cm"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) + 
  coord_equal() 

dev.off()
