
library(dplyr)
library(lubridate)
library(viridis)
library(RColorBrewer)

rm(list=ls())

d3 <- read.csv("~/Documents/Research/WhiteSturgeonOccupancy/white_sturgeon_mm_fixed.csv")
receiverDeployment <- read.csv("~/Documents/Research/WhiteSturgeonOccupancy/receiver_deployments.csv")

receiverRegion <- data.frame(Location=receiverDeployment$Location, 
                             Region=receiverDeployment$Region)
receiverRegion <- unique(receiverRegion)

detects <- left_join(d3, receiverRegion, by=join_by(DetectLocation==Location))
detects$DetectDate <- as.Date(detects$DetectDate)
detects$tagging_date <- as.Date(detects$tagging_date)
detects$DetectDoy <- yday(detects$DetectDate)
detects <- detects[detects$DetectDate >= detects$tagging_date,]

detects<-detects[!detects$StudyID %in% c("DWR_FRP_sturgeon_AS","DWR-FR-Sturgeon"),]

rm(d3)

### Sketchy--histograms of total detections over time in SF Bay and San Pablo Bay -----------------

## san francisco bay
sfbaydetects <- detects[detects$Region=="SF_Bay",]

hist(sfbaydetects$DetectDoy)
abline(v=c(214,244))


## san pablo bay
sanpablodetects <- detects[detects$Region=="San_Pablo_Bay",]

hist(sanpablodetects$DetectDoy)
abline(v=c(214,244))



## Check deployment histories ---------------------------------------------------------------------- 
# 
# receiverDeployment$GeneralLocation <- gsub('[[:digit:]]+', '', receiverDeployment$Location)
# receiverDeployment <- receiverDeployment[receiverDeployment$Region!="",]
# receiverDeployment$RiverKm <- as.numeric(receiverDeployment$RiverKm)
# 
# receiverDeployment$Region[receiverDeployment$Region=="carquinez_Strait"] <- "Carquinez_Strait"
# receiverDeployment$Region[receiverDeployment$Region=="Suisun Bay"] <- "Suisun_Bay"
# receiverDeployment$Region<-factor(receiverDeployment$Region, levels=c("Pacific_Coast","Pt_Reyes","Gulf_Farallones","Golden_Gate","SF_Bay","San_Pablo_Bay",
#                                   "Carquinez_Strait","Suisun_Bay","West_Delta","Central_Delta","North_Delta","Sac_River","Amr_River","Feather_River","Battle_Ck",
#                                   "South Delta","SJ_River","Mok_River","Stan_River","Tuol_River","Merced_River"))
# 
# receiverDeployment$Start <- as.Date(receiverDeployment$Start)
# receiverDeployment$Stop <- as.Date(receiverDeployment$Stop)
# receiverDeployment <- receiverDeployment[order(receiverDeployment$Region),]
# 
# pal <- c(rev(magma(15)),viridis(6))
# 
# 
# genLoc <- unique(receiverDeployment$GeneralLocation)
# nn <- length(unique(receiverDeployment$GeneralLocation))
# 
# layout(matrix(1:2),heights=c(0.25,0.8))
# par(mar=c(9.1,6.1,1.1,1.1))
# image(x=1:21,z=matrix(1:21), col=pal, xaxt="n", yaxt="n", xlab="", ylab="")
# axis(1, at=1:21, labels=levels(receiverDeployment$Region), las=2)
# 
# plot(NA,NA,xlim=range(c(receiverDeployment$Start,receiverDeployment$Stop), na.rm=TRUE),ylim=c(1,nn), xaxt="n")
# axis(1, at=as.Date(c("2009-10-01","2010-10-01","2011-10-01","2012-10-01","2013-10-01","2014-10-01",
#                      "2015-10-01","2016-10-01","2017-10-01")), 
#      labels=c("2009-10-01","2010-10-01","2011-10-01","2012-10-01","2013-10-01","2014-10-01",
#               "2015-10-01","2016-10-01","2017-10-01"), las=2)
# 
# for(ii in 1:nn){
#   tmp <- receiverDeployment[receiverDeployment$GeneralLocation==genLoc[ii],]
#   for(jj in 1:nrow(tmp)){
#     lines(c(tmp$Start[jj],tmp$Stop[jj]), rep(ii,2), lwd=1.5, col=pal[tmp$Region[jj]])
#   }
# }



## Better--proportions of fish detected in a region by 5-day period -------------------------------

## Limit to water years with good presence of receivers in SF Bay and San Pablo Bay
detects <- detects[detects$DetectDate >= as.Date("2009-10-01") & detects$DetectDate < as.Date("2014-10-01"),]

detects_collapsed <- unique(data.frame(FishID=detects$FishID, 
                                       DetectDate=detects$DetectDate, 
                                       Region=detects$Region))
detects_collapsed$DetectDoy <- yday(detects_collapsed$DetectDate)
detects_collapsed$DetectDoy[detects_collapsed$DetectDoy==366] <- 365

dstart <- seq(from=1, to=365, by=5)
dend <- dstart + 4

nfish_sfbay <- rep(NA, length(dstart)-1)
nfish_sanpablo <- rep(NA, length(dstart)-1)
nfish_all <- rep(NA, length(dstart)-1)

for(ii in 1:(length(dstart)-1)){
  
  tmp1 <- detects_collapsed[detects_collapsed$DetectDoy >= dstart[ii] &
                              detects_collapsed$DetectDoy <= dend[ii],]
  nfish_all[ii] <- length(unique(tmp1$FishID))
  nfish_sfbay[ii] <- length(unique(tmp1$FishID[tmp1$Region=="SF_Bay"]))
  nfish_sanpablo[ii] <- length(unique(tmp1$FishID[tmp1$Region=="San_Pablo_Bay"]))
  
}



par(mfrow=c(2,1), mar=c(4.1,4.1,2.1,1.1), tcl=-0.3, mgp=c(2.2,0.8,0))

barplot(nfish_sfbay/nfish_all, width=5,space=0, ylab="Proportion", xlab="Day of year")
axis(1, at=seq(0,365,by=30))
abline(v=c(214,244), lwd=1.5, col="red", lty=2)
mtext("San Francisco Bay", line=0.5)
text(par("usr")[1]+0.05*diff(par("usr")[1:2]),
     par("usr")[4]-0.05*diff(par("usr")[3:4]), "a)")

barplot(nfish_sanpablo/nfish_all, width=5,space=0, ylab="Proportion", xlab="Day of year")
axis(1, at=seq(0,365,by=30))
abline(v=c(214,244), lwd=1.5, col="red", lty=2)
mtext("San Pablo Bay", line=0.5)
text(par("usr")[1]+0.05*diff(par("usr")[1:2]),
     par("usr")[4]-0.05*diff(par("usr")[3:4]), "b)")



## Now break down by size/age class ---------------------------------------------------------------

#convert TL to FL
detects$Length[detects$Length_Type=="TL"] <- -3.55 + 0.93*detects$Length[detects$Length_Type=="TL"] 

age_from_length <- function(ll){
  return((-log(1-(ll/380))/0.027) - 2.36)
}


detects$tag_age <- round(age_from_length(detects$Length/10))
detects$est_age <- floor(detects$tag_age + as.numeric(detects$DetectDate-detects$tagging_date)/365)
detects$age_group <- NA
detects$age_group[detects$est_age<=10] <- "Juvenile"
detects$age_group[detects$est_age >10 & detects$est_age <=15] <- "Transitional"
detects$age_group[detects$est_age>15] <- "Reproductive"

# for below slot fish

detects_blw <- detects[detects$age_group=="Juvenile",]
detects_collapsed_blw <- unique(data.frame(FishID=detects_blw$FishID, 
                                       DetectDate=detects_blw$DetectDate, 
                                       Region=detects_blw$Region))
detects_collapsed_blw$DetectDoy <- yday(detects_collapsed_blw$DetectDate)
detects_collapsed_blw$DetectDoy[detects_collapsed_blw$DetectDoy==366] <- 365

nfish_sfbay_blw <- rep(NA, length(dstart)-1)
nfish_sanpablo_blw <- rep(NA, length(dstart)-1)
nfish_all_blw <- rep(NA, length(dstart)-1)

for(ii in 1:(length(dstart)-1)){
  
  tmp1 <- detects_collapsed_blw[detects_collapsed_blw$DetectDoy >= dstart[ii] &
                              detects_collapsed_blw$DetectDoy <= dend[ii],]
  nfish_all_blw[ii] <- length(unique(tmp1$FishID))
  nfish_sfbay_blw[ii] <- length(unique(tmp1$FishID[tmp1$Region=="SF_Bay"]))
  nfish_sanpablo_blw[ii] <- length(unique(tmp1$FishID[tmp1$Region=="San_Pablo_Bay"]))
  
}

# for within slot fish

detects_slt <- detects[detects$age_group=="Transitional",]
detects_collapsed_slt <- unique(data.frame(FishID=detects_slt$FishID, 
                                           DetectDate=detects_slt$DetectDate, 
                                           Region=detects_slt$Region))
detects_collapsed_slt$DetectDoy <- yday(detects_collapsed_slt$DetectDate)
detects_collapsed_slt$DetectDoy[detects_collapsed_slt$DetectDoy==366] <- 365

nfish_sfbay_slt <- rep(NA, length(dstart)-1)
nfish_sanpablo_slt <- rep(NA, length(dstart)-1)
nfish_all_slt <- rep(NA, length(dstart)-1)

for(ii in 1:(length(dstart)-1)){
  
  tmp1 <- detects_collapsed_slt[detects_collapsed_slt$DetectDoy >= dstart[ii] &
                                  detects_collapsed_slt$DetectDoy <= dend[ii],]
  nfish_all_slt[ii] <- length(unique(tmp1$FishID))
  nfish_sfbay_slt[ii] <- length(unique(tmp1$FishID[tmp1$Region=="SF_Bay"]))
  nfish_sanpablo_slt[ii] <- length(unique(tmp1$FishID[tmp1$Region=="San_Pablo_Bay"]))
  
}

#for above slot fish

detects_abv <- detects[detects$age_group=="Reproductive",]
detects_collapsed_abv <- unique(data.frame(FishID=detects_abv$FishID, 
                                           DetectDate=detects_abv$DetectDate, 
                                           Region=detects_abv$Region))
detects_collapsed_abv$DetectDoy <- yday(detects_collapsed_abv$DetectDate)
detects_collapsed_abv$DetectDoy[detects_collapsed_abv$DetectDoy==366] <- 365

nfish_sfbay_abv <- rep(NA, length(dstart)-1)
nfish_sanpablo_abv <- rep(NA, length(dstart)-1)
nfish_all_abv <- rep(NA, length(dstart)-1)

for(ii in 1:(length(dstart)-1)){
  
  tmp1 <- detects_collapsed_abv[detects_collapsed_abv$DetectDoy >= dstart[ii] &
                                  detects_collapsed_abv$DetectDoy <= dend[ii],]
  nfish_all_abv[ii] <- length(unique(tmp1$FishID))
  nfish_sfbay_abv[ii] <- length(unique(tmp1$FishID[tmp1$Region=="SF_Bay"]))
  nfish_sanpablo_abv[ii] <- length(unique(tmp1$FishID[tmp1$Region=="San_Pablo_Bay"]))
  
}

pal=brewer.pal(3, "Set2")

pdf("fig_hab_exposure_bysize_20250317.pdf", width=6.5, height=6.5)

par(mfcol=c(2,3), mar=c(3.6,3.6,2.1,0.5), tcl=-0.3, mgp=c(2.2,0.8,0))

barplot(nfish_sfbay_blw/nfish_all_blw, width=5,space=0, ylab="Proportion", xlab="Day of year",
        col=pal[1], ylim=c(0,0.225))
axis(1, at=seq(0,365,by=30))
abline(v=c(214,244), lwd=1.5, col="darkgrey", lty=2)
mtext("San Francisco Bay", line=0.5, cex=3/4)
text(par("usr")[1]+0.05*diff(par("usr")[1:2]),
     par("usr")[4]-0.05*diff(par("usr")[3:4]), "a)")
legend("topleft", fill=pal, legend=c("Juvenile","Transitional","Reproductive"), bty="n", inset=c(0,0.08), cex=1.25)

barplot(nfish_sanpablo_blw/nfish_all_blw, width=5,space=0, ylab="Proportion", xlab="Day of year", 
        col=pal[1], ylim=c(0,0.55))
axis(1, at=seq(0,365,by=30))
abline(v=c(214,244), lwd=1.5, col="darkgrey", lty=2)
mtext("San Pablo Bay", line=0.5, cex=3/4)
text(par("usr")[1]+0.05*diff(par("usr")[1:2]),
     par("usr")[4]-0.05*diff(par("usr")[3:4]), "b)")

barplot(nfish_sfbay_slt/nfish_all_slt, width=5,space=0, ylab="Proportion", xlab="Day of year", 
        col=pal[2], ylim=c(0,0.225))
axis(1, at=seq(0,365,by=30))
abline(v=c(214,244), lwd=1.5, col="darkgrey", lty=2)
mtext("San Francisco Bay", line=0.5, cex=3/4)
text(par("usr")[1]+0.05*diff(par("usr")[1:2]),
     par("usr")[4]-0.05*diff(par("usr")[3:4]), "c)")

barplot(nfish_sanpablo_slt/nfish_all_slt, width=5,space=0, ylab="Proportion", xlab="Day of year", 
        col=pal[2], ylim=c(0,0.55))
axis(1, at=seq(0,365,by=30))
abline(v=c(214,244), lwd=1.5, col="darkgrey", lty=2)
mtext("San Pablo Bay", line=0.5, cex=3/4)
text(par("usr")[1]+0.05*diff(par("usr")[1:2]),
     par("usr")[4]-0.05*diff(par("usr")[3:4]), "d)")

barplot(nfish_sfbay_abv/nfish_all_abv, width=5,space=0, ylab="Proportion", xlab="Day of year", 
        col=pal[3], ylim=c(0,0.225))
axis(1, at=seq(0,365,by=30))
abline(v=c(214,244), lwd=1.5, col="darkgrey", lty=2)
mtext("San Francisco Bay", line=0.5, cex=3/4)
text(par("usr")[1]+0.05*diff(par("usr")[1:2]),
     par("usr")[4]-0.05*diff(par("usr")[3:4]), "e)")

barplot(nfish_sanpablo_abv/nfish_all_abv, width=5,space=0, ylab="Proportion", xlab="Day of year", 
        col=pal[3], ylim=c(0,0.55))
axis(1, at=seq(0,365,by=30))
abline(v=c(214,244), lwd=1.5, col="darkgrey", lty=2)
mtext("San Pablo Bay", line=0.5, cex=3/4)
text(par("usr")[1]+0.05*diff(par("usr")[1:2]),
     par("usr")[4]-0.05*diff(par("usr")[3:4]), "f)")

dev.off()