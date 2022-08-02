
rm(list=ls())

library(lubridate)
library(dplyr)
library(RColorBrewer)
library(rgdal)
library(cluster)
library(ecodist)
library(scales)
library(car)
library(vegan)
library(abdiv)
library(ncf)
library(raster)

## ------------------------------------------------------------------------------------------------
## load and format data ---------------------------------------------------------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
dat <- read.csv("WSdetects_filtered_wy.csv")
generalAreaCoords <- read.csv("generalArea_coords.csv")

goodLocs <- generalAreaCoords$GeneralArea
dat$DetectDate <- as.Date(dat$DetectDate)

## try dropping locations in Mok and SJ rivers
dropLocs <- goodLocs[grepl("Mok_", goodLocs) | grepl("SJ_", goodLocs)]
dropLocs <- c(dropLocs, "PotatoSlough", "FranksTractE")

dat$GeneralArea[dat$GeneralArea =="Steam_Marina"] <- "SR_BlwSteam" #merge two very close sites

goodLocs <- goodLocs[!goodLocs %in% c(dropLocs, "Steam_Marina")]
generalAreaCoords <- generalAreaCoords[!generalAreaCoords$GeneralArea %in% c(dropLocs, "Steam_Marina"),]
dat <- dat[!dat$GeneralArea %in% dropLocs,]
generalAreaCoords$displayName <- c("Golden Gate", "Richmond Bridge", "Carquinez Bridge",
                                   "Benicia Bridge", "Decker Island", "Three Mile Slough",
                                   "Rio Vista", "River Mouth", "Steamboat Sl. Down",
                                   "Georgiana Slough", "Delta Cross Channel", "Steamboat Sl. Up",
                                   "Sutter Slough", "Freeport", "Up River")

generalAreaCoords$Lat[generalAreaCoords$GeneralArea=="SR_SteamboatSl"] <- 38.18317
generalAreaCoords$Lon[generalAreaCoords$GeneralArea=="SR_SteamboatSl"] <- -121.6579

generalAreaCoords$Lat[generalAreaCoords$GeneralArea=="UpRiver"] <- 38.78543
generalAreaCoords$Lon[generalAreaCoords$GeneralArea=="SR_SteamboatSl"] <- -121.6221



fishIDs <- unique(dat$FishID)
fishWYIDs <- unique(dat$FishWYID)

dat <- dat[!is.na(dat$Length),]


wyclass <- data.frame(wateryear=as.character(2010:2017),
                      wyclass=factor(c("BN","W","BN","D","C","C","BN","W"))
)


## ------------------------------------------------------------------------------------------------
## Create occupancy matrices ----------------------------------------------------------------------

occupancy_matrix<-matrix(0, nrow=length(fishWYIDs), ncol=length(goodLocs))
rownames(occupancy_matrix)<-fishWYIDs
colnames(occupancy_matrix)<-goodLocs

for(ii in 1:length(fishWYIDs)){
  
  dat.ii<-dat[dat$FishWYID==fishWYIDs[ii],]
  #dat.ii<-dat.ii[order(dat.ii$DetectDate),]
  
  counter=1
  for(tt in 1:nrow(dat.ii)){
    dat.ii.tt<-dat.ii[dat.ii$DetectDate==dat.ii$DetectDate[counter],]
    occupancy_matrix[ii,goodLocs %in% dat.ii.tt$GeneralArea] <- 
      occupancy_matrix[ii,goodLocs %in% dat.ii.tt$GeneralArea] +
      rep(1/length(unique(dat.ii.tt$GeneralArea)), times=length(unique(dat.ii.tt$GeneralArea)))
    counter=counter+nrow(dat.ii.tt)
    if(counter>=nrow(dat.ii)){break}
  }
  
}

## drop fish-water years with too few detections
drop <- which(rowSums(occupancy_matrix)<30)

#sum(rowSums(occupancy_matrix)<30)

fishWYIDs <- fishWYIDs[-drop]
occupancy_matrix <- occupancy_matrix[-drop,]


## normalize to a proportion
occupancy_matrix_norm<-occupancy_matrix
for(ii in 1:nrow(occupancy_matrix)){
  occupancy_matrix_norm[ii,]<-occupancy_matrix[ii,]/sum(occupancy_matrix[ii,])
}


## ------------------------------------------------------------------------------------------------
## summarize information about fish and water years -----------------------------------------------

tagging.date<-rep(NA, length(fishWYIDs))
release.location<-rep(NA, length(fishWYIDs))
sex<-rep(NA, length(fishWYIDs))
fishlength<-rep(NA, length(fishWYIDs))
length.type<-rep(NA, length(fishWYIDs))
detect.length<-rep(NA, length(fishWYIDs))
life.stage<-rep(NA, length(fishWYIDs))
wateryear<-rep(NA, length(fishWYIDs))
studyID<-rep(NA, length(fishWYIDs))

for(ii in 1:length(fishWYIDs)){
  
  tmp.ii<-dat[dat$FishWYID == fishWYIDs[ii],]
  
  tagging.date[ii]<-tmp.ii$tagging_date[1]
  release.location[ii]<-tmp.ii$Release_Location[1]
  sex[ii]<-tmp.ii$Sex[1]
  fishlength[ii]<-tmp.ii$Length[1]
  length.type[ii]<-tmp.ii$Length_Type[1]
  detect.length[ii]<-round(max(tmp.ii$DetectDate,na.rm=T)-min(tmp.ii$DetectDate,na.rm=T))
  life.stage[ii]<-tmp.ii$Life_Stage[1]
  wateryear[ii]<-tmp.ii$waterYear[1]
  studyID[ii]<-tmp.ii$StudyID
}

fish.summary<-cbind(fishWYIDs
                    ,tagging.date
                    ,release.location
                    ,sex
                    ,fishlength
                    ,length.type
                    ,life.stage
                    ,detect.length
                    ,wateryear
                    ,studyID
)


fish.summary<-data.frame(fish.summary)
fish.summary$sex[fish.summary$sex=="f"]<-"F"
fish.summary$release.location[fish.summary$release.location=="Freemont Weir"]<-"Fremont Weir"


fish.summary <- left_join(fish.summary, wyclass)

rm(tagging.date, release.location, sex, fishlength, length.type, life.stage, detect.length, wateryear)

fish.summary$fishlength <- as.numeric(fish.summary$fishlength)
fish.summary$tagging.date <- as.POSIXct(fish.summary$tagging.date)
fish.summary$tagging.year <- year(fish.summary$tagging.date)


fish.summary$fishlength[fish.summary$length.type=="TL"] <- -3.55 + 0.93*fish.summary$fishlength[fish.summary$length.type=="TL"] 

age_from_length <- function(ll){
  return((-log(1-(ll/380))/0.027) - 2.36)
}

fish.summary$tag_age <- floor(age_from_length(fish.summary$fishlength/10))
fish.summary$age_est <- fish.summary$tag_age + (as.numeric(fish.summary$wateryear) - fish.summary$tagging.year)
fish.summary$age_group <- NA
fish.summary$age_group[fish.summary$age_est<=10] <- "Below Slot"
fish.summary$age_group[fish.summary$age_est>10 & fish.summary$age_est <= 15] <- "Fishery"
fish.summary$age_group[fish.summary$age_est>15] <- "Above Slot"
#fish.summary$age_group[fish.summary$age_est>20] <- "Above Slot"
fish.summary$age_group <- as.factor(fish.summary$age_group)
fish.summary$fishID <- substr(fish.summary$fishWYIDs,1,8)

table(fish.summary$fishID)
sum(table(fish.summary$fishID)==1)
sum(table(fish.summary$fishID)>1)
nrow(fish.summary)

fish.summary$n_wy <- NA
for(ii in 1:nrow(fish.summary)){
  fish.summary$n_wy[ii] <- sum(fish.summary$fishID == fish.summary$fishID[ii])
}

# boxplot(fish.summary$fishlength ~ fish.summary$n_wy)
# abline(h=1524, lty=2, col="blue")
# abline(h=1016, lty=2, col="blue")

fish.summary$release.region <- NA
fish.summary$release.region[fish.summary$release.location %in% c("Tisdale Bypass","DWR_YB_Fyke","Fremont Weir")] <- "SR" #sacramento river
fish.summary$release.region[fish.summary$release.location %in% c("Grizzly_Bay","Suisun_Bay","San_Pablo_Bay")] <- "DELTA"


# # subset data
# dat <- dat[dat$FishID %in% fish.summary$fishID[fish.summary$n_wy >1],]
# occupancy_matrix <- occupancy_matrix[fish.summary$n_wy > 1,]
# occupancy_matrix_norm <- occupancy_matrix_norm[fish.summary$n_wy > 1,]
# fish.summary <- fish.summary[fish.summary$n_wy > 1,]


## ------------------------------------------------------------------------------------------------
## Begin dissimilarity analyses -------------------------------------------------------------------

# make distance matrix
occupancy_dist_matrix<-matrix(NA, nrow=nrow(occupancy_matrix), ncol=nrow(occupancy_matrix))
for(ii in 2:nrow(occupancy_dist_matrix)){
  for(jj in 1:(ii-1)){
    occupancy_dist_matrix[ii,jj]<-hellinger(occupancy_matrix[ii,], occupancy_matrix[jj,])
  }
}
occupancy_dist<-as.dist(occupancy_dist_matrix)


## Test for significance of divergence in behavior ------------------------------------------------

occ_medians <- apply(occupancy_matrix_norm, 2, median)
probvect <- cumsum(occ_medians)
ndetects <- rowSums(occupancy_matrix)
nreps <- 1000

occmat_surrogs <- list()

for(nn in 1:nreps){
  
  surrog_occup <- matrix(0, length(fishWYIDs), length(goodLocs))
  
  for(ii in 1:length(fishWYIDs)){
    for(jj in 1:ndetects[ii]){
      p <- runif(1,0,1)
      if(p < probvect[1]){surrog_occup[ii,1] <- surrog_occup[ii,1] + 1}
      if(p >= probvect[1] & p < probvect[2]){surrog_occup[ii,2] <- surrog_occup[ii,2] + 1}
      if(p >= probvect[2] & p < probvect[3]){surrog_occup[ii,3] <- surrog_occup[ii,3] + 1}
      if(p >= probvect[3] & p < probvect[4]){surrog_occup[ii,4] <- surrog_occup[ii,4] + 1}
      if(p >= probvect[4] & p < probvect[5]){surrog_occup[ii,5] <- surrog_occup[ii,5] + 1}
      if(p >= probvect[5] & p < probvect[6]){surrog_occup[ii,6] <- surrog_occup[ii,6] + 1}
      if(p >= probvect[6] & p < probvect[7]){surrog_occup[ii,7] <- surrog_occup[ii,7] + 1}
      if(p >= probvect[7] & p < probvect[8]){surrog_occup[ii,8] <- surrog_occup[ii,8] + 1}
      if(p >= probvect[8] & p < probvect[9]){surrog_occup[ii,9] <- surrog_occup[ii,9] + 1}
      if(p >= probvect[9] & p < probvect[10]){surrog_occup[ii,10] <- surrog_occup[ii,10] + 1}
      if(p >= probvect[10] & p < probvect[11]){surrog_occup[ii,11] <- surrog_occup[ii,11] + 1}
      if(p >= probvect[11] & p < probvect[12]){surrog_occup[ii,12] <- surrog_occup[ii,12] + 1}
      if(p >= probvect[12] & p < probvect[13]){surrog_occup[ii,13] <- surrog_occup[ii,13] + 1}
      if(p >= probvect[13] & p < probvect[14]){surrog_occup[ii,14] <- surrog_occup[ii,14] + 1}
      if(p >= probvect[14] & p < probvect[15]){surrog_occup[ii,15] <- surrog_occup[ii,15] + 1}
    }
  }
  occmat_surrogs[[nn]] <- surrog_occup
  
}

sumDiffs_surr <- rep(NA, length(occmat_surrogs))

for(nn in 1:length(occmat_surrogs)){
  
  occupancy_matrix.nn <- occmat_surrogs[[nn]]
  
  occupancy_dist.nn<-matrix(NA, nrow=nrow(occupancy_matrix.nn), ncol=nrow(occupancy_matrix.nn))
  for(ii in 2:nrow(occupancy_dist.nn)){
    for(jj in 1:(ii-1)){
      occupancy_dist.nn[ii,jj]<-hellinger(occupancy_matrix.nn[ii,], occupancy_matrix.nn[jj,])
    }
  }
  occupancy_dist.nn<-as.dist(occupancy_dist.nn)
  
  sumDiffs_surr[nn] <- sum(occupancy_dist.nn, na.rm=T)
  
}

summary(sumDiffs_surr)
rank(c(sum(occupancy_dist), sumDiffs_surr))[1]/1001 #largest, by a wide margin


# Is there more dispersion among fish or among water years? ---------------------------------------
fishWYIDs <- fish.summary$fishWYIDs

same_fish <- matrix(NA, length(fishWYIDs), length(fishWYIDs))
same_wyear <- matrix(NA, length(fishWYIDs), length(fishWYIDs))

for(ii in 2:length(fishWYIDs)){
  for(jj in 1:(ii-1)){
    same_fish[ii,jj] <- ifelse(substr(fishWYIDs[ii],1,8)==substr(fishWYIDs[jj],1,8), 1, NA)
    same_wyear[ii,jj] <- ifelse(substr(fishWYIDs[ii],10,13)==substr(fishWYIDs[jj],10,13), 1, NA)
  }
}


png("~/GitHub/fishsync/WhiteSturgeonOccupancy/fig_distance_btwfish_btwyr.png", units="in", 
    width=6.5, height=6.5, res=150)
par(mfrow=c(2,1), mar=c(4.1,4.1,2.1,1.1))
hist(occupancy_dist*same_fish[lower.tri(same_fish)], xlim=c(0,1.5), 
     xlab="Hellinger distance", main="same fish, different year")
hist(occupancy_dist*same_wyear[lower.tri(same_wyear)], xlim=c(0,1.5),
     xlab="Hellinger distance", main="same year, different fish")
dev.off()

#test significance using resampling

diff.fishwy <- abs(median(occupancy_dist*same_fish[lower.tri(same_fish)],na.rm=T) - 
                     median(occupancy_dist*same_wyear[lower.tri(same_wyear)], na.rm=T))

nrand = 1000
diff.fishwy.surr <- rep(NA, nrand)

for(ii in 1:nrand){
  
  shuffle <- sample(1:length(fishWYIDs), length(fishWYIDs), replace=TRUE)
  dmat_surr <- occupancy_dist_matrix[shuffle,shuffle]
  diff.fishwy.surr[ii] <- abs(median(dmat_surr*same_fish,na.rm=T) - 
                                median(dmat_surr*same_wyear, na.rm=T))
}

summary(diff.fishwy.surr)
rank(c(diff.fishwy, diff.fishwy.surr))[1]/1001 #largest, by a wide margin



## PERMANOVA statistical model --------------------------------------------------------------------
fish.summary$tagging.year <- factor(fish.summary$tagging.year)

table(fish.summary$release.region, fish.summary$tagging.year)
chisq.test(table(fish.summary$release.region, fish.summary$tagging.year))

table(fish.summary$tagging.year, fish.summary$sex)
chisq.test(table(fish.summary$tagging.year, fish.summary$sex))

table(fish.summary$release.region, fish.summary$sex)
chisq.test(table(fish.summary$release.region, fish.summary$sex))

occupancy_dist_sexed <- as.dist(occupancy_dist_matrix[fish.summary$sex !="U", fish.summary$sex !="U"])
fish.summary.sexed <- fish.summary[fish.summary$sex != "U",]


mod1 <- adonis2(occupancy_dist ~ release.region + tagging.year + wyclass, data=fish.summary, by = "margin", permutations = 9999)
print(mod1)

mod2 <- adonis2(occupancy_dist ~ release.region + tagging.year + sex, data=fish.summary, by = "margin", permutations = 9999)
print(mod2)

mod2b <- adonis2(occupancy_dist_sexed ~ release.region + tagging.year + sex, data=fish.summary.sexed, by = "margin", permutations = 9999)
print(mod2b)

mod3 <- adonis2(occupancy_dist ~ release.region + tagging.year + age_group, data=fish.summary, by = "margin", permutations = 9999)
print(mod3)

## Investigate behavior of important kinds of fish based on significant terms in mod
keytypes <- expand.grid(sex=unique(fish.summary$sex), 
                        age_group=unique(fish.summary$age_group))
keytypes$sex <- as.character(keytypes$sex)
#keytypes$sex <- as.character(keytypes$sex)
keytypes$age_group <- as.character(keytypes$age_group)


typemean <- matrix(NA, nrow=nrow(keytypes), ncol=length(goodLocs))
typeN <- rep(NA, nrow(keytypes))

for(ii in 1:nrow(keytypes)){
  
  typeN[ii] <- sum(fish.summary$sex == keytypes$sex[ii]
                   & fish.summary$age_group == keytypes$age_group[ii])
  if(typeN[ii] > 0){
    if(typeN[ii] == 1){
      typemean[ii,] <- occupancy_matrix_norm[fish.summary$sex == keytypes$sex[ii]
                                             & fish.summary$age_group == keytypes$age_group[ii], ]
    }
    else{
      typemean[ii,] <- apply(occupancy_matrix_norm[fish.summary$sex == keytypes$sex[ii]
                                                   & fish.summary$age_group == keytypes$age_group[ii], ]
                             , MARGIN=2, FUN="mean")
    }
  }
}



pdf("~/GitHub/fishsync/WhiteSturgeonOccupancy/ws_occupancy_barplot_wy_bygroup_allfish.pdf", 
    onefile=TRUE)
par(mar=c(8.5,3.5,3.5,0.5), mgp=c(2.2,0.5,0), tcl=-0.3)

for(ii in 1:nrow(keytypes)){
  
  if(!all(is.na(typemean[ii,]))){
    barplot(typemean[ii,match(generalAreaCoords$GeneralArea, goodLocs)],
            names.arg=generalAreaCoords$displayName, las=2,
            ylab="Proportional occupancy", args.legend=list(bty="n"))
    mtext(paste(as.character(keytypes[ii,]), collapse=" "))
    mtext(paste0("n = ", typeN[ii]), at=0)
    # arrows(bp,
    #        occupancy_group_q75[,match(generalAreaCoords$GeneralArea, goodLocs)],
    #        bp,
    #        occupancy_group_q25[,match(generalAreaCoords$GeneralArea, goodLocs)],
    #        angle=90, code=3, length=0.02)
    # lines(predict(sm1b, x=newx), col=pal2[1], lwd=2, lty=2)
    # lines(predict(sm2b, x=newx), col=pal2[2], lwd=2, lty=2)
    # lines(predict(sm3b, x=newx), col=pal2[3], lwd=2, lty=2)
    #axis(3, at=bp[2,], labels=signif(generalAreaCoords$RiverKm,3), las=2)
    #mtext("River Km", 3, line=2.4)
  }
  
}

dev.off()


## make some maps

aboveSlot_avg <- colMeans(occupancy_matrix_norm[fish.summary$age_group=="Above Slot",])
fishery_avg <- colMeans(occupancy_matrix_norm[fish.summary$age_group=="Fishery",])
belowSlot_avg <- colMeans(occupancy_matrix_norm[fish.summary$age_group=="Below Slot",])


states <- readOGR("statesp020.shp")
states <- states[states$STATE=="California",]
bkgd <- readOGR("backrgound_layer.shp")
water <- readOGR("/Volumes/GoogleDrive-108173464723009405822/My Drive/Synchrony of Movement/Data/GIS/CA_Hydro/CA_hydrology_Polygon.shp")
water <- spTransform(water, CRS(proj4string(states)))


range(generalAreaCoords$Lon)
range(generalAreaCoords$Lat)

ptscale <- function(x, mm, mx, xx=NULL){
  xs <- x - min(x)
  if(is.null(xx)){xs <- xs/max(x)}
  else{xs <- xs/xx}
  xs <- xs*(mx-mm) + mm
  return(xs)
}


png("map_occ_by_size.png", width=6.5, height=3.25, 
    res=300, units="in")

par(mfrow=c(1,3), mar=c(0.5,0.5,1.5,0.5))

sel <- belowSlot_avg>0.01
txtpos <- rep(3,15)
txtpos[5] <- 2
txtpos[6] <- 4
txtpos[7] <- 2
txtpos[8] <- 4
txtpos[10] <- 1
txtpos[11] <- 4
txtpos[12] <- 4

plot(states, col="grey70", border=NA, xlim=c(-122.7,-121.4), ylim=c(37.7,39))
plot(bkgd, col="grey70", border=NA, add=TRUE)
plot(water, col="skyblue3", border="skyblue3", add=TRUE)
points(generalAreaCoords$Lon[sel], generalAreaCoords$Lat[sel], cex=ptscale(belowSlot_avg[sel],1,3,0.5),
       pch=21, col="white", bg="black")
points(generalAreaCoords$Lon[!sel], generalAreaCoords$Lat[!sel], col="grey50")
mtext("Below slot")
legend("topleft",pch=21,pt.cex=c(1,1,2,3),col=c("grey50","white","white","white"),
       pt.bg=c(NA,"black","black","black"),legend=c("<0.01"," 0.01"," 0.25"," 0.50"),
       inset=0.01, bg=NA, title="Proportional\noccupancy", bty="n")
text(x=generalAreaCoords$Lon, generalAreaCoords$Lat, labels=as.character(1:15), pos=txtpos)
text(x= -122.3, y=37.7, "SF Bay", col="darkblue",cex=0.8)
text(x=-122.5, y=38.1, "San Pablo Bay", col="darkblue",cex=0.8)
text(x=-122.15, y=38.25, "Grizzly Bay", col="darkblue",cex=0.8)
segments(x0=-122.15, y0=38.225, x1=-122.05, y1=38.13, col="darkblue", lwd=0.5)
text(x=-122, y=37.95, "Suisun Bay", col="darkblue",cex=0.8)
segments(x0=-122, y0=37.97, x1=-121.95, y1=38.07, col="darkblue", lwd=0.5)
text(x=-121.5, y=38.05, "Delta", col="darkblue", cex=0.8)
text(x=-121.8, y=38.8, "Sacramento River", col="darkblue", cex=0.8, srt=-57)

sel <- fishery_avg>0.01
plot(states, col="grey70", border=NA, xlim=c(-122.7,-121.4), ylim=c(37.7,39))
plot(bkgd, col="grey70", border=NA, add=TRUE)
plot(water, col="skyblue3", border="skyblue3", add=TRUE)
points(generalAreaCoords$Lon[sel], generalAreaCoords$Lat[sel], cex=ptscale(fishery_avg[sel],1,3,0.5),
       pch=21, col="white", bg="black")
points(generalAreaCoords$Lon[!sel], generalAreaCoords$Lat[!sel], col="grey50")
mtext("Fishery")
text(x=generalAreaCoords$Lon, generalAreaCoords$Lat, labels=as.character(1:15), pos=txtpos)
#legend("topleft", pch=as.character(1:7), legend=generalAreaCoords$displayName[1:7], 
#       inset=0.01, bg=NA, trace=TRUE)

text(-122.73, 38.8, pos=4, labels="1 Golden Gate\n2 Richmond Bridge\n3 Carquinez Bridge\n4 Benicia Bridge\n5 Decker Island\n6 Three Mile Slough\n7 Rio Vista")

sel <-aboveSlot_avg>0.01
plot(states, col="grey70", border=NA, xlim=c(-122.7,-121.4), ylim=c(37.7,39))
plot(bkgd, col="grey70", border=NA, add=TRUE)
plot(water, col="skyblue3", border="skyblue3", add=TRUE)
points(generalAreaCoords$Lon[sel], generalAreaCoords$Lat[sel], cex=ptscale(aboveSlot_avg[sel],1,3,0.5),
       pch=21, col="white", bg="black")
points(generalAreaCoords$Lon[!sel], generalAreaCoords$Lat[!sel], col="grey50")
mtext("Above slot")
text(x=generalAreaCoords$Lon, generalAreaCoords$Lat, labels=as.character(1:15), pos=txtpos)
#legend("topleft", pch=as.character(8:15), legend=generalAreaCoords$displayName[8:15], 
#       inset=0.01, bg=NA)
text(-122.73, 38.75, pos=4, labels=" 8 River Mouth\n 9 Steamboat Sl. Down\n10 Georgiana Slough\n11 Delta Cross Channel\n12 Steamboat Sl. Up\n13 Sutter Slough\n14 Freeport \n15 Up River")


dev.off()




# Water diversions --------------------------------------------------------------------------------


divert_locs <- read.csv("/Users/jonathanwalter/Documents/Research/DATA/CaliWaterDiversions/diversion_locations.csv")
divert_data <- read.csv("/Users/jonathanwalter/Documents/Research/DATA/CaliWaterDiversions/diversion_attributes.csv")

divert_locs <- divert_locs[(grepl("DELTA", divert_locs$SOURCE_NAME)
                            & !grepl("MENDOTA", divert_locs$SOURCE_NAME)) |
                             grepl("GRIZZLY BAY", divert_locs$SOURCE_NAME) |
                             (grepl("SACRAMENTO", divert_locs$SOURCE_NAME) &
                                !grepl("SPRING", divert_locs$SOURCE_NAME)) |
                             #grepl("YOLO BYPASS", divert_locs$SOURCE_NAME) |
                             grepl("STEAMBOAT SLOUGH", divert_locs$SOURCE_NAME) |
                             grepl("POTATO SLOUGH", divert_locs$SOURCE_NAME) |
                             grepl("SUTTER SLOUGH", divert_locs$SOURCE_NAME),]

#divert_locs <- divert_locs[divert_locs$LATITUDE > 37.0,]

divert_locs <- divert_locs[!duplicated(divert_locs),]

divert_locs <- divert_locs[,colnames(divert_locs) %in% c("POD_ID","LATITUDE","LONGITUDE","SOURCE_NAME")]

divert_data <- divert_data[,colnames(divert_data) %in% c("CORE_POD_ID", "WR_STATUS", "FACE_VALUE_AMOUNT")]

diversions <- inner_join(divert_locs, divert_data, by=c("POD_ID" = "CORE_POD_ID"))
diversions <- diversions[complete.cases(diversions),]
diversions <- diversions[!diversions$WR_STATUS %in% c("Cancelled","Revoked"),]


plot(diversions$LONGITUDE, diversions$LATITUDE)

dmat <- gcdist(generalAreaCoords$Lon, generalAreaCoords$Lat)

summary(dmat[lower.tri(dmat)])



#assign a diversion to the nearest location <---------------------------------

diversion_to_sensorLoc <- rep(NA, nrow(diversions))
maxdist <- 2

for(ii in 1:nrow(diversions)){
  
  locdist<-rep(NA, nrow(generalAreaCoords))
  for(jj  in 1:nrow(generalAreaCoords)){
    
    #locdist[jj] <- sqrt((diversions$LONGITUDE[ii] - sensorLocations$Lon[jj])^2 + (diversions$LATITUDE[ii] - sensorLocations$Lat[jj])^2)
    locdist[jj] <- gcdist(x = c(diversions$LONGITUDE[ii], generalAreaCoords$Lon[jj]), 
                          y=c(diversions$LATITUDE[ii], generalAreaCoords$Lat[jj]))[2,1]
    
  }
  if(min(locdist) <= maxdist){
    diversion_to_sensorLoc[ii] <- generalAreaCoords$GeneralArea[which.min(locdist)]
  }
  
}

n_diversions<-as.data.frame(table(diversion_to_sensorLoc))

volume_to_sensorLoc <- data.frame(sensorLoc=diversion_to_sensorLoc,
                                  volume=diversions$FACE_VALUE_AMOUNT)

sum_Volume <- aggregate(volume ~ sensorLoc, data=volume_to_sensorLoc, FUN="sum")


#apportion to groups based on proportional occupancy

sensorLocations <- left_join(generalAreaCoords, n_diversions, by=c("GeneralArea" = "diversion_to_sensorLoc"))
sensorLocations$Freq[is.na(sensorLocations$Freq)]<-0
sensorLocations <- left_join(sensorLocations, sum_Volume, by=c("GeneralArea"="sensorLoc"))
sensorLocations$volume[is.na(sensorLocations$volume)]<-0


diversion_risk <- rep(0,length(fishWYIDs))

for(fish in 1:length(fishWYIDs)){
  
  for(ll in 1:ncol(occupancy_matrix)){
    
    diversion_risk[fish] <- (diversion_risk[fish] +
                               log10(sensorLocations$volume[sensorLocations$GeneralArea==goodLocs[ll]] *
                               occupancy_matrix[fish,ll]+1))
  }
  
}


group_diversion_risk <- c(mean(diversion_risk[fish.summary$age_group=="Below Slot"]),
                          mean(diversion_risk[fish.summary$age_group=="Fishery"]),
                          mean(diversion_risk[fish.summary$age_group=="Above Slot"]))
# group_diversion_risk.q25<- c(quantile(diversion_risk[fish.summary$age_group=="Below Slot"], 0.25),
#                              quantile(diversion_risk[fish.summary$age_group=="Fishery"], 0.25),
#                              quantile(diversion_risk[fish.summary$age_group=="Above Slot"], 0.25))
# group_diversion_risk.q75<- c(quantile(diversion_risk[fish.summary$age_group=="Below Slot"], 0.75),
#                              quantile(diversion_risk[fish.summary$age_group=="Fishery"], 0.75),
#                              quantile(diversion_risk[fish.summary$age_group=="Above Slot"], 0.75))
group_diversion_risk.sd <- c(sd(diversion_risk[fish.summary$age_group=="Below Slot"]),
                             sd(diversion_risk[fish.summary$age_group=="Fishery"]),
                             sd(diversion_risk[fish.summary$age_group=="Above Slot"]))


pal=brewer.pal(3,"Set2")

png("~/GitHub/fishsync/WhiteSturgeonOccupancy/diversion_risk.png", units="in", 
    res=200, width=4.5, height=5.5)
par(mar=c(4.1,4.1,1,1))
bp<-barplot(group_diversion_risk, ylab="Diversion exposure index", ylim=c(0,13),
            names.arg=c("Below Slot", "Fishery", "Above Slot"), col=pal)
arrows(x0=c(bp), 
       y0=group_diversion_risk-(group_diversion_risk.sd)/sqrt(c(24,114,119)),
       x1=(bp),
       y1=group_diversion_risk+(group_diversion_risk.sd)/sqrt(c(24,114,119)),
       angle=90, code=3, length=0.1)

dev.off()



## Fishing pressure -------------------------------------------------------------------------------

reportcards.raw<-read.csv("~/GitHub/fishsync/WhiteSturgeonOccupancy/White_Sturgeon_ReportCard_Data_20210414.csv")
loc_to_region<-read.csv("~/GitHub/fishsync/WhiteSturgeonOccupancy/detectLoc_to_reportcardRegion.csv")

#return_rate <- data.frame(year = 2007:2020,
#                          rate = c(18.36, 13.63, 11.80, 11.33, 11.08, 11.20, 20.90,
#                                   24.51, 29.75, 32.92, 33.52, 32.55, 30.84, 31.49)/100)

#return_rate <- return_rate[return_rate$year >=2015,]

#format data for analysis
loc_to_region <- loc_to_region[loc_to_region$detectLocation %in% goodLocs,]
reportcards<-reportcards.raw[reportcards.raw$Year >= 2015 & reportcards.raw$Year <= 2020,]
studyArea <- c("SFBay_south_of_hwy80br","SFBay_north_of_hwy80br","SanPabloBay",
               "CarquinezStrait","SuisunBay","MontezumaSlough","SacR_RioVista_to_ChippsIsland",
               "SacR_DeepW_Ship_Canal","YoloBypass","SacR_KnightsLanding_to_RioVista",
               "SacR_Colusa_to_KnightsLanding","SacR_hwy32br_to_Colusa",
               "SacR_RedBluff_to_hwy32br","SacR_Upstream_of_RedBluff")

reportcards <- reportcards[reportcards$Location %in% studyArea,]



reportcards <- reportcards[!is.na(reportcards$Location),]
kept <- reportcards[reportcards$Fate=="kept",]
# kept_by_region <- data.frame(table(kept$Location))
# caught_by_region <- data.frame(table(reportcards$Location))

ncaught_loc_yr <- as.matrix(table(reportcards$Location, reportcards$Year))
ncaught_loc_yr <- ncaught_loc_yr[match(studyArea, rownames(ncaught_loc_yr)),]
nkept_loc_yr <- as.matrix(table(kept$Location, kept$Year))
nkept_loc_yr <- nkept_loc_yr[match(studyArea, rownames(nkept_loc_yr)),]

nregions <- length(unique(reportcards$Location))

# retrate_mat <- NULL
# for(ii in 1:nregions){
#   retrate_mat <- rbind(retrate_mat, return_rate$rate)
# }
# 
# ncaught_loc_yr.corr <- ncaught_loc_yr/retrate_mat
# nkept_loc_yr.corr <- nkept_loc_yr/retrate_mat

caught_by_region <- data.frame(region = rownames(ncaught_loc_yr),
                               n = rowMeans(ncaught_loc_yr))
kept_by_region <- data.frame(region = rownames(nkept_loc_yr),
                               n = rowMeans(nkept_loc_yr))


detectLocs <- colnames(occupancy_matrix)
detectLocs <- gsub("\\."," ",detectLocs)


## calculate risk index for kept fish
risk_index_kept <- rep(0, length(fishWYIDs))

for(fish in 1:length(fishWYIDs)){
  
  for(ll in 1:nrow(loc_to_region)){
    
    if(loc_to_region$detectLocation[ll] %in% c("Richmond Bridge","Carquinez Bridge",
                                               "Benicia Bridge","SR_RioVista")){
      risk_index_kept[fish] <- (risk_index_kept[fish] + 
                                  0.5* occupancy_matrix[fish,detectLocs==loc_to_region$detectLocation[ll]] *
                                  (kept_by_region$n[kept_by_region$region==loc_to_region$reportcardRegion[ll]]/
                                     sum(loc_to_region$reportcardRegion==loc_to_region$reportcardRegion[ll]))
      )
    }
    
    else{
      risk_index_kept[fish] <- (risk_index_kept[fish] + 
                                  occupancy_matrix[fish,detectLocs==loc_to_region$detectLocation[ll]] *
                                  (kept_by_region$n[kept_by_region$region==loc_to_region$reportcardRegion[ll]]/
                                     sum(loc_to_region$reportcardRegion==loc_to_region$reportcardRegion[ll]))
      )
    }
  }
}

risk_index_kept <- risk_index_kept/sum(kept_by_region$n)#[kept_by_region$region %in% loc_to_region$reportcardRegion])


## calculate risk index for caught fish
risk_index_caught <- rep(0, length(fishWYIDs))

for(fish in 1:length(fishWYIDs)){
  
  for(ll in 1:nrow(loc_to_region)){
    
    if(loc_to_region$detectLocation[ll] %in% c("Richmond Bridge","Carquinez Bridge",
                                               "Benicia Bridge","SR_RioVista")){
      risk_index_caught[fish] <- (risk_index_caught[fish] + 
                                    0.5* occupancy_matrix[fish,detectLocs==loc_to_region$detectLocation[ll]] *
                                    (caught_by_region$n[caught_by_region$region==loc_to_region$reportcardRegion[ll]]/
                                       sum(loc_to_region$reportcardRegion==loc_to_region$reportcardRegion[ll]))
      )
    }
    
    else{
      risk_index_caught[fish] <- (risk_index_caught[fish] + 
                                    occupancy_matrix[fish,detectLocs==loc_to_region$detectLocation[ll]] *
                                    (caught_by_region$n[caught_by_region$region==loc_to_region$reportcardRegion[ll]]/
                                       sum(loc_to_region$reportcardRegion==loc_to_region$reportcardRegion[ll]))
      )
    }
  }
}

risk_index_caught <- risk_index_caught/sum(caught_by_region$n)#$Freq[caught_by_region$Var1 %in% loc_to_region$reportcardRegion])



group_risk_index_kept <- c(mean(risk_index_kept[fish.summary$age_group=="Below Slot"]),
                           mean(risk_index_kept[fish.summary$age_group=="Fishery"]),
                           mean(risk_index_kept[fish.summary$age_group=="Above Slot"]))
# group_risk_index_kept.q25 <- c(quantile(risk_index_kept[fish.summary$age_group=="Below Slot"], 0.25),
#                                quantile(risk_index_kept[fish.summary$age_group=="Fishery"], 0.25),
#                                quantile(risk_index_kept[fish.summary$age_group=="Above Slot"], 0.25))
# group_risk_index_kept.q75 <- c(quantile(risk_index_kept[fish.summary$age_group=="Below Slot"], 0.75),
#                                quantile(risk_index_kept[fish.summary$age_group=="Fishery"], 0.75),
#                                quantile(risk_index_kept[fish.summary$age_group=="Above Slot"], 0.75))
group_risk_index_kept.sd <- c(sd(risk_index_kept[fish.summary$age_group=="Below Slot"]),
                           sd(risk_index_kept[fish.summary$age_group=="Fishery"]),
                           sd(risk_index_kept[fish.summary$age_group=="Above Slot"]))

group_risk_index_caught <- c(mean(risk_index_caught[fish.summary$age_group=="Below Slot"]),
                             mean(risk_index_caught[fish.summary$age_group=="Fishery"]),
                             mean(risk_index_caught[fish.summary$age_group=="Above Slot"]))
# group_risk_index_caught.q25<- c(quantile(risk_index_caught[fish.summary$age_group=="Below Slot"], 0.25),
#                                 quantile(risk_index_caught[fish.summary$age_group=="Fishery"], 0.25),
#                                 quantile(risk_index_caught[fish.summary$age_group=="Above Slot"], 0.25))
# group_risk_index_caught.q75<- c(quantile(risk_index_caught[fish.summary$age_group=="Below Slot"], 0.75),
#                                 quantile(risk_index_caught[fish.summary$age_group=="Fishery"], 0.75),
#                                 quantile(risk_index_caught[fish.summary$age_group=="Above Slot"], 0.75))
group_risk_index_caught.sd <- c(sd(risk_index_caught[fish.summary$age_group=="Below Slot"]),
                             sd(risk_index_caught[fish.summary$age_group=="Fishery"]),
                             sd(risk_index_caught[fish.summary$age_group=="Above Slot"]))

pal1 = brewer.pal(3, "Set2")
pal2 = brewer.pal(3, "Dark2")

png("~/GitHub/fishsync/WhiteSturgeonOccupancy/fishing_risk_wy.png", res=300, units="in", width=4.5, height=5.5)
par(mar=c(4.1,4.1,1,1))
bp<-barplot(rbind(group_risk_index_caught,group_risk_index_kept), beside=T,
            #density=c(0,10),
            col=rep(pal1, each=2), ylim=c(0,10),
            legend.text = NULL, ylab="Fishing capture index",
            names.arg = c("Below slot", "Fishery", "Above slot"),
            args.legend=list(x="topright", inset=0.05))
barplot(rbind(group_risk_index_caught,group_risk_index_kept), beside=T,
        density=c(0,10), add=TRUE,
        col=NA,
        legend.text = c("caught","kept"), ylab="Fishing capture index",
        names.arg = NULL,
        args.legend=list(x="topright", inset=0.05))
arrows(x0=c(bp), 
       y0=c(rbind(group_risk_index_caught,group_risk_index_kept)) -
          c(rbind(group_risk_index_caught.sd/sqrt(c(24,114,119)),
                  group_risk_index_kept.sd/sqrt(c(24,114,119)))),
       x1=c(bp),
       y1=c(rbind(group_risk_index_caught,group_risk_index_kept)) +
          c(rbind(group_risk_index_caught.sd/sqrt(c(24,114,119)),
                  group_risk_index_kept.sd/sqrt(c(24,114,119)))),
       angle=90, code=3, length=0.1)

dev.off()

barplot(colSums(ncaught_loc_yr.corr))
barplot(rowMeans(ncaught_loc_yr.corr), las=2)

png("~/GitHub/fishsync/WhiteSturgeonOccupancy/risk_indices_combined.png", units="in", 
    res=300, width=6.5, height=4)
par(mar=c(3.1,4.1,1,1), mfrow=c(1,2))

#bp<-barplot(rbind(group_risk_index_caught,group_risk_index_kept), beside=T,
bp<-barplot(group_risk_index_caught,
            #density=c(0,10),
            col=pal, ylim=c(0,11),
            legend.text = NULL, ylab="Fishing capture index",
            names.arg = c("Below", "Fishery", "Above"),
            args.legend=list(x="topright", inset=0.05))
# barplot(rbind(group_risk_index_caught,group_risk_index_kept), beside=T,
#         density=c(0,10), add=TRUE,
#         col=NA,
#         legend.text = c("caught","kept"), ylab="Fishing capture index",
#         names.arg = NULL,
#         args.legend=list(x="topright", inset=0.05))
# arrows(x0=c(bp), 
#        y0=c(rbind(group_risk_index_caught,group_risk_index_kept)) -
#          c(rbind(group_risk_index_caught.sd/sqrt(c(24,114,119)),
#                  group_risk_index_kept.sd/sqrt(c(24,114,119)))),
#        x1=c(bp),
#        y1=c(rbind(group_risk_index_caught,group_risk_index_kept)) +
#          c(rbind(group_risk_index_caught.sd/sqrt(c(24,114,119)),
#                  group_risk_index_kept.sd/sqrt(c(24,114,119)))),
#        angle=90, code=3, length=0.1)
arrows(x0=c(bp), 
       y0=group_risk_index_caught -
         group_risk_index_caught.sd/sqrt(c(24,114,119)),
       x1=c(bp),
       y1=group_risk_index_caught +
         group_risk_index_caught.sd/sqrt(c(24,114,119)),
       angle=90, code=3, length=0.1)

bp<-barplot(group_diversion_risk, ylab="Habitat modification index", ylim=c(0,13),
            names.arg=c("Below", "Fishery", "Above"), col=pal)
arrows(x0=c(bp), 
       y0=group_diversion_risk-(group_diversion_risk.sd)/sqrt(c(24,114,119)),
       x1=(bp),
       y1=group_diversion_risk+(group_diversion_risk.sd)/sqrt(c(24,114,119)),
       angle=90, code=3, length=0.1)

dev.off()




## bathymery analyses -----------------------------------------------------------------------------

bathy <- raster("/Users/jonathanwalter/Desktop/SanFranciscoBay/sfbaydeltadem10m2016.asc")

generalAreaPoints <- SpatialPoints(coords=cbind(generalAreaCoords$Lon, 
                                                generalAreaCoords$Lat),
                                   proj4string = CRS("+init=EPSG:4269"))

generalAreaPoints <- spTransform(generalAreaPoints, CRS(proj4string(bathy)))

plot(bathy)
points(generalAreaPoints)

maxdist=2 

generalAreaBuffer <- buffer(generalAreaPoints, width=maxdist*1000, dissolve=FALSE) #convert km to m




water.prj <- spTransform(water, CRS(proj4string(bathy)))

bathy_mask <- mask(bathy, water.prj)

quantile(values(bathy_mask), na.rm=T)

plot(bathy_mask)

median_depth <- rep(NA, length(generalAreaBuffer))

for(ii in 1:length(generalAreaBuffer)){
  tmp <- crop(bathy_mask, generalAreaBuffer[ii])
  tmp <- mask(tmp, generalAreaBuffer[ii])
  median_depth[ii] <- median(values(tmp), na.rm=TRUE)
}


depth_occupancy <- rep(NA, nrow(fish.summary))
depth_occupancy_BAY <- rep(NA, nrow(fish.summary))
depth_occupancy_DELTA <- rep(NA, nrow(fish.summary))

for(ii in 1:length(depth_occupancy)){
  depth_occupancy[ii] <- weighted.mean(median_depth, occupancy_matrix_norm[ii,])
  depth_occupancy_BAY[ii] <- weighted.mean(median_depth[1:4], occupancy_matrix_norm[ii,1:4])
  depth_occupancy_DELTA[ii] <- weighted.mean(median_depth[5:15], occupancy_matrix_norm[ii,5:15])
}

depth_occ_bySize <- c(
  mean(depth_occupancy[fish.summary$age_group=="Below Slot"]),
  mean(depth_occupancy[fish.summary$age_group=="Fishery"]),
  mean(depth_occupancy[fish.summary$age_group=="Above Slot"])
)

depth_occ_bySize.sd <- c(
  sd(depth_occupancy[fish.summary$age_group=="Below Slot"]),
  sd(depth_occupancy[fish.summary$age_group=="Fishery"]),
  sd(depth_occupancy[fish.summary$age_group=="Above Slot"])
)

depth_occ_bySize_BAY <- c(
  mean(depth_occupancy_BAY[fish.summary$age_group=="Below Slot"], na.rm=T),
  mean(depth_occupancy_BAY[fish.summary$age_group=="Fishery"], na.rm=T),
  mean(depth_occupancy_BAY[fish.summary$age_group=="Above Slot"], na.rm=T)
)

depth_occ_bySize_BAY.sd <- c(
  sd(depth_occupancy_BAY[fish.summary$age_group=="Below Slot"], na.rm=T),
  sd(depth_occupancy_BAY[fish.summary$age_group=="Fishery"], na.rm=T),
  sd(depth_occupancy_BAY[fish.summary$age_group=="Above Slot"], na.rm=T)
)

depth_occ_bySize_DELTA <- c(
  mean(depth_occupancy_DELTA[fish.summary$age_group=="Below Slot"], na.rm=T),
  mean(depth_occupancy_DELTA[fish.summary$age_group=="Fishery"], na.rm=T),
  mean(depth_occupancy_DELTA[fish.summary$age_group=="Above Slot"], na.rm=T)
)

depth_occ_bySize_DELTA.sd <- c(
  sd(depth_occupancy_DELTA[fish.summary$age_group=="Below Slot"], na.rm=T),
  sd(depth_occupancy_DELTA[fish.summary$age_group=="Fishery"], na.rm=T),
  sd(depth_occupancy_DELTA[fish.summary$age_group=="Above Slot"], na.rm=T)
)



png("~/GitHub/fishsync/WhiteSturgeonOccupancy/depth.png", units="in", 
    res=300, width=6.5, height=3.5)

par(mfrow=c(1,3), mar=c(2.6,4.1,1.6,1.1))

bp=barplot(depth_occ_bySize, ylim=c(-12,0), col=pal, names.arg=c("Below","Fishery","Above"),
           ylab="Elevation (m)")
arrows(x0=c(bp), 
       y0=depth_occ_bySize-(depth_occ_bySize.sd)/sqrt(c(24,114,119)),
       x1=(bp),
       y1=depth_occ_bySize+(depth_occ_bySize.sd)/sqrt(c(24,114,119)),
       angle=90, code=3, length=0.1)
mtext("All Study Area", line=0.2, cex=2/3)
mtext("a)", line=0.2, cex=2/3, at=0)

bp=barplot(depth_occ_bySize_BAY, ylim=c(-13,0), col=pal, names.arg=c("Below","Fishery","Above"),
           ylab="Elevation (m)")
arrows(x0=c(bp), 
       y0=depth_occ_bySize_BAY-(depth_occ_bySize_BAY.sd)/sqrt(c(24,114,119)),
       x1=(bp),
       y1=depth_occ_bySize_BAY+(depth_occ_bySize_BAY.sd)/sqrt(c(24,114,119)),
       angle=90, code=3, length=0.1)
mtext("Bay", line=0.2, cex=2/3)
mtext("b)", line=0.2, cex=2/3, at=0)

bp=barplot(depth_occ_bySize_DELTA, ylim=c(-9,0), col=pal, names.arg=c("Below","Fishery","Above"),
           ylab="Elevation (m)")
arrows(x0=c(bp), 
       y0=depth_occ_bySize_DELTA-(depth_occ_bySize_DELTA.sd)/sqrt(c(24,114,119)),
       x1=(bp),
       y1=depth_occ_bySize_DELTA+(depth_occ_bySize_DELTA.sd)/sqrt(c(24,114,119)),
       angle=90, code=3, length=0.1)
mtext("Delta & SR", line=0.2, cex=2/3)
mtext("c)", line=0.2, cex=2/3, at=0)

dev.off()
