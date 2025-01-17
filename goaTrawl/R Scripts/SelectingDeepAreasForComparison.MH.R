#Mary modified Ole's code for selecting discrete areas for comparison.
#This code provides projection location for 5 deep areas in the central/western GoA

rm(list=ls())
library(INLA)
library(rgdal)
library(ggplot2)
library(splancs)
library(dplyr)
library(sp)
library(readr)
library(vegan)


#### GO GET THE PROJECTION POINTS
#proj.dir	<- "/Users/ole.shelton/GitHub/exxonValdez_nceas/goaTrawl/"
#setwd(proj.dir)
setwd("~/Desktop/Diversity")
IDinfo<-read.csv("goa_projection_points_temp.csv")
head(IDinfo)

#### GO GET THE PROJECTION POINTS
dat.project	= IDinfo
dat.project$LonUTMAlbers = dat.project$LonUTMAlbers/1000
dat.project$LatUTMAlbers = dat.project$LatUTMAlbers/1000

dat.project$depth = -dat.project$NGDC24_M 
### Go Get the shapefile to provide a shoreline.
#setwd(paste(proj.dir,"/Output plots/_Alaska Shapefile",sep=""))
shp.alaska	 <-	readOGR(dsn=".",layer="Alaska-Albers")
dat.alaska   <- fortify(shp.alaska,"data.frame")
dat.alaska$long.km	<-	dat.alaska$long/1000
dat.alaska$lat.km	<-	dat.alaska$lat/1000

####### 
## Make some preliminary plots
#######

# Some Arguments for making prettier plots
bGrid <-theme(panel.grid =element_blank())
bBack <-theme(panel.background =element_blank())
bAxis <-theme(axis.title.y =element_blank())
bTics <-theme(axis.text =element_blank(), axis.text.y =element_blank(), axis.ticks =element_blank())

## Just a map of depth
z.lim		<-	c(0,max(dat.project$depth)+1)

quartz()
p1	<-	ggplot() +
  scale_size(range = c(1,1))+
  scale_colour_gradientn(limits=z.lim,colours= c("#98F5FF","black"))+
  geom_point(data= dat.project,alpha=0.3,
             mapping=aes(LonUTMAlbers,LatUTMAlbers,colour=depth )) + 
  geom_polygon(data=dat.alaska, fill=grey(0.4),color=NA,aes(long.km,lat.km,group=group)) +
  labs(x = "Eastings",y="Northings",title=paste("Depth"))+
  coord_cartesian(xlim = c(min(dat.project$LonUTMAlbers),max(dat.project$LonUTMAlbers)), ylim = c(min(dat.project$LatUTMAlbers),max(dat.project$LatUTMAlbers)))+
  bGrid  + bBack 
p1


# Subset to look at areas of defined depth
MIN.D	<-	150
MAX.D	<-	300
x.lim	<-	c(min(dat.project$LonUTMAlbers),max(dat.project$LonUTMAlbers))
y.lim	<-	c(min(dat.project$LatUTMAlbers),max(dat.project$LatUTMAlbers))


##### ONLY PLOT SPECIFIED DEPTH RANGE
quartz()
z.lim	<-	c(MIN.D, MAX.D)
p2	<-	ggplot() +
  scale_size(range = c(1,1))+
  scale_colour_gradientn(limits=z.lim,colours= c("red"))+
  geom_point(data= dat.project,alpha=0.3,
             mapping=aes(LonUTMAlbers,LatUTMAlbers,colour=depth )) + 
  geom_polygon(data=dat.alaska, fill=grey(0.4),color=NA,aes(long.km,lat.km,group=group)) +
  labs(x = "Eastings",y="Northings",title=paste("Depth"))+
  coord_cartesian(xlim = x.lim, ylim = y.lim)+
  bGrid  + bBack 
p2

#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
###### Below makes the file "goa_central_gulf(.......).csv"
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
MIN.D	<-	150
MAX.D	<-	300


######## Region 12: East of Shelikof Strait
trim.x	<- c(-25,150)
trim.y	<- c(5967,6035)
x.lim	<-	c(-200,200)
y.lim	<-	c(5850,6150)
#MIN.D<-50
#MAX.D<-300

#trim.x	<- c(-25,150)
#trim.y	<- c(5967,6035)

#6020 108
#5994 -16

dat.select	<-	dat.project[dat.project$depth >= MIN.D & dat.project$depth <= MAX.D &
                            dat.project$LonUTMAlbers >= trim.x[1] &	dat.project$LonUTMAlbers <= trim.x[2] &
                            dat.project$LatUTMAlbers >= trim.y[1] &	dat.project$LatUTMAlbers <= trim.y[2]
                          ,]

dat.select	<-	dat.select[order(dat.select$LatUTMAlbers,dat.select$LonUTMAlbers,decreasing=F),]

dat.select	<- dat.select[dat.select$LatUTMAlbers <= 6020 |
                           dat.select$LonUTMAlbers <=108 ,]
dat.select	<- dat.select[dat.select$LatUTMAlbers >= 5994 |
                           dat.select$LonUTMAlbers >= -16 ,]


quartz()
p3	<-	ggplot() +
  scale_size(range = c(1,1))+
  scale_colour_gradientn(limits=z.lim,colours= c("red"))+
  geom_point(data= dat.select,alpha=0.3,
             mapping=aes(LonUTMAlbers,LatUTMAlbers,colour=depth )) + 
  geom_polygon(data=dat.alaska, fill=grey(0.4),color=NA,aes(long.km,lat.km,group=group)) +
  labs(x = "Eastings",y="Northings",title=paste("Depth"))+
  coord_cartesian(xlim = x.lim, ylim = y.lim)+
  bGrid  + bBack 
p3

dat.region	<- dat.select
dat.region$Area	<- 12

######## Region 13: Shelikof Strait
trim.x	<- c(-285,-160)
trim.y	<- c(5850,5975)
x.lim	<-	c(-400,0)
y.lim	<-	c(5750,6000)
#MIN.D<-50
#MAX.D<-300


dat.select	<-	dat.project[dat.project$depth >= MIN.D & dat.project$depth <= MAX.D &
                            dat.project$LonUTMAlbers >= trim.x[1] &	dat.project$LonUTMAlbers <= trim.x[2] &
                            dat.project$LatUTMAlbers >= trim.y[1] &	dat.project$LatUTMAlbers <= trim.y[2]
                          ,]

dat.select	<-	dat.select[order(dat.select$LatUTMAlbers,dat.select$LonUTMAlbers,decreasing=F),]

dat.select	<- dat.select[dat.select$LatUTMAlbers >= 5885 |
                           dat.select$LonUTMAlbers <= -245 ,]
dat.select	<- dat.select[dat.select$LatUTMAlbers >= 5905 |
                           dat.select$LonUTMAlbers <= -218 ,]
dat.select	<- dat.select[dat.select$LatUTMAlbers >= 5913 |
                           dat.select$LonUTMAlbers <= -202 ,]
dat.select	<- dat.select[dat.select$LatUTMAlbers >= 5935 |
                           dat.select$LonUTMAlbers <= -185 ,]

quartz()
p3	<-	ggplot() +
  scale_size(range = c(1,1))+
  scale_colour_gradientn(limits=z.lim,colours= c("red"))+
  geom_point(data= dat.select,alpha=0.3,
             mapping=aes(LonUTMAlbers,LatUTMAlbers,colour=depth )) + 
  geom_polygon(data=dat.alaska, fill=grey(0.4),color=NA,aes(long.km,lat.km,group=group)) +
  labs(x = "Eastings",y="Northings",title=paste("Depth"))+
  coord_cartesian(xlim = x.lim, ylim = y.lim)+
  bGrid  + bBack 
p3


dat.select$Area <- 13
dat.region <- rbind(dat.region,dat.select)

######## Region 14 - West of Shelikof Strait
trim.x	<- c(-428,-355)
trim.y	<- c(5620,5760)
x.lim	<-	c(-500,-200)
y.lim	<-	c(5500,6000)
#MIN.D<-50
#MAX.D<-300
#trim.x	<- c(-580,-512)
#trim.y	<- c(5580,5700)
#5620
dat.select	<-	dat.project[dat.project$depth >= MIN.D & dat.project$depth <= MAX.D &
                            dat.project$LonUTMAlbers >= trim.x[1] &	dat.project$LonUTMAlbers <= trim.x[2] &
                            dat.project$LatUTMAlbers >= trim.y[1] &	dat.project$LatUTMAlbers <= trim.y[2]
                          ,]

dat.select	<-	dat.select[order(dat.select$LatUTMAlbers,dat.select$LonUTMAlbers,decreasing=F),]
dat.select	<- dat.select[dat.select$LatUTMAlbers >= 5628 |
                           dat.select$LonUTMAlbers <= -375,]

quartz()
p3	<-	ggplot() +
  scale_size(range = c(1,1))+
  scale_colour_gradientn(limits=z.lim,colours= c("red"))+
  geom_point(data= dat.select,alpha=0.3,
             mapping=aes(LonUTMAlbers,LatUTMAlbers,colour=depth )) + 
  geom_polygon(data=dat.alaska, fill=grey(0.4),color=NA,aes(long.km,lat.km,group=group)) +
  labs(x = "Eastings",y="Northings",title=paste("Depth"))+
  coord_cartesian(xlim = x.lim, ylim = y.lim)+
  bGrid  + bBack 
p3


dat.select$Area <- 14
dat.region <- rbind(dat.region,dat.select)

######## Region 15 - West west of Shelikof Strait
trim.x	<- c(-580,-512)
trim.y	<- c(5580,5685)
x.lim	<-	c(-650,-450)
y.lim	<-	c(5500,6000)
#MIN.D<-50
#MAX.D<-300
#trim.x	<- c(-580,-512)
#trim.y	<- c(5580,5700)

dat.select	<-	dat.project[dat.project$depth >= MIN.D & dat.project$depth <= MAX.D &
                            dat.project$LonUTMAlbers >= trim.x[1] &	dat.project$LonUTMAlbers <= trim.x[2] &
                            dat.project$LatUTMAlbers >= trim.y[1] &	dat.project$LatUTMAlbers <= trim.y[2]
                          ,]

dat.select	<-	dat.select[order(dat.select$LatUTMAlbers,dat.select$LonUTMAlbers,decreasing=F),]

dat.select	<- dat.select[dat.select$LatUTMAlbers >= 5630 |
                           dat.select$LonUTMAlbers >= -565 ,]

quartz()
p3	<-	ggplot() +
  scale_size(range = c(1,1))+
  scale_colour_gradientn(limits=z.lim,colours= c("red"))+
  geom_point(data= dat.select,alpha=0.3,
             mapping=aes(LonUTMAlbers,LatUTMAlbers,colour=depth )) + 
  geom_polygon(data=dat.alaska, fill=grey(0.4),color=NA,aes(long.km,lat.km,group=group)) +
  labs(x = "Eastings",y="Northings",title=paste("Depth"))+
  coord_cartesian(xlim = x.lim, ylim = y.lim)+
  bGrid  + bBack 
p3

dat.select$Area <- 15
dat.region <- rbind(dat.region,dat.select)

######## Region 16 - Outershelf off Kodiak
trim.x	<- c(-200,45)
trim.y	<- c(5600,5930)
x.lim	<-	c(-250,250)
y.lim	<-	c(5550,6000)


dat.select	<-	dat.project[dat.project$depth >= MIN.D & dat.project$depth <= MAX.D &
                            dat.project$LonUTMAlbers >= trim.x[1] &	dat.project$LonUTMAlbers <= trim.x[2] &
                            dat.project$LatUTMAlbers >= trim.y[1] &	dat.project$LatUTMAlbers <= trim.y[2]
                          ,]

dat.select	<-	dat.select[order(dat.select$LatUTMAlbers,dat.select$LonUTMAlbers,decreasing=F),]

dat.select	<- dat.select[dat.select$LatUTMAlbers <= 5811 |
                           dat.select$LonUTMAlbers >=-75 ,]
dat.select	<- dat.select[dat.select$LatUTMAlbers <= 5880 |
                           dat.select$LonUTMAlbers >= -41 ,]
dat.select	<- dat.select[dat.select$LatUTMAlbers <= 5750 |
                           dat.select$LonUTMAlbers >= -140 ,]


quartz()
p3	<-	ggplot() +
  scale_size(range = c(1,1))+
  scale_colour_gradientn(limits=z.lim,colours= c("red"))+
  geom_point(data= dat.select,alpha=0.3,
             mapping=aes(LonUTMAlbers,LatUTMAlbers,colour=depth )) + 
  geom_polygon(data=dat.alaska, fill=grey(0.4),color=NA,aes(long.km,lat.km,group=group)) +
  labs(x = "Eastings",y="Northings",title=paste("Depth"))+
  coord_cartesian(xlim = x.lim, ylim = y.lim)+
  bGrid  + bBack 
p3

dat.select$Area <- 16
dat.region <- rbind(dat.region,dat.select)
dat.region<-dat.region[,c(1:9,24,25)]



write.csv(dat.region,"Deep_goa_discrete_areas_for_comparison_150_to_300m.MH.csv")

