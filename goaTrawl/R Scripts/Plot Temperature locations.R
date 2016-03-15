rm(list=ls())
library(INLA)
library(rgdal)
library(ggplot2)
library(splancs)
library(sp)

# Use only species with length-based data, or total biomass
totalBiomass = FALSE
#### CHOOSE A MODEL - "binomial" or "positive"
model = "positive"
single.intercept =TRUE

#Switch for doing cross-validation (TRUE) or not (FALSE)
LEAVE.OUT   <- FALSE
frac.leave.out  <- 0.10

max.depth <- 600 # in meters

#Define directories for the data and for the plots
proj.dir = paste(getwd(),"/goaTrawl",sep="")
plot.dir = paste(proj.dir,"/_Output plots Pres Abs",sep="")
if(model != "binomial") plot.dir = paste(proj.dir,"/_Output plots Pos",sep="")

#### GO GET THE OBSERVED TRAWL DATA
setwd(proj.dir)
df = read.csv("Output Data/goa_trawl_final_albers+temp.csv")
# include switch for whether the length data will be used instead
if(totalBiomass == FALSE) {
  # then use the length data
  df = read.csv("Output Data/goa_trawl_final_size_albers+temp.csv")
}

df = df[order(df$Year,df$Lat),]

### Remove NA entries in BottomDepth
df = df[df$BottomDepth != -9999,]
#df = df[df$BottomTemp != -9999,]

df = df[df$BottomDepth <= max.depth,] 

df$Station = as.character(df$Station)
df$Year = as.numeric(as.character(df$Year))
df$LonUTMAlbers = df$LonUTMAlbers/1000
df$LatUTMAlbers = df$LatUTMAlbers/1000

#### GO GET THE PROJECTION POINTS
dat.project	= read.csv("Output Data/goa_projection_points+temp.csv")

dat.project$LonUTMAlbers = dat.project$LonUTMAlbers/1000
dat.project$LatUTMAlbers = dat.project$LatUTMAlbers/1000

#### Exclude points that end up on land. Check this with Ole, not working -- units may be off
dat.project$NGDC24_M =	-dat.project$NGDC24_M	# depth in m
dat.project$SRTM_M = -dat.project$SRTM_M	# depth in m
dat.project = dat.project[dat.project$NGDC24 > 0,]

THESE	= seq(1,nrow(dat.project),by=15)
dat.trim = dat.project[THESE,]

dat.trim = dat.project[THESE,]
dat.new.trim = dat.trim
# x.lim = c(-200,200)
# x.lim = c(min(dat.project$LonUTMAlbers),max(dat.project$LonUTMAlbers))
# x.lim = c(-800,550)
# y.lim = c(5700,6200)
# y.lim = c(min(dat.project$LatUTMAlbers),max(dat.project$LatUTMAlbers))
# plot(LatUTMAlbers~LonUTMAlbers,data=dat.project[ #&
# 		 dat.project$LonUTMAlbers < 200 & dat.project$LonUTMAlbers > -200
# 		 ,],pch=".",xlim=x.lim,ylim=y.lim)
# par(new=T)
# plot(LatUTMAlbers~LonUTMAlbers,data=dat.trim[dat.trim$NGDC24_M<=500 &
# 	 dat.trim$LonUTMAlbers < 550 & dat.trim$LonUTMAlbers > -800 &
# 	 dat.trim$NGDC24_M > 25
# 	 ,],pch=".",xlim=x.lim,ylim=y.lim,col=2)
# 
# par(new=T)
# plot(LatUTMAlbers~LonUTMAlbers,data=dat.new.trim[dat.new$LonUTMAlbers < 200 & dat.new$LonUTMAlbers > -200,]
#  ,pch=".",xlim=x.lim,ylim=y.lim,col=4)
# dim(dat.new.trim)
# dat.proj.trim	= dat.new.trim

# This is a switch for whether we're using full data, or length-stratified data
maxCol = dim(df)[2]
species = names(df)[17:maxCol]
nCovCol = 16
if(totalBiomass==FALSE) {
  maxCol = 31
  species = names(df)[14:maxCol] 
  temp <- NULL
  for(ZZ in 1:length(species)){
    temp[ZZ] <- substr(species[ZZ],nchar(species[ZZ])-3,nchar(species[ZZ]))
  }
  species <- species[temp =="cpue"]
  nCovCol = 13
}

#########################################################################################################################

### Basic plotting files
NAME	<-	paste(plot.dir,"/Location of Temp data.pdf",sep="")
pdf(NAME,onefile=TRUE,width=15,5)

#### PLOT SPATIAL DISTRIBUTION
# Import the datafile of Alaska shoreline

plot.df <- df
plot.df$BottomTemp[plot.df$BottomTemp >0 ] <- 1
plot.df$BottomTemp[plot.df$BottomTemp <0 ] <- 0
plot.df$BottomTemp <- as.factor(plot.df$BottomTemp)

setwd(paste(plot.dir,"/_Alaska Shapefile",sep=""))
shp.alaska	 <-	readOGR(dsn=".",layer="Alaska-Albers")
dat.alaska   <- fortify(shp.alaska,"data.frame")

# Some Arguments for making prettier plots
bGrid <-theme(panel.grid =element_blank())
bBack <-theme(panel.background =element_blank())
bAxis <-theme(axis.title.y =element_blank())
bTics <-theme(axis.text =element_blank(), axis.text.y =element_blank(), axis.ticks =element_blank())

### Plot presence-absence in space
q1	<-	ggplot() +
  scale_size_manual(values=c(0.8,0.8)) +
  scale_colour_manual(values=c("black","red"),
                      name="")+
  geom_point(data=plot.df,alpha=0.4,
             aes(LonUTMAlbers*1000,LatUTMAlbers*1000,colour=BottomTemp))+
  geom_polygon(data=dat.alaska, fill=grey(0.4),color=NA,aes(long,lat,group=group)) +
  labs(x = "Eastings",y="Northings",title=paste("Bottom Temperature Data")) +
  coord_cartesian(xlim = c(min(df$LonUTMAlbers*1000)*1.01,max(df$LonUTMAlbers*1000)), 
                  ylim = c(min(df$LatUTMAlbers*1000),6300000))+
  bGrid  + bBack

YEARS <- sort(unique(df$Year))
p	<-	list()

for(j in 1:length(YEARS)){
  p[[j]]	<-	ggplot() +
    scale_colour_manual(values=c("black","red"),
                        name="")+
    geom_point(data=plot.df[plot.df$Year==YEARS[j],],alpha=0.4,
               aes(LonUTMAlbers*1000,LatUTMAlbers*1000,colour=BottomTemp))+
    geom_polygon(data=dat.alaska, fill=grey(0.4),color=NA,aes(long,lat,group=group)) +
      labs(x = "Eastings",y="Northings",title=paste(YEARS[j],"Bottom Temperature Data")) +
      coord_cartesian(xlim = c(min(df$LonUTMAlbers*1000)*1.01,max(df$LonUTMAlbers*1000)), 
                      ylim = c(min(df$LatUTMAlbers*1000),6300000))+
      bGrid  + bBack
}

## write plots to file
print(q1)
for(j in 1:length(YEARS)){
  print(p[j])
}

dev.off()

setwd(proj.dir)
