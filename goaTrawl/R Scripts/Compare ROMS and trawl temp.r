library(ncdf)
library(reshape2)
library(rgdal)
library(mgcv)
library(ggplot2)
library(scales)

#import data
nc <- open.ncdf("/Users/ole.shelton/Documents/Science/Active projects/Exxon/ROMS models/shelikof.nc")

#what is inside
print(nc)

#pull out the variables
lon		<-	get.var.ncdf(nc,"LONGITUDE51_134")
lat		<-	get.var.ncdf(nc,"LATITUDE68_234")
Time	<-	get.var.ncdf(nc,"TIME1")
Depth	<-	get.var.ncdf(nc,"ZSALT")
Temp	<-	get.var.ncdf(nc,"TEMP")
Date	<- as.Date(as.POSIXct(nc$dim$TIME1$vals, origin="1900-01-01"))

# Get rid of extra months and years
	MONTHS	<-	which(substr(Date,6,7) == "05" | 
					substr(Date,6,7) == "06" | 
					substr(Date,6,7) == "07" |
					substr(Date,6,7) == "08" )
# 					substr(Date,6,7) == "09" |
# 					substr(Date,6,7) == "10" )

	Date.month	<-	Date[MONTHS]

	YEARS	<-	which(substr(Date,1,4) == "1999" | 
					substr(Date,1,4) == "2001" |
					substr(Date,1,4) == "2003" |
					substr(Date,1,4) == "2005" |
					substr(Date,1,4) == "2007" |
					substr(Date,1,4) == "2009" |
					substr(Date,1,4) == "2011" )

	these.temp		<-	YEARS[match(MONTHS,YEARS)]
	THESE  		<- 	these.temp[is.na(these.temp)==F]

Date.names <- as.character(Date[THESE])

Temp	<-	Temp[,,,c(THESE)]
dimnames(Temp)	<- list(Lon=lon,Lat=lat,Depth=Depth,Date=Date.names)

##########################################################################################
# Subset by year and create a 2D array
##########################################################################################

YEAR	<-	c("1999","2001","2003","2005","2007","2009","2011")
All.years	<- NULL

for (i in 1:length(YEAR)){

	THESE	<-	which(substr(Date.names,1,4)==YEAR[i])
	NAME.temp	<- paste("Y.",YEAR[i],sep="")	
	assign(NAME.temp,apply(Temp[,,,c(THESE)],c(1,2,3),mean,na.rm=T))

	temp.dat	<-	melt(eval(as.name(NAME.temp)),c("lon","lat","depth"))
	temp.dat	<-	temp.dat[is.nan(temp.dat$value)==F,]
	temp.dat	<-	cbind(Year = YEAR[i],temp.dat)

	All.years	<- rbind(All.years,temp.dat)
}

## Cull to only include the maximum depth available at each location
All.years.trim	<-	aggregate(All.years$depth,by=list(lon=All.years$lon,lat=All.years$lat),max)
colnames(All.years.trim)[3]	<-	"depth"

All.temp.dat 		<- merge(All.years,All.years.trim)
All.temp.dat$lon 	<-	 All.temp.dat$lon-360

All.temp.dat	<-	All.temp.dat[All.temp.dat$depth > 20,]
All.temp.dat$log.BottomDepth		<-	log(All.temp.dat$depth)

All.temp.dat	<- All.temp.dat[order(All.temp.dat$Year),]
##########################################################################################
##### Convert to Albers projection
##########################################################################################

aea.proj <- "+proj=aea +lat_1=51 +lat_2=62 +lon_0=-150 +x_0=0 +y_0=0 +datum=WGS84"

new	<-	 cbind(All.temp.dat$lon,All.temp.dat$lat)
new.sp	<- SpatialPoints(new,CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
new.sp.albers	<- spTransform(new.sp, CRS(aea.proj))

new.dat		<- as(new.sp.albers,"data.frame")

All.temp.dat$LonUTMAlbers	<-	new.dat$coords.x1
All.temp.dat$LatUTMAlbers	<-	new.dat$coords.x2

########################################################################################
#  Go get the Trawl Temperature Data 
########################################################################################
## This is a duplicate of the Temperature Map create code

# Read in the data
data.dir	<-	"/Users/ole.shelton/GitHub/exxonValdez_nceas/goaTrawl/"
plot.dir	<-	"/Users/ole.shelton/GitHub/exxonValdez_nceas/goaTrawl/_Temperature Surfaces/"

# Go get the data
setwd(data.dir)
trawl.dat	<-	read.csv(file="goa_trawl_final_albers+temp.csv",header=T)

### GET RID OF NAs
trawl.dat2	<-	trawl.dat[is.na(trawl.dat$BottomDepth)==F,]
trawl.trim	<-	trawl.dat[is.na(trawl.dat$BottomTemp)==F & 
							is.na(trawl.dat$BottomDepth)==F &
							trawl.dat$BottomTemp != -9999
								,1:14]
	
#### Fit a spatial GAM to the observed data for each year
all.dat	<- NULL

for(i in 1:length(YEAR)){
	dat.all.obs	<-	trawl.dat2[trawl.dat2$Year == YEAR[i],]
	dat.all.obs$log.BottomDepth	<- log(dat.all.obs$BottomDepth)

	dat	<-	trawl.trim[trawl.trim$Year == YEAR[i],]
	dat$log.BottomDepth	<- log(dat$BottomDepth)
	
 	out		<-	gam(BottomTemp ~ te(LonUTMAlbers,LatUTMAlbers,k=7)+s(BottomDepth),data=dat)
 	out.2	<-	gam(BottomTemp ~ te(LonUTMAlbers,LatUTMAlbers,k=7)+s(log.BottomDepth),data=dat)

	THESE	<-	c("LonUTMAlbers","LatUTMAlbers","log.BottomDepth")
	new.dat	<-	dat.all.obs[,THESE]
	new.points <- predict.gam(out.2,new.dat)
	
	dat.all.obs$Pred.bot.temp			<-	unlist(new.points)
	dat.all.obs$BottomTemp[dat.all.obs$BottomTemp == -9999]	<- NA
	dat.all.obs$Resid.bot.temp			<- dat.all.obs$BottomTemp - dat.all.obs$Pred.bot.temp

	all.dat	<- rbind(all.dat,dat.all.obs)
	
	### predict for all points in the projection set
	new.dat	<- All.temp.dat[All.temp.dat$Year == YEAR[i],c('LonUTMAlbers','LatUTMAlbers','log.BottomDepth')]
	All.temp.dat[All.temp.dat$Year == YEAR[i],"gam.pred.temp"]		<- predict.gam(out.2,new.dat)
}

All.temp.dat$Resid.temp <- All.temp.dat$value - All.temp.dat$gam.pred.temp



### Plot in Space

setwd(paste(plot.dir,"_Alaska Shapefile",sep=""))
shp.alaska	 <-	readOGR(dsn=".",layer="Alaska-Albers")
dat.alaska   <- fortify(shp.alaska,"data.frame")

# Some Arguments for making prettier plots
bGrid <-theme(panel.grid =element_blank())
bBack <-theme(panel.background =element_blank())
bAxis <-theme(axis.title.y =element_blank())
bTics <-theme(axis.text =element_blank(), axis.text.y =element_blank(), axis.ticks =element_blank())


### ADD trawl locations

## Get rid of deep areas
z.lim	<-	rep(max(abs(min(All.temp.dat$Resid.temp)),abs(max(All.temp.dat$Resid.temp))),2)
z.lim[1]	<- z.lim[1] * -1

z.lim=c(-3,3)

P	<-	list()
for(j in 1:length(YEAR)){
	p1	<-	ggplot() +
		scale_size(range = c(1,1))+
  		scale_colour_gradientn(limits=z.lim,colours= c("blue","white","red"))+
    	geom_point(data= All.temp.dat[All.temp.dat$Year == YEAR[j],] ,alpha=0.3,
    			mapping=aes(LonUTMAlbers,LatUTMAlbers,colour=	Resid.temp )) + 
		geom_polygon(data=dat.alaska, fill=grey(0.4),color=NA,aes(long,lat,group=group)) +
		labs(x = "Eastings",y="Northings",title=paste(YEAR[j],"Difference Bottom Temperature (C)"))+
   		coord_cartesian(xlim = c(min(trawl.dat$LonUTMAlbers),max(trawl.dat$LonUTMAlbers)), ylim = c(min(trawl.dat$LatUTMAlbers),6300000))+
 		bGrid  + bBack 
	P[[j]]	<-	p1
}


## Write plots to file
setwd(plot.dir)
pdf("Compare Trawl and ROMS Temperature.pdf",onefile=TRUE,width=15,5)

par(mfrow=c(2,4))
for(i in 1:length(YEAR)){
	lim	<-	c(min(All.temp.dat$value,All.temp.dat$gam.pred.temp),max(All.temp.dat$value,All.temp.dat$gam.pred.temp))
	plot(value~gam.pred.temp,data=All.temp.dat[All.temp.dat$Year ==YEAR[i],],xlim=lim,ylim=lim )
	abline(0,1,lty=2,col=2,lwd=3)
	title(YEAR[i])
}

## write plots to file
	for(j in 1:length(YEAR)){
		print(P[j])
	}
dev.off()





