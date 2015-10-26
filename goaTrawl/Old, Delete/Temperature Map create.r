rm(list=ls())

# Using GAMs to make temperature surfaces for the Gulf of Alaska
library(mgcv)
library(ggplot2)
library(rgdal)
library(scales)

# Read in the data
data.dir	<-	"/Users/ole.shelton/GitHub/exxonValdez_nceas/goaTrawl/"
plot.dir	<-	"/Users/ole.shelton/GitHub/exxonValdez_nceas/goaTrawl/_Temperature Surfaces/"

# Go get the data
setwd(data.dir)
trawl.dat	<-	read.csv(file="goa_trawl_final_albers.csv",header=T)
project.dat	<-	read.csv(file="goa_projection_points.csv",header=T)

### GET RID OF NAs
trawl.dat2	<-	trawl.dat[is.na(trawl.dat$BottomDepth)==F,]
trawl.trim	<-	trawl.dat[is.na(trawl.dat$BottomTemp)==F & 
							is.na(trawl.dat$BottomDepth)==F &
							trawl.dat$BottomTemp != -9999
							,1:14]
	
project.dat	<- project.dat[project.dat$NGDC24_M< -25,]
project.dat$log.BottomDepth	<-	log(-project.dat$NGDC24_M)
	
#### Fit a spatial GAM to the observed data for each year
YEAR	<- unique(trawl.trim$Year)
all.dat	<- NULL

for(i in 1:length(YEAR)){
	dat.all.obs	<-	trawl.dat2[trawl.dat2$Year == YEAR[i],]
	dat.all.obs$log.BottomDepth	<- log(dat.all.obs$BottomDepth)

	dat	<-	trawl.trim[trawl.trim$Year == YEAR[i],]
	dat$log.BottomDepth	<- log(dat$BottomDepth)
	
 	out		<-	gam(BottomTemp ~ te(LonUTMAlbers,LatUTMAlbers,k=7)+s(BottomDepth),data=dat)
 	out.2	<-	gam(BottomTemp ~ te(LonUTMAlbers,LatUTMAlbers,k=7)+s(log.BottomDepth),data=dat)

# 	summary(out)
# 	summary(out.2)
# 	
# 	plot(out,pages=1,residuals=TRUE)  ## show partial residuals
# 	plot(out.2,pages=1,seWithMean=TRUE) ## `with intercept' CIs
# 	## run some basic model checks, including checking
# 	## smoothing basis dimensions...
# 	gam.check(out.2)
# 
# 	plot(out.2,pages=1) 
# 	plot(out.2,pages=1,scheme=2) ## alternative visualization

	THESE	<-	c("LonUTMAlbers","LatUTMAlbers","log.BottomDepth")
	new.dat	<-	dat.all.obs[,THESE]
	new.points <- predict.gam(out.2,new.dat)
	
	dat.all.obs$Pred.bot.temp			<-	unlist(new.points)
	dat.all.obs$BottomTemp[dat.all.obs$BottomTemp == -9999]	<- NA
	dat.all.obs$Resid.bot.temp			<- dat.all.obs$BottomTemp - dat.all.obs$Pred.bot.temp

	all.dat	<- rbind(all.dat,dat.all.obs)
	
	### predict for all points in the projection set
	new.dat	<- project.dat[,c('LonUTMAlbers','LatUTMAlbers','log.BottomDepth')]
	project.dat[,paste("Bot.Temp.",YEAR[i],sep="")]		<- predict.gam(out.2,new.dat)
}

all.dat2<-all.dat[,c(1:14,72:74,15:71)]

# write.csv(all.dat2,file="goa_trawl_final_albers+temp.csv")
# write.csv(project.dat,file="goa_projection_points+temp.csv")

#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#### PLOT SPATIAL DISTRIBUTION
# Import the datafile of Alaska shoreline

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

project.plot	<-	project.dat[project.dat$NGDC24_M > -600 & project.dat$NGDC24_M < -25,]

z.lim	<-	c(2,12)

ZZZ	<-	list()

NAME	<-	paste("Bot.Temp.",YEAR,sep="")

for(j in 1:length(YEAR)){
	p1 <-	ggplot() +
		scale_size(range = c(1,1))+
  		scale_colour_gradientn(limits=z.lim,colours= c(muted("blue"),"green","yellow","orange","red"))+
    	geom_point(data=project.plot,alpha=0.3,
    			mapping=aes_string("LonUTMAlbers","LatUTMAlbers",colour=	NAME[j])) + 
		geom_polygon(data=dat.alaska, fill=grey(0.4),color=NA,aes(long,lat,group=group)) +
		labs(x = "Eastings",y="Northings",title= paste(YEAR[j],"Bottom Temperature (C)"))+
   		coord_cartesian(xlim = c(min(project.dat$LonUTMAlbers),max(project.dat$LonUTMAlbers)), ylim = c(min(project.dat$LatUTMAlbers),6300000))+
 		bGrid  + bBack  		

	ZZZ[[j]]	<- p1
}


setwd(plot.dir)
pdf("Temperature Plots.pdf",onefile=TRUE,width=15,5)
## write plots to file
	for(j in 1:length(YEAR)){
		print(ZZZ[[j]])
	}
dev.off()








	