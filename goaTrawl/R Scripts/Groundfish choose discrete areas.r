rm(list=ls())
library(INLA)
library(rgdal)
library(ggplot2)
library(splancs)
library(sp)

#### GO GET THE PROJECTION POINTS
proj.dir	<- "/Users/ole.shelton/GitHub/exxonValdez_nceas/goaTrawl/"
setwd(proj.dir)


#### GO GET THE PROJECTION POINTS
dat.project	= read.csv("Output Data/goa_projection_points+temp.csv")
dat.project$LonUTMAlbers = dat.project$LonUTMAlbers/1000
dat.project$LatUTMAlbers = dat.project$LatUTMAlbers/1000

dat.project$depth = -dat.project$NGDC24_M 
### Go Get the shapefile to provide a shoreline.
setwd(paste(proj.dir,"/Output plots/_Alaska Shapefile",sep=""))
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
MIN.D	<-	50
MAX.D	<-	150
x.lim	<-	c(min(dat.project$LonUTMAlbers),max(dat.project$LonUTMAlbers))
y.lim	<-	c(min(dat.project$LatUTMAlbers),max(dat.project$LatUTMAlbers))


##### ONLY PLOT SPECIFIED DEPTH RANGE
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
MIN.D	<-	50
MAX.D	<-	500

trim.x	<- c(-850,335)
trim.y	<- c(0,6170)
x.lim	<-	c(min(dat.project$LonUTMAlbers),max(dat.project$LonUTMAlbers))
y.lim	<-	c(min(dat.project$LatUTMAlbers),max(dat.project$LatUTMAlbers))

dat.select	<-	dat.project[dat.project$depth >= MIN.D & dat.project$depth <= MAX.D &
					dat.project$LonUTMAlbers >= trim.x[1] &	dat.project$LonUTMAlbers <= trim.x[2] &
					dat.project$LatUTMAlbers >= trim.y[1] &	dat.project$LatUTMAlbers <= trim.y[2]
					,]

 dat.select	<- dat.select[dat.select$LatUTMAlbers <= 6040 |
 							dat.select$LonUTMAlbers >= -60 ,]

dat.select	<- dat.select[dat.select$LatUTMAlbers <= 6160 |
 							dat.select$LonUTMAlbers >= 190 ,]
dat.select	<- dat.select[dat.select$LatUTMAlbers <= 6130 |
 							dat.select$LonUTMAlbers >= 150 ,]

dat.select	<-	dat.select[order(dat.select$LatUTMAlbers,dat.select$LonUTMAlbers,decreasing=T),]
dat.region	<- dat.select

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

########## WRITE DATA TO FILE
dat.fin			<-	dat.region[,1:8]
dat.fin$Area	<-	1
dat.fin$depth = -dat.fin$NGDC24_M 

dat.shallow	<-	dat.fin[dat.fin$depth >= 50  & dat.fin$depth <= 150,] 
dat.mid	<-	dat.fin[dat.fin$depth >  150 & dat.fin$depth <= 250,] 
dat.deep	<-	dat.fin[dat.fin$depth >  250,] 

setwd(paste(proj.dir,"Output Data",sep=""))
write.csv(dat.shallow,file=paste("goa_central_gulf(",50,"_to_",150,"m).csv",sep=""),row.names=F)
write.csv(dat.mid,file=paste("goa_central_gulf(",150,"_to_",250,"m).csv",sep=""),row.names=F)
write.csv(dat.deep,file=paste("goa_central_gulf(",250,"_to_deep).csv",sep=""),row.names=F)

### Make four files for particular depth breaks in the central gulf.
dat.fin			<-	dat.region[,1:8]
dat.fin$Area	<-	1

setwd(paste(proj.dir,"Output Data",sep=""))
write.csv(dat.fin,file=paste("goa_central_gulf(",MIN.D,"_to_",MAX.D,"m).csv",sep=""),row.names=F)


#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
###### Below makes the file "goa_discrete_areas_for_comparison(50_to_150m).csv"
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################

# Region 1 (east most)
trim.x	<- c(160,276)
trim.y	<- c(6000,6167)
x.lim	<-	c(160,280)
y.lim	<-	c(6000,6200)


dat.select	<-	dat.project[dat.project$depth >= MIN.D & dat.project$depth <= MAX.D &
					dat.project$LonUTMAlbers >= trim.x[1] &	dat.project$LonUTMAlbers <= trim.x[2] &
					dat.project$LatUTMAlbers >= trim.y[1] &	dat.project$LatUTMAlbers <= trim.y[2]
					,]
dat.select	<-	dat.select[order(dat.select$LatUTMAlbers,dat.select$LonUTMAlbers,decreasing=T),]

dat.select	<- dat.select[dat.select$MASTER_ID != 111753 &
							dat.select$MASTER_ID != 111754 &
							dat.select$MASTER_ID != 111755 &
							dat.select$MASTER_ID != 111532 &
							dat.select$MASTER_ID != 111533 &
							dat.select$MASTER_ID != 111534 &
							dat.select$MASTER_ID != 111541 &
							dat.select$MASTER_ID != 110341 &
							dat.select$MASTER_ID != 110069 &
							dat.select$MASTER_ID != 110601 &
							dat.select$MASTER_ID != 110602 & 
							dat.select$MASTER_ID != 109798 ,]

dat.region	<- dat.select
dat.region$Area	<- 1
# 

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

##########################################################################################
# Region 2 (Right at the sw corner of Montague Is.)
trim.x	<- c(100,170)
trim.y	<- c(6030,6111)
x.lim	<-	c(0,300)
y.lim	<-	c(6000,6200)


dat.select	<-	dat.project[dat.project$depth >= MIN.D & dat.project$depth <= MAX.D &
					dat.project$LonUTMAlbers >= trim.x[1] &	dat.project$LonUTMAlbers <= trim.x[2] &
					dat.project$LatUTMAlbers >= trim.y[1] &	dat.project$LatUTMAlbers <= trim.y[2]
					,]
dat.select	<-	dat.select[order(dat.select$LatUTMAlbers,dat.select$LonUTMAlbers,decreasing=F),]
# 
 dat.select	<- dat.select[dat.select$LatUTMAlbers >= 6050 |
 							dat.select$LonUTMAlbers >= 115 ,]

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


dat.select$Area <- 2
dat.region <- rbind(dat.region,dat.select)
##########################################################################################
# Region 3 (Rise WSW Montague Is.)
trim.x	<- c(63,120)
trim.y	<- c(6000,6090)
x.lim	<-	c(-100,300)
y.lim	<-	c(6000,6200)


dat.select	<-	dat.project[dat.project$depth >= MIN.D & dat.project$depth <= MAX.D &
					dat.project$LonUTMAlbers >= trim.x[1] &	dat.project$LonUTMAlbers <= trim.x[2] &
					dat.project$LatUTMAlbers >= trim.y[1] &	dat.project$LatUTMAlbers <= trim.y[2]
					,]

dat.select	<-	dat.select[order(dat.select$LatUTMAlbers,dat.select$LonUTMAlbers,decreasing=F),]
# 
  dat.select	<- dat.select[dat.select$LatUTMAlbers <= 6050 |
  							dat.select$LonUTMAlbers <= 100 ,]

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

dat.select$Area <- 3
dat.region <- rbind(dat.region,dat.select)

##########################################################################################
# Region 4 (Rise W or region 3 (due south of Kenai Penn))
trim.x	<- c(-100,57)
trim.y	<- c(6000,6083)
x.lim	<-	c(-100,300)
y.lim	<-	c(6000,6200)


dat.select	<-	dat.project[dat.project$depth >= MIN.D & dat.project$depth <= MAX.D &
					dat.project$LonUTMAlbers >= trim.x[1] &	dat.project$LonUTMAlbers <= trim.x[2] &
					dat.project$LatUTMAlbers >= trim.y[1] &	dat.project$LatUTMAlbers <= trim.y[2]
					,]

dat.select	<-	dat.select[order(dat.select$LatUTMAlbers,dat.select$LonUTMAlbers,decreasing=F),]
# 
  dat.select	<- dat.select[dat.select$LatUTMAlbers >= 6050 |
  							dat.select$LonUTMAlbers >= -10 ,]
  dat.select	<- dat.select[dat.select$LatUTMAlbers <= 6045 |
  							dat.select$LonUTMAlbers >= 5 ,]



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

dat.select$Area <- 4
dat.region <- rbind(dat.region,dat.select)

##########################################################################################
# Region 5 (NE corner of Kodiak, SW corner of Kenai Penn))
trim.x	<- c(-145,-95)
trim.y	<- c(5960,6030)
x.lim	<-	c(-250,0)
y.lim	<-	c(5800,6100)

dat.select	<-	dat.project[dat.project$depth >= MIN.D & dat.project$depth <= MAX.D &
					dat.project$LonUTMAlbers >= trim.x[1] &	dat.project$LonUTMAlbers <= trim.x[2] &
					dat.project$LatUTMAlbers >= trim.y[1] &	dat.project$LatUTMAlbers <= trim.y[2]
					,]

dat.select	<-	dat.select[order(dat.select$LatUTMAlbers,dat.select$LonUTMAlbers,decreasing=F),]
# 
#   dat.select	<- dat.select[dat.select$LatUTMAlbers >= 6050 |
#   							dat.select$LonUTMAlbers >= -10 ,]
#   dat.select	<- dat.select[dat.select$LatUTMAlbers <= 6045 |
#   							dat.select$LonUTMAlbers >= 5 ,]



	p3	<-	ggplot() +
		scale_size(range = c(1,1))+
  		scale_colour_gradientn(limits=z.lim,colours= c("red","black"))+
    	geom_point(data= dat.select,alpha=0.3,
    			mapping=aes(LonUTMAlbers,LatUTMAlbers,colour=depth )) + 
		geom_polygon(data=dat.alaska, fill=grey(0.4),color=NA,aes(long.km,lat.km,group=group)) +
		labs(x = "Eastings",y="Northings",title=paste("Depth"))+
   		coord_cartesian(xlim = x.lim, ylim = y.lim)+
 		bGrid  + bBack 
p3

dat.select$Area <- 5
dat.region <- rbind(dat.region,dat.select)

##########################################################################################
# Region 6 (E of Kodiak)
trim.x	<- c(-68,35)
trim.y	<- c(5920,5983)
x.lim	<-	c(-250,100)
y.lim	<-	c(5800,6100)

dat.select	<-	dat.project[dat.project$depth >= MIN.D & dat.project$depth <= MAX.D &
					dat.project$LonUTMAlbers >= trim.x[1] &	dat.project$LonUTMAlbers <= trim.x[2] &
					dat.project$LatUTMAlbers >= trim.y[1] &	dat.project$LatUTMAlbers <= trim.y[2]
					,]

dat.select	<-	dat.select[order(dat.select$LatUTMAlbers,dat.select$LonUTMAlbers,decreasing=F),]
# 
   dat.select	<- dat.select[dat.select$LatUTMAlbers <= 5961 |
   							dat.select$LonUTMAlbers <= 21 ,]
   dat.select	<- dat.select[dat.select$LatUTMAlbers <= 5951 |
   							dat.select$LonUTMAlbers <= 25 ,]
   dat.select	<- dat.select[dat.select$LatUTMAlbers <= 5970 |
   							dat.select$LonUTMAlbers <= 16 ,]
   dat.select	<- dat.select[dat.select$LatUTMAlbers <= 5970 |
   							dat.select$LonUTMAlbers <= 16 ,]
   dat.select	<- dat.select[dat.select$LatUTMAlbers >= 5929 |
   							dat.select$LonUTMAlbers >= -55 ,]
   dat.select	<- dat.select[dat.select$LatUTMAlbers >= 5940|
   							dat.select$LonUTMAlbers >= -62 ,]
   dat.select	<- dat.select[dat.select$LatUTMAlbers >= 5935|
   							dat.select$LonUTMAlbers >= -59 ,]

	p3	<-	ggplot() +
		scale_size(range = c(1,1))+
  		scale_colour_gradientn(limits=z.lim,colours= c("red","black"))+
    	geom_point(data= dat.select,alpha=0.3,
    			mapping=aes(LonUTMAlbers,LatUTMAlbers,colour=depth )) + 
		geom_polygon(data=dat.alaska, fill=grey(0.4),color=NA,aes(long.km,lat.km,group=group)) +
		labs(x = "Eastings",y="Northings",title=paste("Depth"))+
   		coord_cartesian(xlim = x.lim, ylim = y.lim)+
 		bGrid  + bBack 
p3

dat.select$Area <- 6
dat.region <- rbind(dat.region,dat.select)

##########################################################################################
# Region 7 (NE region of 3 that are SE of Kodiak)
trim.x	<-c(-120,0)
trim.y	<- c(5810,5930)
x.lim	<-	c(-250,100)
y.lim	<-	c(5800,6100)

dat.select	<-	dat.project[dat.project$depth >= MIN.D & dat.project$depth <= MAX.D &
					dat.project$LonUTMAlbers >= trim.x[1] &	dat.project$LonUTMAlbers <= trim.x[2] &
					dat.project$LatUTMAlbers >= trim.y[1] &	dat.project$LatUTMAlbers <= trim.y[2]
					,]

dat.select	<-	dat.select[order(dat.select$LatUTMAlbers,dat.select$LonUTMAlbers,decreasing=F),]

	dat.select	<- dat.select[dat.select$LatUTMAlbers <= 5922 |
   							dat.select$LonUTMAlbers <= -77 ,]
	dat.select	<- dat.select[dat.select$LatUTMAlbers <= 5905 |
   							dat.select$LonUTMAlbers <= -38 ,]
	dat.select	<- dat.select[dat.select$LatUTMAlbers <= 5910 |
   							dat.select$LonUTMAlbers <= -55 ,]
	dat.select	<- dat.select[dat.select$LatUTMAlbers >= 5835 |
    							dat.select$LonUTMAlbers >= -80 ,]
    dat.select	<- dat.select[dat.select$LatUTMAlbers <= 5905 |
    							dat.select$LonUTMAlbers >= -100 ,]
     dat.select	<- dat.select[dat.select$LatUTMAlbers <= 5915 |
      							dat.select$LonUTMAlbers >= -86 ,]
    dat.select	<- dat.select[dat.select$LatUTMAlbers >= 5850 |
    							dat.select$LonUTMAlbers >= -90 ,]
    dat.select	<- dat.select[dat.select$LatUTMAlbers >= 5845 |
    							dat.select$LonUTMAlbers >= -83 ,]
    dat.select	<- dat.select[dat.select$LatUTMAlbers >= 5867 |
    							dat.select$LonUTMAlbers >= -110 ,]
    dat.select	<- dat.select[dat.select$LatUTMAlbers >= 5855 |
    							dat.select$LonUTMAlbers >= -100 ,]

	p3	<-	ggplot() +
		scale_size(range = c(1,1))+
  		scale_colour_gradientn(limits=z.lim,colours= c("red","black"))+
    	geom_point(data= dat.select,alpha=0.3,
    			mapping=aes(LonUTMAlbers,LatUTMAlbers,colour=depth )) + 
		geom_polygon(data=dat.alaska, fill=grey(0.4),color=NA,aes(long.km,lat.km,group=group)) +
		labs(x = "Eastings",y="Northings",title=paste("Depth"))+
   		coord_cartesian(xlim = x.lim, ylim = y.lim)+
 		bGrid  + bBack 
p3

dat.select$Area <- 7
dat.region <- rbind(dat.region,dat.select)

##########################################################################################
# Region 8 (Cental region of 3 that are SE of Kodiak)
trim.x	<-c(-150,-82)
trim.y	<- c(5730,5865)
x.lim	<-	c(-250,100)
y.lim	<-	c(5700,6000)

dat.select	<-	dat.project[dat.project$depth >= MIN.D & dat.project$depth <= MAX.D &
					dat.project$LonUTMAlbers >= trim.x[1] &	dat.project$LonUTMAlbers <= trim.x[2] &
					dat.project$LatUTMAlbers >= trim.y[1] &	dat.project$LatUTMAlbers <= trim.y[2]
					,]

dat.select	<-	dat.select[order(dat.select$LatUTMAlbers,dat.select$LonUTMAlbers,decreasing=F),]

    dat.select	<- dat.select[dat.select$LatUTMAlbers <= 5850 |
    							dat.select$LonUTMAlbers <= -95 ,]
    dat.select	<- dat.select[dat.select$LatUTMAlbers <= 5840 |
    							dat.select$LonUTMAlbers <= -85 ,]
    dat.select	<- dat.select[dat.select$LatUTMAlbers <= 5867 |
    							dat.select$LonUTMAlbers <= -110 ,]
    dat.select	<- dat.select[dat.select$LatUTMAlbers <= 5855 |
    							dat.select$LonUTMAlbers <= -112 ,]

# 
	p3	<-	ggplot() +
		scale_size(range = c(1,1))+
  		scale_colour_gradientn(limits=z.lim,colours= c("red","black"))+
    	geom_point(data= dat.select,alpha=0.3,
    			mapping=aes(LonUTMAlbers,LatUTMAlbers,colour=depth )) + 
		geom_polygon(data=dat.alaska, fill=grey(0.4),color=NA,aes(long.km,lat.km,group=group)) +
		labs(x = "Eastings",y="Northings",title=paste("Depth"))+
   		coord_cartesian(xlim = x.lim, ylim = y.lim)+
 		bGrid  + bBack 
p3

dat.select$Area <- 8
dat.region <- rbind(dat.region,dat.select)

##########################################################################################
# Region 9 (SW region of 3 that are SE of Kodiak)
trim.x	<-	c(-197,-155)
trim.y	<- 	c(5700,5810)
x.lim	<-	c(-300,-100)
y.lim	<-	c(5700,5850)

# trim.x	<-	x.lim
# trim.y	<-	y.lim


dat.select	<-	dat.project[dat.project$depth >= MIN.D & dat.project$depth <= MAX.D &
					dat.project$LonUTMAlbers >= trim.x[1] &	dat.project$LonUTMAlbers <= trim.x[2] &
					dat.project$LatUTMAlbers >= trim.y[1] &	dat.project$LatUTMAlbers <= trim.y[2]
					,]

dat.select	<-	dat.select[order(dat.select$LatUTMAlbers,dat.select$LonUTMAlbers,decreasing=F),]

    dat.select	<- dat.select[dat.select$LatUTMAlbers <= 5795 |
    							dat.select$LonUTMAlbers <= -160 ,]
	p3	<-	ggplot() +
		scale_size(range = c(1,1))+
  		scale_colour_gradientn(limits=z.lim,colours= c("red","black"))+
    	geom_point(data= dat.select,alpha=0.3,
    			mapping=aes(LonUTMAlbers,LatUTMAlbers,colour=depth )) + 
		geom_polygon(data=dat.alaska, fill=grey(0.4),color=NA,aes(long.km,lat.km,group=group)) +
		labs(x = "Eastings",y="Northings",title=paste("Depth"))+
   		coord_cartesian(xlim = x.lim, ylim = y.lim)+
 		bGrid  + bBack 
p3

dat.select$Area <- 9
dat.region <- rbind(dat.region,dat.select)


##########################################################################################
# Region 10 (SW of Kodiak)
trim.x	<-	c(-400,-280)
trim.y	<- 	c(5600,5834)
x.lim	<-	c(-500,-100)
y.lim	<-	c(5600,5950)

# trim.x	<-	x.lim
# trim.y	<-	y.lim


dat.select	<-	dat.project[dat.project$depth >= MIN.D & dat.project$depth <= MAX.D &
					dat.project$LonUTMAlbers >= trim.x[1] &	dat.project$LonUTMAlbers <= trim.x[2] &
					dat.project$LatUTMAlbers >= trim.y[1] &	dat.project$LatUTMAlbers <= trim.y[2]
					,]

 dat.select	<-	dat.select[order(dat.select$LatUTMAlbers,dat.select$LonUTMAlbers,decreasing=F),]

    dat.select	<- dat.select[dat.select$LatUTMAlbers <= 5750 |
    							dat.select$LonUTMAlbers >= -350 ,]
    dat.select	<- dat.select[dat.select$LatUTMAlbers >= 5725 |
    							dat.select$LonUTMAlbers <= -320 ,]
    dat.select	<- dat.select[dat.select$LatUTMAlbers >= 5694 |
    							dat.select$LonUTMAlbers <= -356 ,]
    dat.select	<- dat.select[dat.select$LatUTMAlbers >= 5670 |
    							dat.select$LonUTMAlbers <= -364 ,]
    dat.select	<- dat.select[dat.select$LatUTMAlbers >= 5790 |
    							dat.select$LonUTMAlbers <= -289 ,]


	p3	<-	ggplot() +
		scale_size(range = c(1,1))+
  		scale_colour_gradientn(limits=z.lim,colours= c("red","black"))+
    	geom_point(data= dat.select,alpha=0.3,
    			mapping=aes(LonUTMAlbers,LatUTMAlbers,colour=depth )) + 
		geom_polygon(data=dat.alaska, fill=grey(0.4),color=NA,aes(long.km,lat.km,group=group)) +
		labs(x = "Eastings",y="Northings",title=paste("Depth"))+
   		coord_cartesian(xlim = x.lim, ylim = y.lim)+
 		bGrid  + bBack 
p3

dat.select$Area <- 10
dat.region <- rbind(dat.region,dat.select)

##########################################################################################
# Region 11 
trim.x	<-	c(-524,-407.5)
trim.y	<- 	c(5500,5745)
x.lim	<-	c(-600,-400)
y.lim	<-	c(5500,5800)

# trim.x	<-	x.lim
# trim.y	<-	y.lim


dat.select	<-	dat.project[dat.project$depth >= MIN.D & dat.project$depth <= MAX.D &
					dat.project$LonUTMAlbers >= trim.x[1] &	dat.project$LonUTMAlbers <= trim.x[2] &
					dat.project$LatUTMAlbers >= trim.y[1] &	dat.project$LatUTMAlbers <= trim.y[2]
					,]

dat.select	<-	dat.select[order(dat.select$LatUTMAlbers,dat.select$LonUTMAlbers,decreasing=F),]

    dat.select	<- dat.select[dat.select$LatUTMAlbers <= 5695 |
    							dat.select$LonUTMAlbers >= -508 ,]
    dat.select	<- dat.select[dat.select$LatUTMAlbers <= 5700 |
    							dat.select$LonUTMAlbers >= -505 ,]
    dat.select	<- dat.select[dat.select$LatUTMAlbers <= 5705 |
    							dat.select$LonUTMAlbers >= -495 ,]
    dat.select	<- dat.select[dat.select$LatUTMAlbers <= 5709 |
    							dat.select$LonUTMAlbers >= -490 ,]
    dat.select	<- dat.select[dat.select$LatUTMAlbers <= 5720 |
    							dat.select$LonUTMAlbers >= -465 ,]
    dat.select	<- dat.select[dat.select$LatUTMAlbers <= 5727 |
    							dat.select$LonUTMAlbers >= -452 ,]
    dat.select	<- dat.select[dat.select$LatUTMAlbers <= 5732 |
    							dat.select$LonUTMAlbers >= -445 ,]
    dat.select	<- dat.select[dat.select$LatUTMAlbers >= 5650 |
    							dat.select$LonUTMAlbers >= -470 ,]
    dat.select	<- dat.select[dat.select$LatUTMAlbers >= 5645 |
    							dat.select$LonUTMAlbers >= -465 ,]
    dat.select	<- dat.select[dat.select$LatUTMAlbers >= 5675 |
    							dat.select$LonUTMAlbers >= -474 ,]
    dat.select	<- dat.select[dat.select$LatUTMAlbers >= 5605 |
    							dat.select$LonUTMAlbers >= -460 ,]

	p3	<-	ggplot() +
		scale_size(range = c(1,1))+
  		scale_colour_gradientn(limits=z.lim,colours= c("red","black"))+
    	geom_point(data= dat.select,alpha=0.3,
    			mapping=aes(LonUTMAlbers,LatUTMAlbers,colour=depth )) + 
		geom_polygon(data=dat.alaska, fill=grey(0.4),color=NA,aes(long.km,lat.km,group=group)) +
		labs(x = "Eastings",y="Northings",title=paste("Depth"))+
   		coord_cartesian(xlim = x.lim, ylim = y.lim)+
 		bGrid  + bBack 
p3

dat.select$Area <- 11
dat.region <- rbind(dat.region,dat.select)
##########################################################################################
## Plot the 11 Regions
##########################################################################################
x.lim	<-	c(min(dat.project$LonUTMAlbers),max(dat.project$LonUTMAlbers))
y.lim	<-	c(min(dat.project$LatUTMAlbers),max(dat.project$LatUTMAlbers))

x.lim	<-	c(-600,400)
y.lim	<-	c(5500,6300)


COLS	<- brewer.pal(length(unique(dat.region$Area)),"Set3")
dat.region$Area	<- as.factor(dat.region$Area)

	p4	<-	ggplot() +
 		scale_size(range = c(1,1))+
  		scale_colour_manual(values= COLS)+
    	geom_point(data= dat.region,alpha=1,
    			mapping=aes(LonUTMAlbers,LatUTMAlbers,colour=Area )) + 
		geom_polygon(data=dat.alaska, fill=grey(0.4),color=NA,aes(long.km,lat.km,group=group)) +
		labs(x = "Eastings",y="Northings",title=paste("Areas with depths between",MIN.D,"and",MAX.D,"m") )+
   		coord_cartesian(xlim = x.lim, ylim = y.lim)+
 		bGrid  + bBack 
p4


## Write plots to file
setwd(paste(proj.dir,"Areas for comparison",sep=""))
pdf("Shallow Areas for comparison.pdf",onefile=TRUE,width=10,5)
	print(p4)
dev.off()

########## WRITE DATA TO FILE
dat.fin	<-	dat.region[,1:8]
dat.fin$Area	<-	dat.region$Area

setwd(paste(proj.dir,"Output Data",sep=""))
write.csv(dat.fin,file=paste("goa_discrete_areas_for_comparison(",MIN.D,"_to_",MAX.D,"m).csv",sep=""),row.names=F)


