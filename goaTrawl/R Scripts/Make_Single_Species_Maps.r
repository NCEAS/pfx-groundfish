library(rgdal)
library(ggplot2)
library(RColorBrewer)
setwd("/Users/ole.shelton/Documents/GitHub/exxonValdez_nceas/goaTrawl")
# Read in Predicted locations
	dat.project	= read.csv("Output Data/goa_projection_points+temp.csv")
    dat.project$LonUTMAlbers = dat.project$LonUTMAlbers/1000
    dat.project$LatUTMAlbers = dat.project$LatUTMAlbers/1000
    #### Exclude points that end up on land
    dat.project$NGDC24_M =	-dat.project$NGDC24_M	# depth in m
    dat.project$SRTM_M = -dat.project$SRTM_M	# depth in m
    dat.project = dat.project[dat.project$NGDC24 > 0,]

#### PLOT SPATIAL DISTRIBUTION
# Import the datafile of Alaska shoreline
setwd("/Users/ole.shelton/Documents/GitHub/exxonValdez_nceas/goaTrawl/_Output plots Pos/_Alaska Shapefile")
shp.alaska	 <-	readOGR(dsn=".",layer="Alaska-Albers")
#shp.alaska	 <-	readOGR(dsn=".",layer="/_Output plots Pos/_Alaska Shapefile/Alaska-Albers")
dat.alaska   <- fortify(shp.alaska,"data.frame")
dat.alaska$lat	<-	dat.alaska$lat/1000
dat.alaska$long	<-	dat.alaska$long/1000
#######################################################################################
#### PLOT Mean Field for each species in each year (look for nonsense)
#######################################################################################
### Read in Data
dat <-read.csv("/Users/ole.shelton/Dropbox/INLA output/All_species_uncond_pred_goa_shallow.csv")
dat2<-read.csv("/Users/ole.shelton/Dropbox/INLA output/All_species_uncond_pred_goa_mid.csv")
dat3<-read.csv("/Users/ole.shelton/Dropbox/INLA output/All_species_uncond_pred_goa_deep.csv")
 dat <- rbind(dat,dat2)
 dat <- rbind(dat,dat3)

	y.lim	<-	c(5500,6300)
	
	# Some Arguments for making prettier plots
	bGrid <-theme(panel.grid =element_blank())	
	bBack <-theme(panel.background =element_blank())
	bAxis <-theme(axis.title.y =element_blank())
	bTics <-theme(axis.text =element_blank(), axis.text.y =element_blank(), axis.ticks =element_blank())

	species	<-	sort(unique(dat$Species))
	years 	<- 	sort(unique(dat$Year))

# Loop over species
for(i in 1:length(species)){
	temp	<-	dat[dat$Species == species[i],]
	temp	<-	merge(dat.project[,c(1,4,5)],temp)
	max.quant	<-	0.975
	MAX.Z	<-	quantile(temp$Median,max.quant)
	if(MAX.Z < 10){MAX.Z <- 10}
	MIN.Z	<-	1e-3
	z.lim	<-	c(0, MAX.Z)
	temp$Median[temp$Median > MAX.Z]		<- MAX.Z
	temp$Median[temp$Median < MIN.Z]	<- MIN.Z

	Breaks	<- round(10^(seq(log(MIN.Z,10),log(MAX.Z-1,10),length.out=5)),3)
	if(MAX.Z ==10){	Breaks	<- round(10^(seq(log(MIN.Z,10),log(MAX.Z-0.01,10),length.out=5)),3)}

	p	<-	list()
	# Loop over years
	for(j in 1:length(years)){
		temp.plot	<-	temp[temp$Year == years[j],]
	p[[j]]	<-	ggplot() +
 		scale_colour_gradientn(name="CPUE", colours=brewer.pal(9,"OrRd"),trans = "log",limits=c(MIN.Z,MAX.Z),breaks=Breaks)+
    	geom_point(data=temp.plot,alpha=0.4,
    			aes(LonUTMAlbers,LatUTMAlbers,colour=Median))+
#     	geom_point(data=plot.zeros[plot.zeros$Year==YEARS[j],],alpha=0.4,shape="+",colour="black",
#     			aes(LonUTMAlbers*1000,LatUTMAlbers*1000))+
		geom_polygon(data=dat.alaska, fill=grey(0.4),color=NA,aes(long,lat,group=group)) +
		labs(x = "Eastings",y="Northings",title=paste(species[i],years[j],"Uncond Expectation; max =",max.quant,"of max (or 10)")) +
   		coord_cartesian(xlim = c(min(temp.plot$LonUTMAlbers),max(temp.plot$LonUTMAlbers)), 
   				ylim = y.lim)+
   		bGrid  #+ bBack
	}


NAME	<-	paste(species[i],"uncond CPUE all years.pdf")
setwd("/Users/ole.shelton/Dropbox/INLA output/Plots")
pdf(NAME,onefile=TRUE,width=15,height=5)
	## write plots to file
	for(j in 1:length(years)){
		print(p[j])
	}
dev.off()
}

#########################################
# Yearly differences in abundance
#########################################


# Loop over species
# for(i in 1:length(species)){
# 	temp	<-	dat[dat$Species == species[i],]
# 	temp	<-	merge(dat.project[,c(1,4,5)],temp)
# 	max.quant	<-	0.975
# 	MAX.Z	<-	quantile(temp$Median,max.quant)
# 	if(MAX.Z < 10){MAX.Z <- 10}
# 	MIN.Z	<-	1e-3
# 	z.lim	<-	c(0, MAX.Z)
# 	temp$Median[temp$Median > MAX.Z]		<- MAX.Z
# 	temp$Median[temp$Median < MIN.Z]	<- MIN.Z
# 
# 	Breaks	<- round(10^(seq(log(MIN.Z,10),log(MAX.Z-1,10),length.out=5)),3)
# 	if(MAX.Z ==10){	Breaks	<- round(10^(seq(log(MIN.Z,10),log(MAX.Z-0.01,10),length.out=5)),3)}
# 
# 	p	<-	list()
# 	# Loop over years
# 	for(j in 1:length(years)){
# 		temp.plot	<-	temp[temp$Year == years[j],]
# 	p[[j]]	<-	ggplot() +
#  		scale_colour_gradientn(name="CPUE", colours=brewer.pal(9,"OrRd"),trans = "log",limits=c(MIN.Z,MAX.Z),breaks=Breaks)+
#     	geom_point(data=temp.plot,alpha=0.4,
#     			aes(LonUTMAlbers,LatUTMAlbers,colour=Median))+
# #     	geom_point(data=plot.zeros[plot.zeros$Year==YEARS[j],],alpha=0.4,shape="+",colour="black",
# #     			aes(LonUTMAlbers*1000,LatUTMAlbers*1000))+
# 		geom_polygon(data=dat.alaska, fill=grey(0.4),color=NA,aes(long,lat,group=group)) +
# 		labs(x = "Eastings",y="Northings",title=paste(species[i],years[j],"Uncond Expectation; max =",max.quant,"of max (or 10)")) +
#    		coord_cartesian(xlim = c(min(temp.plot$LonUTMAlbers),max(temp.plot$LonUTMAlbers)), 
#    				ylim = y.lim)+
#    		bGrid  #+ bBack
# 	}
# 
# 
# NAME	<-	paste(species[i],"uncond CPUE DIFFERENCE all years.pdf")
# setwd("/Users/ole.shelton/GitHub/exxonValdez_nceas/goaTrawl/_Output plots Pos/_Plots species posteriors")
# pdf(NAME,onefile=TRUE,width=15,height=5)
# 	## write plots to file
# 	for(j in 1:length(years)){
# 		print(p[j])
# 	}
# dev.off()
# }
# 


##########################################################################################
##########################################################################################
#### PLOT PRESENCE-ABSENCE for all species, all years.
##########################################################################################
##########################################################################################
dat<-read.csv("/Users/ole.shelton/Documents/Science/Active projects/Exxon/Groundfish/Projections/All_species_pres_pred_goa_shallow.csv")
dat2<-read.csv("/Users/ole.shelton/Documents/Science/Active projects/Exxon/Groundfish/Projections/All_species_pres_pred_goa_mid.csv")
dat3<-read.csv("/Users/ole.shelton/Documents/Science/Active projects/Exxon/Groundfish/Projections/All_species_pres_pred_goa_deep.csv")
dat <- rbind(dat,dat2)
dat <- rbind(dat,dat3)
	
	y.lim	<-	c(5500,6300)
	
	# Some Arguments for making prettier plots
	bGrid <-theme(panel.grid =element_blank())	
	bBack <-theme(panel.background =element_blank())
	bAxis <-theme(axis.title.y =element_blank())
	bTics <-theme(axis.text =element_blank(), axis.text.y =element_blank(), axis.ticks =element_blank())

	species	<-	sort(unique(dat$Species))
	years 	<- 	sort(unique(dat$Year))

# Loop over species
for(i in 1:length(species)){
	temp	<-	dat[dat$Species == species[i],]
	temp	<-	merge(dat.project[,c(1,4,5)],temp)
# 	MAX.Z	<-	quantile(temp$Median,0.975)
# 	MIN.Z	<-	1e-3
	z.lim	<-	c(0, 1)
	Breaks	<- c(0,0.25,0.5,0.75,1.0)
	
	p	<-	list()
	# Loop over years
	for(j in 1:length(years)){
		temp.plot	<-	temp[temp$Year == years[j],]
	p[[j]]	<-	ggplot() +
 		scale_colour_gradientn(name="Occurrence", colours=brewer.pal(9,"OrRd"),limits=z.lim,breaks=Breaks)+
    	geom_point(data=temp.plot,alpha=0.4,
    			aes(LonUTMAlbers,LatUTMAlbers,colour=Median))+
#     	geom_point(data=plot.zeros[plot.zeros$Year==YEARS[j],],alpha=0.4,shape="+",colour="black",
#     			aes(LonUTMAlbers*1000,LatUTMAlbers*1000))+
		geom_polygon(data=dat.alaska, fill=grey(0.4),color=NA,aes(long,lat,group=group)) +
		labs(x = "Eastings",y="Northings",title=paste(species[i],years[j],"Occurrence")) +
   		coord_cartesian(xlim = c(min(temp.plot$LonUTMAlbers),max(temp.plot$LonUTMAlbers)), 
   				ylim = y.lim)+
   		bGrid  #+ bBack


	}

NAME	<-	paste(species[i],"Pres-Abs all years.pdf")
setwd("/Users/ole.shelton/Dropbox/INLA output/Plots")
pdf(NAME,onefile=TRUE,width=15,height=5)
	## write plots to file
	for(j in 1:length(years)){
		print(p[j])
	}
dev.off()

}




























### Calculate some simple diversity metrics
H	 	<-	aggregate(dat$Median, by=list(MASTER_ID=dat$MASTER_ID,Year=dat$Year),shannon)
D	 	<-	aggregate(dat$Median, by=list(MASTER_ID=dat$MASTER_ID,Year=dat$Year),gini.simp)
BIOMASS	<-	aggregate(dat$Median, by=list(MASTER_ID=dat$MASTER_ID,Year=dat$Year),total.biomass)

colnames(H)[3]			<-	"shannon"
colnames(D)[3]			<-	"gini.simp"
colnames(BIOMASS)[3]	<-	"biomass"






A <- dat[dat$MASTER_ID==THIS,]
B	<-	A[A$Year==2001,]

H	 <-	aggregate(A$Median, by=list(MASTER_ID=A$MASTER_ID,Year=A$Year),shannon)


H	 <-	aggregate(B$Median, by=list(MASTER_ID=B$MASTER_ID,Year=B$Year),shannon)

C	<- B$Median / sum(B$Median)
D	<-	log(C)