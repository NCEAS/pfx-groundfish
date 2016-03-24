
library(knitr)
library(rgdal)
library(ggplot2)
library(RColorBrewer)

proj.dir	<-	"/Users/ole.shelton/Documents/GitHub/pfx-groundfish/goaTrawl"#getwd()# 
data.type 	<-	"abundance" # options are occurrence or abundance
##########################################################################################
## Write functions for calculating various diversity metrics of interest
richness	<-	function(VAL,cutoff){
					VAL[VAL >= cutoff]	<- 1
					VAL[VAL < cutoff]	<- 0					
					rich	<- sum(VAL)
					return(rich)
				}
shannon		<- function(X){
					sum.X	<-	sum(X)
					prop	<-  X/sum.X
					H		<-	- sum(prop * log(prop))
					return(H)
				}
gini.simp	<-	function(X){
					sum.X	<-	sum(X)
					prop	<-  X/sum.X
					D		<-	1- sum(prop^2)
				}
total.biomass	<-	function(X){
					tot			<-	sum(X)
					return(tot)
				}	
##########################################################################################
# Get needed data:
##########################################################################################

# Read in raw trawl data

df = read.csv(paste(proj.dir,"/Output Data/goa_trawl_final_albers+temp.csv",sep=""))

temp	<-	df[df$BottomDepth<150,]

# Read in Projections that are appropriate
	if(data.type == "abundance"){	
		# Choose predictions file
		dat 	<-	read.csv("/Users/ole.shelton/Dropbox/INLA output/Projections/All_species_uncond_pred_goa_shallow.csv")
		dat2	<-  	read.csv("/Users/ole.shelton/Dropbox/INLA output/Projections/All_species_uncond_pred_goa_mid.csv")
		dat3	<-  	read.csv("/Users/ole.shelton/Dropbox/INLA output/Projections/All_species_uncond_pred_goa_deep.csv")
		# Choose points to standardize across
		dat.areas	<- read.csv(paste(proj.dir,"/Output Data/goa_discrete_areas_for_comparison(50_to_150m).csv",sep=""))
	}	
	if(data.type == "occurrence"){	
		# Choose predictions file
		dat <-read.csv("/Users/ole.shelton/Documents/Science/Active projects/Exxon/Groundfish/Projections/All_species_uncond_pred_goa_shallow.csv")
		# Choose points to standardize across
		dat.areas	= read.csv(paste(proj.dir,"/Output Data/goa_discrete_areas_for_comparison(50_to_150m).csv",sep=""))
	}	
	#Summarize the number of locations in each area
		area.summary	<-	aggregate(dat.areas$Area,by=list(Area=dat.areas$Area),length)

	dat	<-	rbind(dat,dat2)
	dat	<-	rbind(dat,dat3)

# Read in Pilot files which determine which species to include or how to group species
	dat.control	= read.csv(paste(proj.dir,"/Output Data/trawl_species_control_file.csv",sep=""))

# Read in Predicted locations for all Gulf of Alaska
	dat.project	= read.csv(paste(proj.dir,"/Output Data/goa_projection_points+temp.csv",sep=""))
    dat.project$LonUTMAlbers = dat.project$LonUTMAlbers/1000
    dat.project$LatUTMAlbers = dat.project$LatUTMAlbers/1000
    #### Exclude points that end up on land
    dat.project$NGDC24_M =	-dat.project$NGDC24_M	# depth in m
    dat.project$SRTM_M = -dat.project$SRTM_M	# depth in m
    dat.project = dat.project[dat.project$NGDC24 > 0,]
	### Exclude deep points
	dat.project	<-	dat.project[dat.project$NGDC24_M < 300,]


#### PLOT SPATIAL DISTRIBUTION
	# Import the datafile of Alaska shoreline for projection purposes
	setwd(paste(proj.dir,"/_Output plots Pos/_Alaska Shapefile",sep=""))
	shp.alaska	 <-	readOGR(dsn=".",layer="Alaska-Albers")
	dat.alaska   <- fortify(shp.alaska,"data.frame")
	dat.alaska$lat	<-	dat.alaska$lat/1000
	dat.alaska$long	<-	dat.alaska$long/1000

#######################################################################################
#### Calculate the Community Metric of interest
#######################################################################################
	# Limit prediction points to the focal areas
	
	# Declare the metric of interest:
		FUNCTION	<-	"total.biomass"
		BY			<-	c("All")
#		BY			<-	c("class","fish.invert","pelagic.benthic","All")
		these		<- c("database.name",BY)
		dat.control.trim	<-	dat.control[,these]

	# Merge up the dataframes
	dat.fin	<-	merge(dat,dat.project[,c("MASTER_ID","LonUTMAlbers","LatUTMAlbers","NGDC24_M")])
	dat.fin$All	<-	"All"
	dat.fin$Var	<-	dat.fin$SE^2
 	
 	# Aggregate by Master_ID
	OUT	<-	list()
	for(Z in 1:length(BY)){
		out.all	<-	aggregate(dat.fin[,c("Mean","Median","Var")],
						by=list(MASTER_ID=dat.fin$MASTER_ID,Year=dat.fin$Year,dat.fin[,BY[Z]]),
						sum)
		colnames(out.all)[3]	<-	c(BY[Z])
		out.all$SD				<-	sqrt(out.all$Var)
		OUT[[Z]]<-out.all
	}
	names(OUT)	<-	BY

##########################################################################################
##########################################################################################
# Make some plots of the spatial patterns
##########################################################################################
##########################################################################################
	y.lim	<-	c(5500,6300)
	
	# Some Arguments for making prettier plots
	bGrid <-theme(panel.grid =element_blank())	
	bBack <-theme(panel.background =element_blank())
	bAxis <-theme(axis.title.y =element_blank())
	bTics <-theme(axis.text =element_blank(), axis.text.y =element_blank(), axis.ticks =element_blank())

	species	<-	sort(unique(dat$Species))
	years 	<- 	sort(unique(dat$Year))

# Loop over Metric of Interest
for(i in 1:length(BY)){
	temp	<-	OUT[[i]]
	temp	<-	merge(dat.project[,c(1,4,5)],temp)
	max.quant	<-	0.95
	min.quant	<-	0.01
	MAX.Z	<-	quantile(temp$Median,max.quant)
	if(MAX.Z < 10){MAX.Z <- 10}
	MIN.Z	<-	quantile(temp$Median,min.quant)
	z.lim	<-	c(0, MAX.Z)
	temp$Median[temp$Median > MAX.Z]	<- MAX.Z
	temp$Median[temp$Median < MIN.Z]	<- MIN.Z

	Breaks	<- round(10^(seq(log(MIN.Z+0.01,10),log(MAX.Z-1,10),length.out=5)),1)
	if(MAX.Z ==10){	Breaks	<- round(10^(seq(log(MIN.Z,10),log(MAX.Z-0.01,10),length.out=5)),3)}

	p	<-	list()
	# Loop over years
	for(j in 1:length(years)){
		temp.plot	<-	temp[temp$Year == years[j],]
	p[[j]]	<-	ggplot() +
 		scale_colour_gradientn(name="CPUE", colours=brewer.pal(9,"OrRd"),trans = "log",limits=c(MIN.Z,MAX.Z),breaks=Breaks)+
    	geom_point(data=temp.plot,alpha=0.4,
    			aes(LonUTMAlbers,LatUTMAlbers,colour=Median))+
		geom_polygon(data=dat.alaska, fill=grey(0.4),color=NA,aes(long,lat,group=group)) +
		labs(x = "Eastings",y="Northings",title=paste(FUNCTION,years[j],"Uncond Expectation; max =",max.quant,"of max (or 10)")) +
   		coord_cartesian(xlim = c(min(temp.plot$LonUTMAlbers),max(temp.plot$LonUTMAlbers)), 
   				ylim = y.lim)+
   		bGrid  #+ bBack
	}


NAME	<-	paste(BY[Z],"uncond CPUE all years.pdf")
setwd("/Users/ole.shelton/Dropbox/INLA output/Plots")
pdf(NAME,onefile=TRUE,width=15,height=5)
	## write plots to file
	for(j in 1:length(years)){
		print(p[j])
	}
dev.off()
}




























for(i in 1:length(BY)){ # loop through the list
	# Start with time series for each Area
	par(mfrow=c(4,3),mar=c(3,3,1.5,1))
	Levels <- unique(temp[,BY[i]])
	for(j in 1:length(AREAS)){ # Loop through the Areas
		
		temp	<-	INDEX[[i]][INDEX[[i]]$Area== j,]	
		y.lim	<-	c(0,max(INDEX[[i]]$Mean.50.))

		for(k in 1:length(Levels)){	
			plot(Mean.50.~Year,data=temp[temp[,1]==Levels[k],],col=k,type="b",
					ylim=y.lim,axes=F,pch=21,lwd=1.5)
			par(new=T)
		}
		par(new=F)
		axis(1);axis(2,las=2)
		box(bty="o")
		mtext(paste("Area", AREAS[j]),side=3,adj=1)
	}
	# add legend in the 12th box
	plot(1:5,1:5,type="n",axes=F,)
	
	legend(x=2,y=4.5,pch=21,col=1:length(Levels),lwd=2,legend=Levels)




















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
setwd("/Users/ole.shelton/GitHub/exxonValdez_nceas/goaTrawl/_Output plots Pos/_Plots species posteriors")
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
# dat<-read.csv("/Users/ole.shelton/Documents/Science/Active projects/Exxon/Groundfish/Projections/All_species_pres_pred_goa_shallow.csv")
# dat2<-read.csv("/Users/ole.shelton/Documents/Science/Active projects/Exxon/Groundfish/Projections/All_species_pres_pred_goa_mid.csv")
# dat3<-read.csv("/Users/ole.shelton/Documents/Science/Active projects/Exxon/Groundfish/Projections/All_species_pres_pred_goa_deep.csv")
# dat <- rbind(dat,dat2)
# dat <- rbind(dat,dat3)
	
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
setwd("/Users/ole.shelton/GitHub/exxonValdez_nceas/goaTrawl/_Output plots Pres Abs/_Plots species posteriors")
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