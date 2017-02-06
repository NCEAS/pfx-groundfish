##############################################################################################
###  Script for Catch and Temperature data                                                 ###
###  Originally written by Ole Shelton for PFX Groundfish paper #1                         ###
###       located pfx-groundfish/Scripts and plots for pubs/Catch and Temperature Script.r ###
###  Modified by Rachael E. Blake in January 2017 for                                      ###
###       groundfish paper #2 fishing pressure analysis                                    ###
##############################################################################################

# load libraries needed by this script
library(rgdal)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(viridis)

### Examine fish catches by stat6 area to include information 
catch <- read.csv("./goaTrawl/_fishing areas gfish/pounds_by_stat6.csv")
# this is the total catch of groundfish
# From Eric
# Yep -- I did it by permit, with only the "B" "C" and "M" permits. So it should be only longline and
# trawl, but would include a number of species as bycatch (including very small amounts of salmon, etc)
# Codes here:  #   https://www.cfec.state.ak.us/misc/FshyDesC.htm

### Read in mapping of stat areas to regions.
map_to_regions <- read.csv("./diversity-data/regions+OLE.csv")

catch1     <- merge(catch, map_to_regions[,c("stat6","final.OLE")])
catch_sum <- catch1 %>% 
             group_by(year,final.OLE) %>% 
             summarize(tot_pound =sum(pounds))
catch_sum$met_ton <-  catch_sum$tot_pound * 0.000453592

### Calculating the area (km2) of each fishery Region:
Area_est <- data.frame(
              matrix(c(
              "Alaska Peninsula",124367,
              "Cook Inlet",39443,
              "Kodiak",147333,
              "PWS",45136),
              4,2,byrow=T)
            )
colnames(Area_est) <- c("final.OLE","km2")
Area_est$km2 <- as.numeric(as.character(Area_est$km2))

catch_sum <- merge(catch_sum,Area_est)
catch_sum$met_ton_km2 <- catch_sum$met_ton / catch_sum$km2
catch_sum$final.OLE <- as.character(catch_sum$final.OLE)
catch_sum$final.OLE[catch_sum$final.OLE == "PWS"] <- "Prince William Sound"

###########################################################################
# Plotting time series of catch by region
COL <- viridis(4,begin=0,end=0.8)
biomass_plot <- ggplot(catch_sum) +
                  geom_point(aes(x=year, y=met_ton_km2, shape=final.OLE), size=2.5) +
                  geom_line(aes(x=year, y=met_ton_km2, group=final.OLE)) +
                  scale_shape(name="Region",solid=F) +
                  coord_cartesian(xlim=c(min(catch_sum$year)-0.75,max(catch_sum$year)+0.75),
                                  ylim=c(0,max(catch_sum$met_ton_km2)*1.05),
                                  expand=c(0)) +
               		labs(x="Year", y=expression("Catch (mt km"^-2*")")) +
                  theme_bw()
#biomass_plot


########################################################################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################
### Analysis for Temperature Time-Series.
########################################################################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################
# 
# library(dplyr)
# library(ncdf)
# library(reshape2)
# library(rgdal)
# library(mgcv)
# library(ggplot2)
# library(scales)
# 
# YEAR	<-	c("1990","1993","1996","1999","2001","2003","2005","2007","2009","2011","2013","2015")
# All.years	<- NULL
# 
# ########################################################################################
# #  Go get the Trawl Temperature Data 
# ########################################################################################
# ## This is a duplicate of the Temperature Map create code
# 
# # Read in the data
# data.dir	<-	"/Users/ole.shelton/GitHub/pfx-groundfish/goaTrawl/"
# plot.dir	<-	"/Users/ole.shelton/GitHub/pfx-groundfish/goaTrawl/_Temperature Surfaces/"
# 
# # Go get the data
# setwd(paste(data.dir,"Output Data/",sep=""))
# trawl.dat	<-	read.csv(file="goa_trawl_final_albers+temp.csv",header=T)
# trawl.dat$LonUTMAlbers <- trawl.dat$LonUTMAlbers/1000
# trawl.dat$LatUTMAlbers <- trawl.dat$LatUTMAlbers/1000
# 
# 
# # Go get the projection points
# setwd(paste(data.dir,"Output Data/",sep=""))
# projection.dat <- read.csv(file="goa_discrete_areas_for_comparison(50_to_150m).csv",header=T)
# projection.dat$log.BottomDepth <- log(-projection.dat$NGDC24_M)
# 
# ### GET RID OF NAs
# trawl.dat2	<-	trawl.dat[is.na(trawl.dat$BottomDepth)==F,]
# trawl.trim	<-	trawl.dat[is.na(trawl.dat$BottomTemp)==F & 
#                           is.na(trawl.dat$BottomDepth)==F &
#                           trawl.dat$BottomTemp != -9999
#                         ,1:14]
# 
# #### Fit a spatial GAM to the observed data for each year
# all.dat	<- NULL
# 
# for(i in 1:length(YEAR)){
#   dat	<-	trawl.trim[trawl.trim$Year == YEAR[i],]
#   dat$log.BottomDepth	<- log(dat$BottomDepth)
#   
#   out		<-	gam(BottomTemp ~ te(LonUTMAlbers,LatUTMAlbers,k=7)+s(BottomDepth),data=dat)
#   out.2	<-	gam(BottomTemp ~ te(LonUTMAlbers,LatUTMAlbers,k=7)+s(log.BottomDepth),data=dat)
#   
#   THESE	<-	c("LonUTMAlbers","LatUTMAlbers","log.BottomDepth")
#   
#   # Resid
#   new.dat	<-	dat[,THESE]
#   new.points <- predict.gam(out.2,new.dat)
#   dat$Pred.bot.temp			<-	unlist(new.points)
#   dat$Resid.bot.temp		<- dat$BottomTemp - dat$Pred.bot.temp
#   
#   # Projection
#   new.dat	<-	projection.dat[,THESE]
#   new.points <- predict.gam(out.2,new.dat)
#   
#   projection.dat$Pred.bot.temp			<-	unlist(new.points)
#   projection.dat$year <- YEAR[i]
#   all.dat	<- rbind(all.dat, projection.dat)
# }
# 
# 
# 
# ### Plot in Space
# 
# setwd(paste(plot.dir,"_Alaska Shapefile",sep=""))
# shp.alaska	 <-	readOGR(dsn=".",layer="Alaska-Albers")
# dat.alaska   <- fortify(shp.alaska,"data.frame")
# dat.alaska$long <- dat.alaska$long /1000 
# dat.alaska$lat  <- dat.alaska$lat /1000 
# 
# # Some Arguments for making prettier plots
# bGrid <-theme(panel.grid =element_blank())
# bBack <-theme(panel.background =element_blank())
# bAxis <-theme(axis.title.y =element_blank())
# bTics <-theme(axis.text =element_blank(), axis.text.y =element_blank(), axis.ticks =element_blank())
# 
# ### ADD trawl locations
# 
# ## Get rid of deep areas
# z.lim	<-	c(min(all.dat$Pred.bot.temp),max(all.dat$Pred.bot.temp))
# #z.lim[1]	<- z.lim[1] * -1
# #z.lim=c(-3,3)
# 
# P	<-	list()
# for(j in 1:length(YEAR)){
#   p1	<-	ggplot() +
#     scale_size(range = c(1,1))+
#     scale_colour_gradientn(limits=z.lim,colours= c("blue","red"))+
#     geom_point(data= all.dat[all.dat$year == YEAR[j],] ,alpha=0.3,
#                mapping=aes(LonUTMAlbers,LatUTMAlbers,colour=	Pred.bot.temp )) + 
#     geom_polygon(data=dat.alaska, fill=grey(0.4),color=NA,aes(long,lat,group=group)) +
#     labs(x = "Eastings",y="Northings",title=paste(YEAR[j],"Difference Bottom Temperature (C)"))+
#     coord_cartesian(xlim = c(min(all.dat$LonUTMAlbers),max(all.dat$LonUTMAlbers)), ylim = c(min(all.dat$LatUTMAlbers),max(all.dat$LatUTMAlbers)))+
#     bGrid  + bBack 
#   P[[j]]	<-	p1
# }
# 
# ####################################################################################
# ####################################################################################
# ####################################################################################
# # Calculate mean values to each area and calulate time=series for each area.
# 
# temp.yr    <- summarise(group_by(all.dat,Area,year),Mean=mean(Pred.bot.temp),SD = sd(Pred.bot.temp))
# grand.mean <- summarise(group_by(temp.yr,Area),grand.mean=mean(Mean))
# 
# temp.dat <- merge(temp.yr,grand.mean)
# temp.dat$cent.temp <- temp.dat$Mean-temp.dat$grand.mean
# 
# #########################################################################################
# 
# time_plot <- function(dat,NAME){
#   AREA <- as.numeric(as.character(unique(dat$Area)))
#   YEAR <- sort(unique(dat$year))
#   y.lim <- c(-1.5,1.5)
#   x.lim <- c(1985,2015)
#   par(mar=c(4.5,4,0.5,1))
#   for(i in 1:length(AREA)){
#     if(AREA[i] !=1){par(new=T)}
#     if(AREA[i] >= 3 & AREA[i] <= 5 ){
#       plot(cent.temp~year,data=dat[dat$Area ==AREA[i],],
#            type="l",axes=F,xlab="",ylim=y.lim,col=2,lwd=1.25,ylab="",xlim=x.lim)
#     }
#     if(AREA[i] < 3 | AREA[i] > 5){
#       plot(cent.temp~year,data=dat[dat$Area ==AREA[i],],
#            type="l",axes=F,xlab="",ylim=y.lim,col=1,lwd=1,ylab="",xlim=x.lim)
#     }
#   }
#   abline(h=0,lty=2)
#   axis(1,at=seq(1985,2015,by=5),las=1,padj=-1,tcl=-0.25,cex.axis=0.8)
#   axis(2,las=2,hadj=0.5,tcl=-0.25,cex.axis=0.8)
#   box(bty="o",lwd=2)
#   title(ylab=NAME,line=2.5)
#   title(xlab="Year",line=2)
#   #abline(v=1989,lty=2,lwd=2)
#   #arrows(x0=1989,x1=1989,y0=y.lim[1],y1=y.lim[1]+c(y.lim[2]-y.lim[1])*0.8,length=0,lty=2,lwd=2)
# }
# 
# ########
# ########################### MAKE PLOT FOR PUB HERE.
# 
# 
# setwd("/Users/ole.shelton/GitHub/pfx-groundfish/goaTrawl/_fishing area gfish")
# 
# quartz(file = "/Users/ole.shelton/GitHub/pfx-groundfish/goaTrawl/_fishing areas gfish/Plots/Catch + Temp by region.pdf",
#                 type="pdf",width=6,height=8,dpi=600)
#   
#   par(mar=c(2,4,3,1),mfrow=c(2,1))
#   
#     NAME <-  unique(catch.sum$final.OLE) 
#     PCH  <-  c(21,24,22,3)
#     catch.sum <- catch.sum[order(catch.sum$year),]
#     
#     y.lim=c(0,1.5)
#     for(i in 1:length(NAME)){
#       plot(met.ton.km2~year,data=catch.sum[catch.sum$final.OLE == NAME[i],],pch=PCH[i],ylim=y.lim,axes=F,xlab="",ylab="") 
#       par(new=T)        
#       plot(met.ton.km2~year,data=catch.sum[catch.sum$final.OLE == NAME[i],],type="l",ylim=y.lim,axes=F,xlab="",ylab="") 
#       par(new=T)
#     }         
#      axis(1,cex.axis=0.8,tcl=-0.25,padj=-1)
#      axis(2,las=2,cex.axis=0.8,tcl=-0.25,hadj=0.5)
#      box(bty="o",lwd=2)
#      title(ylab=expression("Catch (mt km"^-2*")"),line=2.5)
# 
#      legend(x=1985,y=1.5,pch=PCH,legend=NAME,cex=0.8)
# 
#      time_plot(temp.dat,NAME="Temperature (C)")
#          
# dev.off()         
#          
#        
#        
#        
#        
# 
# quartz(file="/Users/ole.shelton/GitHub/pfx-groundfish/goaTrawl/_Temperature Surfaces/Time-series of temperaturesRaw.pdf",type="pdf",dpi=300,height=4,width=5)
# par(mar=c(4,4,1,1))
# 
# dev.off()  
# 
# 
#        
#        
#        
#        
#        
#        
#        
#        
