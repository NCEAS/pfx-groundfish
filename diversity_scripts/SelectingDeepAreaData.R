#Mary modified Ole's code, March 23, 2016

rm(list=ls())
library(INLA)
library(rgdal)
library(ggplot2)
library(splancs)
library(dplyr)
library(sp)
library(readr)
library(vegan)

#load("goaDeepAreas.RData")

#### GO GET THE PROJECTION POINTS
#proj.dir	<- "/Users/ole.shelton/GitHub/exxonValdez_nceas/goaTrawl/"
#setwd(proj.dir)
setwd("~/Desktop/Diversity")
#read.csv<-(")

#### HERE I IMPORT DATA USING READ R - SHOULD BE FASTER
#read_csv(shallowDat,"shallow.uncond.csv",col_names=T)
#write_csv(midDat,"mid.uncond.csv",col_names=T)
#write_csv(deepDat,"deep.uncond.csv",col_names=T)
#write_csv(shallowDat.Occ,"shallow.occ.csv",col_names=T)
#write_csv(midDat.Occ,"mid.occ.csv",col_names=T)
#write_csv(deepDat.Occ,"deep.occ.csv",col_names=T)

## OTHE OPTION IS TO LOAD WORKSPACE - WHERE I COMBINEED SHALLOW and DEEP FILES
load("goaDatForKristin.RData")
allDat.pos<-allDat
head(allDat.pos)
head(allDat.Occ)
allDat.occ<-allDat.Occ
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
MIN.D	<-	50
MAX.D	<-	150
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

######## Region 14 - West west of Shelikof Strait
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


### MERGE DATA THAT MATCH ID NUMBERS IN idDat
goaDat.pos<-merge(dat.region, allDat, by="MASTER_ID")
goaDat.occ<-merge(dat.region, allDat.Occ, by="MASTER_ID")

### CALCULATE AVERAGE CPUE & PROPORTIONS BY SPECIES ACROSS SELECTED REGION
goaGF.pos <- goaDat.pos %>%
  group_by(Area, Species, Year) %>% 
  summarise(meanCPUE = mean(Mean)) %>%
  ungroup() %>%
  print(width = Inf)
write.csv(goaGF.pos,"goaDat.deep.pos.csv")

goaGF.occ <- goaDat.occ %>%
  group_by(Area, Species, Year) %>% 
  summarise(meanOcc = mean(Mean)) %>%
  ungroup() %>%
  print(width = Inf)
write.csv(goaGF.occ,"goaDat.deep.occ.csv")


### Function for later
space_plot <- function(dat,NAME){
  par(mfrow=c(3,2),mar=c(3,4,1,0.5))
  AT <- dat$Area  
  plot(1:5,1:5,type="n",axes=F,xlab="",ylab="")
  text(1.5,4,NAME,cex=3,pos=4)
  y.lim <- c(min(dat$Trend - dat$SE.trend),max(dat$Trend + dat$SE.trend))
  plot(Trend~Area,data=dat,ylim=y.lim,axes=F,pch=21,bg=1)
  abline(h=0,lty=2)
  arrows(x0=dat$Area,x1=dat$Area,
         y0=dat$Trend + dat$SE,y1= dat$Trend - dat$SE,
         length=0,lwd=1.5)
  axis(1,at=AT); box(bty="o")
  axis(2,las=2)
  title(xlab="Area",line=2)
  plot(Grand.w.mean~Area,data=dat,type="b",axes=F,xlab="",pch=21,bg=1)
  axis(1,at=AT); box(bty="o")
  axis(2,las=2)
  title(xlab="Area",line=2)
  plot(SD~Area,data=dat,type="b",axes=F,ylim= c(0,max(dat$SD)),xlab="",pch=21,bg=1)
  axis(1,at=AT); box(bty="o")
  axis(2,las=2)
  title(xlab="Area",line=2)
  plot(CV~Area,data=dat,type="b",axes=F,xlab="",pch=21,bg=1)
  axis(1,at=AT); box(bty="o")
  axis(2,las=2)
  title(xlab="Area",line=2)
  plot(Var~Grand.w.mean,data=dat,axes=F,xlab="",pch=21,bg=1)
  axis(1); box(bty="o")
  axis(2,las=2)
  title(xlab="W.mean",line=2)
}

###### Calculate diversity metrics - 
head(goaGF.pos)
goaGF.pos$Area<-as.factor(goaGF.pos$Area)
DivMetrics.pos <- goaGF.pos %>%
group_by(Area, Year) %>% 
  summarise(SW_Div = diversity(meanCPUE, index="shannon"), Eff_Num_Sp = exp(SW_Div),
            simp_Div = diversity(meanCPUE, index="simpson"), 
            invsimp = diversity(meanCPUE,index="inv")) %>%
  ungroup() %>%
  print(width = Inf)

dat<-DivMetrics.pos
AREA <- as.numeric(as.character(unique(dat$Area)))
YEAR <- sort(unique(dat$Year))

#Plot time series by area -Effective number of species
y.lim <- c(min(dat$Eff_Num_Sp),max(dat$Eff_Num_Sp))
for(i in 1:length(AREA)){
  if(AREA[i] !=12){par(new=T)}
    if(AREA[i] == 12 | AREA[i] == 13 ){
    plot(Eff_Num_Sp~Year,data=dat[dat$Area==AREA[i],],
         type="l",axes=F,xlab="",ylim=y.lim,col=2,lwd=2,ylab="")
  }
  if(AREA[i] == 14 | AREA[i] == 15){
    plot(Eff_Num_Sp~Year,data=dat[dat$Area ==AREA[i],],
         type="l",axes=F,xlab="",ylim=y.lim,col=1,lwd=1.5,ylab="")
  }
}
axis(1,at=YEAR,las=2)
axis(2,las=2)
box(bty="o",lwd=2)
title(ylab="Effective number of species",line=2.5)
#mtext(NAME,adj=0)
abline(v=1989,lty=2,lwd=2)

# Calculate least squares trend in each area and portfolio metrics.
vals <- NULL
for(j in 1:length(AREA)){
  A <- lm(Eff_Num_Sp~Year,data=dat[dat$Area == AREA[j] & dat$Year > 1989,]) #,weights = inv.var )
  B <- lm(Eff_Num_Sp~1,data=dat[dat$Area == AREA[j] & dat$Year > 1989,]) #,weights = inv.var )
  vals <- rbind(vals,c(AREA[j],
                       coef(summary(A))["Year",c("Estimate","Std. Error")],
                       coef(summary(B))[1,c("Estimate","Std. Error")],
                       var(A$residuals),
                       sd(A$residuals)))
  
}
vals <- data.frame(vals)
colnames(vals) <- c("Area","Trend","SE.trend","Grand.w.mean","Grand.w.mean.se","Var","SD")
vals$CV  <- vals$SD / vals$Grand.w.mean 

dat.2 <- vals

space_plot(dat=dat.2,NAME="Effective number of species")


#Plot time series by area - Inverse Simpsons
y.lim <- c(min(dat$invsimp),max(dat$invsimp))
for(i in 1:length(AREA)){
  if(AREA[i] !=12){par(new=T)}
  if(AREA[i] == 12 | AREA[i] == 13 ){
    plot(invsimp~Year,data=dat[dat$Area==AREA[i],],
         type="l",axes=F,xlab="",ylim=y.lim,col=2,lwd=2,ylab="")
  }
  if(AREA[i] == 14 | AREA[i] == 15){
    plot(invsimp~Year,data=dat[dat$Area ==AREA[i],],
         type="l",axes=F,xlab="",ylim=y.lim,col=1,lwd=1.5,ylab="")
  }
}
axis(1,at=YEAR,las=2)
axis(2,las=2)
box(bty="o",lwd=2)
title(ylab="Inverse Simpson",line=2.5)
#mtext(NAME,adj=0)
abline(v=1989,lty=2,lwd=2)

# Calculate least squares trend in each area and portfolio metrics.
vals <- NULL
for(j in 1:length(AREA)){
  A <- lm(invsimp~Year,data=dat[dat$Area == AREA[j] & dat$Year > 1989,]) #,weights = inv.var )
  B <- lm(invsimp~1,data=dat[dat$Area == AREA[j] & dat$Year > 1989,]) #,weights = inv.var )
  vals <- rbind(vals,c(AREA[j],
                       coef(summary(A))["Year",c("Estimate","Std. Error")],
                       coef(summary(B))[1,c("Estimate","Std. Error")],
                       var(A$residuals),
                       sd(A$residuals)))
  
}
vals <- data.frame(vals)
colnames(vals) <- c("Area","Trend","SE.trend","Grand.w.mean","Grand.w.mean.se","Var","SD")
vals$CV  <- vals$SD / vals$Grand.w.mean 

dat.2 <- vals

space_plot(dat=dat.2,NAME="Inverse Simpson")

#### CALCULATE SPECIES RICHNESS
head(goaGF.occ)
goaGF.occ$Area<-as.factor(goaGF.occ$Area)
DivMetrics.Occ <- goaGF.occ %>%
  group_by(Area, Year) %>% 
  summarise(SpRich = sum(meanOcc)) %>%
  ungroup() %>%
  print(width = Inf)

dat<-DivMetrics.Occ
AREA <- as.numeric(as.character(unique(dat$Area)))
YEAR <- sort(unique(dat$Year))

#Plot time series by area -Effective number of species
y.lim <- c(min(dat$SpRich),max(dat$SpRich))
for(i in 1:length(AREA)){
  if(AREA[i] !=12){par(new=T)}
  if(AREA[i] == 12 | AREA[i] == 13 ){
    plot(SpRich~Year,data=dat[dat$Area==AREA[i],],
         type="l",axes=F,xlab="",ylim=y.lim,col=2,lwd=2,ylab="")
  }
  if(AREA[i] == 14 | AREA[i] == 15){
    plot(SpRich~Year,data=dat[dat$Area ==AREA[i],],
         type="l",axes=F,xlab="",ylim=y.lim,col=1,lwd=1.5,ylab="")
  }
}
axis(1,at=YEAR,las=2)
axis(2,las=2)
box(bty="o",lwd=2)
title(ylab="Species Richness",line=2.5)
#mtext(NAME,adj=0)
abline(v=1989,lty=2,lwd=2)

# Calculate least squares trend in each area and portfolio metrics.
vals <- NULL
for(j in 1:length(AREA)){
  A <- lm(SpRich~Year,data=dat[dat$Area == AREA[j] & dat$Year > 1989,]) #,weights = inv.var )
  B <- lm(SpRich~1,data=dat[dat$Area == AREA[j] & dat$Year > 1989,]) #,weights = inv.var )
  vals <- rbind(vals,c(AREA[j],
                       coef(summary(A))["Year",c("Estimate","Std. Error")],
                       coef(summary(B))[1,c("Estimate","Std. Error")],
                       var(A$residuals),
                       sd(A$residuals)))
  
}
vals <- data.frame(vals)
colnames(vals) <- c("Area","Trend","SE.trend","Grand.w.mean","Grand.w.mean.se","Var","SD")
vals$CV  <- vals$SD / vals$Grand.w.mean 

dat.2 <- vals

space_plot(dat=dat.2,NAME="Species richness")


