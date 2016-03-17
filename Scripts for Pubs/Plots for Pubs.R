#### Plots for Pub



## Array of figures for Pub:

library(extrafont)
#font_import() #only do this one time - it takes a while
loadfonts(quiet = T)
#windowsFonts(Times=windowsFont("TT Times New Roman"))

library(ggplot2)
library(maps)
library(mapdata)
library(mapproj)
library(dplyr)
library(tidyr)

# Go get the needed data:

meanCPUE <- read.csv("All_sp_index_meanCPUEByArea.csv") # load data
meanCPUE <- subset(meanCPUE,area!='Total') #remove total field
meanCPUE <- droplevels(meanCPUE) 
meanCPUE$area <- as.numeric(levels(meanCPUE$area))[meanCPUE$area]#use numeric/factor conventions
meanCPUE$Area <- factor(meanCPUE$area)
meanCPUE$vari <- meanCPUE$SD.totalDensity^2 #calculate variance
#meanCPUE$Species <- gsub("[.]","",meanCPUE$Species) #quit changing how names are spelled damnit!

discrete_areas <- read.csv("./goaTrawl/Output Data/goa_discrete_areas_for_comparison(50_to_150m).csv")
#head(discrete_areas)

trawl_species <- read.csv("./goaTrawl/Output Data/trawl_species_control_file.csv")
s <- trawl_species$diet
s <- as.data.frame(s) %>% separate(s, into = paste("diet", 1:2, sep = ""))
trawl_species <- cbind(trawl_species,s)

### Merge in the guild identifiers, fish habit, and other functional categories into the meanCPUE data.frame
temp <- trawl_species[,c("database.name","fish.invert","pelagic.benthic","total.biomass.fish","guild","diet1","diet2")]
colnames(temp)[1] <- "Species"

meanCPUE <- merge(meanCPUE,temp)


time_plot <- function(dat,NAME){
  AREA <- as.numeric(as.character(unique(dat$Area)))
  YEAR <- sort(unique(dat$Year))
  y.lim <- c(min(dat$Mean.totalDensity),max(dat$Mean.totalDensity))
  for(i in 1:length(AREA)){
    if(AREA[i] !=1){par(new=T)}
    if(AREA[i] >= 3 & AREA[i] <= 5 ){
      plot(Mean.totalDensity~Year,data=dat[dat$Area ==AREA[i],],
           type="l",axes=F,xlab="",ylim=y.lim,col=2,lwd=1.25,ylab="")
    }
    if(AREA[i] < 3 | AREA[i] > 5){
      plot(Mean.totalDensity~Year,data=dat[dat$Area ==AREA[i],],
           type="l",axes=F,xlab="",ylim=y.lim,col=1,lwd=1,ylab="")
    }
  }
  axis(1,at=YEAR,las=2,hadj=0.7,tcl=-0.25)
  axis(2,las=2,hadj=0.5,tcl=-0.25)
  box(bty="l",lwd=2)
  #title(ylab=expression("Biomass (kg ha"^"-1"*")"),line=2.5)
  mtext(NAME,adj=0.05,line=-1)
  #abline(v=1989,lty=2,lwd=2)
  arrows(x0=1989,x1=1989,y0=y.lim[1],y1=y.lim[1]+c(y.lim[2]-y.lim[1])*0.8,length=0,lty=2,lwd=2)
}

##############################################################################################
def.par <- par(no.readonly = TRUE)

quartz(file="Time-series of Diversity, Biomass, and Guilds.pdf",type="pdf",dpi=300,height=7,width=7)

MAR <- c(2.5,2,0.5,0.1)
A <- layout(matrix(c(1,2,3,4,1,5,6,7,1,8,9,10,1,11,12,13),4,4,byrow=T),
            widths=c(0.1,1,1,1))

# TOTAL BIOMASS
# Plot shared y axis first
par(mar=c(0,0,0,0))
plot(1,1,type="n",xlab="",ylab="",axes=F)
mtext(expression("Density (kg ha"^"-1"*")"),side=2,line=-2,adj=0.25)
mtext(expression("Density (kg ha"^"-1"*")"),side=2,line=-2,adj=0.98)

dat   <- aggregate(meanCPUE[meanCPUE$fish.invert=="fish",c("Mean.totalDensity","vari")],
                   by=list(Area=meanCPUE$Area[meanCPUE$fish.invert=="fish"],Year=meanCPUE$year[meanCPUE$fish.invert=="fish"]),sum)
dat$inv.var <- dat$vari^(-1)
dat$x.bar <- dat$Mean.totalDensity * dat$inv.var

par(mar=MAR)
time_plot(dat,NAME="Total Biomass")

### DIVERSITY METRIC 
for(i in 1:5){
  par(mar=c(0,0,0,0))
  plot(1,1,type="n",xlab="",ylab="",axes=F)
  mtext(paste("Diversity Metric",i),side=3,line=-5)
}  

# GUILDS
GUILD <- c("A","B","P") #unique(meanCPUE$guild)
for(i in 1:length(GUILD)){
  temp <- meanCPUE[meanCPUE$guild == GUILD[i],]
  dat   <- aggregate(temp[,c("Mean.totalDensity","vari")],by=list(Area=temp$Area,Year=temp$year),sum)
  dat$inv.var <- dat$vari^(-1)
  dat$x.bar <- dat$Mean.totalDensity * dat$inv.var
  
  par(mar=MAR)
  time_plot(dat,NAME=paste("Guild",GUILD[i]))
}

# Lifestyle
DIET <- unique(meanCPUE$diet1)
for(i in 1:length(DIET)){
  temp <- meanCPUE[meanCPUE$diet1 == DIET[i],]
  dat   <- aggregate(temp[,c("Mean.totalDensity","vari")],by=list(Area=temp$Area,Year=temp$year),sum)
  dat$inv.var <- dat$vari^(-1)
  dat$x.bar <- dat$Mean.totalDensity * dat$inv.var
  
  time_plot(dat,NAME=paste("Diet",DIET[i]))
}

dev.off()
####################################################################################################
####################################################################################################
###################### REPEAT FOR TREND ESTIMATTION
####################################################################################################
####################################################################################################
trend_plot <- function(dat,NAME){
  AT <- dat$Area  
  y.lim <- c(min(dat$Trend - dat$SE.trend),max(dat$Trend + dat$SE.trend))
  if(y.lim[1]>0){y.lim[1]=-0.1}
  plot(Trend~Area,data=dat,ylim=y.lim,axes=F,pch=21,bg=1)
  abline(h=0,lty=2)
  abline(h=0,lty=2)
  arrows(x0=dat$Area,x1=dat$Area,
         y0=dat$Trend + dat$SE,y1= dat$Trend - dat$SE,
         length=0,lwd=1.5)
  axis(1,at=AT,las=1,padj=-1,tcl=-0.25,cex.axis=0.85)
  axis(2,las=2,hadj=0.5,tcl=-0.25)
  box(bty="l",lwd=2)
  mtext(NAME,adj=0.05,line=-1)
  
  
  #title(xlab="Area",line=2)
}



quartz(file="Trend of Diversity, Biomass, and Guilds since 1990.pdf",type="pdf",dpi=300,height=7,width=7)

MAR <- c(2.5,2,0.5,0.2)
A <- layout(matrix(c(1,2,3,4,1,5,6,7,1,8,9,10,1,11,12,13,1,14,14,14),5,4,byrow=T),
            widths=c(0.08,1,1,1),heights =c(1,1,1,1,0.1) )
#layout.show(A)
# Plot shared y axis first
par(mar=c(0,0,0,0))
plot(1,1,type="n",xlab="",ylab="",axes=F)
mtext(expression("Linear trend (year"^"-1"*")"),side=2,line=-2)

# TOTAL BIOMASS
dat   <- aggregate(meanCPUE[meanCPUE$fish.invert=="fish",c("Mean.totalDensity","vari")],
                   by=list(Area=meanCPUE$Area[meanCPUE$fish.invert=="fish"],Year=meanCPUE$year[meanCPUE$fish.invert=="fish"]),sum)
dat$inv.var <- dat$vari^(-1)
dat$x.bar <- dat$Mean.totalDensity * dat$inv.var

### Calculate least squares trend in each area and portfolio metrics.
AREA <- unique(dat$Area)
vals <- NULL
for(j in 1:length(AREA)){
  A <- lm(Mean.totalDensity~Year,data=dat[dat$Area == AREA[j] & dat$Year > 1989,],weights = inv.var )
  B <- lm(Mean.totalDensity~1,data=dat[dat$Area == AREA[j] & dat$Year > 1989,],weights = inv.var )
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

par(mar=MAR)
trend_plot(dat=dat.2,NAME=paste("Total Biomass"))


### DIVERSITY METRICS 
for(i in 1:5){
  par(mar=c(0,0,0,0))
  plot(1,1,type="n",xlab="",ylab="",axes=F)
  mtext(paste("Diversity Metric",i),side=3,line=-5)
}  
##

GUILD <- c("A","B","P")#unique(meanCPUE$guild)
for(i in 1:length(GUILD)){
  temp <- meanCPUE[meanCPUE$guild == GUILD[i],]
  dat   <- aggregate(temp[,c("Mean.totalDensity","vari")],by=list(Area=temp$Area,Year=temp$year),sum)
  dat$inv.var <- dat$vari^(-1)
  dat$x.bar <- dat$Mean.totalDensity * dat$inv.var
  
  ### Calculate least squares trend in each area and portfolio metrics.
  AREA <- unique(dat$Area)
  vals <- NULL
  for(j in 1:length(AREA)){
    A <- lm(Mean.totalDensity~Year,data=dat[dat$Area == AREA[j] & dat$Year > 1989,],weights = inv.var )
    B <- lm(Mean.totalDensity~1,data=dat[dat$Area == AREA[j] & dat$Year > 1989,],weights = inv.var )
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
  par(mar=MAR)
  trend_plot(dat=dat.2,NAME=paste("Guild ",GUILD[i]))
}

### BY DIET TYPE.
# Aggregate into groups by guild and make plots
DIET <- unique(meanCPUE$diet1)
for(i in 1:length(DIET)){
  temp <- meanCPUE[meanCPUE$diet1 == DIET[i],]
  dat   <- aggregate(temp[,c("Mean.totalDensity","vari")],by=list(Area=temp$Area,Year=temp$year),sum)
  dat$inv.var <- dat$vari^(-1)
  dat$x.bar <- dat$Mean.totalDensity * dat$inv.var
  
  ### Calculate least squares trend in each area and portfolio metrics.
  AREA <- unique(dat$Area)
  vals <- NULL
  for(j in 1:length(AREA)){
    A <- lm(Mean.totalDensity~Year,data=dat[dat$Area == AREA[j] & dat$Year > 1989,],weights = inv.var )
    B <- lm(Mean.totalDensity~1,data=dat[dat$Area == AREA[j] & dat$Year > 1989,],weights = inv.var )
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
  trend_plot(dat=dat.2,NAME=paste("Diet type" ,DIET[i]))
}

# Shared X-axis label
par(mar=c(0,0,0,0))
plot(1,1,type="n",xlab="",ylab="",axes=F)
mtext(expression("Area"),side=3,line=-1)


dev.off()
