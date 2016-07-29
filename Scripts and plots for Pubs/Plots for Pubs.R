#### Plots for Pub

## Array of figures for Pub:
library(extrafont)
library(RColorBrewer) 
#font_import() #only do this one time - it takes a while
#loadfonts(quiet = T)
#windowsFonts(Times=windowsFont("TT Times New Roman"))

library(ggplot2)
library(maps)
#library(mapdata)
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


# Read in von Bertalanffy parameters for species used to calculate age at recruitment to the trawl survey
L.K <- read.csv("./goaTrawl/Output Data/linf_k.csv")
L.K$trait <- as.character(L.K$trait)
L.K$genus.species <- as.character(L.K$genus.species)

### Write function for deriving age at which fish reach a given length in cm
vonB_age <- function(K,Linf,Size){
    return(- log(1- Size/ Linf) / K)
}

# Add estimated age to recruitment to the trawl_species file

# Linf and K are the median of available data from North Pacific (Ben W. did this)
LENGTH <- 20 # size needed to recruit to survey
trawl_species$age_recruit <- NA
for( i in 1:nrow(trawl_species)){
  print(i)
  print(as.character(trawl_species$database.name[i]))
  Linf <- L.K$value[L.K$genus.species == trawl_species$database.name[i] & 
                        L.K$trait == "Linfinity"]
  K <- L.K$value[L.K$genus.species == trawl_species$database.name[i] & 
                      L.K$trait == "K"]
  if(length(K)>0){
    trawl_species$age_recruit[i] <- vonB_age(K,Linf,LENGTH)  
  }
}

trawl_species$age_recruit[trawl_species$database.name == "Podothecus.accipenserinus" ] <- 1.5
trawl_species$age_recruit[trawl_species$database.name == "Myctophidae" ] <- 1.5

trawl_species$age_class[trawl_species$age_recruit <= 2] <- "Short"
trawl_species$age_class[trawl_species$age_recruit > 2 & trawl_species$age_recruit <= 4 ] <- "Medium"
trawl_species$age_class[trawl_species$age_recruit > 4] <- "Long"

# Testing function to make sure it makes sense
# age <- seq(0,12,by=0.1)
#   Linf <- 100
#   K=0.1
#   Lt <- Linf*(1 - exp(-K*age))
#   plot(Lt~age)
#   vonB_age(K,Linf,20)

####
### Merge in the guild identifiers, fish habit, and other functional categories into the meanCPUE data.frame
temp <- trawl_species[,c("database.name","fish.invert","pelagic.benthic","total.biomass.fish","guild","diet1","diet2","age_recruit","age_class")]
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

quartz(file="Time-series of Biomass and Guilds (20 cm).pdf",type="pdf",dpi=300,height=7,width=7)

MAR <- c(2.5,2,0.5,0.1)
A <- layout(matrix(c(1,2,3,4,1,5,6,7,1,8,9,10,1,11,12,13),4,4,byrow=T),
            widths=c(0.1,1,1,1))

# TOTAL BIOMASS
# Plot shared y axis first
par(mar=c(0,0,0,0))
plot(1,1,type="n",xlab="",ylab="",axes=F)
#mtext(expression("Density (kg ha"^"-1"*")"),side=2,line=-2,adj=0.25)
mtext(expression("Density (kg ha"^"-1"*")"),side=2,line=-2)

dat   <- aggregate(meanCPUE[meanCPUE$fish.invert=="fish",c("Mean.totalDensity","vari")],
                   by=list(Area=meanCPUE$Area[meanCPUE$fish.invert=="fish"],Year=meanCPUE$year[meanCPUE$fish.invert=="fish"]),sum)
dat$inv.var <- dat$vari^(-1)
dat$x.bar <- dat$Mean.totalDensity * dat$inv.var

par(mar=MAR)
time_plot(dat,NAME="Total Biomass")

for(i in 1:2){
  plot(1:1,type="n",axes=F,xlab="",ylab="")
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

### Age Lagged METRIC 

AGE <- c("Short","Medium","Long") #unique(meanCPUE$guild)
name.age <- c("Short","Medium","Long")
for(i in 1:length(AGE)){
  temp <- meanCPUE[meanCPUE$age_class == AGE[i],]
  dat   <- aggregate(temp[,c("Mean.totalDensity","vari")],by=list(Area=temp$Area,Year=temp$year),sum)
  dat$inv.var <- dat$vari^(-1)
  dat$x.bar <- dat$Mean.totalDensity * dat$inv.var
  
  par(mar=MAR)
  time_plot(dat,NAME=paste(name.age[i]))
}

dev.off()

####################################################################################################
####################################################################################################
### Make plot standardized to mean.
####################################################################################################
####################################################################################################
####
### Merge in the guild identifiers, fish habit, and other functional categories into the meanCPUE data.frame
temp <- trawl_species[,c("database.name","fish.invert","pelagic.benthic","total.biomass.fish","guild","diet1","diet2","age_recruit","age_class")]
colnames(temp)[1] <- "Species"

meanCPUE <- merge(meanCPUE,temp)

cent_time_plot <- function(dat,NAME){
  AREA <- as.numeric(as.character(unique(dat$Area)))
  YEAR <- sort(unique(dat$Year))
  y.lim <- c(min(dat$cent_density),max(dat$cent_density))
  for(i in 1:length(AREA)){
    if(AREA[i] !=1){par(new=T)}
    if(AREA[i] >= 3 & AREA[i] <= 5 ){
      plot(cent_density~Year,data=dat[dat$Area ==AREA[i],],
           type="l",axes=F,xlab="",ylim=y.lim,col=2,lwd=1.25,ylab="")
    }
    if(AREA[i] < 3 | AREA[i] > 5){
      plot(cent_density~Year,data=dat[dat$Area ==AREA[i],],
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

quartz(file="Centered Time-series of Biomass and Guilds (20 cm).pdf",type="pdf",dpi=300,height=7,width=7)

MAR <- c(2.5,2,0.5,0.1)
A <- layout(matrix(c(1,2,3,4,1,5,6,7,1,8,9,10,1,11,12,13),4,4,byrow=T),
            widths=c(0.1,1,1,1))

# TOTAL BIOMASS
# Plot shared y axis first
par(mar=c(0,0,0,0))
plot(1,1,type="n",xlab="",ylab="",axes=F)
#mtext(expression("Density (kg ha"^"-1"*")"),side=2,line=-2,adj=0.25)
mtext(expression("Centered density (kg ha"^"-1"*")"),side=2,line=-2)

dat   <- aggregate(meanCPUE[meanCPUE$fish.invert=="fish",c("Mean.totalDensity","vari")],
                   by=list(Area=meanCPUE$Area[meanCPUE$fish.invert=="fish"],Year=meanCPUE$year[meanCPUE$fish.invert=="fish"]),sum)
dat$inv.var <- dat$vari^(-1)
dat$x.bar <- dat$Mean.totalDensity * dat$inv.var
grand.mean <- aggregate(dat$Mean.totalDensity,by=list(Area=dat$Area),mean)
dat$cent_density <-0
for( j in 1:11){
  dat$cent_density[dat$Area == j] <- dat$Mean.totalDensity[dat$Area == j] - grand.mean$x[grand.mean$Area == j] 
}

par(mar=MAR)
cent_time_plot(dat,NAME="Total Biomass")
abline(h=0, lty=2,lwd=1.5)

for(i in 1:2){
  plot(1:1,type="n",axes=F,xlab="",ylab="")
}

# GUILDS
GUILD <- c("A","B","P") #unique(meanCPUE$guild)
for(i in 1:length(GUILD)){
  temp <- meanCPUE[meanCPUE$guild == GUILD[i],]
  dat   <- aggregate(temp[,c("Mean.totalDensity","vari")],by=list(Area=temp$Area,Year=temp$year),sum)
  dat$inv.var <- dat$vari^(-1)
  dat$x.bar <- dat$Mean.totalDensity * dat$inv.var
  grand.mean <- aggregate(dat$Mean.totalDensity,by=list(Area=dat$Area),mean)
  dat$cent_density <-0
  for( j in 1:11){
    dat$cent_density[dat$Area == j] <- dat$Mean.totalDensity[dat$Area == j] - grand.mean$x[grand.mean$Area == j] 
  }
  
  par(mar=MAR)
  cent_time_plot(dat,NAME=paste("Guild",GUILD[i]))
  abline(h=0, lty=2,lwd=1.5)
}

# Lifestyle
DIET <- unique(meanCPUE$diet1)
for(i in 1:length(DIET)){
  temp <- meanCPUE[meanCPUE$diet1 == DIET[i],]
  dat   <- aggregate(temp[,c("Mean.totalDensity","vari")],by=list(Area=temp$Area,Year=temp$year),sum)
  dat$inv.var <- dat$vari^(-1)
  dat$x.bar <- dat$Mean.totalDensity * dat$inv.var
  grand.mean <- aggregate(dat$Mean.totalDensity,by=list(Area=dat$Area),mean)
  dat$cent_density <-0
  for( j in 1:11){
    dat$cent_density[dat$Area == j] <- dat$Mean.totalDensity[dat$Area == j] - grand.mean$x[grand.mean$Area == j] 
  }
  
  cent_time_plot(dat,NAME=paste("Diet",DIET[i]))
  abline(h=0, lty=2,lwd=1.5)
}

### Age Lagged METRIC 

AGE <- c("Short","Medium","Long") #unique(meanCPUE$guild)
name.age <- c("Short","Medium","Long")
for(i in 1:length(AGE)){
  temp <- meanCPUE[meanCPUE$age_class == AGE[i],]
  dat   <- aggregate(temp[,c("Mean.totalDensity","vari")],by=list(Area=temp$Area,Year=temp$year),sum)
  dat$inv.var <- dat$vari^(-1)
  dat$x.bar <- dat$Mean.totalDensity * dat$inv.var
  grand.mean <- aggregate(dat$Mean.totalDensity,by=list(Area=dat$Area),mean)
  dat$cent_density <-0
  for( j in 1:11){
    dat$cent_density[dat$Area == j] <- dat$Mean.totalDensity[dat$Area == j] - grand.mean$x[grand.mean$Area == j] 
  }
  
  par(mar=MAR)
  cent_time_plot(dat,NAME=paste(name.age[i]))
  abline(h=0, lty=2,lwd=1.5)
  
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
  plot(Trend~Area,data=dat,ylim=y.lim,axes=F,pch=21,bg=1,type="n")
  polygon(y=c(y.lim[1],0.95*(y.lim[2]-y.lim[1])+y.lim[1],0.95*(y.lim[2]-y.lim[1])+y.lim[1],y.lim[1]),
              x=c(1.5,1.5,6.5,6.5),border=NA,col=grey(0.9))
  polygon(y=c(y.lim[1],0.95*(y.lim[2]-y.lim[1])+y.lim[1],0.95*(y.lim[2]-y.lim[1])+y.lim[1],y.lim[1]),
              x=c(2.5,2.5,5.5,5.5),border=NA,col=grey(0.8))
  par(new=T)
  plot(Trend~Area,data=dat,ylim=y.lim,axes=F,pch=21,bg=1)
  
  abline(h=0,lty=2)
  abline(h=0,lty=2)
  arrows(x0=dat$Area,x1=dat$Area,
         y0=dat$Trend + dat$SE,y1= dat$Trend - dat$SE,
         length=0,lwd=1.5)
  axis(1,at=AT,las=1,padj=-1,tcl=-0.25,cex.axis=0.85)
  axis(2,las=2,hadj=0.5,tcl=-0.25)
  box(bty="l",lwd=2)
  mtext(NAME,adj=0.05,line=-0.75)
  
  
  #title(xlab="Area",line=2)
}

##############################################################
##############################################################
##############################################################

quartz(file="Trend of Biomass and Guilds since 1990.pdf",type="pdf",dpi=300,height=7,width=7)

MAR <- c(2.5,2,0.5,0.2)
A <- layout(matrix(c(1,2,3,4,1,5,6,7,1,8,9,10,1,11,12,13,1,14,14,14),5,4,byrow=T),
            widths=c(0.08,1,1,1),heights =c(1,1,1,1,0.1) )
#layout.show(A)
# Plot shared y axis first
par(mar=c(0,0,0,0))
plot(1,1,type="n",xlab="",ylab="",axes=F)
mtext(expression("Linear trend (kg ha"^"-1"*"year"^"-1"*")"),side=2,line=-2)

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
for(i in 1:2){
  par(mar=c(0,0,0,0))
  plot(1,1,type="n",xlab="",ylab="",axes=F)
  #mtext(paste("Diversity Metric",i),side=3,line=-5)
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

### BY DIET TYPE.
# Aggregate into groups by guild and make plots
AGE <- unique(meanCPUE$age_class[is.na(meanCPUE$age_class)==F])
name.age <- c("Short","Medium","Long")
for(i in 1:length(AGE)){
  temp <- meanCPUE[meanCPUE$age_class == AGE[i],]
  dat   <- aggregate(temp[,c("Mean.totalDensity","vari")],by=list(Area=temp$Area,Year=temp$year),sum)
  dat$inv.var <- dat$vari^(-1)
  dat$x.bar   <- dat$Mean.totalDensity * dat$inv.var
  
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
  trend_plot(dat=dat.2,NAME=paste(name.age[i]))
}

# Shared X-axis label
par(mar=c(0,0,0,0))
plot(1,1,type="n",xlab="",ylab="",axes=F)
mtext(expression("Area"),side=3,line=-1)

dev.off()

######################################################################################
#### CV PLOT
######################################################################################
cv_plot <- function(dat,NAME){
  AT <- dat$Area  
  y.lim <- c(0,round(max(dat$CV.boot+dat$CV.boot.se)*1.15,1))

  if(y.lim[1]>0){y.lim[1]=-0.1}
  plot(CV.boot~Area,data=dat,ylim=y.lim,axes=F,pch=21,bg=1,type="n")
  polygon(y=c(y.lim[1],0.95*(y.lim[2]-y.lim[1])+y.lim[1],0.95*(y.lim[2]-y.lim[1])+y.lim[1],y.lim[1]),
          x=c(1.5,1.5,6.5,6.5),border=NA,col=grey(0.9))
  polygon(y=c(y.lim[1],0.95*(y.lim[2]-y.lim[1])+y.lim[1],0.95*(y.lim[2]-y.lim[1])+y.lim[1],y.lim[1]),
          x=c(2.5,2.5,5.5,5.5),border=NA,col=grey(0.8))
  par(new=T)
  plot(CV.boot~Area,data=dat,ylim=y.lim,axes=F,pch=21,bg=1)
  
  arrows(x0=dat$Area,x1=dat$Area,
         y0=dat$CV.boot + dat$CV.boot.se,y1= dat$CV.boot - dat$CV.boot.se,
         length=0,lwd=1.5)
  axis(1,at=AT,las=1,padj=-1,tcl=-0.25,cex.axis=0.85)
  axis(2,las=2,hadj=0.7,tcl=-0.25)
  box(bty="l",lwd=2)
  mtext(NAME,adj=0.05,line=-0.75)
}

quartz(file="CV of Biomass and Guilds since 1990.pdf",type="pdf",dpi=300,height=7,width=7)

MAR <- c(2.5,3,0.5,0.2)
A <- layout(matrix(c(1,2,3,4,1,5,6,7,1,8,9,10,1,11,12,13,1,14,14,14),5,4,byrow=T),
            widths=c(0.08,1,1,1),heights =c(1,1,1,1,0.1) )
#layout.show(A)
# Plot shared y axis first
par(mar=c(0,0,0,0))
plot(1,1,type="n",xlab="",ylab="",axes=F)
mtext(expression("Coefficient of Variation"),side=2,line=-2)

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
  store <- NULL 
  for(i in 1:1000){
    store[i] <- sd(sample(A$residuals,11,replace=T))
  }
  vals <- rbind(vals,c(AREA[j],
                       coef(summary(A))["Year",c("Estimate","Std. Error")],
                       coef(summary(B))[1,c("Estimate","Std. Error")],
                       var(A$residuals),
                       sd(A$residuals),
                       mean(store / coef(summary(B))[1,c("Estimate")]),
                       sd(store / coef(summary(B))[1,c("Estimate")])))

}
vals <- data.frame(vals)
colnames(vals) <- c("Area","Trend","SE.trend",
                      "Grand.w.mean","Grand.w.mean.se",
                      "Var","SD","CV.boot","CV.boot.se")
vals$CV.point <- vals$SD / vals$Grand.w.mean

dat.2 <- vals

par(mar=MAR)
cv_plot(dat=dat.2,NAME=paste("Total Biomass"))

### DIVERSITY METRICS 
for(i in 1:2){
  par(mar=c(0,0,0,0))
  plot(1,1,type="n",xlab="",ylab="",axes=F)
  #mtext(paste("Diversity Metric",i),side=3,line=-5)
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
    store <- NULL 
      for(k in 1:1000){
        store[k] <- sd(sample(A$residuals,11,replace=T))
      }
      vals <- rbind(vals,c(AREA[j],
                           coef(summary(A))["Year",c("Estimate","Std. Error")],
                           coef(summary(B))[1,c("Estimate","Std. Error")],
                           var(A$residuals),
                           sd(A$residuals),
                           mean(store / coef(summary(B))[1,c("Estimate")]),
                           sd(store / coef(summary(B))[1,c("Estimate")])))
      
    }
    vals <- data.frame(vals)
    colnames(vals) <- c("Area","Trend","SE.trend",
                        "Grand.w.mean","Grand.w.mean.se",
                        "Var","SD","CV.boot","CV.boot.se")
    vals$CV.point <- vals$SD / vals$Grand.w.mean
    
    dat.2 <- vals
    par(mar=MAR)
    cv_plot(dat=dat.2,NAME=paste("Guild ",GUILD[i]))
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
    store <- NULL 
    for(k in 1:1000){
      store[k] <- sd(sample(A$residuals,11,replace=T))
    }
    vals <- rbind(vals,c(AREA[j],
                         coef(summary(A))["Year",c("Estimate","Std. Error")],
                         coef(summary(B))[1,c("Estimate","Std. Error")],
                         var(A$residuals),
                         sd(A$residuals),
                         mean(store / coef(summary(B))[1,c("Estimate")]),
                         sd(store / coef(summary(B))[1,c("Estimate")])))
    
  }
  vals <- data.frame(vals)
  colnames(vals) <- c("Area","Trend","SE.trend",
                      "Grand.w.mean","Grand.w.mean.se",
                      "Var","SD","CV.boot","CV.boot.se")
  vals$CV.point <- vals$SD / vals$Grand.w.mean
  dat.2 <- vals
  cv_plot(dat=dat.2,NAME=paste("Diet type" ,DIET[i]))
}

### BY AGE CLASS TYPE.
# Aggregate into groups by guild and make plots
AGE <- unique(meanCPUE$age_class[is.na(meanCPUE$age_class)==F])
for(i in 1:length(AGE)){
  temp        <- meanCPUE[meanCPUE$age_class == AGE[i],]
  dat         <- aggregate(temp[,c("Mean.totalDensity","vari")],by=list(Area=temp$Area,Year=temp$year),sum)
  dat$inv.var <- dat$vari^(-1)
  dat$x.bar   <- dat$Mean.totalDensity * dat$inv.var
  
  ### Calculate least squares trend in each area and portfolio metrics.
  AREA <- unique(dat$Area)
  vals <- NULL
  for(j in 1:length(AREA)){
    A <- lm(Mean.totalDensity~Year,data=dat[dat$Area == AREA[j] & dat$Year > 1989,],weights = inv.var )
    B <- lm(Mean.totalDensity~1,data=dat[dat$Area == AREA[j] & dat$Year > 1989,],weights = inv.var )
    store <- NULL 
    for(k in 1:1000){
      store[k] <- sd(sample(A$residuals,11,replace=T))
    }
    vals <- rbind(vals,c(AREA[j],
                         coef(summary(A))["Year",c("Estimate","Std. Error")],
                         coef(summary(B))[1,c("Estimate","Std. Error")],
                         var(A$residuals),
                         sd(A$residuals),
                         mean(store / coef(summary(B))[1,c("Estimate")]),
                         sd(store / coef(summary(B))[1,c("Estimate")])))
    
  }
  vals <- data.frame(vals)
  colnames(vals) <- c("Area","Trend","SE.trend",
                      "Grand.w.mean","Grand.w.mean.se",
                      "Var","SD","CV.boot","CV.boot.se")
  vals$CV.point <- vals$SD / vals$Grand.w.mean
  dat.2 <- vals
  cv_plot(dat=dat.2,NAME=paste(name.age[i]))
}

# Shared X-axis label
par(mar=c(0,0,0,0))
plot(1,1,type="n",xlab="",ylab="",axes=F)
mtext(expression("Area"),side=3,line=-1)

dev.off()

######################################################################################
#### SD PLOT
######################################################################################
sd_plot <- function(dat,NAME){
  AT <- dat$Area  
  y.lim <- c(0,round(max(dat$SD)*1.1,1))
  plot(SD ~ Area,data=dat,type="b",axes=F,xlab="",pch=21,bg=1,ylim=y.lim)
  axis(1,at=AT,las=2,hadj=0.5,tcl=-0.25)
  axis(2,las=2,hadj=0.75,tcl=-0.25)
  box(bty="l",lwd=2)
  mtext(NAME,adj=0.05,line=-1)
}

quartz(file="SD of Diversity, Biomass, and Guilds since 1990.pdf",type="pdf",dpi=300,height=7,width=7)

MAR <- c(2.5,3,0.5,0.2)
A <- layout(matrix(c(1,2,3,4,1,5,6,7,1,8,9,10,1,11,12,13,1,14,14,14),5,4,byrow=T),
            widths=c(0.08,1,1,1),heights =c(1,1,1,1,0.1) )
#layout.show(A)
# Plot shared y axis first
par(mar=c(0,0,0,0))
plot(1,1,type="n",xlab="",ylab="",axes=F)
mtext(expression("Standard Deviation"),side=2,line=-2)

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
sd_plot(dat=dat.2,NAME=paste("Total Biomass"))


### DIVERSITY METRICS 
for(i in 1:2){
  par(mar=c(0,0,0,0))
  plot(1,1,type="n",xlab="",ylab="",axes=F)
  #mtext(paste("Diversity Metric",i),side=3,line=-5)
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
  sd_plot(dat=dat.2,NAME=paste("Guild ",GUILD[i]))
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
  sd_plot(dat=dat.2,NAME=paste("Diet type" ,DIET[i]))
}

### BY AGE CLASS TYPE.
# Aggregate into groups by guild and make plots
AGE <- unique(meanCPUE$age_class[is.na(meanCPUE$age_class)==F])
for(i in 1:length(AGE)){
  temp        <- meanCPUE[meanCPUE$age_class == AGE[i],]
  dat         <- aggregate(temp[,c("Mean.totalDensity","vari")],by=list(Area=temp$Area,Year=temp$year),sum)
  dat$inv.var <- dat$vari^(-1)
  dat$x.bar   <- dat$Mean.totalDensity * dat$inv.var
  
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
  sd_plot(dat=dat.2,NAME=paste(AGE[i]))
}

# Shared X-axis label
par(mar=c(0,0,0,0))
plot(1,1,type="n",xlab="",ylab="",axes=F)
mtext(expression("Area"),side=3,line=-1)

dev.off()

############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
#### TERNARY PLOTS
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################

library(ggplot2)
theme_set(theme_bw(base_size=12)+ 
            theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank()))
library(dplyr)
library(tidyr)
library(reshape2)
library(ggtern)

#Pull data for plotting
meanCPUE %>% 
  select(Species,Median.totalDensity, year, Area, diet1, guild,age_class) %>% 
  mutate(Spill = ifelse(year<1990, 'Pre', 'Post'),
         Year = factor(year))-> dat

dat$age_class <- as.factor(dat$age_class)

#Diet data
dat %>% 
  group_by(Year, diet1, Area, Spill) %>% 
  summarize(density = sum(Median.totalDensity)) -> dat1

dat1 <- spread(dat1, diet1, density)

#Guild data
dat %>% 
  group_by(Year, guild, Area, Spill) %>% 
  summarize(density = sum(Median.totalDensity)) -> dat2

dat2 <- spread(dat2, guild, density)

dat2 %>%
  group_by(A,B,P,Area) %>% 
  summarize(D = mean(density)) -> dat2a

dat2a <- dat2
dat2a$A <- dat2$A / (dat2$A+ dat2$B+ dat2$P)
dat2a$B <- dat2$B / (dat2$A+ dat2$B+ dat2$P)
dat2a$P <- dat2$P / (dat2$A+ dat2$B+ dat2$P)

dat2b <- aggregate(dat2a[,c("A","B","P")],by=list(Area=dat2a$Area),mean)
dat2b$Area <- factor(dat2b$Area)

#Guild data
dat %>% 
  group_by(Year, age_class, Area, Spill) %>% 
  summarize(density = sum(Median.totalDensity)) %>% 
  data.frame() -> dat3

dat3 <-   spread(dat3, age_class, density)

############### PLOT
quartz(file="Ternary Plots Guilds, Diet, Age.pdf",type="pdf",dpi=300,height=7,width=7)

#Diet by area
A <- ggtern(data = dat1, aes(x = F, y = G, z = I)) + 
  geom_point(aes(fill = Area),
             size = 4, 
             shape = 21, 
             color = "black",
             alpha=.35) + 
  facet_wrap(~Year) +
  theme_hidegrid() + theme_hidelabels() + 
  theme(legend.justification=c(1,0), legend.position=c(1,0)) + 
  guides(fill=guide_legend(ncol=3))
print(A)

#Diet by year
A <- ggtern(data = dat1, aes(x = F, y = G, z = I)) + 
  geom_point(aes(fill = Year),
             size = 4, 
             shape = 21, 
             color = "black",
             alpha=.35) +facet_wrap(~Area)+theme_hidegrid()+theme_hidelabels()+
  theme(legend.justification=c(1,0), legend.position=c(1,0)) + 
  guides(fill=guide_legend(ncol=3))
print(A)
#Diet pre/post spill

A <- ggtern(data = dat1, aes(x = F, y = G, z = I)) + 
      geom_point(aes(fill = Spill),
             size = 4, 
             shape = 21, 
             color = "black",
             alpha=.35) +facet_wrap(~Area)+
             theme_hidegrid()+theme_hidelabels() +
             theme(legend.justification=c(1,0), legend.position=c(1,0))
print(A)
#Guild by area
A <- ggtern(data = dat2, aes(x = A, y = B, z = P)) + 
  geom_point(aes(fill = Area),
             size = 4, 
             shape = 21, 
             color = "black",
             alpha=.5) +facet_wrap(~Year)+
  theme_hidegrid()+theme_hidelabels()+ 
  theme(legend.justification=c(1,0), legend.position=c(1,0)) + 
  guides(fill=guide_legend(ncol=3))
print(A)

#Guild by year
A <- ggtern(data = dat2, aes(x = A, y = B, z = P)) + 
  geom_point(aes(fill = Year),
             size = 4, 
             shape = 21, 
             color = "black",
             alpha=.5) +facet_wrap(~Area)+
  theme_hidegrid()+theme_hidelabels()+ 
  theme(legend.justification=c(1,0), legend.position=c(1,0)) + 
  guides(fill=guide_legend(ncol=3))
print(A)

#Guild by pre/post
A <- ggtern(data = dat2, aes(x = A, y = B, z = P)) + 
  geom_point(aes(fill = Spill),
             size = 4, 
             shape = 21, 
             color = "black",
             alpha=.5) +facet_wrap(~Area) +theme_hidegrid()+theme_hidelabels()
print(A)

###### Age groups

#AGE by area
A <- ggtern(data = dat3, aes(x = Long, y = Medium, z = Short)) + 
  geom_point(aes(fill = Area),
             size = 4, 
             shape = 21, 
             color = "black",
             alpha=.5) +facet_wrap(~Year)+
  theme_hidegrid()+theme_hidelabels()+ 
  theme(legend.justification=c(1,0), legend.position=c(1,0)) + 
  guides(fill=guide_legend(ncol=3))
print(A)

#Age by year
A <- ggtern(data = dat3, aes(x = Long, y = Medium, z = Short)) + 
  geom_point(aes(fill = Year),
             size = 4, 
             shape = 21, 
             color = "black",
             alpha=.5) +facet_wrap(~Area)+
  theme_hidegrid()+theme_hidelabels()+ 
  theme(legend.justification=c(1,0), legend.position=c(1,0)) + 
  guides(fill=guide_legend(ncol=3))
print(A)

#Age by pre/post
A <- ggtern(data = dat3, aes(x = Long, y = Medium, z = Short)) + 
  geom_point(aes(fill = Spill),
             size = 4, 
             shape = 21, 
             color = "black",
             alpha=.5) +facet_wrap(~Area) +theme_hidegrid()+theme_hidelabels()
print(A)

dev.off()



quartz(file="Ternary Plots Guilds ONLY.pdf",type="pdf",dpi=300,height=7,width=7)

dat2b$area.col <- as.numeric(as.character(dat2b$Area))

# Guild by year
A <- ggtern() + 
  geom_point(data= dat2, aes(x = A, y = B, z = P),
#              size = 4, 
#              shape = 21, 
#              color = "black",
             alpha=.15) +#facet_wrap(~Area)+
  geom_point(data= dat2b, aes(x = A, y = B, z = P, fill = area.col), 
                           size = 5,
                           shape =21,
                           color="black",alpha=0.8
                           ) +
  scale_fill_gradient2(low = "red",high="blue",breaks=1:11,name="Area")+
  #theme_hidegrid()+theme_hidelabels()+ 
  theme(legend.position=c(0.88,0.5)) +
  guides(fill=guide_legend(ncol=1)) 
  print(A)

dev.off()



