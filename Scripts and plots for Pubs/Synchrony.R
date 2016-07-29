### Calculations of synchrony 
library(synchrony)
library(ggplot2)
library(maps)
#library(mapdata)
library(mapproj)
library(dplyr)
library(tidyr)
library(reshape2)

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
meanCPUE <- merge(meanCPUE,temp)############### THIS IS THE BASIC IDEX-STANDARDIZED DATA FROM PROJECTION

### Start Analysis on the total biomass

# TOTAL BIOMASS
# Plot shared y axis first

dat   <- aggregate(meanCPUE[meanCPUE$fish.invert=="fish",c("Mean.totalDensity","vari")],
                   by=list(Area=meanCPUE$Area[meanCPUE$fish.invert=="fish"],Year=meanCPUE$year[meanCPUE$fish.invert=="fish"]),sum)

dat.wide <- reshape(dat,timevar=c("Area"),idvar="Year",direction="wide",drop="vari")

#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
# all pairwise correlations (pearson)
n.Area <- 11
cor.mat <- NULL
for(i in 1:n.Area){
  for(j in 1:n.Area){
    A <- cor.test(dat.wide[,i+1],dat.wide[,j+1],method="pearson")
    cor.mat   <- rbind(cor.mat,data.frame(Area.1 = i, Area.2=j, Pearson = A$estimate))
    
  }
}
cor.mat <- cor.mat[cor.mat$Area.1 <= cor.mat$Area.2,]
################################
################################ WINDOW!!
################################
################################
################################
################################
################################
################################

# Moving window approach
# Use nine year windows (will incorporate varying numbers of observations due to spacing of surveys)
WINDOW <- 9

start.year <- sort(unique(dat$Year))
start.year <- start.year[start.year <= (max(start.year)-WINDOW)]
end.year   <- start.year + WINDOW


mean.year <- NULL
for(k in 1:length(start.year) ){
  mean.year[k] <- min(dat.wide$Year[dat.wide$Year >= start.year[k] & dat.wide$Year <= end.year[k]]) + WINDOW / 2
}

cor.window  <- NULL
for(k in 1:length(start.year) ){
  temp        <- matrix(0,n.Area,n.Area)
  for(i in 1:n.Area){
    for(j in 1:n.Area){
      A <-  cor.test(dat.wide[dat.wide$Year >= start.year[k] & dat.wide$Year <= end.year[k] ,i+1],
                        dat.wide[dat.wide$Year >= start.year[k] & dat.wide$Year <= end.year[k],j+1],
                        method="pearson")
     cor.window   <- rbind(cor.window,data.frame(Area.1 = i, Area.2=j, Year = mean.year[k], Pearson = A$estimate))
    }
  }
  # colnames(temp) <- 
  # temp <- data.frame(Area=1:n.Area, year= rep(mean.year[k],n.Area),temp)
  # temp$Area <- (1:n.Area)
  # temp$year <- mean.year[k]
  # cor.window <- rbind(cor.window,temp)  
}
cor.window.trim <- cor.window[cor.window$Area.1 <= cor.window$Area.2,]

#####################################
#####################################
#####################################
## MAKE SOME PLOTS OF THE CORRELATION MATRIX
#####################################
#####################################
#####################################
cor.mat.plot <- cor.mat
cor.mat.plot$Area.1 <- factor(cor.mat.plot$Area.1)
cor.mat.plot$Area.2 <- factor(cor.mat.plot$Area.2)
cor.mat.plot$Pearson[cor.mat.plot$Pearson > 0.99999999] <- -99
### MAKE A PLOT OF 4 grids (true, mean est, sd, difference.)              
p.all <- ggplot(cor.mat.plot,
                 aes(x = Area.1, y = Area.2, fill = Pearson)) +
  geom_tile() +
  scale_fill_gradientn(colours = c("blue",grey(0.8),"red"),limits=c(-1,1),na.value="black") +
    coord_equal() +
  scale_x_discrete(breaks=1:11)+
  scale_y_discrete(breaks=1:11,limits = rev(levels(cor.mat.plot$Area.1))) +
  labs(x= "Area", y= "Area") +
  ggtitle("1984-2015") +
  theme_bw()

print(p.all)
##### YEAR BY YEAR

cor.wind.plot <- cor.window.trim
cor.wind.plot$Area.1 <- factor(cor.wind.plot$Area.1)
cor.wind.plot$Area.2 <- factor(cor.wind.plot$Area.2)
cor.wind.plot$Pearson[cor.wind.plot$Pearson > 0.999999] <- -99

P <- list()

for(i in 1:length(mean.year)){
  
  p.temp <- ggplot(cor.wind.plot[cor.wind.plot$Year == mean.year[i],],
                  aes(x = Area.1, y = Area.2, fill = Pearson)) +
    geom_tile() +
    scale_fill_gradientn(colours = c("blue",grey(0.8),"red"),limits=c(-1,1),na.value="black") +
    coord_equal() +
    scale_x_discrete(breaks=1:11)+
    scale_y_discrete(breaks=1:11,limits = rev(levels(cor.mat.plot$Area.1))) +
    labs(x= "Area", y= "Area") +
    ggtitle(paste(start.year[i],"to",end.year[i])) +
    theme_bw()
  
  P[[i]] <- p.temp
}


quartz(file=paste("./Scripts and plots for Pubs/Pearson cor matrices,window=",WINDOW,".pdf"),dpi=300,height=7,width=7,type="pdf")
  print(p.all)
  for(i in 1:length(mean.year)){
    print(P[[i]])
  }
dev.off()

##### 
# Add indicator variable for EVOS - EVOS, EVOS-non, and non-non combinations of 


cor.window$tag.1[cor.window$Area.1 >=2 & cor.window$Area.2 <=6 ]  <- "EVOS"
cor.window$tag.1[cor.window$Area.1 ==1 & cor.window$Area.2 >=2 &  cor.window$Area.2 <=6 ] <- "EVOS-Control"
cor.window$tag.1[cor.window$Area.2 >=7 & cor.window$Area.1 >=2 &  cor.window$Area.1 <=6 ] <- "EVOS-Control"
cor.window$tag.1[cor.window$Area.1 ==1 & cor.window$Area.2 > 6 ]  <- "Control"
cor.window$tag.1[cor.window$Area.1 >=7 & cor.window$Area.2 >= 7 ] <- "Control"

cor.window$tag.2[cor.window$Area.1 >=3 & cor.window$Area.2 <=5 ]  <- "EVOS"
cor.window$tag.2[cor.window$Area.1 <=2 & cor.window$Area.2 >=3 &  cor.window$Area.2 <=5 ] <- "EVOS-Control"
cor.window$tag.2[cor.window$Area.2 >=6 & cor.window$Area.1 >=3 &  cor.window$Area.1 <=5 ] <- "EVOS-Control"
cor.window$tag.2[cor.window$Area.1 <=2 & cor.window$Area.2 > 5 ]  <- "Control"
cor.window$tag.2[cor.window$Area.1 <=2 & cor.window$Area.2 <= 2 ]  <- "Control"
cor.window$tag.2[cor.window$Area.1 >=6 & cor.window$Area.2 >= 6 ] <- "Control"
cor.wind.plot <- cor.window[cor.window$Area.1 != cor.window$Area.2, ]

cor.wind.agg <- aggregate(cor.wind.plot$Pearson, by=list(Year=cor.wind.plot$Year,tag.1=cor.wind.plot$tag.1),median)
colnames(cor.wind.agg)[3] <- "Pearson.median"
cor.wind.agg2 <- aggregate(cor.wind.plot$Pearson, by=list(Year=cor.wind.plot$Year,tag.1=cor.wind.plot$tag.1),mean)
colnames(cor.wind.agg2)[3] <- "Pearson.mean"
cor.tag.1 <- merge(cor.wind.agg,cor.wind.agg2)

cor.wind.agg <- aggregate(cor.wind.plot$Pearson, by=list(Year=cor.wind.plot$Year,tag.2=cor.wind.plot$tag.2),median)
colnames(cor.wind.agg)[3] <- "Pearson.median"
cor.wind.agg2 <- aggregate(cor.wind.plot$Pearson, by=list(Year=cor.wind.plot$Year,tag.2=cor.wind.plot$tag.2),mean)
colnames(cor.wind.agg2)[3] <- "Pearson.mean"
cor.tag.2 <- merge(cor.wind.agg,cor.wind.agg2)

###
p.temp <- ggplot(cor.wind.plot,
                 aes(x = Year, y = Pearson, colour=factor(tag.1))) +
  geom_jitter() 

p.temp

  # scale_fill_gradientn(colours = c("blue",grey(0.8),"red"),limits=c(-1,1),na.value="black") +
  # coord_equal() +
  # scale_x_discrete(breaks=1:11)+
  # scale_y_discrete(breaks=1:11,limits = rev(levels(cor.mat.plot$Area.1))) +
  # labs(x= "Area", y= "Area") +
  # ggtitle(paste(start.year[i],"to",end.year[i])) +
  # theme_bw()

  p.line <- ggplot(cor.tag.1,
                   aes(x = Year, y = Pearson.median, colour=factor(tag.1))) +
     geom_jitter(data=cor.wind.plot,
                 aes(x = Year, y = Pearson, colour=factor(tag.1)),width=0.5,height=0,alpha=0.3) +
    geom_line(size=1.2) +
    geom_point(size=3) +
    geom_hline(yintercept=0,linetype=3)+
    scale_y_continuous(limits=c(-0.75,1)) +
    theme_bw()
  p.line
  
  p.line.mean <- ggplot(cor.tag.1,
                   aes(x = Year, y = Pearson.mean, colour=factor(tag.1))) +
     geom_jitter(data=cor.wind.plot,
                 aes(x = Year, y = Pearson, colour=factor(tag.1)),width=0.5,height=0,alpha=0.3) +
    geom_line(size=1.2) +
    geom_point(size=3) +
    geom_hline(yintercept=0,linetype=3)+
    scale_y_continuous(limits=c(-0.75,1)) +
    theme_bw()
  
  
  quartz(file=paste("./Scripts and plots for Pubs/Pearson time-series,Biomass,window=",WINDOW,".pdf"),dpi=300,height=4,width=6,type="pdf")
      print(p.line)
      print(p.line.mean)
  dev.off()
  
  #### TAG.2 (more stringent EVOS classification)
  
  
  p.line <- ggplot(cor.tag.2,
                   aes(x = Year, y = Pearson.median, colour=factor(tag.2))) +
     geom_jitter(data=cor.wind.plot,
                 aes(x = Year, y = Pearson, colour=factor(tag.1)),width=0.5,height=0,alpha=0.3) +
    geom_line(size=1.2) +
    geom_point(size=3) +
    geom_hline(yintercept=0,linetype=3)+
    scale_y_continuous(limits=c(-0.75,1)) +
    theme_bw()
  p.line
  
  p.line.mean <- ggplot(cor.tag.2,
                        aes(x = Year, y = Pearson.mean, colour=factor(tag.2))) +
    geom_jitter(data=cor.wind.plot,
                 aes(x = Year, y = Pearson, colour=factor(tag.1)),width=0.5,height=0,alpha=0.3) +
    geom_line(size=1.2) +
    geom_point(size=3) +
    geom_hline(yintercept=0,linetype=3)+
    scale_y_continuous(limits=c(-0.75,1)) +
    theme_bw()
  
  setwd("./Scripts and plots for Pubs/")
  quartz(file=paste("Pearson time series, Biomass (Tag 2),window=",WINDOW,".pdf"),dpi=300,height=4,width=6,type="pdf")
    print(p.line)
    print(p.line.mean)
  dev.off()
  
  ###################################################################################################
  ###################################################################################################
  ###################################################################################################
  ###################################################################################################
  ###################################################################################################
  #### REPEAT FOR EACH OF THE FUNCITONAL GROUPS
  ###################################################################################################
  ###################################################################################################
  ###################################################################################################
  ###################################################################################################
  ###################################################################################################
  
  #Guilds  
  GUILD <- unique(meanCPUE$guild)
  
  for(QQ in 1:3){
    dat   <- aggregate(meanCPUE[meanCPUE$guild==GUILD[QQ],c("Mean.totalDensity","vari")],
                     by=list(Area=meanCPUE$Area[meanCPUE$guild==GUILD[QQ]],Year=meanCPUE$year[meanCPUE$guild==GUILD[QQ]]),sum)
  
    dat.wide <- reshape(dat,timevar=c("Area"),idvar="Year",direction="wide",drop="vari")
  #########################################################
  #########################################################
  # all pairwise correlations (pearson)
  n.Area <- 11
  cor.mat <- NULL
  for(i in 1:n.Area){
    for(j in 1:n.Area){
      A <- cor.test(dat.wide[,i+1],dat.wide[,j+1],method="pearson")
      cor.mat   <- rbind(cor.mat,data.frame(Area.1 = i, Area.2=j, Year = mean.year[k], Pearson = A$estimate))
      
    }
  }
  cor.mat <- cor.mat[cor.mat$Area.1 <= cor.mat$Area.2,]

  # Moving window approach
  # Use six year windows (will incorporate varying numbers of observations due to )
  
  start.year <- sort(unique(dat$Year))
  start.year <- start.year[start.year <= (max(start.year)-WINDOW)]
  end.year   <- start.year + WINDOW

    mean.year <- NULL
  for(k in 1:length(start.year) ){
    mean.year[k] <- dat.wide$Year[dat.wide$Year >= start.year[k] & dat.wide$Year <= end.year[k]] + WINDOW / 2
  }
  
  cor.window  <- NULL
  for(k in 1:length(start.year) ){
    temp        <- matrix(0,n.Area,n.Area)
    for(i in 1:n.Area){
      for(j in 1:n.Area){
        A <-  cor.test(dat.wide[dat.wide$Year >= start.year[k] & dat.wide$Year <= end.year[k] ,i+1],
                       dat.wide[dat.wide$Year >= start.year[k] & dat.wide$Year <= end.year[k],j+1],
                       method="pearson")
        cor.window   <- rbind(cor.window,data.frame(Area.1 = i, Area.2=j, Year = mean.year[k], Pearson = A$estimate))
      }
    }
    # colnames(temp) <- 
    # temp <- data.frame(Area=1:n.Area, year= rep(mean.year[k],n.Area),temp)
    # temp$Area <- (1:n.Area)
    # temp$year <- mean.year[k]
    # cor.window <- rbind(cor.window,temp)  
  }
  cor.window.trim <- cor.window[cor.window$Area.1 <= cor.window$Area.2,]
  
  #####################################
  #####################################
  #####################################
  ## MAKE SOME PLOTS OF THE CORRELATION MATRIX
  #####################################
  #####################################
  #####################################
  cor.mat.plot <- cor.mat
  cor.mat.plot$Area.1 <- factor(cor.mat.plot$Area.1)
  cor.mat.plot$Area.2 <- factor(cor.mat.plot$Area.2)
  cor.mat.plot$Pearson[cor.mat.plot$Pearson > 0.99999999] <- -99
  ### MAKE A PLOT OF 4 grids (true, mean est, sd, difference.)              
  p.all <- ggplot(cor.mat.plot,
                  aes(x = Area.1, y = Area.2, fill = Pearson)) +
    geom_tile() +
    scale_fill_gradientn(colours = c("blue",grey(0.8),"red"),limits=c(-1,1),na.value="black") +
    coord_equal() +
    scale_x_discrete(breaks=1:11)+
    scale_y_discrete(breaks=1:11,limits = rev(levels(cor.mat.plot$Area.1))) +
    labs(x= "Area", y= "Area") +
    ggtitle(paste("Guild",GUILD[QQ],"1984-2015")) +
    theme_bw()
  
  p.all
  ##### YEAR BY YEAR
  
  cor.wind.plot <- cor.window.trim
  cor.wind.plot$Area.1 <- factor(cor.wind.plot$Area.1)
  cor.wind.plot$Area.2 <- factor(cor.wind.plot$Area.2)
  cor.wind.plot$Pearson[cor.wind.plot$Pearson > 0.999999] <- -99
  
  P <- list()
  
  for(i in 1:length(mean.year)){
    
    p.temp <- ggplot(cor.wind.plot[cor.wind.plot$Year == mean.year[i],],
                     aes(x = Area.1, y = Area.2, fill = Pearson)) +
      geom_tile() +
      scale_fill_gradientn(colours = c("blue",grey(0.8),"red"),limits=c(-1,1),na.value="black") +
      coord_equal() +
      scale_x_discrete(breaks=1:11)+
      scale_y_discrete(breaks=1:11,limits = rev(levels(cor.mat.plot$Area.1))) +
      labs(x= "Area", y= "Area") +
      ggtitle(paste("Guild",GUILD[QQ],start.year[i],"to",end.year[i])) +
      theme_bw()
    
    P[[i]] <- p.temp
  }
  
  
  quartz(file=paste("./Scripts and plots for Pubs/Pearson cor matrices, Guild=",GUILD[QQ],", window=",WINDOW,".pdf"),dpi=300,height=7,width=7,type="pdf")
  print(p.all)
  for(i in 1:length(mean.year)){
    print(P[[i]])
  }
  dev.off()
  
  ##### 
  # Add indicator variable for EVOS - EVOS, EVOS-non, and non-non combinations of 
  
  
  cor.window$tag.1[cor.window$Area.1 >=2 & cor.window$Area.2 <=6 ]  <- "EVOS"
  cor.window$tag.1[cor.window$Area.1 ==1 & cor.window$Area.2 >=2 &  cor.window$Area.2 <=6 ] <- "EVOS-Control"
  cor.window$tag.1[cor.window$Area.2 >=7 & cor.window$Area.1 >=2 &  cor.window$Area.1 <=6 ] <- "EVOS-Control"
  cor.window$tag.1[cor.window$Area.1 ==1 & cor.window$Area.2 > 6 ]  <- "Control"
  cor.window$tag.1[cor.window$Area.1 >=7 & cor.window$Area.2 >= 7 ] <- "Control"
  
  cor.window$tag.2[cor.window$Area.1 >=3 & cor.window$Area.2 <=5 ]  <- "EVOS"
  cor.window$tag.2[cor.window$Area.1 <=2 & cor.window$Area.2 >=3 &  cor.window$Area.2 <=5 ] <- "EVOS-Control"
  cor.window$tag.2[cor.window$Area.2 >=6 & cor.window$Area.1 >=3 &  cor.window$Area.1 <=5 ] <- "EVOS-Control"
  cor.window$tag.2[cor.window$Area.1 <=2 & cor.window$Area.2 > 5 ]  <- "Control"
  cor.window$tag.2[cor.window$Area.1 <=2 & cor.window$Area.2 <= 2 ]  <- "Control"
  cor.window$tag.2[cor.window$Area.1 >=6 & cor.window$Area.2 >= 6 ] <- "Control"
  cor.wind.plot <- cor.window[cor.window$Area.1 != cor.window$Area.2, ]
  
  cor.wind.agg <- aggregate(cor.wind.plot$Pearson, by=list(Year=cor.wind.plot$Year,tag.1=cor.wind.plot$tag.1),median)
  colnames(cor.wind.agg)[3] <- "Pearson.median"
  cor.wind.agg2 <- aggregate(cor.wind.plot$Pearson, by=list(Year=cor.wind.plot$Year,tag.1=cor.wind.plot$tag.1),mean)
  colnames(cor.wind.agg2)[3] <- "Pearson.mean"
  cor.tag.1 <- merge(cor.wind.agg,cor.wind.agg2)
  
  cor.wind.agg <- aggregate(cor.wind.plot$Pearson, by=list(Year=cor.wind.plot$Year,tag.2=cor.wind.plot$tag.2),median)
  colnames(cor.wind.agg)[3] <- "Pearson.median"
  cor.wind.agg2 <- aggregate(cor.wind.plot$Pearson, by=list(Year=cor.wind.plot$Year,tag.2=cor.wind.plot$tag.2),mean)
  colnames(cor.wind.agg2)[3] <- "Pearson.mean"
  cor.tag.2 <- merge(cor.wind.agg,cor.wind.agg2)
  
  ###
  p.temp <- ggplot(cor.wind.plot,
                   aes(x = Year, y = Pearson, colour=factor(tag.1))) +
    geom_jitter() 
  
  p.temp
  
  # scale_fill_gradientn(colours = c("blue",grey(0.8),"red"),limits=c(-1,1),na.value="black") +
  # coord_equal() +
  # scale_x_discrete(breaks=1:11)+
  # scale_y_discrete(breaks=1:11,limits = rev(levels(cor.mat.plot$Area.1))) +
  # labs(x= "Area", y= "Area") +
  # ggtitle(paste(start.year[i],"to",end.year[i])) +
  # theme_bw()
  
  p.line <- ggplot(cor.tag.1,
                   aes(x = Year, y = Pearson.median, colour=factor(tag.1))) +
    geom_jitter(data=cor.wind.plot,
                aes(x = Year, y = Pearson, colour=factor(tag.1)),width=0.5,height=0,alpha=0.3) +
    geom_line(size=1.2) +
    geom_point(size=3) +
    geom_hline(yintercept=0,linetype=3)+
    scale_y_continuous(limits=c(-0.75,1)) +
    theme_bw()
  p.line
  
  p.line.mean <- ggplot(cor.tag.1,
                        aes(x = Year, y = Pearson.mean, colour=factor(tag.1))) +
    geom_jitter(data=cor.wind.plot,
                aes(x = Year, y = Pearson, colour=factor(tag.1)),width=0.5,height=0,alpha=0.3) +
    geom_line(size=1.2) +
    geom_point(size=3) +
    geom_hline(yintercept=0,linetype=3)+
    scale_y_continuous(limits=c(-0.75,1)) +
    theme_bw()
  
  
  quartz(file=paste("./Scripts and plots for Pubs/Pearson time series, Guild=",GUILD[QQ],", window=",WINDOW,".pdf"),dpi=300,height=4,width=6,type="pdf")
  print(p.line)
  print(p.line.mean)
  dev.off()
  
  #### TAG.2 (more stringent EVOS classification)
  
  
  p.line <- ggplot(cor.tag.2,
                   aes(x = Year, y = Pearson.median, colour=factor(tag.2))) +
    geom_jitter(data=cor.wind.plot,
                aes(x = Year, y = Pearson, colour=factor(tag.1)),width=0.5,height=0,alpha=0.3) +
    geom_line(size=1.2) +
    geom_point(size=3) +
    geom_hline(yintercept=0,linetype=3)+
    scale_y_continuous(limits=c(-0.75,1)) +
    theme_bw()
  p.line
  
  p.line.mean <- ggplot(cor.tag.2,
                        aes(x = Year, y = Pearson.mean, colour=factor(tag.2))) +
    geom_jitter(data=cor.wind.plot,
                aes(x = Year, y = Pearson, colour=factor(tag.1)),width=0.5,height=0,alpha=0.3) +
    geom_line(size=1.2) +
    geom_point(size=3) +
    geom_hline(yintercept=0,linetype=3)+
    scale_y_continuous(limits=c(-0.75,1)) +
    theme_bw()
  
  
  quartz(file=paste("./Scripts and plots for Pubs/Pearson time series (tag 2), Guild=",GUILD[QQ],", window=",WINDOW,".pdf"),dpi=300,height=4,width=6,type="pdf")
  print(p.line)
  print(p.line.mean)
  dev.off()
  
  }
  
  ###################################################################################################
  ###################################################################################################
  ###################################################################################################
  
  #DIET  
  DIET <- unique(meanCPUE$diet1)
  
  for(QQ in 1:3){
    dat   <- aggregate(meanCPUE[meanCPUE$diet1==DIET[QQ],c("Mean.totalDensity","vari")],
                       by=list(Area=meanCPUE$Area[meanCPUE$diet1==DIET[QQ]],Year=meanCPUE$year[meanCPUE$diet1==DIET[QQ]]),sum)
    
    dat.wide <- reshape(dat,timevar=c("Area"),idvar="Year",direction="wide",drop="vari")
    #########################################################
    #########################################################
    # all pairwise correlations (pearson)
    n.Area <- 11
    cor.mat <- NULL
    for(i in 1:n.Area){
      for(j in 1:n.Area){
        A <- cor.test(dat.wide[,i+1],dat.wide[,j+1],method="pearson")
        cor.mat   <- rbind(cor.mat,data.frame(Area.1 = i, Area.2=j, Year = mean.year[k], Pearson = A$estimate))
        
      }
    }
    cor.mat <- cor.mat[cor.mat$Area.1 <= cor.mat$Area.2,]
    
    # Moving window approach
    # Use six year windows (will incorporate varying numbers of observations due to )
    
    start.year <- sort(unique(dat$Year))
    start.year <- start.year[start.year <= (max(start.year)-WINDOW)]
    end.year   <- start.year + WINDOW
    
    mean.year <- NULL
    for(k in 1:length(start.year) ){
      mean.year[k] <- dat.wide$Year[dat.wide$Year >= start.year[k] & dat.wide$Year <= end.year[k]] + WINDOW / 2
    }
    
    cor.window  <- NULL
    for(k in 1:length(start.year) ){
      temp        <- matrix(0,n.Area,n.Area)
      for(i in 1:n.Area){
        for(j in 1:n.Area){
          A <-  cor.test(dat.wide[dat.wide$Year >= start.year[k] & dat.wide$Year <= end.year[k] ,i+1],
                         dat.wide[dat.wide$Year >= start.year[k] & dat.wide$Year <= end.year[k],j+1],
                         method="pearson")
          cor.window   <- rbind(cor.window,data.frame(Area.1 = i, Area.2=j, Year = mean.year[k], Pearson = A$estimate))
        }
      }
      # colnames(temp) <- 
      # temp <- data.frame(Area=1:n.Area, year= rep(mean.year[k],n.Area),temp)
      # temp$Area <- (1:n.Area)
      # temp$year <- mean.year[k]
      # cor.window <- rbind(cor.window,temp)  
    }
    cor.window.trim <- cor.window[cor.window$Area.1 <= cor.window$Area.2,]
    
    #####################################
    #####################################
    #####################################
    ## MAKE SOME PLOTS OF THE CORRELATION MATRIX
    #####################################
    #####################################
    #####################################
    cor.mat.plot <- cor.mat
    cor.mat.plot$Area.1 <- factor(cor.mat.plot$Area.1)
    cor.mat.plot$Area.2 <- factor(cor.mat.plot$Area.2)
    cor.mat.plot$Pearson[cor.mat.plot$Pearson > 0.99999999] <- -99
    ### MAKE A PLOT OF 4 grids (true, mean est, sd, difference.)              
    p.all <- ggplot(cor.mat.plot,
                    aes(x = Area.1, y = Area.2, fill = Pearson)) +
      geom_tile() +
      scale_fill_gradientn(colours = c("blue",grey(0.8),"red"),limits=c(-1,1),na.value="black") +
      coord_equal() +
      scale_x_discrete(breaks=1:11)+
      scale_y_discrete(breaks=1:11,limits = rev(levels(cor.mat.plot$Area.1))) +
      labs(x= "Area", y= "Area") +
      ggtitle(paste("Diet",DIET[QQ],"1984-2015")) +
      theme_bw()
    
    p.all
    ##### YEAR BY YEAR
    
    cor.wind.plot <- cor.window.trim
    cor.wind.plot$Area.1 <- factor(cor.wind.plot$Area.1)
    cor.wind.plot$Area.2 <- factor(cor.wind.plot$Area.2)
    cor.wind.plot$Pearson[cor.wind.plot$Pearson > 0.999999] <- -99
    
    P <- list()
    
    for(i in 1:length(mean.year)){
      
      p.temp <- ggplot(cor.wind.plot[cor.wind.plot$Year == mean.year[i],],
                       aes(x = Area.1, y = Area.2, fill = Pearson)) +
        geom_tile() +
        scale_fill_gradientn(colours = c("blue",grey(0.8),"red"),limits=c(-1,1),na.value="black") +
        coord_equal() +
        scale_x_discrete(breaks=1:11)+
        scale_y_discrete(breaks=1:11,limits = rev(levels(cor.mat.plot$Area.1))) +
        labs(x= "Area", y= "Area") +
        ggtitle(paste("Diet",DIET[QQ],start.year[i],"to",end.year[i])) +
        theme_bw()
      
      P[[i]] <- p.temp
    }
    
    
    quartz(file=paste("./Scripts and plots for Pubs/Pearson cor matrices, Diet=",DIET[QQ],", window=",WINDOW,".pdf"),dpi=300,height=7,width=7,type="pdf")
    print(p.all)
    for(i in 1:length(mean.year)){
      print(P[[i]])
    }
    dev.off()
    
    ##### 
    # Add indicator variable for EVOS - EVOS, EVOS-non, and non-non combinations of 
    
    
    cor.window$tag.1[cor.window$Area.1 >=2 & cor.window$Area.2 <=6 ]  <- "EVOS"
    cor.window$tag.1[cor.window$Area.1 ==1 & cor.window$Area.2 >=2 &  cor.window$Area.2 <=6 ] <- "EVOS-Control"
    cor.window$tag.1[cor.window$Area.2 >=7 & cor.window$Area.1 >=2 &  cor.window$Area.1 <=6 ] <- "EVOS-Control"
    cor.window$tag.1[cor.window$Area.1 ==1 & cor.window$Area.2 > 6 ]  <- "Control"
    cor.window$tag.1[cor.window$Area.1 >=7 & cor.window$Area.2 >= 7 ] <- "Control"
    
    cor.window$tag.2[cor.window$Area.1 >=3 & cor.window$Area.2 <=5 ]  <- "EVOS"
    cor.window$tag.2[cor.window$Area.1 <=2 & cor.window$Area.2 >=3 &  cor.window$Area.2 <=5 ] <- "EVOS-Control"
    cor.window$tag.2[cor.window$Area.2 >=6 & cor.window$Area.1 >=3 &  cor.window$Area.1 <=5 ] <- "EVOS-Control"
    cor.window$tag.2[cor.window$Area.1 <=2 & cor.window$Area.2 > 5 ]  <- "Control"
    cor.window$tag.2[cor.window$Area.1 <=2 & cor.window$Area.2 <= 2 ]  <- "Control"
    cor.window$tag.2[cor.window$Area.1 >=6 & cor.window$Area.2 >= 6 ] <- "Control"
    cor.wind.plot <- cor.window[cor.window$Area.1 != cor.window$Area.2, ]
    
    cor.wind.agg <- aggregate(cor.wind.plot$Pearson, by=list(Year=cor.wind.plot$Year,tag.1=cor.wind.plot$tag.1),median)
    colnames(cor.wind.agg)[3] <- "Pearson.median"
    cor.wind.agg2 <- aggregate(cor.wind.plot$Pearson, by=list(Year=cor.wind.plot$Year,tag.1=cor.wind.plot$tag.1),mean)
    colnames(cor.wind.agg2)[3] <- "Pearson.mean"
    cor.tag.1 <- merge(cor.wind.agg,cor.wind.agg2)
    
    cor.wind.agg <- aggregate(cor.wind.plot$Pearson, by=list(Year=cor.wind.plot$Year,tag.2=cor.wind.plot$tag.2),median)
    colnames(cor.wind.agg)[3] <- "Pearson.median"
    cor.wind.agg2 <- aggregate(cor.wind.plot$Pearson, by=list(Year=cor.wind.plot$Year,tag.2=cor.wind.plot$tag.2),mean)
    colnames(cor.wind.agg2)[3] <- "Pearson.mean"
    cor.tag.2 <- merge(cor.wind.agg,cor.wind.agg2)
    
    ###
    p.temp <- ggplot(cor.wind.plot,
                     aes(x = Year, y = Pearson, colour=factor(tag.1))) +
      geom_jitter() 
    
    p.temp
    
    # scale_fill_gradientn(colours = c("blue",grey(0.8),"red"),limits=c(-1,1),na.value="black") +
    # coord_equal() +
    # scale_x_discrete(breaks=1:11)+
    # scale_y_discrete(breaks=1:11,limits = rev(levels(cor.mat.plot$Area.1))) +
    # labs(x= "Area", y= "Area") +
    # ggtitle(paste(start.year[i],"to",end.year[i])) +
    # theme_bw()
    
    p.line <- ggplot(cor.tag.1,
                     aes(x = Year, y = Pearson.median, colour=factor(tag.1))) +
      geom_jitter(data=cor.wind.plot,
                  aes(x = Year, y = Pearson, colour=factor(tag.1)),width=0.5,height=0,alpha=0.3) +
      geom_line(size=1.2) +
      geom_point(size=3) +
      geom_hline(yintercept=0,linetype=3)+
      scale_y_continuous(limits=c(-0.75,1)) +
      theme_bw()
    p.line
    
    p.line.mean <- ggplot(cor.tag.1,
                          aes(x = Year, y = Pearson.mean, colour=factor(tag.1))) +
      geom_jitter(data=cor.wind.plot,
                  aes(x = Year, y = Pearson, colour=factor(tag.1)),width=0.5,height=0,alpha=0.3) +
      geom_line(size=1.2) +
      geom_point(size=3) +
      geom_hline(yintercept=0,linetype=3)+
      scale_y_continuous(limits=c(-0.75,1)) +
      theme_bw()
    
    
    quartz(file=paste("./Scripts and plots for Pubs/Pearson time series, Diet=",DIET[QQ],", window=",WINDOW,".pdf"),dpi=300,height=4,width=6,type="pdf")
    print(p.line)
    print(p.line.mean)
    dev.off()
    
    #### TAG.2 (more stringent EVOS classification)
    
    
    p.line <- ggplot(cor.tag.2,
                     aes(x = Year, y = Pearson.median, colour=factor(tag.2))) +
      geom_jitter(data=cor.wind.plot,
                  aes(x = Year, y = Pearson, colour=factor(tag.1)),width=0.5,height=0,alpha=0.3) +
      geom_line(size=1.2) +
      geom_point(size=3) +
      geom_hline(yintercept=0,linetype=3)+
      scale_y_continuous(limits=c(-0.75,1)) +
      theme_bw()
    p.line
    
    p.line.mean <- ggplot(cor.tag.2,
                          aes(x = Year, y = Pearson.mean, colour=factor(tag.2))) +
      geom_jitter(data=cor.wind.plot,
                  aes(x = Year, y = Pearson, colour=factor(tag.1)),width=0.5,height=0,alpha=0.3) +
      geom_line(size=1.2) +
      geom_point(size=3) +
      geom_hline(yintercept=0,linetype=3)+
      scale_y_continuous(limits=c(-0.75,1)) +
      theme_bw()
    
    
    quartz(file=paste("./Scripts and plots for Pubs/Pearson time series (tag 2), Diet=",DIET[QQ],", window=",WINDOW,".pdf"),dpi=300,height=4,width=6,type="pdf")
    print(p.line)
    print(p.line.mean)
    dev.off()
  }
  
  
  ###################################################################################################
  
  #RECRUIT LAG
  RECRUIT <- unique(meanCPUE$age_class)
  
  for(QQ in 1:3){
    dat   <- aggregate(meanCPUE[meanCPUE$age_class==RECUIT[QQ],c("Mean.totalDensity","vari")],
                       by=list(Area=meanCPUE$Area[meanCPUE$age_class==RECRUIT[QQ]],Year=meanCPUE$year[meanCPUE$age_class==RECRUIT[QQ]]),sum)
    
    dat.wide <- reshape(dat,timevar=c("Area"),idvar="Year",direction="wide",drop="vari")
    #########################################################
    #########################################################
    # all pairwise correlations (pearson)
    n.Area <- 11
    cor.mat <- NULL
    for(i in 1:n.Area){
      for(j in 1:n.Area){
        A <- cor.test(dat.wide[,i+1],dat.wide[,j+1],method="pearson")
        cor.mat   <- rbind(cor.mat,data.frame(Area.1 = i, Area.2=j, Year = mean.year[k], Pearson = A$estimate))
        
      }
    }
    cor.mat <- cor.mat[cor.mat$Area.1 <= cor.mat$Area.2,]
    
    # Moving window approach
    # Use six year windows (will incorporate varying numbers of observations due to )
    
    start.year <- sort(unique(dat$Year))
    start.year <- start.year[start.year <= (max(start.year)-WINDOW)]
    end.year   <- start.year + WINDOW
    
    mean.year <- NULL
    for(k in 1:length(start.year) ){
      mean.year[k] <- dat.wide$Year[dat.wide$Year >= start.year[k] & dat.wide$Year <= end.year[k]] + WINDOW / 2
    }
    
    cor.window  <- NULL
    for(k in 1:length(start.year) ){
      temp        <- matrix(0,n.Area,n.Area)
      for(i in 1:n.Area){
        for(j in 1:n.Area){
          A <-  cor.test(dat.wide[dat.wide$Year >= start.year[k] & dat.wide$Year <= end.year[k] ,i+1],
                         dat.wide[dat.wide$Year >= start.year[k] & dat.wide$Year <= end.year[k],j+1],
                         method="pearson")
          cor.window   <- rbind(cor.window,data.frame(Area.1 = i, Area.2=j, Year = mean.year[k], Pearson = A$estimate))
        }
      }
      # colnames(temp) <- 
      # temp <- data.frame(Area=1:n.Area, year= rep(mean.year[k],n.Area),temp)
      # temp$Area <- (1:n.Area)
      # temp$year <- mean.year[k]
      # cor.window <- rbind(cor.window,temp)  
    }
    cor.window.trim <- cor.window[cor.window$Area.1 <= cor.window$Area.2,]
    
    #####################################
    #####################################
    #####################################
    ## MAKE SOME PLOTS OF THE CORRELATION MATRIX
    #####################################
    #####################################
    #####################################
    cor.mat.plot <- cor.mat
    cor.mat.plot$Area.1 <- factor(cor.mat.plot$Area.1)
    cor.mat.plot$Area.2 <- factor(cor.mat.plot$Area.2)
    cor.mat.plot$Pearson[cor.mat.plot$Pearson > 0.99999999] <- -99
    ### MAKE A PLOT OF 4 grids (true, mean est, sd, difference.)              
    p.all <- ggplot(cor.mat.plot,
                    aes(x = Area.1, y = Area.2, fill = Pearson)) +
      geom_tile() +
      scale_fill_gradientn(colours = c("blue",grey(0.8),"red"),limits=c(-1,1),na.value="black") +
      coord_equal() +
      scale_x_discrete(breaks=1:11)+
      scale_y_discrete(breaks=1:11,limits = rev(levels(cor.mat.plot$Area.1))) +
      labs(x= "Area", y= "Area") +
      ggtitle(paste("Recruit",RECRUIT[QQ],"1984-2015")) +
      theme_bw()
    
    p.all
    ##### YEAR BY YEAR
    
    cor.wind.plot <- cor.window.trim
    cor.wind.plot$Area.1 <- factor(cor.wind.plot$Area.1)
    cor.wind.plot$Area.2 <- factor(cor.wind.plot$Area.2)
    cor.wind.plot$Pearson[cor.wind.plot$Pearson > 0.999999] <- -99
    
    P <- list()
    
    for(i in 1:length(mean.year)){
      
      p.temp <- ggplot(cor.wind.plot[cor.wind.plot$Year == mean.year[i],],
                       aes(x = Area.1, y = Area.2, fill = Pearson)) +
        geom_tile() +
        scale_fill_gradientn(colours = c("blue",grey(0.8),"red"),limits=c(-1,1),na.value="black") +
        coord_equal() +
        scale_x_discrete(breaks=1:11)+
        scale_y_discrete(breaks=1:11,limits = rev(levels(cor.mat.plot$Area.1))) +
        labs(x= "Area", y= "Area") +
        ggtitle(paste("Recruit",RECRUIT[QQ],start.year[i],"to",end.year[i])) +
        theme_bw()
      
      P[[i]] <- p.temp
    }
    
    
    quartz(file=paste("./Scripts and plots for Pubs/Pearson cor matrices, Diet=",DIET[QQ],", window=",WINDOW,".pdf"),dpi=300,height=7,width=7,type="pdf")
    print(p.all)
    for(i in 1:length(mean.year)){
      print(P[[i]])
    }
    dev.off()
    
    ##### 
    # Add indicator variable for EVOS - EVOS, EVOS-non, and non-non combinations of 
    
    
    cor.window$tag.1[cor.window$Area.1 >=2 & cor.window$Area.2 <=6 ]  <- "EVOS"
    cor.window$tag.1[cor.window$Area.1 ==1 & cor.window$Area.2 >=2 &  cor.window$Area.2 <=6 ] <- "EVOS-Control"
    cor.window$tag.1[cor.window$Area.2 >=7 & cor.window$Area.1 >=2 &  cor.window$Area.1 <=6 ] <- "EVOS-Control"
    cor.window$tag.1[cor.window$Area.1 ==1 & cor.window$Area.2 > 6 ]  <- "Control"
    cor.window$tag.1[cor.window$Area.1 >=7 & cor.window$Area.2 >= 7 ] <- "Control"
    
    cor.window$tag.2[cor.window$Area.1 >=3 & cor.window$Area.2 <=5 ]  <- "EVOS"
    cor.window$tag.2[cor.window$Area.1 <=2 & cor.window$Area.2 >=3 &  cor.window$Area.2 <=5 ] <- "EVOS-Control"
    cor.window$tag.2[cor.window$Area.2 >=6 & cor.window$Area.1 >=3 &  cor.window$Area.1 <=5 ] <- "EVOS-Control"
    cor.window$tag.2[cor.window$Area.1 <=2 & cor.window$Area.2 > 5 ]  <- "Control"
    cor.window$tag.2[cor.window$Area.1 <=2 & cor.window$Area.2 <= 2 ]  <- "Control"
    cor.window$tag.2[cor.window$Area.1 >=6 & cor.window$Area.2 >= 6 ] <- "Control"
    cor.wind.plot <- cor.window[cor.window$Area.1 != cor.window$Area.2, ]
    
    cor.wind.agg <- aggregate(cor.wind.plot$Pearson, by=list(Year=cor.wind.plot$Year,tag.1=cor.wind.plot$tag.1),median)
    colnames(cor.wind.agg)[3] <- "Pearson.median"
    cor.wind.agg2 <- aggregate(cor.wind.plot$Pearson, by=list(Year=cor.wind.plot$Year,tag.1=cor.wind.plot$tag.1),mean)
    colnames(cor.wind.agg2)[3] <- "Pearson.mean"
    cor.tag.1 <- merge(cor.wind.agg,cor.wind.agg2)
    
    cor.wind.agg <- aggregate(cor.wind.plot$Pearson, by=list(Year=cor.wind.plot$Year,tag.2=cor.wind.plot$tag.2),median)
    colnames(cor.wind.agg)[3] <- "Pearson.median"
    cor.wind.agg2 <- aggregate(cor.wind.plot$Pearson, by=list(Year=cor.wind.plot$Year,tag.2=cor.wind.plot$tag.2),mean)
    colnames(cor.wind.agg2)[3] <- "Pearson.mean"
    cor.tag.2 <- merge(cor.wind.agg,cor.wind.agg2)
    
    ###
    p.temp <- ggplot(cor.wind.plot,
                     aes(x = Year, y = Pearson, colour=factor(tag.1))) +
      geom_jitter() 
    
    p.temp
    
    # scale_fill_gradientn(colours = c("blue",grey(0.8),"red"),limits=c(-1,1),na.value="black") +
    # coord_equal() +
    # scale_x_discrete(breaks=1:11)+
    # scale_y_discrete(breaks=1:11,limits = rev(levels(cor.mat.plot$Area.1))) +
    # labs(x= "Area", y= "Area") +
    # ggtitle(paste(start.year[i],"to",end.year[i])) +
    # theme_bw()
    
    p.line <- ggplot(cor.tag.1,
                     aes(x = Year, y = Pearson.median, colour=factor(tag.1))) +
      geom_jitter(data=cor.wind.plot,
                  aes(x = Year, y = Pearson, colour=factor(tag.1)),width=0.5,height=0,alpha=0.3) +
      geom_line(size=1.2) +
      geom_point(size=3) +
      geom_hline(yintercept=0,linetype=3)+
      scale_y_continuous(limits=c(-0.75,1)) +
      theme_bw()
    p.line
    
    p.line.mean <- ggplot(cor.tag.1,
                          aes(x = Year, y = Pearson.mean, colour=factor(tag.1))) +
      geom_jitter(data=cor.wind.plot,
                  aes(x = Year, y = Pearson, colour=factor(tag.1)),width=0.5,height=0,alpha=0.3) +
      geom_line(size=1.2) +
      geom_point(size=3) +
      geom_hline(yintercept=0,linetype=3)+
      scale_y_continuous(limits=c(-0.75,1)) +
      theme_bw()
    
    
    quartz(file=paste("./Scripts and plots for Pubs/Pearson time series, Recruit=",RECRUIT[QQ],", window=",WINDOW,".pdf"),dpi=300,height=4,width=6,type="pdf")
    print(p.line)
    print(p.line.mean)
    dev.off()
    
    #### TAG.2 (more stringent EVOS classification)
    
    
    p.line <- ggplot(cor.tag.2,
                     aes(x = Year, y = Pearson.median, colour=factor(tag.2))) +
      geom_jitter(data=cor.wind.plot,
                  aes(x = Year, y = Pearson, colour=factor(tag.1)),width=0.5,height=0,alpha=0.3) +
      geom_line(size=1.2) +
      geom_point(size=3) +
      geom_hline(yintercept=0,linetype=3)+
      scale_y_continuous(limits=c(-0.75,1)) +
      theme_bw()
    p.line
    
    p.line.mean <- ggplot(cor.tag.2,
                          aes(x = Year, y = Pearson.mean, colour=factor(tag.2))) +
      geom_jitter(data=cor.wind.plot,
                  aes(x = Year, y = Pearson, colour=factor(tag.1)),width=0.5,height=0,alpha=0.3) +
      geom_line(size=1.2) +
      geom_point(size=3) +
      geom_hline(yintercept=0,linetype=3)+
      scale_y_continuous(limits=c(-0.75,1)) +
      theme_bw()
    
    
    quartz(file=paste("./Scripts and plots for Pubs/Pearson time series (tag 2), Recruit=",RECRUIT[QQ],", window=",WINDOW,".pdf"),dpi=300,height=4,width=6,type="pdf")
    print(p.line)
    print(p.line.mean)
    dev.off()
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
#   
#   
#   
#   
#   p.line <- ggplot(cor.tag.2,
#                    aes(x = Year, y = Pearson.median, colour=factor(tag.2))) +
#     geom_line() +
#     scale_y_continuous(limits=c(-1,1)) 
#   p.line
#   
# #### Trim to only include pairwise correlations up to 3 steps away
#   cor.wind.plot <- cor.wind.plot[cor.wind.plot$Area.1 <= (cor.wind.plot$Area.2 + 3), ]
#   
#   
#   cor.wind.agg <- aggregate(cor.wind.plot$Pearson, by=list(Year=cor.wind.plot$Year,tag.1=cor.wind.plot$tag.1),median)
#   colnames(cor.wind.agg)[3] <- "Pearson.median"
#   cor.wind.agg2 <- aggregate(cor.wind.plot$Pearson, by=list(Year=cor.wind.plot$Year,tag.1=cor.wind.plot$tag.1),mean)
#   colnames(cor.wind.agg2)[3] <- "Pearson.mean"
#   cor.tag.1.simp <- merge(cor.wind.agg,cor.wind.agg2)
#   
#   cor.wind.agg <- aggregate(cor.wind.plot$Pearson, by=list(Year=cor.wind.plot$Year,tag.2=cor.wind.plot$tag.2),median)
#   colnames(cor.wind.agg)[3] <- "Pearson.median"
#   cor.wind.agg2 <- aggregate(cor.wind.plot$Pearson, by=list(Year=cor.wind.plot$Year,tag.2=cor.wind.plot$tag.2),mean)
#   colnames(cor.wind.agg2)[3] <- "Pearson.mean"
#   cor.tag.2.simp <- merge(cor.wind.agg,cor.wind.agg2)
#   
#   
#   p.temp <- ggplot(cor.wind.plot,
#                    aes(x = Year, y = Pearson, colour=factor(tag.1))) +
#     geom_jitter() 
#   p.temp
#  
#   p.temp <- ggplot(cor.wind.plot,
#                    aes(x = Year, y = Pearson, colour=factor(tag.2))) +
#     geom_jitter() 
#   p.temp
#   
#    
#   
#   p.line <- ggplot(cor.tag.1.simp,
#                    aes(x = Year, y = Pearson.mean, colour=factor(tag.1))) +
#     geom_line() +
#     geom_point() +
#     geom_jitter(cor.wind.plot,
#                aes(x = Year, y = Pearson, colour=factor(tag.1))) +
#     geom_hline(yintercept=0,linetype=3)+
#     scale_y_continuous(limits=c(-0.5,1)) +
#     theme_bw()
#   p.line
#   
#   p.line <- ggplot(cor.tag.2.simp,
#                    aes(x = Year, y = Pearson.mean, colour=factor(tag.2))) +
#     geom_line() +
#     geom_hline(yintercept=0,linetype=3)+
#     scale_y_continuous(limits=c(-0.5,1)) +
#     theme_bw()
#   p.line
##############################################################
##############################################################
##############################################################
### Calculate community synchrony.. (across Area )
##############################################################
##############################################################
##############################################################
# 
#   # FIRST Calculate the community synchrony for the entire time series
#   community.sync(data = dat.wide[,2:ncol(dat.wide)], nrands = 999, type=1, quiet= TRUE)
# 
#   community.sync(data = dat.wide[,3:7], nrands = 999, type=1, quiet= TRUE)
#   
#   community.sync(data = dat.wide[,c(1,8:12)], nrands = 999, type=1, quiet= TRUE)
#   Int <- 2
#   start.loc <- 1:9
#   Out <- NULL
#   for(i in 1:length(start.loc)){
#     temp <- dat.wide[,(1+start.loc[i]):(1+start.loc[i]+Int)]
#     A <- community.sync(data = temp, nrands = 999, type=1, quiet= TRUE)
#     Out  <- rbind(Out, data.frame(start.loc = start.loc[i],stop.loc=start.loc[i]+Int,A$obs,A$meancorr,A$pval))
#   }
#   
#   ## Calculate the same community synchrony metrics using a moving window.
#   
#   WIN <- 9
#   Int <- 2
#   start.loc <- 1:9
#   Out.window <- NULL
#   
#   start.year <- sort(unique(dat$Year))
#   start.year <- start.year[start.year <= (max(start.year)-WINDOW)]
#   end.year   <- start.year + WINDOW
#   
#   for(j in 1:length(start.year)){
#     for(i in 1:length(start.loc)){
#       temp        <- dat.wide[dat.wide$Year>=start.year[j] & dat.wide$Year<=end.year[j],(1+start.loc[i]):(1+start.loc[i]+Int)]
#       A           <- community.sync(data = temp, nrands = 999, type=1, quiet= TRUE)
#       Out.window  <- rbind(Out.window, data.frame(start.year=start.year[j],end.year=end.year[j],start.loc = start.loc[i],stop.loc=start.loc[i]+Int,A$obs,A$meancorr,A$pval))
#     }
#   }
#   
#   ############# Synchronity for community in each area
#   
#   trim <- meanCPUE[,c("Species","area","Mean.totalDensity","year")]
#   trim <- spread(data= trim,key=Species,value=Mean.totalDensity)
# 
#   Out.comm <- NULL
#   for(i in 1:n.Area){
#     temp <- trim[trim$area == i,3:ncol(trim)]
#     A <-    community.sync(data = temp, nrands = 999, type=1, quiet= TRUE)
#     Out.comm  <- rbind(Out.comm, data.frame(Area = i,A$obs,A$meancorr,A$pval))
#   }


