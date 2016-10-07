rm(list=ls())

## Plot for publication is about line 265

### Calculations of synchrony 
library(synchrony)
library(ggplot2)
library(maps)
library(mapproj)
library(dplyr)
library(tidyr)
library(reshape2)

setwd("/Users/ole.shelton/GitHub/pfx-groundfish/")

#### GO GET multiplot script.  #### STOPPED HERE.
source("/Users/ole.shelton/GitHub/pfx-groundfish/Scripts and plots for Pubs/multiplot.r")

# Go get the needed data:
meanCPUE <- read.csv("/Users/ole.shelton/GitHub/pfx-groundfish/All_sp_index_meanCPUEByArea.csv") # load data
meanCPUE <- subset(meanCPUE,area!='Total') #remove total field
meanCPUE <- droplevels(meanCPUE) 
meanCPUE$area <- as.numeric(levels(meanCPUE$area))[meanCPUE$area]#use numeric/factor conventions
meanCPUE$Area <- factor(meanCPUE$area)
meanCPUE$vari <- meanCPUE$SD.totalDensity^2 #calculate variance
#meanCPUE$Species <- gsub("[.]","",meanCPUE$Species) #quit changing how names are spelled damnit!

# Read in the area identification file
discrete_areas <- read.csv("./goaTrawl/Output Data/goa_discrete_areas_for_comparison(50_to_150m).csv")
#head(discrete_areas)

# Read in the control file for species that identifies functional categories for each species
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
# Moving window approach
# Use nine year windows (will incorporate varying numbers of observations due to spacing of surveys)
WINDOW <- 9

start.year <- sort(unique(dat$Year))
start.year <- start.year[start.year <= (max(start.year)-WINDOW)]
start.year <- c(start.year,2007)
end.year   <- start.year + WINDOW


mean.year <- NULL
for(k in 1:length(start.year) ){
  mean.year[k] <- c(dat.wide$Year[dat.wide$Year >= start.year[k] & dat.wide$Year <= end.year[k]] + WINDOW / 2)[1]
}

cor.window  <- NULL
for(k in 1:length(start.year) ){
  temp        <- matrix(0,n.Area,n.Area)
  for(i in 1:n.Area){
    for(j in 1:n.Area){
      A <-  cor.test(dat.wide[dat.wide$Year >= start.year[k] & dat.wide$Year <= end.year[k] ,i+1],
                     dat.wide[dat.wide$Year >= start.year[k] & dat.wide$Year <= end.year[k],j+1],
                     method="pearson")
      cor.window   <- rbind(cor.window,data.frame(Area.1 = i, Area.2=j, Year = mean.year[k], Pearson = A$estimate,p.val=A$p.value))
    }
  }
  # colnames(temp) <- 
  # temp <- data.frame(Area=1:n.Area, year= rep(mean.year[k],n.Area),temp)
  # temp$Area <- (1:n.Area)
  # temp$year <- mean.year[k]
  # cor.window <- rbind(cor.window,temp)  
}
cor.window.trim <- cor.window[cor.window$Area.1 <= cor.window$Area.2,]


#######################################################
#######################################################
### STOP DATA PROCESSING
#######################################################
#######################################################



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
cor.wind.plot.tile <- cor.wind.plot
cor.wind.plot <- cor.window[cor.window$Area.1 != cor.window$Area.2, ]
cor.wind.plot <- cor.wind.plot[cor.wind.plot$Area.1 < cor.wind.plot$Area.2,]

cor.wind.agg <- aggregate(cor.wind.plot$Pearson, by=list(Year=cor.wind.plot$Year,tag.1=cor.wind.plot$tag.1),median)
colnames(cor.wind.agg)[3] <- "Pearson.median"
cor.wind.agg2 <- aggregate(cor.wind.plot$Pearson, by=list(Year=cor.wind.plot$Year,tag.1=cor.wind.plot$tag.1),mean)
colnames(cor.wind.agg2)[3] <- "Pearson.mean"
cor.wind.agg3 <- aggregate(cor.wind.plot$Pearson, by=list(Year=cor.wind.plot$Year,tag.1=cor.wind.plot$tag.1),sd)
colnames(cor.wind.agg3)[3] <- "Pearson.sd"

A <- aggregate(cor.wind.plot$tag.1,by=list(tag.1 = cor.wind.plot$tag.1,year=cor.wind.plot$Year),length)
A <- A[1:3,c("tag.1","x")]

cor.tag.1 <- merge(cor.wind.agg,cor.wind.agg2)
cor.tag.1 <- merge(cor.tag.1,cor.wind.agg3)
cor.tag.1$Pearson.N <- A$x[match(cor.tag.1$tag.1,A$tag.1)]
cor.tag.1$SE <- cor.tag.1$Pearson.sd / sqrt(cor.tag.1$Pearson.N)

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

# Define some nice colors
COLS <- c(  "EVOS"="#1b9e77",
            "Control"="#d95f02",
            "EVOS-Control"  = "#7570b3")

#scale_fill_manual(values=COLS,name="Origin") 
###############################################################################################
###############################################################################################
###############################################################################################
#####     ##########################################################################################
###############################################################################################
############################# HERE IS A PLOT FOR PUBLICATION.
###############################################################################################  
cor.tag.1$tag.1 <- as.factor(cor.tag.1$tag.1)

  p.line <- ggplot(data=cor.tag.1,aes(x = Year, y = Pearson.mean)) + 
            geom_ribbon(aes(x=Year, ymin = Pearson.mean - SE, ymax = Pearson.mean + SE, fill=tag.1),alpha=0.3) +
            geom_line(  aes(x = Year, y = Pearson.mean, colour=tag.1),size=1) +
            geom_point( aes(x = Year, y = Pearson.mean, colour=tag.1,shape = tag.1),size=3) +
            geom_hline(yintercept=0,linetype=3)+
            labs(x = "Year", y = "Pearson Correlation") +
            scale_y_continuous(limits=c(-0.75,1)) +
            scale_color_manual(values=COLS,name="") +
            scale_fill_manual(values=COLS,name="") +
            scale_shape_manual(values=c(19,17,15),name="")+
            ggtitle("a)")+
            theme_bw() +
            theme(legend.position = c(0.645, 1.035),
                      legend.justification = c(0, 1),
                      legend.text = element_text(size = 8),
                      legend.key.size=unit(0.5,"cm"),
                      legend.title = element_text(size=0),
                      legend.background=element_rect(colour = "black"),
                      panel.border=element_rect(size=1.5,color="black"),
                      plot.title=element_text(hjust=0,size=12))
  p.line
  
  # Pick the Peak year of the correlation matrix to match with the time-series. 
  P.1990.99 <- ggplot(cor.wind.plot.tile[cor.wind.plot.tile$Year == 1994.5,],
                     aes(x = Area.1, y = Area.2, fill = Pearson)) +
      geom_tile() +
      scale_fill_gradientn(colours = c("blue",grey(0.8),"red"),limits=c(-1,1),na.value="black") +
      coord_equal() +
      scale_x_discrete(breaks=1:11)+
      scale_y_discrete(breaks=1:11,limits = rev(levels(cor.mat.plot$Area.1))) +
      labs(x= "Area", y= "Area") +
      ggtitle("b)")+
      theme_bw() +
        theme(legend.position = c(0.7, 1.0),
          legend.justification = c(0, 1),
          legend.text = element_text(size = 8),
          legend.background=element_rect(colour = "black"),
          legend.title = element_text(size=0),
          panel.border=element_rect(size=1.5,color="black"),
          plot.title=element_text(hjust=0,size=12))
  
  P.all <- ggplot(cor.mat.plot,
                      aes(x = Area.1, y = Area.2, fill = Pearson)) +
    geom_tile() +
    scale_fill_gradientn(colours = c("blue",grey(0.8),"red"),limits=c(-1,1),na.value="black") +
    coord_equal() +
    scale_x_discrete(breaks=1:11)+
    scale_y_discrete(breaks=1:11,limits = rev(levels(cor.mat.plot$Area.1))) +
    labs(x= "Area", y= "Area") +
    ggtitle("b)")+
    theme_bw() +
    theme(legend.position = c(0.7, 1.0),
          legend.justification = c(0, 1),
          legend.text = element_text(size = 8),
          legend.background=element_rect(colour = "black"),
          legend.title = element_text(size=0),
          panel.border=element_rect(size=1.5,color="black"),
          plot.title=element_text(hjust=0,size=12))
  
  
  #P.1990.99
  if(WINDOW==9){
    quartz(file=paste("./Scripts and plots for Pubs/Pearson Matrix + Time-series,Biomass,window=",WINDOW,".jpeg"),dpi=300,height=8,width=4,type="jpeg")
        Layout= matrix(c(1,2),nrow=2,ncol=1,byrow=T)
        QQ <- list(p.line,P.1990.99)
        multiplot(plotlist=QQ ,layout= Layout)
    dev.off()
  }
  #P.all
  if(WINDOW==9){
    quartz(file=paste("./Scripts and plots for Pubs/Pearson Matrix + Time-series,Biomass,window=",WINDOW,"full time-series matrix.jpeg"),dpi=300,height=8,width=4,type="jpeg")
    Layout= matrix(c(1,2),nrow=2,ncol=1,byrow=T)
    QQ <- list(p.line,P.all)
    multiplot(plotlist=QQ ,layout= Layout)
    dev.off()
  }
  
  
    ###############################################################################################
  ###############################################################################################
  ###############################################################################################
  ###############################################################################################
  ###############################################################################################
  ###############################################################################################
  ###############################################################################################
  #######################################

    
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
  dev.off()
  
  #### TAG.2 (more stringent EVOS classification)
  
  # 
  # p.line <- ggplot(cor.tag.2,
  #                  aes(x = Year, y = Pearson.median, colour=factor(tag.2))) +
  #    geom_jitter(data=cor.wind.plot,
  #                aes(x = Year, y = Pearson, colour=factor(tag.1)),width=0.5,height=0,alpha=0.3) +
  #   geom_line(size=1.2) +
  #   geom_point(size=3) +
  #   geom_hline(yintercept=0,linetype=3)+
  #   scale_y_continuous(limits=c(-0.75,1)) +
  #   theme_bw()
  # p.line
  # 
  # p.line.mean <- ggplot(cor.tag.2,
  #                       aes(x = Year, y = Pearson.mean, colour=factor(tag.2))) +
  #   geom_jitter(data=cor.wind.plot,
  #                aes(x = Year, y = Pearson, colour=factor(tag.1)),width=0.5,height=0,alpha=0.3) +
  #   geom_line(size=1.2) +
  #   geom_point(size=3) +
  #   geom_hline(yintercept=0,linetype=3)+
  #   scale_y_continuous(limits=c(-0.75,1)) +
  #   theme_bw()
  # 
  # setwd("./Scripts and plots for Pubs/")
  # quartz(file=paste("Pearson time series, Biomass (Tag 2),window=",WINDOW,".pdf"),dpi=300,height=4,width=6,type="pdf")
  #   print(p.line)
  #   print(p.line.mean)
  # dev.off()
  
  ###################################################################################################
  ###################################################################################################
  ###################################################################################################
  ###################################################################################################
  ###################################################################################################
  #### REPEAT FOR EACH OF THE FUNCITONAL GROUPS (Guild, diet, Recruit)
  ###################################################################################################
  ###################################################################################################
  ###################################################################################################
  ###################################################################################################
  ###################################################################################################
  
  #Guilds  
  GUILD <- unique(meanCPUE$guild)
  
  for(QQ in 1:3){
    dat   <- aggregate(meanCPUE[meanCPUE$fish.invert=="fish" & meanCPUE$guild==GUILD[QQ],c("Mean.totalDensity","vari")],
                     by=list(Area=meanCPUE$Area[meanCPUE$fish.invert=="fish" & meanCPUE$guild==GUILD[QQ]],
                             Year=meanCPUE$year[meanCPUE$fish.invert=="fish" & meanCPUE$guild==GUILD[QQ]]),sum)
  
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
    mean.year[k] <- c(dat.wide$Year[dat.wide$Year >= start.year[k] & dat.wide$Year <= end.year[k]] + WINDOW / 2)[1]
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
  
  cor.wind.plot <- cor.window[cor.window$Area.1 != cor.window$Area.2, ]
  cor.wind.plot <- cor.wind.plot[cor.wind.plot$Area.1 < cor.wind.plot$Area.2,]
  
  cor.wind.agg <- aggregate(cor.wind.plot$Pearson, by=list(Year=cor.wind.plot$Year,tag.1=cor.wind.plot$tag.1),median)
  colnames(cor.wind.agg)[3] <- "Pearson.median"
  cor.wind.agg2 <- aggregate(cor.wind.plot$Pearson, by=list(Year=cor.wind.plot$Year,tag.1=cor.wind.plot$tag.1),mean)
  colnames(cor.wind.agg2)[3] <- "Pearson.mean"
  cor.wind.agg3 <- aggregate(cor.wind.plot$Pearson, by=list(Year=cor.wind.plot$Year,tag.1=cor.wind.plot$tag.1),sd)
  colnames(cor.wind.agg3)[3] <- "Pearson.sd"
  
  A <- aggregate(cor.wind.plot$tag.1,by=list(tag.1 = cor.wind.plot$tag.1,year=cor.wind.plot$Year),length)
  A <- A[1:3,c("tag.1","x")]
  
  cor.tag.1 <- merge(cor.wind.agg,cor.wind.agg2)
  cor.tag.1 <- merge(cor.tag.1,cor.wind.agg3)
  cor.tag.1$Pearson.N <- A$x[match(cor.tag.1$tag.1,A$tag.1)]
  cor.tag.1$SE <- cor.tag.1$Pearson.sd / sqrt(cor.tag.1$Pearson.N)
  
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
  
  p.line<- ggplot(data=cor.tag.1,aes(x = Year, y = Pearson.mean)) + 
            geom_ribbon(aes(x=Year, ymin = Pearson.mean - SE, ymax = Pearson.mean + SE, fill=tag.1),alpha=0.3) +
            geom_line(  aes(x = Year, y = Pearson.mean, colour=tag.1),size=1.1) +
            geom_point( aes(x = Year, y = Pearson.mean, colour=tag.1,shape = tag.1),size=3) +
            geom_hline(yintercept=0,linetype=3)+
            labs(x = "Year", y = "Pearson Correlation") +
            scale_y_continuous(limits=c(-0.75,1)) +
            scale_color_manual(values=COLS,name="") +
            scale_fill_manual(values=COLS,name="") +
            scale_shape_manual(values=c(19,17,15),name="")+
            theme_bw() +
            theme(legend.position = c(0.63, 1.02),
              legend.justification = c(0, 1),
              legend.text = element_text(size = 8),
              legend.key.size=unit(0.5,"cm"),
              legend.title = element_text(size=0),
              legend.background=element_rect(colour = "black"),
              panel.border=element_rect(size=1.5,color="black"))
  p.line
  
  quartz(file=paste("./Scripts and plots for Pubs/Pearson time series, Guild=",GUILD[QQ],", window=",WINDOW,".pdf"),dpi=300,height=4,width=6,type="pdf")
    print(p.line)
  dev.off()
  
  #### TAG.2 (more stringent EVOS classification)
  
  # 
  # p.line <- ggplot(cor.tag.2,
  #                  aes(x = Year, y = Pearson.median, colour=factor(tag.2))) +
  #   geom_jitter(data=cor.wind.plot,
  #               aes(x = Year, y = Pearson, colour=factor(tag.1)),width=0.5,height=0,alpha=0.3) +
  #   geom_line(size=1.2) +
  #   geom_point(size=3) +
  #   geom_hline(yintercept=0,linetype=3)+
  #   scale_y_continuous(limits=c(-0.75,1)) +
  #   theme_bw()
  # p.line
  # 
  # p.line.mean <- ggplot(cor.tag.2,
  #                       aes(x = Year, y = Pearson.mean, colour=factor(tag.2))) +
  #   geom_jitter(data=cor.wind.plot,
  #               aes(x = Year, y = Pearson, colour=factor(tag.1)),width=0.5,height=0,alpha=0.3) +
  #   geom_line(size=1.2) +
  #   geom_point(size=3) +
  #   geom_hline(yintercept=0,linetype=3)+
  #   scale_y_continuous(limits=c(-0.75,1)) +
  #   theme_bw()
  # 
  # 
  # quartz(file=paste("./Scripts and plots for Pubs/Pearson time series (tag 2), Guild=",GUILD[QQ],", window=",WINDOW,".pdf"),dpi=300,height=4,width=6,type="pdf")
  # print(p.line)
  # print(p.line.mean)
  # dev.off()
  # 
  }
  
  ###################################################################################################
  ###################################################################################################
  ###################################################################################################
  
  #DIET  
  DIET <- unique(meanCPUE$diet1)
  
  for(QQ in 1:3){
    dat   <- aggregate(meanCPUE[meanCPUE$fish.invert=="fish" & meanCPUE$diet1==DIET[QQ],c("Mean.totalDensity","vari")],
                       by=list(Area=meanCPUE$Area[meanCPUE$fish.invert=="fish" & meanCPUE$diet1==DIET[QQ]],
                               Year=meanCPUE$year[meanCPUE$fish.invert=="fish" & meanCPUE$diet1==DIET[QQ]]),sum)
    
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
      mean.year[k] <- c(dat.wide$Year[dat.wide$Year >= start.year[k] & dat.wide$Year <= end.year[k]] + WINDOW / 2)[1]
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
    cor.wind.plot <- cor.wind.plot[cor.wind.plot$Area.1 < cor.wind.plot$Area.2,]
    
    cor.wind.agg <- aggregate(cor.wind.plot$Pearson, by=list(Year=cor.wind.plot$Year,tag.1=cor.wind.plot$tag.1),median)
    colnames(cor.wind.agg)[3] <- "Pearson.median"
    cor.wind.agg2 <- aggregate(cor.wind.plot$Pearson, by=list(Year=cor.wind.plot$Year,tag.1=cor.wind.plot$tag.1),mean)
    colnames(cor.wind.agg2)[3] <- "Pearson.mean"
    cor.wind.agg3 <- aggregate(cor.wind.plot$Pearson, by=list(Year=cor.wind.plot$Year,tag.1=cor.wind.plot$tag.1),sd)
    colnames(cor.wind.agg3)[3] <- "Pearson.sd"
    
    A <- aggregate(cor.wind.plot$tag.1,by=list(tag.1 = cor.wind.plot$tag.1,year=cor.wind.plot$Year),length)
    A <- A[1:3,c("tag.1","x")]
    
    cor.tag.1 <- merge(cor.wind.agg,cor.wind.agg2)
    cor.tag.1 <- merge(cor.tag.1,cor.wind.agg3)
    cor.tag.1$Pearson.N <- A$x[match(cor.tag.1$tag.1,A$tag.1)]
    cor.tag.1$SE <- cor.tag.1$Pearson.sd / sqrt(cor.tag.1$Pearson.N)
    
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
    
    p.line<- ggplot(data=cor.tag.1,aes(x = Year, y = Pearson.mean)) + 
      geom_ribbon(aes(x=Year, ymin = Pearson.mean - SE, ymax = Pearson.mean + SE, fill=tag.1),alpha=0.3) +
      geom_line(  aes(x = Year, y = Pearson.mean, colour=tag.1),size=1.1) +
      geom_point( aes(x = Year, y = Pearson.mean, colour=tag.1,shape = tag.1),size=3) +
      geom_hline(yintercept=0,linetype=3)+
      labs(x = "Year", y = "Pearson Correlation") +
      scale_y_continuous(limits=c(-0.75,1)) +
      scale_color_manual(values=COLS,name="") +
      scale_fill_manual(values=COLS,name="") +
      scale_shape_manual(values=c(19,17,15),name="")+
      theme_bw() +
      theme(legend.position = c(0.63, 1.02),
            legend.justification = c(0, 1),
            legend.text = element_text(size = 8),
            legend.key.size=unit(0.5,"cm"),
            legend.title = element_text(size=0),
            legend.background=element_rect(colour = "black"),
            panel.border=element_rect(size=1.5,color="black"))
    
    
    quartz(file=paste("./Scripts and plots for Pubs/Pearson time series, Diet=",DIET[QQ],", window=",WINDOW,".pdf"),dpi=300,height=4,width=6,type="pdf")
      print(p.line)
    dev.off()
    
    #### TAG.2 (more stringent EVOS classification)
    
  #   
  #   p.line <- ggplot(cor.tag.2,
  #                    aes(x = Year, y = Pearson.median, colour=factor(tag.2))) +
  #     geom_jitter(data=cor.wind.plot,
  #                 aes(x = Year, y = Pearson, colour=factor(tag.1)),width=0.5,height=0,alpha=0.3) +
  #     geom_line(size=1.2) +
  #     geom_point(size=3) +
  #     geom_hline(yintercept=0,linetype=3)+
  #     scale_y_continuous(limits=c(-0.75,1)) +
  #     theme_bw()
  #   p.line
  #   
  #   p.line.mean <- ggplot(cor.tag.2,
  #                         aes(x = Year, y = Pearson.mean, colour=factor(tag.2))) +
  #     geom_jitter(data=cor.wind.plot,
  #                 aes(x = Year, y = Pearson, colour=factor(tag.1)),width=0.5,height=0,alpha=0.3) +
  #     geom_line(size=1.2) +
  #     geom_point(size=3) +
  #     geom_hline(yintercept=0,linetype=3)+
  #     scale_y_continuous(limits=c(-0.75,1)) +
  #     theme_bw()
  #   
  #   
  #   quartz(file=paste("./Scripts and plots for Pubs/Pearson time series (tag 2), Diet=",DIET[QQ],", window=",WINDOW,".pdf"),dpi=300,height=4,width=6,type="pdf")
  #   print(p.line)
  #   print(p.line.mean)
  #   dev.off()
   }
   
  
  ###################################################################################################
  
  #RECRUIT LAG
  RECRUIT <- unique(meanCPUE$age_class)
  RECRUIT <- RECRUIT[is.na(RECRUIT)==F]
  for(QQ in 1:3){
    dat   <- aggregate(meanCPUE[meanCPUE$fish.invert=="fish" & meanCPUE$age_class==RECRUIT[QQ],c("Mean.totalDensity","vari")],
                       by=list(Area=meanCPUE$Area[meanCPUE$fish.invert=="fish" & meanCPUE$age_class==RECRUIT[QQ]],
                               Year=meanCPUE$year[meanCPUE$fish.invert=="fish" & meanCPUE$age_class==RECRUIT[QQ]]),sum)
    
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
      mean.year[k] <- c(dat.wide$Year[dat.wide$Year >= start.year[k] & dat.wide$Year <= end.year[k]] + WINDOW / 2)[1]
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
    
    
    quartz(file=paste("./Scripts and plots for Pubs/Pearson cor matrices, Recruit=",RECRUIT[QQ],", window=",WINDOW,".pdf"),dpi=300,height=7,width=7,type="pdf")
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
    cor.wind.plot <- cor.wind.plot[cor.wind.plot$Area.1 < cor.wind.plot$Area.2,]
    
    cor.wind.agg <- aggregate(cor.wind.plot$Pearson, by=list(Year=cor.wind.plot$Year,tag.1=cor.wind.plot$tag.1),median)
    colnames(cor.wind.agg)[3] <- "Pearson.median"
    cor.wind.agg2 <- aggregate(cor.wind.plot$Pearson, by=list(Year=cor.wind.plot$Year,tag.1=cor.wind.plot$tag.1),mean)
    colnames(cor.wind.agg2)[3] <- "Pearson.mean"
    cor.wind.agg3 <- aggregate(cor.wind.plot$Pearson, by=list(Year=cor.wind.plot$Year,tag.1=cor.wind.plot$tag.1),sd)
    colnames(cor.wind.agg3)[3] <- "Pearson.sd"
    
    A <- aggregate(cor.wind.plot$tag.1,by=list(tag.1 = cor.wind.plot$tag.1,year=cor.wind.plot$Year),length)
    A <- A[1:3,c("tag.1","x")]
    
    cor.tag.1 <- merge(cor.wind.agg,cor.wind.agg2)
    cor.tag.1 <- merge(cor.tag.1,cor.wind.agg3)
    cor.tag.1$Pearson.N <- A$x[match(cor.tag.1$tag.1,A$tag.1)]
    cor.tag.1$SE <- cor.tag.1$Pearson.sd / sqrt(cor.tag.1$Pearson.N)
    
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
    
    p.line<- ggplot(data=cor.tag.1,aes(x = Year, y = Pearson.mean)) + 
      geom_ribbon(aes(x=Year, ymin = Pearson.mean - SE, ymax = Pearson.mean + SE, fill=tag.1),alpha=0.3) +
      geom_line(  aes(x = Year, y = Pearson.mean, colour=tag.1),size=1.1) +
      geom_point( aes(x = Year, y = Pearson.mean, colour=tag.1,shape = tag.1),size=3) +
      geom_hline(yintercept=0,linetype=3)+
      labs(x = "Year", y = "Pearson Correlation") +
      scale_y_continuous(limits=c(-0.75,1)) +
      scale_color_manual(values=COLS,name="") +
      scale_fill_manual(values=COLS,name="") +
      scale_shape_manual(values=c(19,17,15),name="")+
      theme_bw() +
      theme(legend.position = c(0.63, 1.02),
            legend.justification = c(0, 1),
            legend.text = element_text(size = 8),
            legend.key.size=unit(0.5,"cm"),
            legend.title = element_text(size=0),
            legend.background=element_rect(colour = "black"),
            panel.border=element_rect(size=1.5,color="black"))
    
    quartz(file=paste("./Scripts and plots for Pubs/Pearson time series, Recruit=",RECRUIT[QQ],", window=",WINDOW,".pdf"),dpi=300,height=4,width=6,type="pdf")
    print(p.line)
    dev.off()
  }
    #########################################################################################
    #########################################################################################
    #########################################################################################
    #########################################################################################
    #########################################################################################
    #### LOREAU style synchrony calculations.
    #########################################################################################
    #########################################################################################
    #########################################################################################
    #########################################################################################
    #########################################################################################
    
########################################################
##############################################################
##############################################################
### Calculate community synchrony.. (across Areas )
##############################################################
##############################################################
##############################################################

    dat   <- aggregate(meanCPUE[meanCPUE$fish.invert=="fish" ,c("Mean.totalDensity","vari")],
                       by=list(Area=meanCPUE$Area[meanCPUE$fish.invert=="fish" ],
                               Year=meanCPUE$year[meanCPUE$fish.invert=="fish" ]),sum)
    
    Synch <- NULL
    Synch.wind <-NULL
    WIND <- 9
    
    for(i in 1:length(unique(dat$Area))){
      temp <- meanCPUE[meanCPUE$area == i ,c("Species","year","Mean.totalDensity")]
      dat.wide <- reshape(temp,timevar=c("Species"),idvar=c("year"),direction="wide")
      A <- community.sync(data = dat.wide[,2:ncol(dat.wide)], nrands = 1000, type=2, quiet= TRUE)
      print(i)
      print(A)
      Synch <- rbind(Synch,c(i,A$obs,A$meancorr,A$pval))
      for(j in 1:length(start.year)){
        dat.wide.trim <- dat.wide[dat.wide$year>=start.year[j] & dat.wide$year <= end.year[j],]
        B <- community.sync(data = dat.wide.trim[,2:ncol(dat.wide.trim)], nrands = 1000, type=1, quiet= TRUE)  
        Synch.wind <- rbind(Synch.wind,c(i,mean.year[j],B$obs,B$meancorr,B$pval))
      }
    }
  colnames(Synch) = c("Area","Synch","MeanCorr","pval")
  colnames(Synch.wind) = c("Area","Mean.Year","Synch","MeanCorr","pval")
  Synch <- data.frame(Synch)
  Synch.wind <- data.frame(Synch.wind)
  
  Synch.wind$EVOS <- "Control"
  Synch.wind$EVOS[Synch.wind$Area >= 2 & Synch.wind$Area <=6] <- "EVOS"
  
  Synch.wind.sum<- summarise(group_by(Synch.wind,EVOS,Mean.Year),M=mean(Synch),SD=sd(Synch),N=length(EVOS))
  Synch.wind.sum$SE <- Synch.wind.sum$SD / sqrt(Synch.wind.sum$N)
  Synch.wind.sum <- data.frame(Synch.wind.sum)
  
  # Define some nice colors
  COLS <- c(  "EVOS"=grey(0.5),
              "Control"=grey(0.5),
              "EVOS-Control"  = "#7570b3")
  
  #scale_fill_manual(values=COLS,name="Origin") 
  ###############################################################################################
  ###############################################################################################
  ###############################################################################################
  #####     ##########################################################################################
  ###############################################################################################
  ############################# HERE IS A PLOT FOR moving window of Synchrony (all Species)
  ###############################################################################################  


  p.line <- ggplot() + 
    geom_line(data=Synch.wind,  aes(x = Mean.Year, y = Synch, colour=EVOS,group=Area,linetype=EVOS),size=0.5,alpha=0.4) +
    geom_line(    data=Synch.wind.sum,  aes(x=Mean.Year, y = M ,fill=EVOS,color=EVOS,linetype=EVOS),size=1,color=1) +
    geom_ribbon(  data=Synch.wind.sum,  aes(x=Mean.Year, ymin = M-SE , ymax = M+SE ,fill=EVOS),alpha=0.3) +
   #     geom_hline(yintercept=0,linetype=3)+
    labs(x = "Year", y = "Synchrony") +
    expand_limits(x = 1988.5, y = 0) +
    scale_y_continuous(limits=c(0.0,0.45),expand = c(0,0)) +
    scale_color_manual(values=COLS,name="") +
    scale_fill_manual(values=COLS,name="") +
    scale_linetype_manual(values=c(1,2),name="") +
    #ggtitle("a)")+
    theme_bw() +
     theme(legend.position = "none",
           panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.border=element_rect(color=1))
  p.line
  
  quartz(file=paste("./Scripts and plots for Pubs/Synchrony,window=",WINDOW,".jpeg"),dpi=600,height=3,width=4.5,type="jpeg")
    # Layout= matrix(c(1,2),nrow=2,ncol=1,byrow=T)
    # QQ <- list(p.line,P.1990.99)
    # multiplot(plotlist=QQ ,layout= Layout)
    print(p.line)
  dev.off()
  
  ###############################################################################################
  ###############################################################################################
  ###############################################################################################
  #####     ##########################################################################################
  ###############################################################################################
  ############################# Repeat for Guilds
  ###############################################################################################  
  

  #Guilds  
  GUILD <- unique(meanCPUE$guild)
  
  Synch <- NULL
  Synch.wind <-NULL
  
  Synch.all  <- NULL
  Synch.wind.all <- NULL
  Synch.wind.sum.all <-NULL
  
  
    WIND <- 9
  
  for(k in 1:length(GUILD)){
  
  for(i in 1:length(unique(dat$Area))){
    temp <- meanCPUE[meanCPUE$area == i & meanCPUE$guild == GUILD[k] ,c("Species","year","Mean.totalDensity")]
    dat.wide <- reshape(temp,timevar=c("Species"),idvar=c("year"),direction="wide")
    A <- community.sync(data = dat.wide[,2:ncol(dat.wide)], nrands = 1000, type=2, quiet= TRUE)
    print(i)
    print(A)
    Synch <- rbind(Synch,c(i,A$obs,A$meancorr,A$pval))
    for(j in 1:length(start.year)){
      dat.wide.trim <- dat.wide[dat.wide$year>=start.year[j] & dat.wide$year <= end.year[j],]
      B <- community.sync(data = dat.wide.trim[,2:ncol(dat.wide.trim)], nrands = 1000, type=2, quiet= TRUE)  
      Synch.wind <- rbind(Synch.wind,c(i,mean.year[j],B$obs,B$meancorr,B$pval))
    }
  }
  colnames(Synch) = c("Area","Synch","MeanCorr","pval")
  colnames(Synch.wind) = c("Area","Mean.Year","Synch","MeanCorr","pval")
  Synch <- data.frame(Synch)
  Synch.wind <- data.frame(Synch.wind)
  
  Synch.wind$EVOS <- "Control"
  Synch.wind$EVOS[Synch.wind$Area >= 2 & Synch.wind$Area <=6] <- "EVOS"
  
  Synch.wind.sum<- summarise(group_by(Synch.wind,EVOS,Mean.Year),M=mean(Synch),SD=sd(Synch),N=length(EVOS))
  Synch.wind.sum$SE <- Synch.wind.sum$SD / sqrt(Synch.wind.sum$N)
  Synch.wind.sum <- data.frame(Synch.wind.sum)
  Synch.wind.sum$Guild <- GUILD[k]

  Synch$Guild        <- GUILD[k]
  Synch.wind$Guild   <- GUILD[k]

  
  Synch.all          <- rbind(Synch.all, Synch)
  Synch.wind.all     <- rbind(Synch.wind.all, Synch.wind)
  Synch.wind.sum.all <- rbind(Synch.wind.sum.all, Synch.wind.sum)
  
  
}  
  
    
  # Define some nice colors
  COLS <- c(  "EVOS"=grey(0.5),
              "Control"=grey(0.5),
              "EVOS-Control"  = "#7570b3")

  
  
    
  # Plots for APEX Pred

  p.line <- ggplot() + 
    geom_line(data=Synch.wind,  aes(x = Mean.Year, y = Synch, colour=EVOS,group=Area,linetype=EVOS),size=0.5,alpha=0.4) +
    geom_line(    data=Synch.wind.sum,  aes(x=Mean.Year, y = M ,fill=EVOS,color=EVOS,linetype=EVOS),size=1,color=1) +
    geom_ribbon(  data=Synch.wind.sum,  aes(x=Mean.Year, ymin = M-SE , ymax = M+SE ,fill=EVOS),alpha=0.3) +
    #     geom_hline(yintercept=0,linetype=3)+
    labs(x = "Year", y = "Synchrony") +
    expand_limits(x = 1988.5, y = 0) +
    scale_y_continuous(limits=c(0.0,1.0),expand = c(0,0)) +
    scale_color_manual(values=COLS,name="") +
    scale_fill_manual(values=COLS,name="") +
    scale_linetype_manual(values=c(1,2),name="") +
    #ggtitle("a)")+
    theme_bw() +
    theme(legend.position = "none",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border=element_rect(color=1))
  p.line
  
  quartz(file=paste("./Scripts and plots for Pubs/Synchrony, APEX, window=",WINDOW,".jpeg"),dpi=600,height=3,width=4.5,type="jpeg")
  # Layout= matrix(c(1,2),nrow=2,ncol=1,byrow=T)
  # QQ <- list(p.line,P.1990.99)
  # multiplot(plotlist=QQ ,layout= Layout)
  print(p.line)
  dev.off()
  
  
  
  
  
  
  
  
  
  
  

