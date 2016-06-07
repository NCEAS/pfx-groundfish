##Boostrapping CPUE and Occurrence data to to identify variance around diversity metrics
#Created by M. Hunsicker, O. Shelton & J. Watson

rm(list=ls())

# Load all necessary packages
library(vegan); library(plyr) ; library(dplyr); library(tidyr) ; library(stringr) ;library(rgdal) ; 
library(sp) ; library(maptools) ; library(raster) ; library(rgeos) ; library(ggplot2) ; 
library(maps) ; library(mapdata) ; library(mapproj) ; library(httr) ; library(extrafont)

### Beta Method of Moments - for occurrence data
beta.mom <- function(MEAN,SD){
  VAR   <- SD^2
  alpha <- MEAN * ((MEAN *(1 - MEAN)) / VAR - 1)
  beta  <- (1-MEAN) * ((MEAN *(1 - MEAN)) / VAR - 1)
  OUT   <- list(alpha=alpha,beta=beta) 
  #return(OUT)
  return(list(alpha=alpha,beta=beta))
}

#MEAN <- 0.7
#SD   <- 0.05
#dat <- beta.mom(MEAN,SD)

#X <- rbeta(1e5,dat$alpha,dat$beta)

#MEAN; mean(X)
#SD; sd(X)

#hist(X,breaks=seq(0,1,by=0.01))

### GAMMA method of Moment -  for unconditional data
#  (alpha,beta) parameterization

gamma.mom <- function(MEAN,SD){
  VAR   <- SD^2
  alpha <- MEAN^2 / VAR
  beta  <- MEAN / VAR
  #OUT   <- list(alpha=alpha,beta=beta) 
  #return(OUT)
  return(list(alpha=alpha,beta=beta)) #################### Can save a few pieces of processing/memory. 
  ####### And save directly to df so you don't have to unlist later.
}

#MEAN <- 1.87e-07
#SD   <- 2.47e-06
#dat <- gamma.mom(MEAN,SD)

#X <- rgamma(1,dat$alpha,dat$beta)

#MEAN; mean(X)
#SD; sd(X)

#hist(X)

setwd("~/Documents/git_repos/pfx-groundfish")
setwd("//nmfs.local/AKC-ABL/Users2/jordan.watson/Desktop/AFSC/GOA/pfx_Groundfish")
#uncondDat<-read.csv("All_sp_index_meanCPUEByArea.Shallow.combined.MH.csv")
#uncondDat<-read.csv("All_sp_index_meanCPUEByArea.Deep.MH.csv")
#occDat<-read.csv("All_sp_index_occurrenceByArea.Shallow.combined.MH.csv")
occDat<-read.csv("All_sp_index_occurrenceByArea.Deep.MH.csv")


#-----------------------------------------------
dat<-uncondDat
dat<-dat[dat$area!="Total",]

Div.Metrics<-data.frame()  
BootstrapLength<-1000

years <- unique(dat$year)

system.time({
  for (i in 1:length(years)){
    
    dati <- dat[dat$year == unique(dat$year)[i],]
    
    areas <- unique(dati$area) #  Create a vector of areas per year
    
    for (j in 1:length(areas)){  
      
      datj <- dati[dati$area == areas[j],]
      
      ab <- gamma.mom(datj$Mean.totalDensity,datj$SD.totalDensity)
      
      # ----- boot loop
      
      for(b in 1:BootstrapLength){
        
        #  generate a vector of abundance / species in ab
        obs.sp <- rgamma(1:length(ab$alpha),ab$alpha,ab$beta) ##################
        
        #calculate diversity metrcs
        year=years[i] 
        area=areas[j] 
        SW_Div <- diversity(obs.sp, index="shannon")
        Simp_Div <- diversity(obs.sp, index="simpson")
        Invsimp <-diversity(obs.sp,index="invsimpson")
        
        Div.Metrics <- rbind(Div.Metrics,data.frame(cbind(year,area,SW_Div,Simp_Div,Invsimp)))  
        
      }# END OF BOOT
    } # end of j
    print(i)
  } # end of i
})

colnames(Div.Metrics)<-c("YEAR","AREA" ,"SW", "Simp","Invsimp")

Div.Metrics$Eff_Num_Sp <- exp(Div.Metrics$SW)

#write.csv(Div.Metrics,"DiversityMetrics.Deep.Bootstrapped.csv")
#write.csv(Div.Metrics,"DiversityMetrics.Shallow.Bootstrapped.csv")


##### To calculate Species richness
rm("dat")

dat<-occDat
dat<-dat[dat$area!="Total",]

Sp.Richness<-data.frame()

years <- unique(dat$year)

system.time({
  for (i in 1:length(years)){
    
    dati <- dat[dat$year == unique(dat$year)[i],]
    
    areas <- unique(dati$area) #  Create a vector of areas per year
    
    for (j in 1:length(areas)){  
      
      datj <- dati[dati$area == areas[j],]
      
      ab <- beta.mom(datj$Mean.avgPresence,datj$SD.totalPresence)
      
      # ----- boot loop
      
      for(b in 1:BootstrapLength){
        
        #  generate a vector of abundance / species in ab
        obs.sp <- rbeta(1:length(ab$alpha),ab$alpha,ab$beta) ##################
        
        #calculate diversity metrcs
        year=years[i] 
        area=areas[j] 
        Sp_rich <- sum(obs.sp)
        
        Sp.Richness <- rbind(Sp.Richness,data.frame(cbind(year,area,Sp_rich)))  
        
      }# END OF BOOT
    } # end of j
    print(i)
  } # end of i
})

colnames(Sp.Richness)<-c("YEAR","AREA" ,"Sp_rich")

#Div.Metrics<-cbind(Div.Metrics, Sp.Richness)

#write.csv(Sp.Richness,"Sp.Richness.Shallow.Bootstrapped.csv")
#write.csv(Sp.Richness,"Sp.Richness.Deep.Bootstrapped.csv")




###Plotting  -- This code is also in the Portfolio Diversity Plots R markdown folder

# Make custom theme for plotting
theme_boxplot <- function(base_size = 12){
  theme_bw(base_size)%+replace%
    theme(legend.key.size=unit(15,"points"),
          legend.text=element_text(size=14),
          legend.key=element_blank(),
          legend.title=element_blank(),
          legend.background=element_rect(colour="white", fill="transparent"),
          plot.margin=unit(c(0.5,1,0.5,1), "lines"),  
          panel.border=element_blank(),
          panel.margin=unit(0,"lines"),
          panel.background=element_rect(fill=NA, colour=NA),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.ticks.length=unit(1,"mm"),
          axis.text.x = element_text(margin=margin(5,0,0,0), size=15), 
          axis.text.y = element_text(margin=margin(0,5,0,0), size=15),
          axis.title.x=element_text(size=15, margin=margin(15,0,0,0)),
          axis.title.y=element_text(size=15, angle=90, margin=margin(0,15,0,0)),
          strip.text.x=element_text(size=14),
          strip.background=element_rect(colour="black", fill='white'))
}

## Data for 9 deep areas
Div.Metrics<-read.csv("DiversityMetrics.Shallow.Bootstrapped.csv")
Div.Metrics<-Div.Metrics[Div.Metrics$AREA!=10,]
Sp.Richness<-read.csv("Sp.Richness.Shallow.Bootstrapped.csv")
Sp.Richness<-Sp.Richness[Sp.Richness$AREA!=10,]
shallowDat<-cbind(Div.Metrics,Sp.Richness)
shallowDat<-shallowDat[,c(1:6,10)]

## Data for 5 deep areas
Div.Metrics<-read.csv("DiversityMetrics.Deep.Bootstrapped.csv")
Div.Metrics<-Div.Metrics[Div.Metrics$AREA!=6,]
Sp.Richness<-read.csv("Sp.Richness.Deep.Bootstrapped.csv")
Sp.Richness<-Sp.Richness[Sp.Richness$AREA!=6,]
deepDat<-cbind(Div.Metrics,Sp.Richness)
deepwDat<-deepDat[,c(1:6,10)]

shallowDat$AREA<-as.factor(shallowDat$AREA)
deepDat$AREA<-as.factor(deepDat$AREA)

##Boxplots Alpha Diversity Shallow Areas (exp H'):
sp_box1 <- ggplot(data=shallowDat, aes(AREA, y=Eff_Num_Sp)) + 
  geom_boxplot() + theme_boxplot() + xlab("Area (West <-> East)") +
  ylab("Alpha Diversity (exp H') - Shallow areas") +
  xlim("9","8","7","6","5","4","3","2","1") + #shallow areas
  #  scale_fill_manual(values=barcolor) + ylim(0,0.5) +
  theme(legend.position="none", plot.background=element_blank(),
        axis.text.x=element_text(size=15),
        plot.title=element_text(colour="black", size=15,
                                hjust=0.04, vjust=0.5, face="bold"))
sp_box1

##Boxplots Alpha Diversity Deep Areas (exp H'):
sp_box2 <- ggplot(data=deepDat, aes(AREA, y=Eff_Num_Sp)) + 
  geom_boxplot() + theme_boxplot() + xlab("Area (West <-> East)") +
  ylab("Alpha Diversity (exp H') - Deep Areas") +
  xlim("5","4","3","2","1") 
  #  scale_fill_manual(values=barcolor) + ylim(0,0.5) +
  theme(legend.position="none", plot.background=element_blank(),
        axis.text.x=element_text(size=15),
        plot.title=element_text(colour="black", size=15,
                                hjust=0.04, vjust=0.5, face="bold"))
sp_box2

##Boxplots Species Richness shallow areas:
sp_box3 <- ggplot(data=Sp.Richness, aes(AREA, y=Sp_rich)) + 
  geom_boxplot() + theme_boxplot() + xlab("Area (West <-> East)") +
  ylab("Species Richness - Shallow Areas") +
  xlim("9","8","7","6","5","4","3","2","1") + 
  #  scale_fill_manual(values=barcolor) + ylim(0,0.5) +
  theme(legend.position="none", plot.background=element_blank(),
        axis.text.x=element_text(size=15),
        plot.title=element_text(colour="black", size=15,
                                hjust=0.04, vjust=0.5, face="bold"))
sp_box3


##Boxplots Species Richness Deep areas:
sp_box4 <- ggplot(data=deepDat, aes(AREA, y=Sp_rich)) + 
  geom_boxplot() + theme_boxplot() + xlab("Area (West <-> East)") +
  ylab("Species Richness - Deep Areas") +
  xlim("5","4","3","2","1") + 
  #  scale_fill_manual(values=barcolor) + ylim(0,0.5) +
  theme(legend.position="none", plot.background=element_blank(),
        axis.text.x=element_text(size=15),
        plot.title=element_text(colour="black", size=15,
                                hjust=0.04, vjust=0.5, face="bold"))
sp_box4


##Boxplots Inverse Simpson shallow areas:
sp_box3 <- ggplot(data=shallowDat, aes(AREA, y=Invsimp)) + 
  geom_boxplot() + theme_boxplot() + xlab("Area (West <-> East)") +
  ylab("Inverse Simpson - Shallow Areas") +
  xlim("9","8","7","6","5","4","3","2","1") + 
  #  scale_fill_manual(values=barcolor) + ylim(0,0.5) +
  theme(legend.position="none", plot.background=element_blank(),
        axis.text.x=element_text(size=15),
        plot.title=element_text(colour="black", size=15,
                                hjust=0.04, vjust=0.5, face="bold"))
sp_box3


##Boxplots Inverse Simpson Deep areas:
sp_box4 <- ggplot(data=deepDat, aes(AREA, y=Invsimp)) + 
  geom_boxplot() + theme_boxplot() + xlab("Area (West <-> East)") +
  ylab("Inverse Simpson - Deep Areas") +
  xlim("5","4","3","2","1") + 
  #  scale_fill_manual(values=barcolor) + ylim(0,0.5) +
  theme(legend.position="none", plot.background=element_blank(),
        axis.text.x=element_text(size=15),
        plot.title=element_text(colour="black", size=15,
                                hjust=0.04, vjust=0.5, face="bold"))
sp_box4




######### To create plots with mean and CIs or SEs
sumStat <- Div.Metrics %>%
  filter(AREA!=10) %>%
  group_by(AREA) %>% 
  summarise(Mean_Eff_Num_Sp=mean(Eff_Num_Sp), 
            SD = sd(Eff_Num_Sp),
            N = length(Eff_Num_Sp)) %>%
  ungroup() %>%
  print(width = Inf)

sumStat$se <- sumStat$SD / sqrt(sumStat$N)  # Calculate standard error of the mean

# Confidence interval multiplier for standard error
# e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
conf.interval<-0.95
ciMult <- qt(conf.interval/2 + .5, sumStat$N-1)
sumStat$ci <- sumStat$se * ciMult

## Add CIs to boxplots
sp_box2<-ggplot(sumStat, aes(x=AREA, y=Mean_Eff_Num_Sp, colour="black"))+
  geom_errorbar(aes(ymin=Mean_Eff_Num_Sp-ci, ymax=Mean_Eff_Num_Sp+ci), width=0.5)+
      xlab("Area (West <-> East)") +
      ylab("Mean Alpha Diversity (exp H')") +
      xlim("9","8","7","6","5","4","3","2","1") +
      ylim(0,15) +
      #  scale_fill_manual(values=barcolor) + ylim(0,0.5) +
      theme(legend.position="none", plot.background=element_blank(),
      axis.text.x=element_text(size=15),
      plot.title=element_text(colour="black", size=15,
      hjust=0.04, vjust=0.5, face="bold"))

sp_box2

