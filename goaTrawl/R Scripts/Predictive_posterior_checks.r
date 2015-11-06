## This is a posterior predictive check on the projection of the fitted model objects

rm(list=ls())
library(INLA)
library(rgdal)
library(ggplot2)
#library(splancs)
library(sp)

# Users of this script can pass in other locations to project to. Search for the read.csv() line
# that stores the info in the object 'myAreas' below

runFromDB  = FALSE
keep.MCMC  = FALSE
projection = "discrete_areas" # options: "goa_shallow","goa_mid","goa_deep","discrete_areas" "goa_to_500"
size.data  = FALSE

proj.dir <- getwd() #"/Users/ole.shelton/Documents/GitHub/exxonValdez_nceas/goaTrawl/" 
setwd(proj.dir)
sppList = as.character(read.csv("speciesNames.csv",header=F)[,1])

data.dir  <- "/Users/ole.shelton/Dropbox/INLA output"
setwd(data.dir)
dat.pres    <- read.csv("All_species_pres_pred_goa_trawl_data.csv")
dat.uncond  <- read.csv("All_species_uncond_pred_goa_trawl_data.csv")

dat.trawl = read.csv("Output Data/trawl_data_with_MASTER_ID.csv")
dat.trawl$LonUTMAlbers <- dat.trawl$LonUTMAlbers / 1000
dat.trawl$LatUTMAlbers <- dat.trawl$LatUTMAlbers / 1000

setwd(proj.dir)
# Read in projections from MCMC draws to generate predictive comparisons for each observed locations.

if(size.data == TRUE){
  maxCol = 31
  df = read.csv("Output Data/goa_trawl_final_size_albers+temp.csv")
  species = names(df)[14:maxCol] 
  temp <- NULL
  for(ZZ in 1:length(species)){
    temp[ZZ] <- substr(species[ZZ],nchar(species[ZZ])-3,nchar(species[ZZ]))
  }
  sppList <- species[temp =="cpue"]
}

all.sp.pos	<-	NULL
all.sp.pres	<-	NULL

for(i in 1:length(sppList)) {
  
  spp = sppList[i] # focal species
  this.long = nchar(spp)
  for(mod in c("pos","binomial")) {
    # load in respective fitted model workspace
    if(runFromDB == FALSE) {
      #if(mod=="binomial") load(paste("/users/eric.ward/downloads/",spp,"_","binomial",".Rdata",sep="")) # load(paste("fittedModels/",spp,"_","binomial",".Rdata",sep=""))
      #if(mod!="binomial") load(paste("/users/eric.ward/downloads/",spp,"_","pos",".Rdata",sep=""))
      #if(mod=="binomial") load(paste("/Users/ole.shelton/Documents/Science/Active projects/Exxon/Groundfish/Binomial Output/",spp,"_","binomial",".Rdata",sep="")) 
      #if(mod!="binomial") load(paste("/Users/ole.shelton/Documents/Science/Active projects/Exxon/Groundfish/Positive Output/",spp,"_","pos",".Rdata",sep=""))
      if(size.data ==FALSE){
        if(mod=="binomial"){
          #setwd("/Users/ole.shelton/Documents/Science/Active projects/Exxon/Groundfish/Binomial Output/_RData/_Best")
          setwd("/Users/ole.shelton/Dropbox/INLA output/pres")
          DIR <- dir()
          load(DIR[which(substr(DIR,1,this.long)==spp)])
        } 
        if(mod!="binomial"){        
          #setwd("/Users/ole.shelton/Documents/Science/Active projects/Exxon/Groundfish/Positive Output/_RData/_Best")
          setwd("/Users/ole.shelton/Dropbox/INLA output/pos")
          DIR <- dir()
          load(DIR[which(substr(DIR,1,this.long)==spp)])
        } 
      }
      if(size.data ==TRUE){
        if(mod=="binomial"){
          setwd("/Users/ole.shelton/Documents/Science/Active projects/Exxon/Groundfish/Binomial Output/_RData/_Best Size")
          DIR <- dir()
          load(DIR[which(substr(DIR,1,this.long)==spp)])
        } 
        if(mod!="binomial"){        
          setwd("/Users/ole.shelton/Documents/Science/Active projects/Exxon/Groundfish/Positive Output/_RData/_Best Size")
          DIR <- dir()
          load(DIR[which(substr(DIR,1,this.long)==spp)])
        } 
      }      
    }  # end run from DB if statement    # Summarize the model
    summary(Output$INLA.mod)
    YEARS	<-	sort(unique(Output$Data$Year))
    # re-make spde (for projections)
    spde=inla.spde2.matern(Output$Mesh, alpha=3/2)  
    
    # Things we're interested in: 
    # 1. Diagnostics, maybe plots of observed/ expected values for pos, and AUC for binomial?
    # 2. MCMC draws of biomass / presence-absence at grid cell locations
    # 3. Look at hyperparameter for spatial matern function (range / scale, etc)
    #hyper = inla.spde2.result(Output$INLA.mod, "i2D", spde)
    
    # Merge the predictions from within the model, the observed values, and the posterior predictions from the MCMC sampels
        # first merge trawl id with the raw trawl data
    
    DAT.out	<-	data.frame(cbind(Output$Data,
                                Mean =  Output$INLA.mod$summary.fitted.values$mean[1:nrow(Output$Data)],
                                Median =  Output$INLA.mod$summary.fitted.values$"0.5quant"[1:nrow(Output$Data)],
                                #SD =  Output$INLA.mod$summary.fitted.values$sd[1:nrow(Output$Data)],
                                q.025 = Output$INLA.mod$summary.fitted.values$"0.025quant"[1:nrow(Output$Data)],
                                q.975 = Output$INLA.mod$summary.fitted.values$"0.975quant"[1:nrow(Output$Data)]))
    #DAT.out$y.jit <- DAT.out[,SPECIES] +runif(nrow(DAT.out),-JIT,JIT)
    DAT.out$Mean[is.na(DAT.out[,spp])==T]  <- exp(DAT.out$Mean[is.na(DAT.out[,spp])==T])
    DAT.out$Median[is.na(DAT.out[,spp])==T]  <- exp(DAT.out$Median[is.na(DAT.out[,spp])==T])
    DAT.out$q.025[is.na(DAT.out[,spp])==T]  <- exp(DAT.out$q.025[is.na(DAT.out[,spp])==T])
    DAT.out$q.975[is.na(DAT.out[,spp])==T]  <- exp(DAT.out$q.975[is.na(DAT.out[,spp])==T])
    
    dat.trawl.all <- merge(dat.trawl,DAT.out)
    
    temp.dat  <- dat.uncond[dat.uncond$Species == spp,]
    colnames(temp.dat)[4:6] <- c("Mean.mcmc","Median.mcmc","SE.mcmc")
    
    combine <- merge(dat.trawl.all,temp.dat)
                     
                     
                     