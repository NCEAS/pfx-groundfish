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
projection = "discrete_areas" # options: 
                              #"goa_shallow","goa_mid","goa_deep","discrete_areas" "goa_to_500"
                              #goa_trawl_data 
size.data  = FALSE

proj.dir <-"/Users/ole.shelton/Documents/GitHub/pfx-groundfish/goaTrawl/" #getwd()
setwd(proj.dir)
sppList = as.character(read.csv("speciesNames.csv",header=F)[,1])

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

 all.depth.param <- NULL
 
 
for(i in 1:length(sppList)) {
spp = sppList[i] # focal species
this.long = nchar(spp)
if(spp != "Merluccius.productus" & spp!="Hydrolagus.colliei"){# Exclude a couple of species that have missing model components.
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
          setwd("/Users/ole.shelton/Dropbox/INLA output/pres/Best")
          DIR <- dir()
          load(DIR[which(substr(DIR,1,this.long)==spp)])
        } 
        if(mod!="binomial"){        
          #setwd("/Users/ole.shelton/Documents/Science/Active projects/Exxon/Groundfish/Positive Output/_RData/_Best")
          setwd("/Users/ole.shelton/Dropbox/INLA output/pos/Best")
          DIR <- dir()
          load(DIR[which(substr(DIR,1,this.long)==spp)])
        } 
    }
    if(size.data ==TRUE){
        if(mod=="binomial"){
          #setwd("/Users/ole.shelton/Documents/Science/Active projects/Exxon/Groundfish/Binomial Output/_RData/_Best Size")
          setwd("/Users/ole.shelton/Dropbox/INLA output/Size Data/Pres")
          DIR <- dir()
          load(DIR[which(substr(DIR,1,this.long)==spp)])
        } 
        if(mod!="binomial"){        
          #setwd("/Users/ole.shelton/Documents/Science/Active projects/Exxon/Groundfish/Positive Output/_RData/_Best Size")
          setwd("/Users/ole.shelton/Dropbox/INLA output/Size Data/Pos")
          DIR <- dir()
          load(DIR[which(substr(DIR,1,this.long)==spp)])
        } 
    }      
    }  # end run from DB if statement    # Summarize the model
    summary(Output$INLA.mod)
    dat   <- Output$INLA.mod
    log.center <- mean(Output$Covar$log.depth)
    
    temp <- data.frame(Species=spp,Model=mod,Coef.name=rownames(dat$summary.fixed),dat$summary.fixed,log.center=log.center)
    all.depth.param <- rbind(all.depth.param,temp)     
    } # end mod loop
  }  # end if catch for merluccius.productus, etc.
}    # end sppList
    
  	setwd(proj.dir)
if(size.data==FALSE){
  	 write.csv(all.depth.param, file = paste("All_sp_depth_coef.csv"), row.names = F)
}
if(size.data==TRUE){
    write.csv(all.depth.param, file = paste("All_sp_depth_coef_size.csv"), row.names = F)
}

