rm(list=ls())
library(INLA)
library(rgdal)
library(ggplot2)
library(splancs)
library(sp)

# Users of this script can pass in other locations to project to. Search for the read.csv() line
# that stores the info in the object 'myAreas' below

runFromDB  = FALSE
keep.MCMC  = FALSE
projection = "goa_mid" # options: "goa_shallow","goa_mid","goa_deep","discrete_areas" "goa_to_500"
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

All.uncond<-NULL
All.pres<-NULL
All.pos <- NULL

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
          setwd("Users/ole.shelton/Dropbox/INLA output/Size Data/Pres")
          DIR <- dir()
          load(DIR[which(substr(DIR,1,this.long)==spp)])
        } 
        if(mod!="binomial"){        
          #setwd("/Users/ole.shelton/Documents/Science/Active projects/Exxon/Groundfish/Positive Output/_RData/_Best Size")
          setwd("Users/ole.shelton/Dropbox/INLA output/Size Data/Pos")
          DIR <- dir()
          load(DIR[which(substr(DIR,1,this.long)==spp)])
        } 
    }      
    }  # end run from DB if statement
      
    if(runFromDB==TRUE) {
      library(rdrop2)
      if(mod=="binomial") {
        presFile = paste(spp,"_","binomial.Rdata",sep="")  
        drop_get(paste("/GoA INLA output/pres/",presFile,sep=""),overwrite=TRUE)     
        load(presFile)
      }
      if(mod!="binomial") {
        posFile = paste(spp,"_","pos.Rdata",sep="")  
        drop_get(paste("/GoA INLA output/pos/",posFile,sep=""),overwrite=TRUE)     
        load(posFile)        
      }
    }
    # Summarize the model
    summary(Output$INLA.mod)
    YEARS	<-	sort(unique(Output$Data$Year))
    nYears <- length(YEARS)
    
    # re-make spde (for projections)
    spde=inla.spde2.matern(Output$Mesh, alpha=3/2)  

    # Things we're interested in: 
    # 1. Diagnostics, maybe plots of observed/ expected values for pos, and AUC for binomial?
    # 2. MCMC draws of biomass / presence-absence at grid cell locations
    # 3. Look at hyperparameter for spatial matern function (range / scale, etc)
    # hyper = inla.spde2.result(Output$INLA.mod, "i2D", spde)

    # Generate predictions on the spatial grid including fixed/random effects
    nMCMC = 1000 # takes about 6 minutes to do the sampling
    inla.mcmc = inla.posterior.sample(nMCMC, Output$INLA.mod)
    #inla.mcmc[[1]]$latent # latent predictions of spatial random field / 

    # get indices of parameters in MCMC output
    re.indx = grep("i2D",rownames(inla.mcmc[[1]]$latent))# Get indx of random effects
    depth.indx1 = grep("cent.log.depth.1",rownames(inla.mcmc[[1]]$latent))
    depth.indx2 = grep("cent.log.depth2.1",rownames(inla.mcmc[[1]]$latent))
    fe.indx = grep("Y1984",rownames(inla.mcmc[[1]]$latent))# Fixed effects stacked after random

    # Read in locations and knots, form projection matrix
    knotLocs = Output$Mesh$loc[,c(1,2)] # locations of knots

	# Read in giant file of all prediction points from the Gulf of Alaska
    setwd(proj.dir)
    dat.project	= read.csv("Output Data/goa_projection_points+temp.csv")
    dat.project$LonUTMAlbers = dat.project$LonUTMAlbers/1000
    dat.project$LatUTMAlbers = dat.project$LatUTMAlbers/1000
    #### Exclude points that end up on land
    dat.project$NGDC24_M =	-dat.project$NGDC24_M	# depth in m
    dat.project$SRTM_M = -dat.project$SRTM_M	# depth in m
    dat.project = dat.project[dat.project$NGDC24 > 0,]

	  #remember to use log.BottomDepth and to center based on the data used in estimation
	  dat.project$cent.log.depth <- dat.project$log.BottomDepth - mean(Output$Covar$log.BottomDepth,na.rm=T)
	  dat.project$cent.log.depth2 <- dat.project$cent.log.depth^2

    # Subset these points to only include the points we're interested 
    # in projecting too, not the entire GoA.  See above for the options.
	  
	  if(projection=="goa_trawl_data") myAreas = read.csv("Output Data/goa_trawl_final_albers+temp.csv")
	  if(projection=="goa_shallow") myAreas = read.csv("Output Data/goa_central_gulf(50_to_150m).csv")
	  if(projection=="goa_mid")  myAreas = read.csv("Output Data/goa_central_gulf(150_to_250m).csv")
	  if(projection=="goa_deep") myAreas = read.csv("Output Data/goa_central_gulf(250_to_deep).csv")
	  if(projection=="goa_to_500")   myAreas = read.csv("Output Data/goa_central_gulf(50_to_500m).csv")
	  if(projection=="discrete_areas") myAreas = read.csv("Output Data/goa_discrete_areas_for_comparison(50_to_150m).csv")
	  
	  if(projection!="goa_trawl_data"){ 
	    dat.project = dat.project[which(is.na(match(dat.project$MASTER_ID, myAreas$MASTER_ID))==F), ]
	  }
	  if(projection=="goa_trawl_data"){
	    dat.project = Output$Data
	    #remember to use log.BottomDepth and to center based on the data used in estimation
	    Output$Data$log.BD	<-  log(Output$Data$BottomDepth)
	    dat.project$cent.log.depth  <- dat.project$log.BottomDepth - mean(Output$Data$log.BD,na.rm=T)
	    dat.project$cent.log.depth2 <- dat.project$cent.log.depth^2
	    dat.project$MASTER_ID <- paste("X",1:nrow(dat.project),sep=".")
	  }
	
    gridLocs = dat.project[,c("LonUTMAlbers","LatUTMAlbers")]
    # Projections will be stored as an array
    projectedLatentGrid = array(0, dim = c(dim(gridLocs)[1], nMCMC, nYears))

    #######Make projections:
    projMatrix <- inla.spde.make.A(Output$Mesh, loc=as.matrix(gridLocs))
    
    if(Output$single.intercept==FALSE){
    # Loop over first year, and do MCMC projections for that year
    for(yr in 1:1) {
      #Grab random effects from this year
      indx = which(Output$iset$i2D.group==yr)

      # Multiply this projection matrix x 
      for(n in 1:nMCMC) {
        # combine random effects, fixed effects
        projectedLatentGrid[,n,yr] = as.numeric(projMatrix%*%inla.mcmc[[n]]$latent[re.indx][indx]) + 
        dat.project[,"cent.log.depth"]*inla.mcmc[[n]]$latent[depth.indx1] + (dat.project[,"cent.log.depth2"])*inla.mcmc[[n]]$latent[depth.indx2] + 
        inla.mcmc[[n]]$latent[fe.indx]
      } # end mcmc loop
    } # end year loop

    # Loop over all subsequent years, and do MCMC projections for that year
    for(yr in 2:nYears) {
      #Grab random effects from this year
      indx = which(Output$iset$i2D.group==yr)
  
      # Multiply this projection matrix x 
      for(n in 1:nMCMC) {
        # combine random effects, fixed effects
        projectedLatentGrid[,n,yr] = as.numeric(projMatrix%*%inla.mcmc[[n]]$latent[re.indx][indx]) + 
        dat.project[,"cent.log.depth"]*inla.mcmc[[n]]$latent[depth.indx1] + (dat.project[,"cent.log.depth2"])*inla.mcmc[[n]]$latent[depth.indx2] + 
        inla.mcmc[[n]]$latent[fe.indx] +   inla.mcmc[[n]]$latent[fe.indx - 1 + yr] 
      } # end mcmc loop
    } # end year loop
    } # end if loop
    
    if(Output$single.intercept==TRUE){
      for(yr in 1:nYears) {
        #Grab random effects from this year
        indx = which(Output$iset$i2D.group==yr)
        
        # Multiply this projection matrix x 
        for(n in 1:nMCMC) {
          # combine random effects, fixed effects
          projectedLatentGrid[,n,yr] = as.numeric(projMatrix%*%inla.mcmc[[n]]$latent[re.indx][indx]) + 
            dat.project[,"cent.log.depth"]*inla.mcmc[[n]]$latent[depth.indx1] + (dat.project[,"cent.log.depth2"])*inla.mcmc[[n]]$latent[depth.indx2] + 
            inla.mcmc[[n]]$latent[fe.indx ]
        } # end mcmc loop
      } # end year loop
    }
     
    # Save files to workspaces for future use -- only if folks need raw MCMC draws by cell
  if(keep.MCMC==FALSE){
    if(mod=="binomial"){logit.projGrid = projectedLatentGrid}
    if(mod=="pos"){log.projGrid = projectedLatentGrid}
  }
    
  if(keep.MCMC==TRUE){
      if(mod == "binomial") {
      logit.projGrid = projectedLatentGrid
      projectionOutput = list("logit.projGrid"=logit.projGrid, "MASTER_ID"=dat.project[,"MASTER_ID"], "gridLocs"=gridLocs, "knotLocs"=knotLocs)
      save(projectionOutput, file=paste(spp,"_","binomial.MCMC",".Rdata",sep=""))
    }
    if(mod != "binomial") {
      log.projGrid = projectedLatentGrid
      projectionOutput = list("log.projGrid"=log.projGrid, "MASTER_ID"=dat.project[,"MASTER_ID"], "gridLocs"=gridLocs, "knotLocs"=knotLocs)
      save(projectionOutput, file=paste(spp,"_","pos.MCMC",".Rdata",sep=""))
    }
  }
  } # end loop over models

  ##############################################################################
  # Now we can do the exciting processing stuff, we we generate indices by area, etc.
  ##############################################################################
  #nYears = dim(logit.projGrid)[3]
  presOutput = plogis(logit.projGrid)# convert to normal space
  posOutput = exp(log.projGrid) # convert to normal space
  # Calculate average density by grid cell and MCMC draw. This is not drawing separate MCMC draws by binomial / pos
  # models, but they're already independent and this was done separately above. 
  totalDensityByCell = presOutput * posOutput

  # Calculate expectation, median, etc. for each cell in each year. [Unconditional expectation]
    DensityByCell = apply(totalDensityByCell, c(1,3), mean)
    DensityMEDIANByCell = apply(totalDensityByCell, c(1,3), median)
    DensitySEByCell = apply(totalDensityByCell, c(1,3), sd)
  
    loc<-dat.project$MASTER_ID
    years<-sort(unique(Output$Data$Year))
    SpeciesPred<-data.frame(Species=spp,MASTER_ID=rep(loc,length(years)),Year=sort(rep(years,length(loc))),Mean=c(DensityByCell),Median = c(DensityMEDIANByCell), SE=c(DensitySEByCell))

  # Calculate just the occurrence component of the model projection
    DensityByCell = apply(posOutput, c(1,3), mean)
    DensityMEDIANByCell = apply(posOutput, c(1,3), median)
    DensitySEByCell = apply(posOutput, c(1,3), sd)
    
    loc<-dat.project$MASTER_ID
    years<-sort(unique(Output$Data$Year))
    SpeciesPred.pos<-data.frame(Species=spp,MASTER_ID=rep(loc,length(years)),Year=sort(rep(years,length(loc))),Mean=c(DensityByCell),Median = c(DensityMEDIANByCell), SE=c(DensitySEByCell))
    
  # Calculate just the occurrence component of the model projection
    DensityByCell = apply(presOutput, c(1,3), mean)
    DensityMEDIANByCell = apply(presOutput, c(1,3), median)
    DensitySEByCell = apply(presOutput, c(1,3), sd)
    
    loc<-dat.project$MASTER_ID
    years<-sort(unique(Output$Data$Year))
    SpeciesPred.pres<-data.frame(Species=spp,MASTER_ID=rep(loc,length(years)),Year=sort(rep(years,length(loc))),Mean=c(DensityByCell),Median = c(DensityMEDIANByCell), SE=c(DensitySEByCell))
    
  ##########################################
  ### Combine all the species into a single file for writing to disk.
  ##########################################  
    
    All.uncond <- rbind(All.uncond,SpeciesPred)
    All.pos    <- rbind(All.pos,SpeciesPred.pos)
    All.pres   <- rbind(All.pres,SpeciesPred.pres)

    write.csv(SpeciesPred,file=paste(spp,"_species_uncond_pred_",projection,".csv",sep=""),row.names=F)
    write.csv(SpeciesPred.pres,file=paste(spp,"_species_pres_pred_",projection,".csv",sep=""),row.names=F)
    write.csv(SpeciesPred.pos,file=paste(spp,"_species_pos_pred_",projection,".csv",sep=""),row.names=F)
    
} #End Species if statement
print(spp)
} # End loop over species 

if(size.data==FALSE){
  write.csv(All.uncond,file=paste("All_species_uncond_pred_",projection,".csv",sep=""),row.names=F)
  write.csv(All.pos,file=paste("All_species_pos_pred_",projection,".csv",sep=""),row.names=F)
  write.csv(All.pres,file=paste("All_species_pres_pred_",projection,".csv",sep=""),row.names=F)
}
if(size.data==TRUE){
  write.csv(All.uncond,file=paste("All_size_species_uncond_pred_",projection,".csv",sep=""),row.names=F)
  write.csv(All.pres,file=paste("All_size_species_pres_pred_",projection,".csv",sep=""),row.names=F)
}