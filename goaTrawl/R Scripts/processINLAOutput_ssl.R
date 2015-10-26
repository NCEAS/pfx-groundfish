rm(list=ls())
library(INLA)
library(rgdal)
library(ggplot2)
library(splancs)
library(sp)

# Users of this script can pass in other locations to project to. Search for the read.csv() line
# that stores the info in the object 'myAreas' below

runFromDB = FALSE
keep.MCMC = FALSE
projection = "goa_deep" # options: "goa_shallow","goa_mid","goa_deep","discrete_areas" "goa_to_500"
size.data = FALSE

proj.dir <-getwd()#"/Users/ole.shelton/GitHub/exxonValdez_nceas/goaTrawl/"
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

for(i in 1:length(sppList)) {

spp = sppList[i] # focal species

  for(mod in c("pos","binomial")) {
    # load in respective fitted model workspace
    if(runFromDB == FALSE) {
      if(mod=="binomial") load(paste("/users/eric.ward/downloads/",spp,"_","binomial",".Rdata",sep="")) # load(paste("fittedModels/",spp,"_","binomial",".Rdata",sep=""))
      if(mod!="binomial") load(paste("/users/eric.ward/downloads/",spp,"_","pos",".Rdata",sep=""))
    }
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
    # re-make spde (for projections)
    spde=inla.spde2.matern(Output$Mesh, alpha=3/2)  

    # Things we're interested in: 
    # 1. Diagnostics, maybe plots of observed/ expected values for pos, and AUC for binomial?
    # 2. MCMC draws of biomass / presence-absence at grid cell locations
    # 3. Look at hyperparameter for spatial matern function (range / scale, etc)
    #hyper = inla.spde2.result(Output$INLA.mod, "i2D", spde)

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

    dat.project	= read.csv("Output Data/goa_projection_points+temp.csv")
    dat.project$LonUTMAlbers = dat.project$LonUTMAlbers/1000
    dat.project$LatUTMAlbers = dat.project$LatUTMAlbers/1000
    #### Exclude points that end up on land
    dat.project$NGDC24_M =	-dat.project$NGDC24_M	# depth in m
    dat.project$SRTM_M = -dat.project$SRTM_M	# depth in m
    dat.project = dat.project[dat.project$NGDC24 > 0,]

	  #remember to use log.BottomDepth and to center based on the data used in estimation
	  dat.project$cent.log.depth <- dat.project$log.BottomDepth - mean(Output$Data$log.BottomDepth,na.rm=T)
	  dat.project$cent.log.depth2 <- dat.project$cent.log.depth^2

    # Subset these points to only include the points we're interested 
    # in projecting too, not the entire GoA. Start with Ole's depth areas
	  myAreas = read.table("Output Data/SSLprojectionLocations.csv",header=T,sep=" ")
	  myAreas$MASTER_ID = as.numeric(as.character(myAreas$MASTER_ID))
    dat.project = dat.project[which(is.na(match(dat.project$MASTER_ID, myAreas$MASTER_ID))==F), ]

    gridLocs = dat.project[,c("LonUTMAlbers","LatUTMAlbers")]
    # Projections will be stored as an array
    projectedLatentGrid = array(0, dim = c(dim(gridLocs)[1], nMCMC, 12))

    #######Make projections:
    if(Output$single.intercept==FALSE){
      # Loop over first year, and do MCMC projections for that year
      for(yr in 1:1) {
        #Grab random effects from this year
        indx = which(Output$iset$i2D.group==yr)
        
        projMatrix <- inla.spde.make.A(Output$Mesh, loc=as.matrix(gridLocs))
        
        # Multiply this projection matrix x 
        for(n in 1:nMCMC) {
          # combine random effects, fixed effects
          projectedLatentGrid[,n,yr] = as.numeric(projMatrix%*%inla.mcmc[[n]]$latent[re.indx][indx]) + 
            dat.project[,"cent.log.depth"]*inla.mcmc[[n]]$latent[depth.indx1] + (dat.project[,"cent.log.depth2"])*inla.mcmc[[n]]$latent[depth.indx2] + 
            inla.mcmc[[n]]$latent[fe.indx - 1 + 1]
        } # end mcmc loop
      } # end year loop
      
      # Loop over all subsequent years, and do MCMC projections for that year
      for(yr in 2:12) {
        #Grab random effects from this year
        indx = which(Output$iset$i2D.group==yr)
        
        projMatrix <- inla.spde.make.A(Output$Mesh, loc=as.matrix(gridLocs))
        
        # Multiply this projection matrix x 
        for(n in 1:nMCMC) {
          # combine random effects, fixed effects
          projectedLatentGrid[,n,yr] = as.numeric(projMatrix%*%inla.mcmc[[n]]$latent[re.indx][indx]) + 
            dat.project[,"cent.log.depth"]*inla.mcmc[[n]]$latent[depth.indx1] + (dat.project[,"cent.log.depth2"])*inla.mcmc[[n]]$latent[depth.indx2] + 
            inla.mcmc[[n]]$latent[fe.indx - 1 + 1] +   inla.mcmc[[n]]$latent[fe.indx - 1 + yr] 
        } # end mcmc loop
      } # end year loop
    } # end if loop
    
    if(Output$single.intercept==TRUE){
      for(yr in 1:12) {
        #Grab random effects from this year
        indx = which(Output$iset$i2D.group==yr)
        
        projMatrix <- inla.spde.make.A(Output$Mesh, loc=as.matrix(gridLocs))
        
        # Multiply this projection matrix x 
        for(n in 1:nMCMC) {
          # combine random effects, fixed effects
          projectedLatentGrid[,n,yr] = as.numeric(projMatrix%*%inla.mcmc[[n]]$latent[re.indx][indx]) + 
            dat.project[,"cent.log.depth"]*inla.mcmc[[n]]$latent[depth.indx1] + (dat.project[,"cent.log.depth2"])*inla.mcmc[[n]]$latent[depth.indx2] + 
            inla.mcmc[[n]]$latent[fe.indx - 1 + 1 ]
        } # end mcmc loop
      } # end year loop
    }
    

    # Save files to workspaces for future use -- only if folks need raw MCMC draws by cell
    if(mod == "binomial") {
      logit.projGrid = projectedLatentGrid
      #projectionOutput = list("logit.projGrid"=logit.projGrid, "MASTER_ID"=dat.project[,"MASTER_ID"], "gridLocs"=gridLocs, "knotLocs"=knotLocs)
      #save(projectionOutput, file=paste(spp,"_","binomial.MCMC",".Rdata",sep=""))
    }
    if(mod != "binomial") {
      log.projGrid = projectedLatentGrid
      #projectionOutput = list("log.projGrid"=log.projGrid, "MASTER_ID"=dat.project[,"MASTER_ID"], "gridLocs"=gridLocs, "knotLocs"=knotLocs)
      #save(projectionOutput, file=paste(spp,"_","pos.MCMC",".Rdata",sep=""))
    }
    
  } # end loop over models

  ##############################################################################
  # Now we can do the exciting processing stuff, we we generate indices by area, etc.
  ##############################################################################
  nYears = dim(logit.projGrid)[3]
  presOutput = plogis(logit.projGrid)# convert to normal space
  posOutput = exp(log.projGrid) # convert to normal space
  # Calculate average density by grid cell and MCMC draw. This is not drawing separate MCMC draws by binomial / pos
  # models, but they're already independent and this was done separately above. 
  totalDensityByCell = presOutput * posOutput

  # apply each MCMC draw to SSL areas, and sum across all cells in the area
  sslAreas = as.matrix(myAreas[,-c(1:21)])
  nAreas = dim(sslAreas)[2]

  output = t(totalDensityByCell[,,1]) %*% sslAreas
  for(y in 2:12) {
    output = rbind(output, t(totalDensityByCell[,,y]) %*% sslAreas)
  }
  colnames(output) = names(myAreas[,-c(1:21)])
  output = as.data.frame(output)
  output$Year = sort(rep(seq(1,12), dim(totalDensityByCell)[2]))
  
  write.table(output, paste(spp,"_SSL_MCMC.csv",sep=""),sep=",",row.names=F,col.names=T)
  
  if(runFromDB==TRUE) {
    # upload results to DB
      upFile = paste(spp,"_","binomial_SSL.MCMC.Rdata",sep="")  
      #drop_upload(upFile, paste("/GoA INLA output/MCMC/",upFile,sep=""),overwrite=TRUE)     
      # delete / clean up downloaded file
      #file.remove(presFile)

      upFile = paste(spp,"_","pos_SSL.MCMC.Rdata",sep="")  
      #drop_upload(upFile, paste("/GoA INLA output/MCMC/",upFile,sep=""),overwrite=TRUE) 
      # delete / clean up downloaded file
      #file.remove(posFile)

  }  
}