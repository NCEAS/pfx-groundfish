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

proj.dir <-"/Users/ole.shelton/Documents/GitHub/exxonValdez_nceas/goaTrawl/" #getwd()
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
	  Output$Data$log.BD	<-log(Output$Data$BottomDepth)
	  dat.project$cent.log.depth <- dat.project$log.BottomDepth - mean(Output$Data$log.BD,na.rm=T)
	  dat.project$cent.log.depth2 <- dat.project$cent.log.depth^2

    # Subset these points to only include the points we're interested 
    # in projecting too, not the entire GoA. Start with Ole's depth areas
	
		# This one is for a series of areas that are between 50 and 150 m deep
	  if(projection=="goa_shallow") myAreas = read.csv("Output Data/goa_central_gulf(50_to_150m).csv")
	  if(projection=="goa_mid")  myAreas = read.csv("Output Data/goa_central_gulf(150_to_250m).csv")
	  if(projection=="goa_deep") myAreas = read.csv("Output Data/goa_central_gulf(250_to_deep).csv")
	  if(projection=="goa_to_500")   myAreas = read.csv("Output Data/goa_central_gulf(50_to_500m).csv")
	  if(projection=="discrete_areas") myAreas = read.csv("Output Data/goa_discrete_areas_for_comparison(50_to_150m).csv")

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
	if(keep.MCMC==FALSE){
    	logit.projGrid = projectedLatentGrid
    	log.projGrid = projectedLatentGrid
  	}
    if(keep.MCMC ==TRUE){
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
  # Now we can do the exciting processing stuff, we generate indices by area, etc.
  ##############################################################################
  # Calculate the number of spatial locations in each defined area.
 	n.per.area	<-	aggregate(myAreas$Area,by=list(area=myAreas$Area),length)
 	colnames(n.per.area)[2]	<-	"N.obs.Area"
 	n.per.area	<-	rbind(n.per.area,c("Total",sum(n.per.area$N.obs.Area)))
 	
  nYears = dim(logit.projGrid)[3]
  presOutput = plogis(logit.projGrid)# convert to normal space
  posOutput = exp(log.projGrid) # convert to normal space
  # Calculate average density by grid cell and MCMC draw. This is not drawing separate MCMC draws by binomial / pos
  # models, but they're already independent and this was done separately above. 
  totalDensityByCell = presOutput * posOutput

  # Start by calculating total density estimate for these areas (summed spatially)
  totalDensity = apply(totalDensityByCell, c(2,3), mean)
  summaryStats = cbind(apply(totalDensity,2,mean), apply(totalDensity,2,median), apply(totalDensity,2,sd), apply(totalDensity,2,quantile,0.025), apply(totalDensity,2,quantile,0.975))
  colnames(summaryStats) = c("Mean.totalDensity","Median.totalDensity","SD.totalDensity","lower95","upper95")
  # Do the same thing with presence-absence
  totalDensity.pa = apply(presOutput, c(2,3), mean)
  summaryStats.pa = cbind(apply(totalDensity.pa,2,mean), apply(totalDensity.pa,2,median), apply(totalDensity.pa,2,sd), apply(totalDensity.pa,2,quantile,0.025), apply(totalDensity.pa,2,quantile,0.975))
  colnames(summaryStats.pa) = c("Mean.avgPresence","Median.avgPresence","SD.totalPresence","lower95","upper95")
  
  summaryStatsByArea = list()
  summaryStatsByArea.pa = list()  
  for(a in 1:max(myAreas$Area)) {
    # for each area, subset MCMC draws and create time series
    indx = which(myAreas$Area==a)
    # calculate summary stats using just these cells in the focal area
    totalDensity = (apply(totalDensityByCell[indx,,],c(2,3),mean))
    summaryStatsByArea[[a]] = cbind(apply(totalDensity,2,mean), apply(totalDensity,2,median), apply(totalDensity,2,sd), apply(totalDensity,2,quantile,0.025), apply(totalDensity,2,quantile,0.975))
    colnames(summaryStatsByArea[[a]]) = c("Mean.totalDensity","Median.totalDensity","SD.totalDensity","lower95","upper95")    
    # Do same for presence - absence
    totalDensity.pa = (apply(presOutput[indx,,],c(2,3),mean))
    summaryStatsByArea.pa[[a]] = cbind(apply(totalDensity.pa,2,mean), apply(totalDensity.pa,2,median), apply(totalDensity.pa,2,sd), apply(totalDensity.pa,2,quantile,0.025), apply(totalDensity.pa,2,quantile,0.975))
    colnames(summaryStatsByArea.pa[[a]]) = c("Mean.avgPresence","Median.avgPresence","SD.avgPresence","lower95","upper95")
  }
  

  # Create an output data frame of total density (goA) and densities by depth areas
  outputDF = summaryStats
  outputDF.pa = summaryStats.pa
  for(a in 1:max(myAreas$Area)) {
    outputDF = rbind(outputDF, summaryStatsByArea[[a]])
    outputDF.pa = rbind(outputDF.pa, summaryStatsByArea.pa[[a]])    
  }
  
  outputDF = as.data.frame(outputDF)
  outputDF$year.numb = rep(seq(1,nYears), c(max(myAreas$Area)+1))
  outputDF$year = rep(YEARS, c(max(myAreas$Area)+1))
  outputDF$area = c(rep("Total",nYears), sort(rep(seq(1,max(myAreas$Area)),nYears)))

  outputDF.pa = as.data.frame(outputDF.pa)
  outputDF.pa$year.numb = rep(seq(1,nYears), c(max(myAreas$Area)+1))
  outputDF$year = rep(YEARS, c(max(myAreas$Area)+1))
  outputDF.pa$area = c(rep("Total",nYears), sort(rep(seq(1,max(myAreas$Area)),nYears)))

	  outputDF		<- merge(outputDF,n.per.area)
	  outputDF.pa	<- merge(outputDF.pa,n.per.area)
 
 	outputDF$Species <- spp
 	outputDF.pa$Species <- spp 
 	#combine species into on dataframe
 	
 	all.sp.pos	<-	rbind(all.sp.pos,outputDF)
 	all.sp.pres	<-	rbind(all.sp.pres,outputDF.pa)	

 
  # Write outputDF to file
	#   write.table(outputDF, file = paste(spp,"_summedCPUEByArea.csv",sep=""), row.names = F, col.names=T, sep=",")
	#   write.table(outputDF.pa, file = paste(spp,"_occurrenceByArea.csv",sep=""), row.names = F, col.names=T, sep=",")  
  # Make matplot of it
  #matplot(matrix(log(outputDF$Mean.totalDensity[-c(1:12)]), nrow = 12), type = "l", ylab="log(mean total density)")
  
  if(runFromDB==TRUE) {
    # upload results to DB
      upFile = paste(spp,"_","binomial.MCMC.Rdata",sep="")  
      drop_upload(upFile, paste("/GoA INLA output/MCMC/",upFile,sep=""),overwrite=TRUE)     
      # delete / clean up downloaded file
      file.remove(presFile)

      upFile = paste(spp,"_","pos.MCMC.Rdata",sep="")  
      drop_upload(upFile, paste("/GoA INLA output/MCMC/",upFile,sep=""),overwrite=TRUE) 
      # delete / clean up downloaded file
      file.remove(posFile)

  } 
  print(spp) 
} # End loop over species

  # Write outputDF to file
	setwd(proj.dir)
if(size.data==FALSE){
	   write.csv(all.sp.pos, file = paste("All_sp_index_meanCPUEByArea.csv"), row.names = F, col.names=T, sep=",")
	   write.csv(all.sp.pres, file = paste("All_sp_index_occurrenceByArea.csv"), row.names = F, col.names=T, sep=",")  
}
if(size.data==FALSE){
	   write.csv(all.sp.pos, file = paste("All_size_sp_index_meanCPUEByArea.csv"), row.names = F, col.names=T, sep=",")
	   write.csv(all.sp.pres, file = paste("All_size_sp_index_occurrenceByArea.csv"), row.names = F, col.names=T, sep=",")  
}

