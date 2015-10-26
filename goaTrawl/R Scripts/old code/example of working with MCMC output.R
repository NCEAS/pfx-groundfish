

spec = as.character(read.csv("speciesNames.csv",header=F)[,1])


for(i in 3:3) {
  # read in positive model
  load(paste(spec[i],"_pos.MCMC.Rdata",sep=""))
  posOutput = projectionOutput$log.projGrid
  # convert positive model from log -> normal space
  posOutput = lapply(posOutput,exp)
  
  nCells = dim(projectionOutput$log.projGrid[[1]])[1]
  nYears = 12

  # read in binomial model
  load(paste(spec[i],"_pres.MCMC.Rdata",sep=""))
  presOutput = projectionOutput$logit.projGrid
  # convert binomial model from logit -> normal space
  presOutput = lapply(presOutput,plogis)
    
  # Calculate average density by cell 
  avgDensity = array(NA, dim = c(dim(presOutput[[1]])[1],dim(presOutput[[1]])[2],  12))
  for(y in 1:nYears) {
    avgDensity[,,y] = presOutput[[y]] * posOutput[[y]]
  }
  
  # Calculate total density in GoA by year and MCMC draw combination - this gives us variance
  totalDensity = (apply(avgDensity,c(2,3),sum))
  
  summaryStats = cbind(apply(totalDensity,2,mean), apply(totalDensity,2,median), apply(totalDensity,2,quantile,0.025), apply(totalDensity,2,quantile,0.975))
  colnames(summaryStats) = c("Mean.totalDensity","Median.totalDensity","lower95","upper95")
  
  # coordinates of 
  coords = projectionOutput$gridLocs
  masterID = projectionOutput$MASTER_ID
  knotLocs = projectionOutput$knotLocs
  
  # Create time series of total density by Ole's area
  depthAreas = read.csv("Output Data/goa_discrete_areas_for_comparison(50_to_150m).csv")
  
  summaryStatsByArea = list()
  for(a in 1:max(depthAreas$Area)) {
    # for each area, subset MCMC draws and create time series
    ids = depthAreas$MASTER_ID[which(depthAreas$Area==a)]
    indx = which( masterID %in%ids)
    # calculate summary stats using just these values
    totalDensity = (apply(avgDensity[indx,,],c(2,3),sum))
    
    summaryStatsByArea[[a]] = cbind(apply(totalDensity,2,mean), apply(totalDensity,2,median), apply(totalDensity,2,quantile,0.025), apply(totalDensity,2,quantile,0.975))
    colnames(summaryStatsByArea[[a]]) = c("Mean.totalDensity","Median.totalDensity","lower95","upper95")    
    
  }
  
  # Create a data frame of total density (goA) and densities by depth areas
  outputDF = summaryStats
  for(a in 1:max(depthAreas$Area)) {
    outputDF = rbind(outputDF, summaryStatsByArea[[a]])
  }
  outputDF = as.data.frame(outputDF)
  outputDF$year = rep(seq(1,nYears), c(max(depthAreas$Area)+1))
  outputDF$area = c(rep("Total",12), sort(rep(seq(1,max(depthAreas$Area)),nYears)))
  
  # Write outputDF to file
  
  # Make matplot of it
  matplot(matrix(log(outputDF$Mean.totalDensity[-c(1:12)]), nrow = 12), type = "l", ylab="log(mean total density)")
  
  
}