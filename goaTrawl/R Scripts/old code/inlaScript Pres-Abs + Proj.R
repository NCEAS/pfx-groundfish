rm(list=ls())
library(INLA)
library(rgdal)
library(ggplot2)
library(splancs)
library(sp)

# Use only species with length-based data, or total biomass
totalBiomass = TRUE
#### CHOOSE A MODEL - "binomial" or "positive"
model = "positive"

#Define directories for the data and for the plots
proj.dir = getwd()
plot.dir = paste(proj.dir,"/_Output plots Pres Abs",sep="")
if(model != "binomial") plot.dir = paste(proj.dir,"/_Output plots Pos",sep="")

#### GO GET THE OBSERVED TRAWL DATA
setwd(proj.dir)
df = read.csv("Output Data/goa_trawl_final_albers+temp.csv")
# include switch for whether the length data will be used instead
if(totalBiomass == FALSE) {
  # then use the length data
  df = read.csv("Output Data/goa_trawl_final_size_albers+temp.csv")
}

df = df[order(df$Year,df$Lat),]

### Remove NA entries in BottomDepth and Bottom Temp
df = df[df$BottomDepth != -9999,]
df = df[df$BottomTemp != -9999,]

df = df[df$BottomDepth <= 500,] # Check with Ole. 500 cutoff ok?

df$Station = as.character(df$Station)
df$Year = as.numeric(df$Year)
df$LonUTMAlbers = df$LonUTMAlbers/1000
df$LatUTMAlbers = df$LatUTMAlbers/1000

#### GO GET THE PROJECTION POINTS
dat.project	= read.csv("Output Data/goa_projection_points+temp.csv")

dat.project$LonUTMAlbers = dat.project$LonUTMAlbers/1000
dat.project$LatUTMAlbers = dat.project$LatUTMAlbers/1000

#### Exclude points that end up on land. Check this with Ole, not working -- units may be off
dat.project$NGDC24_M =	-dat.project$NGDC24_M	# depth in m
dat.project$SRTM_M = -dat.project$SRTM_M	# depth in m
dat.project = dat.project[dat.project$NGDC24 > 0,]

THESE	= seq(1,nrow(dat.project),by=15)
dat.trim = dat.project[THESE,]

dat.trim = dat.project[THESE,]
dat.new.trim = dat.trim
x.lim = c(-200,200)
x.lim = c(min(dat.project$LonUTMAlbers),max(dat.project$LonUTMAlbers))
x.lim = c(-800,550)
y.lim = c(5700,6200)
y.lim = c(min(dat.project$LatUTMAlbers),max(dat.project$LatUTMAlbers))
	#plot(LatUTMAlbers~LonUTMAlbers,data=dat.project[ #&
	#		 dat.project$LonUTMAlbers < 200 & dat.project$LonUTMAlbers > -200
	#		 ,],pch=".",xlim=x.lim,ylim=y.lim)
	#par(new=T)
plot(LatUTMAlbers~LonUTMAlbers,data=dat.trim[dat.trim$NGDC24_M<=500 &
	 dat.trim$LonUTMAlbers < 550 & dat.trim$LonUTMAlbers > -800 &
	 dat.trim$NGDC24_M > 25
	 ,],pch=".",xlim=x.lim,ylim=y.lim,col=2)

par(new=T)
plot(LatUTMAlbers~LonUTMAlbers,data=dat.new.trim#[
#		 dat.new$LonUTMAlbers < 200 & dat.new$LonUTMAlbers > -200
		 #,]
 ,pch=".",xlim=x.lim,ylim=y.lim,col=4)
dim(dat.new.trim)
dat.proj.trim	= dat.new.trim

# This is a switch for whether we're using full data, or length-stratified data
maxCol = dim(df)[2]
species = names(df)[17:maxCol]
nCovCol = 16
if(totalBiomass==FALSE) {
  maxCol = 31
  species = names(df)[14:maxCol] 
  nCovCol = 13
}

#################################### START INLA LOOP
for(i in 1:length(species)) {
  # Fit the model for species XX
   subdat = df[,c(1:nCovCol,which(names(df)==species[i]))]
   rows2drop = which(is.na(subdat$LonUTMAlbers+subdat$LatUTMAlbers))
   if(length(rows2drop) > 0)subdat = subdat[-rows2drop,]

 	#get rid of NAs in subdat response
 	subdat[is.na(subdat[,species[i]])==T,species[i]] =	0
 	#if binomial replace positive values with 1
 	if(model =="binomial"){
 	  # any positive value gets transformed -> 1
 		subdat[subdat[,species[i]]>0,species[i]] = 1
 	}
  if(model != "binomial") {
    # drop out zeros, and only model positive response. note: this code chunk got moved below
    #subdat = subdat[subdat[,species[i]]>0,]
  }
 	
  #Center the covariates 
  Covar	<- subdat[,1:nCovCol]
  Covar$log.depth = log(Covar$BottomDepth)
  Covar$cent.log.depth = Covar$log.depth - mean(Covar$log.depth,na.rm=T)
  Covar$cent.log.depth.2 = Covar$cent.log.depth^2 
  Covar$cent.temp = Covar$BottomTemp - mean(Covar$BottomTemp)
  Covar$cent.temp.2 = Covar$cent.temp^2 
  
	#call basic plotting routine for raw data
	  if(model == "binomial"){
	  	source("R Scripts/trawl_plot_binom.r") 
	  }
	#call basic plotting routine for raw data
	  if(model == "positive"){
	  	source("R Scripts/trawl_plot_positive.r")
	  }

  # Grab X-Y coords in UTM space
  subcoords = cbind(subdat$LonUTMAlbers[match(unique(subdat$Station),subdat$Station)],subdat$LatUTMAlbers[match(unique(subdat$Station),subdat$Station)])

  bnd = inla.nonconvex.hull(subcoords, convex=80)
  # increase cutoff to ~ 150 to create much coarser mesh
  mesh1 = inla.mesh.2d(boundary=bnd,max.edge=c(60,1500),cutoff=61,offset=c(120,180))
  plot(mesh1)
  summary(mesh1)

  # Make SPDE based on mesh
  spde=inla.spde2.matern(mesh1, alpha=3/2)
  n= max(as.numeric(as.factor(as.character(subdat$Station)))) # unique stations  
  
  ## To think about: treat each row as unique location?  
  subdat$yearID = match(subdat$Year,as.numeric(names(table(subdat$Year))))
  subdat$stationID = as.numeric(as.factor(subdat$Station))

  k = max(subdat$yearID)
  ymat = matrix(NA,n,k)
  ymat01 = matrix(NA,n,k)
  cent.log.depth = matrix(NA,n,k)
  cent.temp = matrix(NA,n,k)
  yearF = matrix(0,n,k)
  for(j in 1:dim(subdat)[1]) {
    # ymat is for use with positive model
    if(subdat[j,"sum.CPUE"]>0) ymat[subdat$stationID[j],subdat$yearID[j]] = subdat[j,"sum.CPUE"] # species always in col 10
    # ymat01 is used by presence-absence model 
    ymat01[subdat$stationID[j],subdat$yearID[j]] = ceiling(subdat[j,"sum.CPUE"]/1.0e10) # convert to 0/1
    
    # include depth, depth2, temp, temp2 as fixed effects
    cent.temp[subdat$stationID[j],subdat$yearID[j]] = Covar$cent.temp[j]
    cent.log.depth[subdat$stationID[j],subdat$yearID[j]] = Covar$cent.log.depth[j]
    # Include year as fixed effect design matrix
    yearF[subdat$stationID[j],subdat$yearID[j]] = 1
  }
  
  # make depth and temperature values for unobserved locations 
  # (calculate mean for each station and replace missing values with average covariate value)
  temp.depth = rowMeans(cent.log.depth,na.rm=T)
  temp.temp	= rowMeans(cent.temp,na.rm=T)
  
  for(j in 1:n){
    cent.log.depth[j,][is.na(cent.log.depth[j,])==T] = temp.depth[j]
    cent.temp[j,][is.na(cent.temp[j,])==T] = temp.temp[j]
  }
  
  # Specify which dataset to use
  if(model == "binomial") {
    z = ymat01
  }
  if(model != "binomial") {
    z = ymat
  }
  
  dat = data.frame(y=as.vector((z)), time=rep(1:k, each=n), xcoo=rep(subcoords[,1], k),
  		ycoo=rep(subcoords[,2], k), #cent.temp = as.vector(cent.temp), cent.temp2 = as.vector(cent.temp^2),
 		 cent.log.depth=as.vector(cent.log.depth),cent.log.depth2=as.vector(cent.log.depth^2))

	dat.proj = data.frame(y=rep(NA,nrow(dat.proj.trim)*k),
		time=rep(1:k, each=nrow(dat.proj.trim)), 
		xcoo=rep(dat.proj.trim$LonUTMAlbers,k),
		ycoo=rep(dat.proj.trim$LatUTMAlbers,k), #cent.temp = as.vector(cent.temp), cent.temp2 = as.vector(cent.temp^2),
  		cent.log.depth=as.vector(log(dat.proj.trim$NGDC24_M)-mean(Covar$log.depth,na.rm=T)),
  		cent.log.depth2=as.vector(log(dat.proj.trim$NGDC24_M)-mean(Covar$log.depth,na.rm=T))^2)

  # tack on year effects
  YEARS = paste("Y",names(table(subdat$Year)),sep="")
  dat[YEARS] = 0
  dat.proj[YEARS] = 0

	# Make a design matrix where the first year is the intercept
  dat[,YEARS[1]] = 1
  dat.proj[,YEARS[1]] = 1

  for(j in 1:length(YEARS)){
		dat[dat$time == j,YEARS[j]]	<-	1 
 		dat.proj[dat.proj$time == j,YEARS[j]]	<-	1
  }
      
  iset = inla.spde.make.index("i2D", n.spde=mesh1$n, n.group = k) 
  
  # Make the covariates
  X.1 = dat[,-c(1:4)]
  Covar.names <- colnames(X.1)
  XX.list <- as.list(X.1)
  effect.list <- list()    				
  #   effect.list[[1]] <- c(iset, list(Intercept=1))
  effect.list[[1]] <- c(iset)
  for (Z in 1:ncol(X.1)) effect.list[[Z+1]] <- XX.list[[Z]]
  names(effect.list) <- c("1", Covar.names)

  # Make the covariates for the prediction matrix
  X.1.proj  = dat.proj[,-c(1:4)]
  Covar.names <- colnames(X.1.proj)
  XX.list <- as.list(X.1.proj)
  effect.list.proj <- list()    				
  effect.list.proj[[1]] <- c(iset)
  for (Z in 1:ncol(X.1.proj)) effect.list.proj[[Z+1]] <- XX.list[[Z]]
  names(effect.list.proj) <- c("1", Covar.names)
  
  ### Make projection points stack.
  A <- inla.spde.make.A(mesh=mesh1, loc=cbind(dat$xcoo, dat$ycoo), group = dat$time)
  A.list = list()
  A.list[[1]] = A
  for (Z in 1:ncol(X.1)) A.list[[Z+1]] <- 1

  ### Make projection points stack.
  A.proj <- inla.spde.make.A(mesh=mesh1, loc=cbind(dat.proj$xcoo, dat.proj$ycoo), group = dat.proj$time)
  A.list.proj = list()
  A.list.proj[[1]] = A.proj
  for (Z in 1:ncol(X.1.proj)) A.list.proj[[Z+1]] <- 1
  
  ### Make projection points stack.
  Ntrials <- dat$y
  Ntrials[is.na(Ntrials)==F]<-1

  Ntrials.proj <- dat.proj$y

  sdat 		<- inla.stack(tag='stdata', data=list(y=dat$y,link=1,Ntrials=Ntrials), A=A.list, effects=effect.list)
  sdat.proj <- inla.stack(tag='pred', data=list(y=dat.proj$y,link=1,Ntrials=NA), A=A.list.proj, effects=effect.list.proj)
 
 	sdat.all	<-	inla.stack(sdat,sdat.proj) 
  
  formula = as.formula(paste0("y ~ -1 +",  paste(Covar.names, collapse="+"), "+ f(i2D, model=spde, group = i2D.group, control.group = list(model='ar1'))"))		# field evolves with AR1 by year

    if(model=="binomial") {

	     inlaModel <- inla(formula, family = "binomial", data=inla.stack.data(sdat),
    				control.predictor=list(compute=TRUE, A=inla.stack.A(sdat)), 
    				verbose = TRUE, debug=TRUE, keep=FALSE,
    				control.compute = list(dic=TRUE, cpo=TRUE), 
    				control.fixed = list(correlation.matrix=TRUE))
	     Output	<-	 list( Data=df,INLA.mod=inlaModel,Covar=Covar,Mesh=mesh1,iset=iset)
	     save(Output,file=paste(species[i],"_binomial.Rdata",sep=""))# Save this fitted thing to a workspace	     
    }
    if(model!="binomial") {
    
       inlaModel <- inla(formula, family = "gamma", data=inla.stack.data(sdat),
            control.predictor=list(compute=TRUE, A=inla.stack.A(sdat)), 
            verbose = TRUE, debug=TRUE, keep=FALSE,
            control.compute = list(dic=TRUE, cpo=TRUE), 
            control.fixed = list(correlation.matrix=TRUE))
       Output	<-	 list( Data=df,INLA.mod=inlaModel,Covar=Covar,Mesh=mesh1,iset=iset)
       save(Output,file=paste(species[i],"_pos.Rdata",sep=""))# Save this fitted thing to a workspace
    
    }    
}