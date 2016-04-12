rm(list=ls())
library(INLA)
library(rgdal)
library(ggplot2)
library(splancs)
library(sp)

# Use only species with length-based data, or total biomass
totalBiomass = FALSE
#### CHOOSE A MODEL - "binomial" or "positive"
model = "binomial"
single.intercept =TRUE

#Switch for doing cross-validation (TRUE) or not (FALSE)
LEAVE.OUT   <- FALSE
frac.leave.out  <- 0.10

max.depth <- 600 # in meters

#Define directories for the data and for the plots
proj.dir = paste(getwd(),"/goaTrawl",sep="")
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

### Remove NA entries in BottomDepth
df = df[df$BottomDepth != -9999,]
#df = df[df$BottomTemp != -9999,]

df = df[df$BottomDepth <= max.depth,] 

df$Station = as.character(df$Station)
df$Year = as.numeric(as.character(df$Year))
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
# x.lim = c(-200,200)
# x.lim = c(min(dat.project$LonUTMAlbers),max(dat.project$LonUTMAlbers))
# x.lim = c(-800,550)
# y.lim = c(5700,6200)
# y.lim = c(min(dat.project$LatUTMAlbers),max(dat.project$LatUTMAlbers))
# plot(LatUTMAlbers~LonUTMAlbers,data=dat.project[ #&
# 		 dat.project$LonUTMAlbers < 200 & dat.project$LonUTMAlbers > -200
# 		 ,],pch=".",xlim=x.lim,ylim=y.lim)
# par(new=T)
# plot(LatUTMAlbers~LonUTMAlbers,data=dat.trim[dat.trim$NGDC24_M<=500 &
# 	 dat.trim$LonUTMAlbers < 550 & dat.trim$LonUTMAlbers > -800 &
# 	 dat.trim$NGDC24_M > 25
# 	 ,],pch=".",xlim=x.lim,ylim=y.lim,col=2)
# 
# par(new=T)
# plot(LatUTMAlbers~LonUTMAlbers,data=dat.new.trim[dat.new$LonUTMAlbers < 200 & dat.new$LonUTMAlbers > -200,]
#  ,pch=".",xlim=x.lim,ylim=y.lim,col=4)
# dim(dat.new.trim)
# dat.proj.trim	= dat.new.trim

# This is a switch for whether we're using full data, or length-stratified data
maxCol = dim(df)[2]
species = names(df)[17:maxCol]
nCovCol = 16
if(totalBiomass==FALSE) {
  maxCol = 31
  species = names(df)[14:maxCol] 
  temp <- NULL
  for(ZZ in 1:length(species)){
      temp[ZZ] <- substr(species[ZZ],nchar(species[ZZ])-3,nchar(species[ZZ]))
  }
  species <- species[temp =="cpue"]
  nCovCol = 13
}

######################################################################################################################
#################################### START INLA LOOP ACROSS SPECIES
######################################################################################################################
for(i in 9:9) {
  SPECIES <- species[i]
  # Fit the model for species XX
  subdat = df[,c(1:nCovCol,which(names(df)==species[i]))]
  
    if(totalBiomass == FALSE){
      subdat[subdat[,SPECIES] == -9999,SPECIES ] <- NA
    }
    rows2drop = which(is.na(subdat$LonUTMAlbers+subdat$LatUTMAlbers))
   if(length(rows2drop) > 0)subdat = subdat[-rows2drop,]

    if(totalBiomass==FALSE){
      subdat <- subdat[is.na(subdat[,species[i]])==F,] 
    }
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
    subdat[which(subdat[,species[i]]==0),species[i]] = NA
  }

  #Center the covariates 
  Covar	<- subdat[,1:nCovCol]
  Covar$log.depth = log(Covar$BottomDepth)
  Covar$cent.log.depth = Covar$log.depth - mean(Covar$log.depth,na.rm=T)
  Covar$cent.log.depth.2 = Covar$cent.log.depth^2 
  Covar$cent.temp = Covar$BottomTemp - mean(Covar$BottomTemp)
  Covar$cent.temp.2 = Covar$cent.temp^2 
  
  # Add covariate columns to subdat
  subdat$cent.log.depth   <- Covar$cent.log.depth
  subdat$cent.log.depth.2 <- Covar$cent.log.depth.2
  ## 
  subdat$yearID = match(subdat$Year,as.numeric(names(table(subdat$Year))))

  #call basic plotting routine for raw data
	  if(model == "binomial"){
	  	source("R Scripts/trawl_plot_binom.r") 
	  }
	#call basic plotting routine for raw data
	  if(model == "positive"){
	  	source("R Scripts/trawl_plot_positive.r")
	  }

  # Grab X-Y coords in UTM space
  subcoords = cbind(subdat$LonUTMAlbers,subdat$LatUTMAlbers)

  bnd = inla.nonconvex.hull(subcoords, convex=80)
  #bnd = inla.nonconvex.hull(subcoords, convex=150)
    # increase cutoff to ~ 150 to create much coarser mesh
  #mesh1 = inla.mesh.2d(boundary=bnd,max.edge=c(60,1500),cutoff=150,offset=c(120,180))
  #mesh1 = inla.mesh.2d(boundary=bnd,max.edge=c(60,1500),cutoff=65,offset=c(120,180))
  mesh1 = inla.mesh.2d(boundary=bnd,max.edge=c(60,1500),cutoff=59,offset=c(110,180))
  plot(mesh1)
  summary(mesh1)

  if(model=="positive"){
    subdat.all<- subdat
    subdat   <-   subdat[is.na(subdat[,SPECIES])==F,]
  }
  
  # Pull leave out for predictions
  if(LEAVE.OUT == TRUE){
    ### Split subdat into training and prediction sets
    # predict to frac.leave.out from each year
  
    count.year <- aggregate(subdat$Year,by=list(Year=subdat$Year),length)
    count.year$numb.leave.out <- round(frac.leave.out*count.year$x)
  
    subdat.pred   <- NULL
    subdat.train  <- NULL
    for(QQQ in 1:nrow(count.year)){
      temp  <- subdat[subdat$Year == count.year$Year[QQQ],]
      THESE  <- sort(sample(1:count.year$x[QQQ],count.year$numb.leave.out[QQQ],replace=F))
      temp.pred <- temp[THESE,]
      temp.train   <- temp[-THESE,]
      
      subdat.pred <- rbind(subdat.pred,temp.pred)
      subdat.train <- rbind(subdat.train,temp.train)
    }
  subdat<-subdat.train
 } # End of leave out if statement
  
  # Make SPDE based on mesh
  spde=inla.spde2.matern(mesh1, alpha=3/2)
  n=nrow(subdat)

    # Specify which dataset to use
  if(model == "binomial") {
    z = subdat[,species[i]]
  }
  if(model != "binomial") {
    z = subdat[,species[i]]
  }
  
  YEARS<-unique(subdat$Year)
  k=length(unique(subdat$Year))
  YEARS.lab = paste("Y",names(table(subdat$Year)),sep="")
  
    dat = data.frame(y=as.vector((z)), time=subdat$yearID, xcoo=subdat$LonUTMAlbers,
  	  	ycoo=subdat$LatUTMAlbers, 
 		    cent.log.depth=subdat$cent.log.depth,cent.log.depth2=subdat$cent.log.depth.2)
  
  # Make a design matrix where the first year is the intercept
  if(single.intercept==FALSE){
    # tack on year effects
    dat[YEARS.lab] = 0
    dat[,YEARS.lab[1]] = 1
    for(j in 1:length(YEARS.lab)){
      dat[dat$time == j,YEARS.lab[j]]	<-	1 
    }
  }
  if(single.intercept==TRUE){
      # tack on year effects
      YEARS.lab = paste("Y",names(table(subdat$Year)),sep="")
      dat[,YEARS.lab[1]] = 1
  }
    
  if(LEAVE.OUT == TRUE){  ### LEAVE OUT LOOP
      dat.pred = data.frame(y=rep(NA,nrow(subdat.pred)), time=subdat.pred$yearID, xcoo=subdat.pred$LonUTMAlbers,
                   ycoo=subdat.pred$LatUTMAlbers, 
                   cent.log.depth=subdat.pred$cent.log.depth,cent.log.depth2=subdat.pred$cent.log.depth.2)

      if(single.intercept==FALSE){
        # tack on year effects
        dat.pred[YEARS.lab] = 0
        dat.pred[,YEARS.lab[1]] = 1
        for(j in 1:length(YEARS.lab)){
          dat.pred[dat.pred$time == j,YEARS.lab[j]]	<-	1 
        }
      }
      if(single.intercept==TRUE){
        # tack on year effects
        YEARS.lab = paste("Y",names(table(subdat$Year)),sep="")
        dat.pred[,YEARS.lab[1]] = 1
      }
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

  ### Make data stack.
  A <- inla.spde.make.A(mesh=mesh1, loc=cbind(dat$xcoo, dat$ycoo), group = dat$time)
  A.list = list()
  A.list[[1]] = A
  for (Z in 1:ncol(X.1)) A.list[[Z+1]] <- 1

  ### Make projection points stack.
  Ntrials <- rep(1,length(dat$y))

  sdat 		<- inla.stack(tag='stdata', data=list(y=dat$y,link=1,Ntrials=Ntrials), A=A.list, effects=effect.list)

  #################################################################################################################################  
  if(LEAVE.OUT == TRUE){  ### LEAVE OUT LOOP
    
  # Make the covariates
  X.1.pred = dat.pred[,-c(1:4)]
  Covar.names.pred <- colnames(X.1.pred)
  XX.list.pred <- as.list(X.1.pred)
  effect.list.pred <- list()    				
  #   effect.list[[1]] <- c(iset, list(Intercept=1))
  effect.list.pred[[1]] <- c(iset)
  for (Z in 1:ncol(X.1.pred)) effect.list.pred[[Z+1]] <- XX.list.pred[[Z]]
  names(effect.list.pred) <- c("1", Covar.names.pred)
  
  ### Make data stack.
  A.pred <- inla.spde.make.A(mesh=mesh1, loc=cbind(dat.pred$xcoo, dat.pred$ycoo), group = dat.pred$time)
  A.list.pred = list()
  A.list.pred[[1]] = A.pred
  for (Z in 1:ncol(X.1)) A.list.pred[[Z+1]] <- 1
  
  ### Make projection points stack.
  Ntrials.pred <- rep(1,length(dat.pred$y))
  
  sdat.pred 		<- inla.stack(tag='stdata.pred', data=list(y=dat.pred$y,link=1,Ntrials=Ntrials.pred), A=A.list.pred, effects=effect.list.pred)
  sdat <- inla.stack(sdat,sdat.pred)
  } 
  
 formula = as.formula(paste0("y ~ -1 +",  paste(Covar.names, collapse="+"), 
                "+ f(i2D, model=spde, group = i2D.group, control.group = list(model='ar1'))"))		# field evolves with AR1 by year

    if(model=="binomial") {
       # config = TRUE included so that inla.posterior.sample() can be applied
       # lincomb.derived.correlation.matrix included to give Var-Cov matrix of random effects
	     inlaModel <- inla(formula, family = "binomial", data=inla.stack.data(sdat),
    				control.predictor=list(compute=TRUE, A=inla.stack.A(sdat)), 
    				verbose = TRUE, debug=TRUE, keep=FALSE,
    				control.compute = list(dic=TRUE, cpo=TRUE, config=TRUE), 
    				control.fixed = list(correlation.matrix=TRUE),
    				control.inla = list(lincomb.derived.correlation.matrix=TRUE))
	     if(LEAVE.OUT == FALSE){Output	<-	 list( Data=subdat,INLA.mod=inlaModel,Covar=Covar,Mesh=mesh1,iset=iset,single.intercept=single.intercept)}
	     if(LEAVE.OUT == TRUE){Output	<-	 list( Data=subdat,Pred=subdat.pred,INLA.mod=inlaModel,Covar=Covar,Mesh=mesh1,iset=iset,single.intercept=single.intercept)}
	     save(Output,file=paste(species[i],"_binomial_;leave.out=",LEAVE.OUT,"; sing_int=",single.intercept,".Rdata",sep=""))# Save this fitted thing to a workspace	     
	     source(paste(proj.dir,"/R Scripts/Plot model diagnostics pres.R",sep=""))
	      }

  if(model!="binomial") {
       # config = TRUE included so that inla.posterior.sample() can be applied
       # lincomb.derived.correlation.matrix included to give Var-Cov matrix of random effects
       inlaModel <- inla(formula, family = "gamma", data=inla.stack.data(sdat),
            control.predictor=list(compute=TRUE, A=inla.stack.A(sdat)), 
            verbose = TRUE, debug=TRUE, keep=FALSE,
            control.compute = list(dic=TRUE, cpo=TRUE, config=TRUE), 
            control.fixed = list(correlation.matrix=TRUE),
            control.inla = list(lincomb.derived.correlation.matrix=TRUE))
       if(LEAVE.OUT == FALSE){Output	<-	 list( Data=subdat,All.Data=subdat.all,INLA.mod=inlaModel,Covar=Covar,Mesh=mesh1,iset=iset,single.intercept=single.intercept)}
       if(LEAVE.OUT == TRUE){Output	<-	 list( Data=subdat,All.Data=subdat.all,Pred=subdat.pred,INLA.mod=inlaModel,Covar=Covar,Mesh=mesh1,iset=iset,single.intercept=single.intercept)}
       save(Output,file=paste(species[i],"_pos_;leave.out=",LEAVE.OUT,"; sing_int=",single.intercept,".Rdata",sep=""))# Save this fitted thing to a workspace
       source(paste(proj.dir,"/R Scripts/Plot model diagnostics pos.R",sep=""))
     }    
  setwd(proj.dir)
}