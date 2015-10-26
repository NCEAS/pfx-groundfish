rm(list=ls())
library(INLA)


#Define directories for the data and for the plots
plot.dir	<-	"/Users/ole.shelton/GitHub/exxonValdez_nceas/goaTrawl/Output plots/"
data.dir	<-	"/Users/ole.shelton/GitHub/exxonValdez_nceas/goaTrawl"

#### GO GET THE OBSERVED TRAWL DATA
setwd(data.dir)
df = read.csv("goa_500trawls_albers.csv")
df	<-	df[order(df$Year,df$Lat),]

### Remove NA entries in BottomDepth and Bottom Temp for now
df	<-	df[is.na(df$BottomDepth)==F,]
df	<-	df[is.na(df$BottomTemp)==F,]

df$Station = as.character(df$Station)
df$Year = as.numeric(df$Year)
df$LonUTMAlbers <- df$LonUTMAlbers/1000
df$LatUTMAlbers <- df$LatUTMAlbers/1000

#### GO GET THE PROJECTION POINTS
dat.project	<- read.csv("goa_projection_points.csv")

	#### Exclude points that end up on land.
	dat.project$NGDC24_M	<-	-dat.project$NGDC24_M	
	dat.project$SRTM_M		<-	-dat.project$SRTM_M	
	dat.project				<-	dat.project[dat.project$NGDC24 > 0,]

#### CHOOSE A MODEL
model = "binomial"
species = names(df)[10:dim(df)[2]]

#################################### START INLA LOOP
for(i in 1:length(species)) {
  # Fit the model for species XX
   subdat = df[,c(1:9,which(names(df)==species[i]))]
 	#get rid of NAs in subdat
 	subdat[is.na(subdat[,species[i]])==T,species[i]]	<-	0
 	#if binomial replace positive values with 1
 	if(model =="binomial"){
 		subdat[subdat[,species[i]]>0,species[i]] <- 1
 	}
 
  #Center the covariates 
  Covar	<- subdat[,1:9]
  Covar$log.depth			<-	log(Covar$BottomDepth)
  Covar$cent.log.depth		<-	Covar$log.depth - mean(Covar$log.depth,na.rm=T)
  Covar$cent.log.depth.2	<-	Covar$cent.log.depth^2 
  Covar$cent.temp			<-	Covar$BottomTemp - mean(Covar$BottomTemp)
  Covar$cent.temp.2			<-	Covar$cent.temp^2 
  
	#call basic plotting routine for raw data
	  if(model == "binomial"){
	  	source("trawl_plot_binom.r")
	  }
	#call basic plotting routine for raw data
	  if(model == "positive"){
	  	source("trawl_plot_positive.r")
	  }

  # Grab X-Y coords in UTM space
  subcoords = cbind(subdat$LonUTMAlbers[match(unique(subdat$Station),subdat$Station)],subdat$LatUTMAlbers[match(unique(subdat$Station),subdat$Station)])
  if(model == "binomial"){
  	# Make a nonconvex hull for the spatial stuff
	bnd = inla.nonconvex.hull(subcoords, convex=80)
	mesh1 = inla.mesh.2d(boundary=bnd,max.edge=c(100,1200),cutoff=90)
	# "cutoff" parameter is used to avoid building many small triangles around clustered input locations, 
	# "offset" species the size of the inner and outer extensions around the data locations,
	# "max.edge" species the maximum allowed triangle edge lengths in the inner domain and in the outer extension.
  }

   plot(mesh1)
   summary(mesh1)
#    points(subcoords,col="red")
  
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
    if(subdat[j,10]>0) ymat[subdat$stationID[j],subdat$yearID[j]] = subdat[j,10] # species always in col 10
    
    ymat01[subdat$stationID[j],subdat$yearID[j]] = subdat[j,10] # convert to 0/1
    # include depth, depth2, temp, temp2 as fixed effects
    cent.temp[subdat$stationID[j],subdat$yearID[j]] = Covar$cent.temp[j]
    cent.log.depth[subdat$stationID[j],subdat$yearID[j]] = Covar$cent.log.depth[j]
    # Include year as fixed effect design matrix
    yearF[subdat$stationID[j],subdat$yearID[j]] = 1
  }

	# make depth and temperature values for unobserved locations 
	# (calculate mean for each station and replace missing values with average covariate value)
	
	temp.depth	<- rowMeans(cent.log.depth,na.rm=T)
	temp.temp	<- rowMeans(cent.temp,na.rm=T)

	for(j in 1:n){
		cent.log.depth[j,][is.na(cent.log.depth[j,])==T]	<-	 temp.depth[j]
		cent.temp[j,][is.na(cent.temp[j,])==T]	<-	 temp.temp[j]
	}
  
  if(model == "binomial") {
    z = ymat01
  }
  if(model != "binomial") {
    z = ymat
  }
  dat <- data.frame(y=as.vector((z)), time=rep(1:k, each=n), xcoo=rep(subcoords[,1], k),
  ycoo=rep(subcoords[,2], k), cent.temp = as.vector(cent.temp), cent.temp2 = as.vector(cent.temp^2),
  cent.log.depth=as.vector(cent.log.depth),cent.log.depth2=as.vector(cent.log.depth^2))

  # tack on year effects
  YEARS <- paste("Y",names(table(subdat$Year)),sep="")
  dat[YEARS] = 0
#   dat[paste("Y",names(table(subdat$Year)),sep="")] = NA

	# Make a design matrix where the first year is the intercept
  dat[,YEARS[1]]	<-	1
  for(j in 1:length(YEARS)){
		dat[dat$time == j,YEARS[j]]	<-	1
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
  
  A <- inla.spde.make.A(mesh=mesh1, loc=cbind(dat$xcoo, dat$ycoo), group = dat$time)
  A.list = list()
  A.list[[1]] = A
  for (Z in 1:ncol(X.1)) A.list[[Z+1]] <- 1
  
  Ntrials <- dat$y
  Ntrials[is.na(Ntrials)==F]<-1
  
  sdat <- inla.stack(tag='stdata', data=list(y=dat$y,link=1,Ntrials=Ntrials), A=A.list, effects=effect.list)
  
  formula = as.formula(paste0("y ~ -1 +",  paste(Covar.names, collapse="+"), "+ f(i2D, model=spde, group = i2D.group, control.group = list(model='ar1'))"))		# field evolves with AR1 by year

  if(model=="binomial") {
    inlaModel <- inla(formula, family = "binomial", data=inla.stack.data(sdat),
    				control.predictor=list(compute=TRUE, A=inla.stack.A(sdat)), 
    				verbose = TRUE, debug=TRUE, keep=FALSE,
    				control.compute = list(dic=TRUE, cpo=TRUE), 
    				control.fixed = list(correlation.matrix=TRUE))
   save.image(paste(species[i],"_binomial.Rdata",sep=""))# Save this fitted thing to a workspace
  }
  if(model!="binomial") {
    inlaModel <- inla(formula, family = "gamma", data=inla.stack.data(sdat),control.predictor=list(compute=TRUE, A=inla.stack.A(sdat)), verbose = TRUE, debug=TRUE, keep=FALSE,control.compute = list(dic=TRUE, cpo=TRUE), control.fixed = list(correlation.matrix=TRUE))
    save.image(paste(species[i],"_pos.Rdata",sep=""))# Save this fitted thing to a workspace
  }  
  
}