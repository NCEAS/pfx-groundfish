library(INLA)
library(rgdal)
library(sp)
library(ggplot2)

#### GO GET THE OBSERVED TRAWL DATA

myWD = "/Users/eric.ward/documents/exxonValdez_nceas/goaTrawl/"

setwd(myWD)
df = read.csv("goa_500trawls_albers.csv")
df<-df[order(df$Year,df$Lat),]

### Remove NA entries in BottomDepth and Bottom Temp for now
df<-df[is.na(df$BottomDepth)==F,]
df<-df[is.na(df$BottomTemp)==F,]

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
model = "positive"
species = names(df)[10:dim(df)[2]]
speciesList = read.csv("species names.csv")
species = species[which(speciesList[,2]==1)]
comm = speciesList$Common.Name[which(speciesList[,2]==1)]

#par(mfrow =c(6,5),mai=c(0.1,0.2,0.1,0.1))
#for(i in 61:90) {#length(species)) {
#	subdat = df[,c(1:9,which(names(df)==species[i]))]
#	subdat[which(is.na(subdat[,10])),10] = 0
#	agg = aggregate(ceiling(subdat[,10]/1.0e10),by=list(subdat$Year),mean,na.rm=T)
#	plot(agg[,1],agg[,2],main=species[i],xlab="",ylab="",type="l",cex.main=0.6)
#}
fitModel = TRUE

for(i in 1:length(species)) {
  # Fit the model for species XX
  
  subdat = df[,c(1:9,which(names(df)==species[i]))]
  #Center the covariates 
  Covar	<- subdat[,1:9]
  Covar$log.depth			<-	log(Covar$BottomDepth)
  Covar$cent.log.depth	<-	Covar$log.depth - mean(Covar$log.depth,na.rm=T)
  Covar$cent.log.depth.2	<-	Covar$log.depth^2 
  Covar$cent.temp			<-	Covar$BottomTemp - mean(Covar$BottomTemp)
  Covar$cent.temp.2		<-	Covar$cent.temp^2 
  setwd(myWD)
	#call basic plotting routine for raw data
	  #if(model == "binomial"){
	  #	source("trawl_plot_binom.r")
	  #}
	#call basic plotting routine for raw data
	  if(model == "positive"){
	  	source("trawl_plot_positive.r")
	  }

  # Grab X-Y coords in UTM space
  subcoords = cbind(subdat$LonUTMAlbers[match(unique(subdat$Station),subdat$Station)],subdat$LatUTMAlbers[match(unique(subdat$Station),subdat$Station)])
  bnd = inla.nonconvex.hull(subcoords, convex=80)
  mesh1 = inla.mesh.2d(boundary=bnd,max.edge=c(100,1200),cutoff=90)
  # "cutoff" parameter is used to avoid building many small triangles around clustered input locations, 
  # "offset" species the size of the inner and outer extensions around the data locations,
  # "max.edge" species the maximum allowed triangle edge lengths in the inner domain and in the outer extension.
	
  #mesh1 = inla.mesh.2d(subcoords,max.edge=c(150,150),cutoff=50)
  #plot(mesh1)
  #points(subcoords,col="red")
  
  # Make SPDE based on mesh
  spde=inla.spde2.matern(mesh1, alpha=3/2)
  n= max(as.numeric(as.factor(as.character(subdat$Station)))) # unique stations
  subdat$yearID = match(subdat$Year,as.numeric(names(table(subdat$Year))))
  subdat$stationID = as.numeric(as.factor(subdat$Station))
  
  k = max(subdat$yearID)
  ymat = matrix(NA,n,k)
  ymat01 = matrix(0,n,k)
  cent.log.depth = matrix(NA,n,k)
  cent.temp = matrix(NA,n,k)
  #yearF = matrix(0,n,k)
  for(j in 1:dim(subdat)[1]) {
    if(is.na(subdat[j,10])==F & subdat[j,10]>0) ymat[subdat$stationID[j],subdat$yearID[j]] = subdat[j,10] # species always in col 10
    if(is.na(subdat[j,10])==F) ymat01[subdat$stationID[j],subdat$yearID[j]] = ceiling(subdat[i,10]/1.0e10) # convert to 0/1
    # include depth, depth2, temp, temp2 as fixed effects
    cent.temp[subdat$stationID[j],subdat$yearID[j]] = Covar$cent.temp[j]
    cent.log.depth[subdat$stationID[j],subdat$yearID[j]] = Covar$cent.log.depth[j]
    # Include year as fixed effect design matrix
    #yearF[subdat$stationID[j],subdat$yearID[j]] = 1
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
    z = log(ymat)
  }
  #dat <- data.frame(y=as.vector((z)), time=rep(1:k, each=n), xcoo=rep(subcoords[,1], k),ycoo=rep(subcoords[,2], k), cent.temp = as.vector(cent.temp), cent.temp2 = as.vector(cent.temp^2),
  #cent.log.depth=as.vector(cent.log.depth),cent.log.depth2=as.vector(cent.log.depth^2))
  dat <- data.frame(y=as.vector((z)), time=rep(1:k, each=n), xcoo=rep(subcoords[,1], k),ycoo=rep(subcoords[,2], k),cent.log.depth=as.vector(cent.log.depth),cent.log.depth2=as.vector(cent.log.depth^2))
  if(fitModel==TRUE) {  
  # Make a design matrix where the first year is the intercept, tack on year effects
  YEARS <- paste("Y",names(table(subdat$Year)),sep="")
  dat[YEARS] = 0	
  dat[,YEARS[1]]	 <- 1
  for(j in 1:length(YEARS)){
		dat[dat$time == j,YEARS[j]]	<-	1
  }
    
  iset = inla.spde.make.index("i2D", n.spde=mesh1$n, n.group = k)  
  
  # Make the covariates
  X.1 = dat[,-c(1:4)]
  Covar.names <- colnames(X.1)
   XX.list <- as.list(X.1)
   effect.list <- list()						
   effect.list[[1]] <- c(iset, list(Intercept=1))
   for (j in 1:ncol(X.1)) effect.list[[j+1]] <- XX.list[[j]]
   names(effect.list) <- c("1", Covar.names)

   A <- inla.spde.make.A(mesh=mesh1, loc=cbind(dat$xcoo, dat$ycoo),group = dat$time)
   A.list = list()
   A.list[[1]] = A
   for (j in 1:ncol(X.1)) A.list[[j+1]] <- 1
   sdat <- inla.stack(tag='stdata', data=list(y=dat$y), A=A.list, effects=effect.list)

  formula = as.formula(paste0("y ~ -1 +",  paste(Covar.names, collapse="+"), "+ f(i2D, model=spde, group = i2D.group, control.group = list(model='ar1'))"))		# field evolves with AR1 by year

  if(model=="binomial") {
    inlaModel <- inla(formula, family = "binomial", data=inla.stack.data(sdat),control.predictor=list(compute=TRUE, A=inla.stack.A(sdat)), verbose = FALSE, debug=FALSE, keep=FALSE,control.compute = list(dic=TRUE, cpo=TRUE), control.fixed = list(correlation.matrix=TRUE))
    save.image(paste(species[i],"_binomial.Rdata",sep=""))# Save this fitted thing to a workspace
  }
  if(model!="binomial") {
    inlaModel <- inla(formula, family = "gaussian", data=inla.stack.data(sdat),control.predictor=list(compute=TRUE, A=inla.stack.A(sdat)), verbose = TRUE, debug=TRUE, keep=FALSE,control.compute = list(dic=TRUE, cpo=TRUE), control.fixed = list(correlation.matrix=TRUE),control.results=list(return.marginals.random=F))
    #save.image(paste(species[i],"_pos.Rdata",sep=""))# Save this fitted thing to a workspace
        #inlaModel <- inla(formula, family = "gamma", data=inla.stack.data(sdat),control.predictor=list(compute=TRUE, A=inla.stack.A(sdat)), verbose = TRUE, debug=TRUE, keep=FALSE,control.compute = list(dic=TRUE, cpo=TRUE), control.fixed = list(correlation.matrix=TRUE), control.results=list(return.marginals.random=F,return.marginals.fixed=F),control.inla(lincomb.derived.correlation.matrix = TRUE))
Output<-list( Data=df,INLA.mod = inlaModel,Covar=Covar,Mesh=mesh1,Mesh=mesh1,iset=iset)
save(Output,file=paste(species[i],"_pos.Rdata",sep=""))# Save this fitted thing to a workspace
   }
  }  
  
}

pdf("Annual_fixed_effects_pos_allSpecies.pdf")
par(mfrow = c(8,4),mai=c(0.2,0.2,0.08,0.02),mgp=c(2,0.5,0))
for(i in 1:length(species)) {
  if(file.exists(paste(species[i],"_pos.Rdata",sep=""))) {
    # Then load the file in 
    load(paste("output/",species[i],"_pos.Rdata",sep=""))
    meansd = Output[[2]]$summary.fixed[-c(1:2),1:2]
    Xs = c(1984,1987,1990,1993,1996,1999,2001,2003,2005,2007,2009,2011)
    lower = meansd[,1] - 2*meansd[,2]
    upper = meansd[,1] + 2*meansd[,2]
    plot(Xs,meansd[,1],ylim=c(min(lower),max(upper)),type="b",lwd=3,xlab="",ylab="log(CPUE)",main=comm[i],cex.main=0.6,cex.axis=0.8,col="white")
    polygon(c(Xs,rev(Xs)),c(lower,rev(upper)),border=NA,col="tomato1")
    lines(Xs,meansd[,1])
    points(Xs,meansd[,1])
    lines(c(1989,1989),c(-10,10),lty=3)
    #lines(Xs,lower,col="grey")
    #lines(Xs,upper,col="grey")
  }
}
dev.off()


xmean <- list()
projgrid <- inla.mesh.projector(mesh1, xlim = c(-1350,1150), ylim=c(5450,6200),dims = c(1150--1350+1,6200-5450+1))
   coefs.pres = res1$summary.fixed[,1]
  
   # Loop over k years. Use inla.mesh.project to predict response to new locations on grid, and because covars
   # don't exist, hold them constant at their means 
   for (j in 1:k) {
      xmean[[j]] <- inla.mesh.project(projgrid, Output[[2]]$summary.random$i2D$mean[Output[[6]]$i2D.group==j]) 
   }

# look at 1990-1987
xmean2 = list()
xmean2[[1]] = xmean[[3]]-xmean[[2]]
 require(gridExtra)
 require(splancs)
 require(lattice)
 do.call(function(...) grid.arrange(..., nrow=1),
 lapply(xmean2, levelplot, xlab='', ylab='',
 col.regions=topo.colors(16), scale=list(draw=FALSE)))

z = melt(xmean2[[1]])
ggplot(z) + stat_contour()
