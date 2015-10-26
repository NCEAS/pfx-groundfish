rm(list=ls())
library(INLA)
library(ggplot2)
library(rgdal)
library(MASS)

ANTI.LOGIT	<- function(log.mu){
					A<-	1/(1+exp(-log.mu))
					return(A)
			}


#Define directories for the data and for the plots
plot.dir	<-	"/Users/ole.shelton/GitHub/exxonValdez_nceas/goaTrawl/_Output plots Pres Abs/"
data.dir	<-	"/Users/ole.shelton/GitHub/exxonValdez_nceas/goaTrawl/_Results Pres Abs/"

#### GO GET THE MODEL OBJECT OF INTEREST
setwd(data.dir)
SPECIES	<-	"Theragra chalcogramma"
sp.name <-	"Theragrachalcogramma"
NAME	<-	"Theragrachalcogramma_slow_binomial.Rdata" #"Anoplopomafimbria_slow_binomial.Rdata"

# SPECIES	<-	"Anoplopoma fimbria"
# sp.name <-	"Anoplopomafimbria"
# NAME	<-	"Anoplopomafimbria_slow_binomial.Rdata"

# SPECIES	<-	"Atheresthes stomias"
# sp.name <-	"Atheresthesstomias"
# NAME	<-	"Atheresthesstomias_slow_binomial.Rdata"

 SPECIES	<-	"Gadus macrocephalus"
 sp.name 	<-	"Gadusmacrocephalus"
 NAME		<-	"Gadusmacrocephalus_slow_binomial.Rdata"

 SPECIES	<-	"Sebastescalutus"
 sp.name 	<-	"Sebastesalutus"
 NAME		<-	"Sebastesalutus_slow_binomial.Rdata"



load(NAME)

#### FIX ME 
Data 	 	<- Output$Data
inla.mod 	<- Output$INLA.mod
Covar		<- Output$Covar
Mesh		<- Output$Mesh
iset		<- Output$iset

#### GO GET THE PROJECTION POINTS
setwd(data.dir)
dat.project	<- read.csv("goa_projection_points.csv")

	#### Exclude points that end up on land.
	dat.project$NGDC24_M	<-	-dat.project$NGDC24_M	
	dat.project$SRTM_M		<-	-dat.project$SRTM_M	
	dat.project				<-	dat.project[dat.project$NGDC24 > 25,]

	dat.project.trim		<- dat.project[dat.project$NGDC24_M < 1200,]
	dim(dat.project.trim)

## MAKE THE COVARIATE VALUES FOR THE PROJECTION LOCATIONS
cent.log.depth	<-	log(dat.project.trim$NGDC24_M) - mean(Covar$log.depth)
cent.temp		<-	0

### GO GET THE ALASKA MAP
setwd(paste(plot.dir,"_Alaska Shapefile",sep=""))
shp.alaska	 <-	readOGR(dsn=".",layer="Alaska-Albers")
dat.alaska   <- fortify(shp.alaska,"data.frame")

THESE	<-	seq(1,nrow(dat.project.trim),by=3)
dat.project.trim	<-	dat.project.trim[THESE,]

#### Summarize Covariance for field and fixed effects
    COV.MAT	<- inla.mod$misc$lincomb.derived.covariance.matrix
    FIXED	<-	inla.mod$summary.fixed

##########################################################################################
##########################################################################################
##########################################################################################
## Random Field Plots
##########################################################################################
##########################################################################################
##########################################################################################
setwd(plot.dir)

 # make the projected grid
	LOC	<-	cbind(dat.project.trim[,"LonUTMAlbers"], dat.project.trim[,"LatUTMAlbers"]) / 1000
   	projgrid <- inla.mesh.projector(Mesh, loc=LOC)

   # Loop over years. Use inla.mesh.project to predict response to new locations on grid,   
   xmean	<- list()
   xsd		<- list()
	for(Z in 1:max(iset$i2D.group)){
			xmean[[Z]] 	<- inla.mesh.project(projgrid, inla.mod$summary.random$i2D$mean[iset$i2D.group==Z])
    		xsd[[Z]]	<-	inla.mesh.project(projgrid, inla.mod$summary.random$i2D$sd[iset$i2D.group==Z])
    }
    
	## TEST PLOTS
    # Some Arguments for making prettier plots
	bGrid <-theme(panel.grid =element_blank())
	bBack <-theme(panel.background =element_blank())
	bAxis <-theme(axis.title.y =element_blank())
	bTics <-theme(axis.text =element_blank(), axis.text.y =element_blank(), axis.ticks =element_blank())

### Plot Random Field for each year
q	<- list()
	z.lim	<-	max(c(abs(c(min(unlist(xmean)),max(unlist(xmean))))))
YEARS	<- sort(unique(Covar$Year))
for(Z in 1: length(YEARS)){
    PROJ	<-	data.frame(LOC,xmean[[Z]],xsd[[Z]])
    colnames(PROJ)	<- c("LonUTMAlbers","LatUTMAlbers","xmean","xsd")

	q[[Z]]	<-	ggplot() +
#  		scale_size_manual(values=c(0.8,0.8)) +
  		scale_colour_gradientn(limits=c(-z.lim, z.lim),"Random Field", colours = c("blue",grey(0.8),"red") )+
    	geom_point(data=PROJ,alpha=0.4,
    			aes(LonUTMAlbers*1000,LatUTMAlbers*1000,colour=xmean))+
		geom_polygon(data=dat.alaska, fill=grey(0.4),color=NA,aes(long,lat,group=group)) +
		labs(x = "Eastings",y="Northings",title=paste(SPECIES,YEARS[Z])) +
   		coord_cartesian(xlim = c(min(Data$LonUTMAlbers*1000)*1.01,max(Data$LonUTMAlbers*1000)), 
   				ylim = c(min(Data$LatUTMAlbers*1000),6300000))+
   		bGrid  + bBack
} 

     ## write random field plot to file
	pdf(paste(SPECIES,"Random Field.pdf"),onefile=TRUE,width=15,5)
	for(j in 1:length(YEARS)){
		print(q[[j]])
	}
	dev.off()
    
##########################################################################################
##########################################################################################
##########################################################################################
## Make Difference in Random field Plot
##########################################################################################
##########################################################################################
##########################################################################################
  
Diff	<-	list()

for(Z in 1: (length(YEARS)-1)){
    Diff[[Z]]	<-	xmean[[Z+1]] - xmean[[Z]]
}   

w	<-list()  
	z.lim	<-	max(c(abs(c(min(unlist(Diff)),max(unlist(Diff))))))
for(Z in 1: (length(YEARS)-1)){    
    PROJ.diff	<-	data.frame(LOC,Diff[[Z]])
    colnames(PROJ.diff)	<- c("LonUTMAlbers","LatUTMAlbers","Diff")

	w[[Z]]	<-	ggplot() +
  		scale_colour_gradientn(limits=c(-z.lim, z.lim),"Random Field", colours = c("red",grey(0.9),"green") )+
    	geom_point(data=PROJ.diff,alpha=0.4,
    			aes(LonUTMAlbers*1000,LatUTMAlbers*1000,colour=Diff))+
		geom_polygon(data=dat.alaska, fill=grey(0.4),color=NA,aes(long,lat,group=group)) +
		labs(x = "Eastings",y="Northings",title=paste(SPECIES,"Diff",YEARS[Z+1],"-",YEARS[Z])) +
   		coord_cartesian(xlim = c(min(Data$LonUTMAlbers*1000)*1.01,max(Data$LonUTMAlbers*1000)), 
   				ylim = c(min(Data$LatUTMAlbers*1000),6300000))+
   		bGrid  + bBack
} 

     ## write random field difference plot to file
	pdf(paste(SPECIES,"Difference Random Field.pdf"),onefile=TRUE,width=15,5)
	for(j in 1:(length(YEARS)-1)){
		print(w[[j]])
	}
	dev.off()
    
##########################################################################################
##########################################################################################
# Plot Predicted Density including covariates + difference map 
# (Based on the point estimates of the Covariates and Random Field.
##########################################################################################
##########################################################################################
    
    # make the projected grid
	LOC	<-	cbind(dat.project.trim[,"LonUTMAlbers"], dat.project.trim[,"LatUTMAlbers"]) / 1000
   	projgrid <- inla.mesh.projector(Mesh, loc=LOC)

	cent.log.depth	<-	log(dat.project.trim[,"NGDC24_M"]) - mean(Covar$log.depth)
	
	YEARS	<- sort(unique(Covar$Year))
   # Loop over years. Use inla.mesh.project to predict response to new locations on grid,   
   xmean	<- list()
   xsd		<- list()
	for(Z in 1:max(iset$i2D.group)){
			THIS		<-	which(substr(rownames(inla.mod$summary.fixed),2,5)==YEARS[Z])
			xmean[[Z]] 	<- inla.mesh.project(projgrid, inla.mod$summary.random$i2D$mean[iset$i2D.group==Z]) +
								cent.log.depth * inla.mod$summary.fixed["cent.log.depth","mean"] +
								cent.log.depth^2 * inla.mod$summary.fixed["cent.log.depth2","mean"] +
								inla.mod$summary.fixed["Y1984","mean"] +
								inla.mod$summary.fixed[THIS,"mean"]
			if(Z==1){xmean[[Z]] = xmean[[Z]] - inla.mod$summary.fixed["Y1984","mean"]}								
    }

	Fixed.pred	<- list()
   	for(Z in 1:length(YEARS)){
	    Fixed.pred[[Z]]	<-	ANTI.LOGIT(xmean[[Z]])
	}    
    
    Diff.pred	<-	list()
    for(Z in 1: (length(YEARS)-1)){
    	Diff.pred[[Z]]	<-	Fixed.pred[[Z+1]] - Fixed.pred[[Z]]
	}   

 
    
    # Some Arguments for making prettier plots
	bGrid <-theme(panel.grid =element_blank())
	bBack <-theme(panel.background =element_blank())
	bAxis <-theme(axis.title.y =element_blank())
	bTics <-theme(axis.text =element_blank(), axis.text.y =element_blank(), axis.ticks =element_blank())
#############################
### Plot OCCURRENCE
#############################

q	<- list()
	z.lim	<-	c(0,1)
YEARS	<- sort(unique(Covar$Year))
for(Z in 1: length(YEARS)){
    PROJ	<-	data.frame(LOC,Fixed.pred[[Z]])
    colnames(PROJ)	<- c("LonUTMAlbers","LatUTMAlbers","Fixed.pred")

	q[[Z]]	<-	ggplot() +
   		scale_colour_gradientn(limits=z.lim,"Occurrence", colours = c("blue","green","yellow","orange","red") )+
#   		scale_colour_gradientn(limits=z.lim,"Occurrence", colours = topo.colors(10) )+
    	geom_point(data=PROJ,alpha=0.4,
    			aes(LonUTMAlbers*1000,LatUTMAlbers*1000,colour=Fixed.pred))+
		geom_polygon(data=dat.alaska, fill=grey(0.4),color=NA,aes(long,lat,group=group)) +
		labs(x = "Eastings",y="Northings",title=paste(SPECIES,YEARS[Z])) +
   		coord_cartesian(xlim = c(min(Data$LonUTMAlbers*1000)*1.01,max(Data$LonUTMAlbers*1000)), 
   				ylim = c(min(Data$LatUTMAlbers*1000),6300000))+
   		bGrid  + bBack
} 
    ## write random field plot to file
	pdf(paste(SPECIES,"Predicted Occur.pdf"),onefile=TRUE,width=15,5)
	for(j in 1:length(YEARS)){
		print(q[[j]])
	}
	dev.off()

#############################
### Plot Change in OCCURRENCE
#############################

q	<- list()
	z.lim	<-	max(c(abs(c(min(unlist(Diff.pred)),max(unlist(Diff.pred))))))

YEARS	<- sort(unique(Covar$Year))
for(Z in 1: (length(YEARS)-1)){
    PROJ	<-	data.frame(LOC,Diff.pred[[Z]])
    colnames(PROJ)	<- c("LonUTMAlbers","LatUTMAlbers","Diff.pred")

	q[[Z]]	<-	ggplot() +
  		scale_colour_gradientn(limits=c(-z.lim, z.lim),"Random Field", colours = c("blue",grey(0.8),"red") )+
    	geom_point(data=PROJ,alpha=0.4,
    			aes(LonUTMAlbers*1000,LatUTMAlbers*1000,colour=Diff.pred))+
		geom_polygon(data=dat.alaska, fill=grey(0.4),color=NA,aes(long,lat,group=group)) +
		labs(x = "Eastings",y="Northings",title=paste(SPECIES,"Diff",YEARS[Z+1],"-",YEARS[Z])) +
   		coord_cartesian(xlim = c(min(Data$LonUTMAlbers*1000)*1.01,max(Data$LonUTMAlbers*1000)), 
   				ylim = c(min(Data$LatUTMAlbers*1000),6300000))+
   		bGrid  + bBack
} 
    ## write random field plot to file
	pdf(paste(SPECIES,"Change in Predicted Occur.pdf"),onefile=TRUE,width=15,5)
	for(j in 1:length(YEARS)){
		print(q[[j]])
	}
	dev.off()
##########################################################################################
### Time Series of Fixed Effects:    
##########################################################################################
    
    #### Summarize Covariance for field and fixed effects
    COV.MAT	<- inla.mod$misc$lincomb.derived.covariance.matrix
    MEAN	<-	inla.mod$summary.fixed$mean

    # Generate simulation machine for the fixed effects
    sim.dat	<-	mvrnorm(500,MEAN,COV.MAT)
		# get rid of temp data if necessary
# 	sim.dat	<-	sim.dat[,-c(1,2)]

 
    cent.log.depth	<-	seq(3-mean(Covar$log.depth),max(Covar$cent.log.depth),length.out=1000)
    Fixed.mean	<-	data.frame(cent.log.depth=cent.log.depth,
    								log.depth=cent.log.depth + mean(Covar$log.depth))
    Fixed.quant	<-	data.frame(cent.log.depth=cent.log.depth,
    								log.depth=cent.log.depth + mean(Covar$log.depth))
    
    for(Q in 1: length(YEARS)){
    	X.mat	<- cbind(cent.log.depth,cent.log.depth^2)
    	X.year	<- matrix(0,nrow(X.mat),length(YEARS))
    	X.year[,1]	<-	1
    	X.year[,Q]	<-	1
    	X.mat	<-	cbind(X.mat,X.year)
    
    	Pred	<- X.mat %*% t(sim.dat)
    
 		pred.mean	<-  rowMeans(Pred)   
    	pred.quant	<-	as.matrix(t(apply(Pred,1,quantile,probs=c(0.05,0.95))))
    	Fixed.mean[,paste("Y",YEARS[Q],".mean",sep="")]	<-	pred.mean
    	Fixed.quant		<-	data.frame(Fixed.quant, pred.quant[,1], pred.quant[,2])
    	colnames(Fixed.quant)[(ncol(Fixed.quant)-1):ncol(Fixed.quant)]	<- c(paste("Y",YEARS[Q],"q.05",sep=""),paste("Y",YEARS[Q],"q.95",sep=""))
    }
    	Fixed.mean.exp	<-	data.frame(Fixed.mean[,c(1,2)],ANTI.LOGIT(Fixed.mean[,-c(1,2)]))
    	Fixed.quant.exp	<-	data.frame(Fixed.quant[,c(1,2)],ANTI.LOGIT(Fixed.quant[,-c(1,2)]))

	#  Start.Plot of Depth
	pres.abs	<- data.frame(Data[,c(1:9)],sp= Data[,sp.name])
	pres.abs$sp[is.na(	pres.abs$sp == T)]	<- 0
	pres.abs$sp[pres.abs$sp > 0]	<-	1
	
	JIT= 0.15
	
	#### PLOT STARTS HERE
	pdf(paste(SPECIES,"Marginals.pdf"),onefile=TRUE)
	
	par(mfrow=c(4,3),mar=c(3,3,1,0.5))
	
	for(Z in 1:length(YEARS)){
		x.lim = c(3,7)
		y.lim= c(-0.15,1.15)
		plot(y=pres.abs$sp[pres.abs$Year == YEARS[Z]]+runif(nrow(pres.abs[pres.abs$Year == YEARS[Z],]),-JIT,JIT),
			x=log(pres.abs$BottomDepth[pres.abs$Year == YEARS[Z]]),
			ylab="Occurrence",axes=F,xlim=x.lim,ylim=y.lim)
		axis(1)
		axis(2,las=2,at=seq(0,1,by=0.2))
	# 	axis(2,las=2,at=c(0,1),labels=c("Abs","Pres"))
		box(bty="o",lwd=2)
		title(paste(SPECIES,YEARS[Z]),line=0.1)
		title(xlab="log(Depth (m))",line=1.75)
	
		par(new=T)
		THIS.ONE	<-	which(substr(colnames(Fixed.mean.exp),2,5)==YEARS[Z])
		plot(y=Fixed.mean.exp[,THIS.ONE],x=Fixed.mean.exp[,"log.depth"],xlim=x.lim,ylim=y.lim,
			axes=F,type="l",xlab="",ylab="",col=2,lty=1,lwd=2)
    	THESE	<-	which(substr(colnames(Fixed.quant.exp),2,5)==YEARS[Z])
		par(new=T)
		plot(y=Fixed.quant.exp[,THESE[1]],x=Fixed.quant.exp[,"log.depth"],xlim=x.lim,ylim=y.lim,
			axes=F,type="l",xlab="",ylab="",col=2,lty=2)
		par(new=T)
		plot(y=Fixed.quant.exp[,THESE[2]],x=Fixed.quant.exp[,"log.depth"],xlim=x.lim,ylim=y.lim,
			axes=F,type="l",xlab="",ylab="",col=2,lty=2)
	}

	#### Plot intercept at mean(log(Depth))

	  DEPTH	=	200
  	  cent.log.depth<-	rep(log(DEPTH)-mean(Covar$log.depth),length.out=1000)
  	  Fixed.mean	<-	data.frame(cent.log.depth=cent.log.depth,
  	  								log.depth=cent.log.depth + mean(Covar$log.depth))
  	  Fixed.quant	<-	data.frame(cent.log.depth=cent.log.depth,
    								log.depth=cent.log.depth + mean(Covar$log.depth))

 		sim.dat.2	<-	sim.dat[,-c(1,2)]
		temp	<- ANTI.LOGIT(sim.dat.2[,1])
		Int		<- NULL
		Int		<-	rbind(Int,c(mean(temp),as.matrix(quantile(temp,probs=c(0.05,0.25,0.75,0.95)))))

		for(Q in 2:length(YEARS)){
			temp	<- ANTI.LOGIT(sim.dat.2[,1] + sim.dat.2[,Q])
			Int		<- rbind(Int,c(mean(temp),as.matrix(quantile(temp,probs=c(0.05,0.25,0.75,0.95)))))
		}
		colnames(Int)	<-	c("Mean","q.05","q.25","q.75","q.95")
		Int				<- data.frame(Year=YEARS,Int)

		par(mfrow=c(1,1),mar=c(4,4,1,1))	
		x.lim = c(min(Int$Year),max(Int$Year))
		y.lim= c(0,1)
		plot(y=Int$Mean,x=Int$Year,
			ylab="Occurrence",axes=F,xlim=x.lim,ylim=y.lim,xlab="",type="b")
		arrows(x0=Int$Year,x1=Int$Year,y0=Int$q.25,y1=Int$q.75,lwd=3,length=0)
		arrows(x0=Int$Year,x1=Int$Year,y0=Int$q.05,y1=Int$q.95,lwd=1,length=0)
		par(new=T)
		plot(y=Int$Mean,x=Int$Year,
			ylab="Occurrence",axes=F,xlim=x.lim,ylim=y.lim,pch=21,bg="white",xlab="")
		
		axis(1,at=Int$Year)
		axis(2,las=2,at=seq(0,1,by=0.2))
	# 	axis(2,las=2,at=c(0,1),labels=c("Abs","Pres"))
		box(bty="o",lwd=2)
		title(paste(SPECIES,";",DEPTH,"m"),line=0.1)

    dev.off()	
    	
    
