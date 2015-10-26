#rm(list=ls())
library(ggplot2)
library(rgdal)
library(MASS)

#Define directories for the data and for the plots
setwd(paste(proj.dir,"/_Output plots Pos",sep=""))
#data.dir	<-	"/Users/ole.shelton/GitHub/exxonValdez_nceas/goaTrawl/_Results Pres Abs/"

Data 	 	  <- Output$Data
inla.mod 	  <- Output$INLA.mod
Covar		  <- Output$Covar
Mesh		  <- Output$Mesh
iset		  <- Output$iset

# Make an index to help with the inla objects
obs.index = inla.stack.index(sdat, "stdata")$data

DAT.out	<-	data.frame(cbind(Data,
                            Mean =  inla.mod$summary.fitted.values$mean[obs.index],
                            Median =  inla.mod$summary.fitted.values$"0.5quant"[obs.index],
                            SD =  inla.mod$summary.fitted.values$sd[obs.index],
                            q.025 = inla.mod$summary.fitted.values$"0.025quant"[obs.index],
                            q.975 = inla.mod$summary.fitted.values$"0.975quant"[obs.index]))

DAT.out$y     <- DAT.out[,SPECIES]
##### MAKE DIAGNOSTIC PLOTS 
pdf(paste(SPECIES,"; Pos Diagnostics; Single Int =",single.intercept,".pdf"),onefile=T)

par(mfrow=c(1,1),mar=c(1,1,1,1))
plot(1:5,1:5,type="n",xlab="",ylab="",axes=F)
text(SPECIES,x=1,y=4.5,pos=4,cex=1.5)
text(paste("Covariates:  ", paste(Covar.names[1:7], collapse=" ")),x=1,y=3.8,pos=4, cex= 0.7)
if(length(Covar.names)>7){
  text(paste(Covar.names[8:14], collapse=" "),x=1.4,y=3.6,pos=4, cex= 0.7)
}
if(length(Covar.names)>14){
  text(paste(Covar.names[15:22], collapse=" "),x=1.4,y=3.4,pos=4, cex= 0.7)
}
if(length(Covar.names)>22){
  text(paste(Covar.names[23:length(Covar.names)], collapse=" "),x=1.4,y=3.2,pos=4, cex= 0.7)
}
text(paste("N = ",sum(Ntrials),"; Pos = ",sum(is.na(Data[,SPECIES])==F)),x=1,y=3,pos=4)
text(paste("DIC:",round(inla.mod$dic$dic,2),", p.eff = ",round(inla.mod$dic$p.eff,2)),x=1,y=2.8,pos=4)
text(paste("Run Time =",round(inla.mod$cpu.used[4]/ (60*60),2),"hours"),x=1,y=2.6,pos=4)
text(paste("N.2D knots = ",spde$n.spde),x=1,y=2.3,pos=4)

##########################################################################################
### OBSERVED VERSUS PREDICTED (Regular Space)
##########################################################################################
		x.lim=c(0,max(DAT.out$Mean[is.na(DAT.out$y)==F],na.rm=T))
		y.lim=c(0,max(DAT.out$y,na.rm=T))

COL="#7F7F7F50"								# polygon color (grey 50)

par(mfrow=c(1,1),mar=c(4,4,1,0.5))
plot(y=DAT.out$y,x=DAT.out$Mean,xlab="Mean Predict",ylab="Observed",main="",axes=F,xlim=x.lim,ylim=y.lim,col=1,cex=0.8)
axis(1)
axis(2,las=2)
box(bty="o",lwd=2)
title("Observed vs. Predicted")

# Make a cheesy loess trend line
A		<-	loess(y~Mean,DAT.out,span=0.4)
XXX		<-	seq(0.01,max(DAT.out$Mean) ,length.out=1000)
B		<-	predict(A, data.frame(Mean = XXX), se = TRUE)

# NWFSC model  
par(new=T)
poly.1	<-	data.frame(XXX,c(B$fit + 1.96*B$se.fit))
poly.2	<-	data.frame(XXX,c(B$fit - 1.96*B$se.fit))
poly.2	<-	poly.2[order(poly.2[,1],decreasing=T),]
poly.all	<- rbind(poly.1,poly.2)

polygon(poly.all[,1],poly.all[,2],border=COL,col=COL)
par(new=T)
plot(XXX,B$fit,axes=F,xlab="",ylab="",col="white",type="l",xlim=x.lim,ylim=y.lim,lwd=3.5)
par(new=T)
plot(XXX,B$fit,axes=F,xlab="",ylab="",col=1,type="l",xlim=x.lim,ylim=y.lim,lwd=1.5)
abline(0,1,lty=2)

##########################################################################################
### OBSERVED VERSUS PREDICTED (Log Space)
##########################################################################################
		x.lim=c(log(min(DAT.out$Mean[is.na(DAT.out$y)==F])),log(max(DAT.out$Mean[is.na(DAT.out$y)==F])))
		y.lim=c(log(min(DAT.out$y,na.rm=T)),log(max(DAT.out$y,na.rm=T)))

COL="#7F7F7F50"								# polygon color (grey 50)

par(mfrow=c(1,1),mar=c(4,4,1,0.5))
plot(y=log(DAT.out$y[is.na(DAT.out$y)==F]),x=log(DAT.out$Mean[is.na(DAT.out$y)==F]),xlab="log(Mean Predict)",ylab="log(Observed)",main="",axes=F,xlim=x.lim,ylim=y.lim,col=1,cex=0.8)
axis(1)
axis(2,las=2)
box(bty="o",lwd=2)
title("Observed vs. Predicted")

# Make a cheesy loess trend line
temp	<-	data.frame(log.y= log(DAT.out$y[is.na(DAT.out$y)==F]),
						log.Mean =log(DAT.out$Mean[is.na(DAT.out$y)==F]))

A		<-	loess(log.y~log.Mean,data=temp,span=0.4)
XXX		<-	seq(x.lim[1],x.lim[2],length.out=1000)
B		<-	predict(A, data.frame(log.Mean = XXX), se = TRUE)

# NWFSC model  
par(new=T)
poly.1	<-	data.frame(XXX,c(B$fit + 1.96*B$se.fit))
poly.2	<-	data.frame(XXX,c(B$fit - 1.96*B$se.fit))
poly.2	<-	poly.2[order(poly.2[,1],decreasing=T),]
poly.all	<- rbind(poly.1,poly.2)

polygon(poly.all[,1],poly.all[,2],border=COL,col=COL)
par(new=T)
plot(XXX,B$fit,axes=F,xlab="",ylab="",col="white",type="l",xlim=x.lim,ylim=y.lim,lwd=3.5)
par(new=T)
plot(XXX,B$fit,axes=F,xlab="",ylab="",col=1,type="l",xlim=x.lim,ylim=y.lim,lwd=1.5)
abline(0,1,lty=2)



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
if(Output$single.intercept ==FALSE){
	for(Q in 1: length(YEARS)){
	  X.mat	<- cbind(cent.log.depth,cent.log.depth^2)
	  X.year	<- matrix(0,nrow(X.mat),length(YEARS))
	  X.year[,1]	<-	1
	  X.year[,Q]	<-	1
	  X.mat	<-	cbind(X.mat,X.year)
  
	  Pred	<- X.mat %*% t(sim.dat)
  
	  pred.mean		<-  rowMeans(Pred)   
	  pred.quant	<-	as.matrix(t(apply(Pred,1,quantile,probs=c(0.05,0.95))))
	  Fixed.mean[,paste("Y",YEARS[Q],".mean",sep="")]	<-	pred.mean
	  Fixed.quant		<-	data.frame(Fixed.quant, pred.quant[,1], pred.quant[,2])
	  colnames(Fixed.quant)[(ncol(Fixed.quant)-1):ncol(Fixed.quant)]	<- c(paste("Y",YEARS[Q],"q.05",sep=""),paste("Y",YEARS[Q],"q.95",sep=""))
	}
	Fixed.mean.exp	<-	data.frame(Fixed.mean[,c(1,2)],exp(Fixed.mean[,-c(1,2)]))
	Fixed.quant.exp	<-	data.frame(Fixed.quant[,c(1,2)],exp(Fixed.quant[,-c(1,2)]))
}
if(Output$single.intercept ==TRUE){
	  X.mat	<- cbind(cent.log.depth,cent.log.depth^2)
	  X.year	<- matrix(0,nrow(X.mat),1)
	  X.year[,1]	<-	1
	  X.mat	<-	cbind(X.mat,X.year)
  
	  Pred	<- X.mat %*% t(sim.dat)
  
	  pred.mean	<-  rowMeans(Pred)   
	  pred.quant	<-	as.matrix(t(apply(Pred,1,quantile,probs=c(0.05,0.95))))
	  Fixed.mean[,paste("Y",YEARS[1],".mean",sep="")]	<-	pred.mean
	  Fixed.quant		<-	data.frame(Fixed.quant, pred.quant[,1], pred.quant[,2])
	  colnames(Fixed.quant)[(ncol(Fixed.quant)-1):ncol(Fixed.quant)]	<- c(paste("Y",YEARS[1],"q.05",sep=""),paste("Y",YEARS[1],"q.95",sep=""))

	Fixed.mean.exp	<-	data.frame(Fixed.mean[,c(1,2)],exp(Fixed.mean[,-c(1,2)]))
		colnames(Fixed.mean.exp)[3]	<- paste("Y",YEARS[1],".mean",sep="")
	Fixed.quant.exp	<-	data.frame(Fixed.quant[,c(1,2)],exp(Fixed.quant[,-c(1,2)]))
}

####  Start.Plot of Depth
pos	<- data.frame(Data)

#### PLOT STARTS HERE
if(Output$single.intercept ==FALSE){
	par(mfrow=c(4,3),mar=c(3,3,1,0.5))
	for(Z in 1:length(YEARS)){
	  x.lim = c(3,7)
	  y.lim=c(0,quantile(DAT.out$y,probs=0.95,na.rm=T))
	  plot(y=pos[pos$Year == YEARS[Z],SPECIES],
	       x=log(pos$BottomDepth[pos$Year == YEARS[Z]]),
	       ylab="CPUE",axes=F,xlim=x.lim,ylim=y.lim)
	  axis(1)
	  axis(2,las=2)
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
}
if(Output$single.intercept ==TRUE){
	par(mfrow=c(1,1),mar=c(3,3,1,0.5))
	  x.lim = c(3,7)
	  y.lim=c(0,quantile(DAT.out$y,probs=0.95,na.rm=T))
	  plot(y=pos[pos$Year == YEARS[1],SPECIES],
	       x=log(pos$BottomDepth[pos$Year == YEARS[1]]),
	       ylab="CPUE",axes=F,xlim=x.lim,ylim=y.lim)
	  axis(1)
	  axis(2,las=2)
	  # 	axis(2,las=2,at=c(0,1),labels=c("Abs","Pres"))
	  box(bty="o",lwd=2)
	  title(paste(SPECIES,YEARS[1]),line=0.1)
	  title(xlab="log(Depth (m))",line=1.75)
  
	  par(new=T)
	  THIS.ONE	<-	which(substr(colnames(Fixed.mean.exp),2,5)==YEARS[1])
	  plot(y=Fixed.mean.exp[,THIS.ONE],x=Fixed.mean.exp[,"log.depth"],xlim=x.lim,ylim=y.lim,
	       axes=F,type="l",xlab="",ylab="",col=2,lty=1,lwd=2)
	  THESE	<-	which(substr(colnames(Fixed.quant.exp),2,5)==YEARS[1])
	  par(new=T)
	  plot(y=Fixed.quant.exp[,THESE[1]],x=Fixed.quant.exp[,"log.depth"],xlim=x.lim,ylim=y.lim,
	       axes=F,type="l",xlab="",ylab="",col=2,lty=2)
	  par(new=T)
	  plot(y=Fixed.quant.exp[,THESE[2]],x=Fixed.quant.exp[,"log.depth"],xlim=x.lim,ylim=y.lim,
	       axes=F,type="l",xlab="",ylab="",col=2,lty=2)

}

#### Plot intercept at specified Depth.

DEPTH			=	100
cent.log.depth<-	rep(log(DEPTH)-mean(Covar$log.depth),length.out=1000)
Fixed.mean	<-	data.frame(cent.log.depth=cent.log.depth,
                         log.depth=cent.log.depth + mean(Covar$log.depth))
Fixed.quant	<-	data.frame(cent.log.depth=cent.log.depth,
                          log.depth=cent.log.depth + mean(Covar$log.depth))
if(single.intercept==FALSE){
	sim.dat.2	<-	sim.dat[,-c(1,2)]
	temp	<- exp(sim.dat.2[,1])
	Int		<- NULL
	Int		<-	rbind(Int,c(mean(temp),as.matrix(quantile(temp,probs=c(0.05,0.25,0.75,0.95)))))

	for(Q in 2:length(YEARS)){
	  temp	<- exp(sim.dat.2[,1] + sim.dat.2[,Q])
	  Int	<- rbind(Int,c(mean(temp),as.matrix(quantile(temp,probs=c(0.05,0.25,0.75,0.95)))))
	}
	colnames(Int)	<-	c("Mean","q.05","q.25","q.75","q.95")
	Int				<- data.frame(Year=YEARS,Int)

	par(mfrow=c(1,1),mar=c(4,4,1,1))	
	x.lim = c(min(Int$Year),max(Int$Year))
	  y.lim=c(0,quantile(DAT.out$y,probs=0.95,na.rm=T))
	plot(y=Int$Mean,x=Int$Year,
	     ylab="CPUE",axes=F,xlim=x.lim,ylim=y.lim,xlab="",type="b")
	arrows(x0=Int$Year,x1=Int$Year,y0=Int$q.25,y1=Int$q.75,lwd=3,length=0)
	arrows(x0=Int$Year,x1=Int$Year,y0=Int$q.05,y1=Int$q.95,lwd=1,length=0)
	par(new=T)
	plot(y=Int$Mean,x=Int$Year,
	     ylab="CPUE",axes=F,xlim=x.lim,ylim=y.lim,pch=21,bg="white",xlab="")
}

if(single.intercept==TRUE){
	sim.dat.2	<-	sim.dat[,-c(1,2)]
	temp	<- exp(sim.dat.2)
	Int		<- NULL
	Int		<-	rbind(Int,c(mean(temp),as.matrix(quantile(temp,probs=c(0.05,0.25,0.75,0.95)))))

	colnames(Int)	<-	c("Mean","q.05","q.25","q.75","q.95")
	Int				<- data.frame(Year=YEARS,Int)

	par(mfrow=c(1,1),mar=c(4,4,1,1))	
	x.lim = c(min(Int$Year),max(Int$Year))
	  y.lim=c(0,quantile(DAT.out$y,probs=0.95,na.rm=T))
	plot(y=Int$Mean,x=Int$Year,
	     ylab="CPUE",axes=F,xlim=x.lim,ylim=y.lim,xlab="",type="b")
	arrows(x0=Int$Year,x1=Int$Year,y0=Int$q.25,y1=Int$q.75,lwd=3,length=0)
	arrows(x0=Int$Year,x1=Int$Year,y0=Int$q.05,y1=Int$q.95,lwd=1,length=0)
	par(new=T)
	plot(y=Int$Mean,x=Int$Year,
	     ylab="CPUE",axes=F,xlim=x.lim,ylim=y.lim,pch=21,bg="white",xlab="")
}
	axis(1,at=Int$Year,las=2)
	axis(2,las=2)
	# 	axis(2,las=2,at=c(0,1),labels=c("Abs","Pres"))
	box(bty="o",lwd=2)
	title(paste(SPECIES,";",DEPTH,"m"),line=0.1)


dev.off()	
