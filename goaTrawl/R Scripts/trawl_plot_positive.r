
NAME	<-	paste(plot.dir,"/",species[i], "Positive Preliminary plots.pdf",sep="")

pdf(NAME,onefile=TRUE,width=15,5)

## Plot against the covariates 
ZZ	<- subdat[,species[i]]
# ZZ[is.na(ZZ)==T] <- 0

oldpar = par()

par(mfrow=c(1,2))
plot(y=ZZ,x=subdat$BottomTemp,xlab="Bottom Temperature (C)",ylab="CPUE",axes=F)
axis(1)
axis(2,las=2)
box(bty="o",lwd=2)
title(species[i])

plot(y=ZZ,x=log(subdat$BottomDepth),xlab="log(Depth (m))",ylab="CPUE",axes=F)
axis(1)
axis(2,las=2)
box(bty="o",lwd=2)
title(species[i])

par(mfrow = c(1,1))#dev.off()
par(mai = oldpar$mai)

#### PLOT SPATIAL DISTRIBUTION
# Import the datafile of Alaska shoreline
setwd(paste(plot.dir,"/_Alaska Shapefile",sep=""))
shp.alaska	 <-	readOGR(dsn=".",layer="Alaska-Albers")
dat.alaska   <- fortify(shp.alaska,"data.frame")

THESE		<-	is.na(subdat[,species[i]])==F
plot.subdat	<-	subdat[THESE,]

q95			<- round(quantile(plot.subdat[,species[i]],probs=0.95),-1)
THESE		<-	which(plot.subdat[,species[i]]>q95)
plot.subdat[THESE,species[i]]	<-	q95

THESE		<-	is.na(subdat[,species[i]])==T
plot.zeros = subdat
if(length(which(THESE==TRUE)) > 0) plot.zeros	<-	subdat[THESE,]
plot.zeros[,species[i]]	<- 0

# THESE	<-	which(plot.subdat[,species[i]]>0)
# plot.subdat[THESE,species[i]]	<- 1
# THESE	<-	is.na(plot.subdat[,species[i]])==T
# plot.subdat[THESE,species[i]]	<- 0


# Some Arguments for making prettier plots
bGrid <-theme(panel.grid =element_blank())
bBack <-theme(panel.background =element_blank())
bAxis <-theme(axis.title.y =element_blank())
bTics <-theme(axis.text =element_blank(), axis.text.y =element_blank(), axis.ticks =element_blank())


### Plot presence-absence in space
q1	<-	ggplot() +
 		scale_size_manual(values=c(0.8,0.8)) +
 		scale_colour_gradient(name="CPUE (max is 95th quantile)", low="blue", high="red",limits=c(0,max(plot.subdat[,species[i]])))+
    	geom_point(data=plot.subdat,alpha=0.4,
    			aes(LonUTMAlbers*1000,LatUTMAlbers*1000,colour=eval(as.name(species[i]))))+
    	geom_point(data=plot.zeros,alpha=0.4,shape="+",colour="black",
    			aes(LonUTMAlbers*1000,LatUTMAlbers*1000))+
		geom_polygon(data=dat.alaska, fill=grey(0.4),color=NA,aes(long,lat,group=group)) +
		labs(x = "Eastings",y="Northings",title=paste(species[i],"All Years Positive")) +
   		coord_cartesian(xlim = c(min(subdat$LonUTMAlbers*1000)*1.01,max(subdat$LonUTMAlbers*1000)), 
   				ylim = c(min(subdat$LatUTMAlbers*1000),6300000))+
   		bGrid  + bBack
q1

YEARS <- sort(unique(subdat$Year))
p	<-	list()

for(j in 1:length(YEARS)){
p[[j]]	<-	ggplot() +
# 		scale_size(range = c(1,1))+
 		scale_colour_gradient(name="CPUE (max is 95th quantile)", low="blue", high="red",limits=c(0,max(plot.subdat[,species[i]])))+
    	geom_point(data=plot.subdat[plot.subdat$Year==YEARS[j],],alpha=0.4,
    			aes(LonUTMAlbers*1000,LatUTMAlbers*1000,colour=eval(as.name(species[i]))))+
    	geom_point(data=plot.zeros[plot.zeros$Year==YEARS[j],],alpha=0.4,shape="+",colour="black",
    			aes(LonUTMAlbers*1000,LatUTMAlbers*1000))+
		geom_polygon(data=dat.alaska, fill=grey(0.4),color=NA,aes(long,lat,group=group)) +
		labs(x = "Eastings",y="Northings",title=paste(species[i],YEARS[j],"Positive")) +
   		coord_cartesian(xlim = c(min(subdat$LonUTMAlbers*1000)*1.01,max(subdat$LonUTMAlbers*1000)), 
   				ylim = c(min(subdat$LatUTMAlbers*1000),6300000))+
   		bGrid  + bBack
}

## write plots to file
print(q1)
for(j in 1:length(YEARS)){
	print(p[j])
}

dev.off()

setwd(proj.dir)