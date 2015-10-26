### Basic plotting files

NAME	<-	paste(plot.dir,"/",species[i], "Pres-Abs Preliminary plots.pdf",sep="")
pdf(NAME,onefile=TRUE,width=15,5)

## Plot against the covariates 
ZZ	<-	subdat[,species[i]]

JIT = 0.15
par(mfrow=c(1,2))
plot(y=ZZ+runif(nrow(subdat),-JIT,JIT),x=subdat$BottomTemp,xlab="Bottom Temperature (C)",ylab="jittered occurrence",axes=F)
axis(1)
axis(2,las=2,at=c(0,1),labels=c("Abs","Pres"))
box(bty="o",lwd=2)
title(species[i])

plot(y=ZZ+runif(nrow(subdat),-JIT,JIT),x=log(subdat$BottomDepth),xlab="log(Depth (m))",ylab="jittered occurrence",axes=F)
axis(1)
axis(2,las=2,at=c(0,1),labels=c("Abs","Pres"))
box(bty="o",lwd=2)
title(species[i])

#### PLOT SPATIAL DISTRIBUTION
# Import the datafile of Alaska shoreline

setwd(paste(plot.dir,"/_Alaska Shapefile",sep=""))
shp.alaska	 <-	readOGR(dsn=".",layer="Alaska-Albers")
dat.alaska   <- fortify(shp.alaska,"data.frame")

plot.subdat	<-	subdat
# THESE	<-	which(plot.subdat[,species[i]]>0)
# plot.subdat[THESE,species[i]]	<- 1
# THESE	<-	is.na(plot.subdat[,species[i]])==T
# plot.subdat[THESE,species[i]]	<- 0
plot.subdat[,species[i]]	<- as.factor(plot.subdat[,species[i]])

# Some Arguments for making prettier plots
bGrid <-theme(panel.grid =element_blank())
bBack <-theme(panel.background =element_blank())
bAxis <-theme(axis.title.y =element_blank())
bTics <-theme(axis.text =element_blank(), axis.text.y =element_blank(), axis.ticks =element_blank())


### Plot presence-absence in space
q1	<-	ggplot() +
 		scale_size_manual(values=c(0.8,0.8)) +
 		scale_colour_manual(values=c("black","red"),
 						 name="",
                         labels=c("Absent", "Present"))+
    	geom_point(data=plot.subdat,alpha=0.4,
    			aes(LonUTMAlbers*1000,LatUTMAlbers*1000,colour=eval(as.name(species[i]))))+
		geom_polygon(data=dat.alaska, fill=grey(0.4),color=NA,aes(long,lat,group=group)) +
		labs(x = "Eastings",y="Northings",title=paste(species[i],"All Years Occurrence")) +
   		coord_cartesian(xlim = c(min(subdat$LonUTMAlbers*1000)*1.01,max(subdat$LonUTMAlbers*1000)), 
   				ylim = c(min(subdat$LatUTMAlbers*1000),6300000))+
   		bGrid  + bBack

YEARS <- sort(unique(subdat$Year))
p	<-	list()

for(j in 1:length(YEARS)){
p[[j]]	<-	ggplot() +
 		scale_colour_manual(values=c("black","red"),
 						 name="", labels=c("Absent", "Present"))+
    geom_point(data=plot.subdat[which(plot.subdat$Year==YEARS[j] & is.na(subdat$LonUTMAlbers)==FALSE & is.na(subdat$LatUTMAlbers)==FALSE),],alpha=0.4,
             aes(LonUTMAlbers*1000,LatUTMAlbers*1000,colour=eval(as.name(species[i]))))+
		geom_polygon(data=dat.alaska, fill=grey(0.4),color=NA,aes(long,lat,group=group)) +
		labs(x = "Eastings",y="Northings",title=paste(species[i],YEARS[j],"Occurrence")) +
   		coord_cartesian(xlim = c(min(subdat$LonUTMAlbers*1000,na.rm=T)*1.01,max(subdat$LonUTMAlbers*1000,na.rm=T)), 
   				ylim = c(min(subdat$LatUTMAlbers*1000,na.rm=T),6300000))+
   		bGrid  + bBack
}

## write plots to file
print(q1)
for(j in 1:length(YEARS)){
	print(p[j])
}

dev.off()

setwd(proj.dir)
