
## Mary's code for calculating species richness

SpByArea<-read.csv("All_sp_index_occurrenceByArea.csv",head=TRUE)
head(SpByArea)
length(SpByArea$area)


-------------------------------------------------------
##We want to look at diversity by area by year.  

SpByArea <- SpByArea %>%
  mutate(year = ifelse((year.numb=="1"),'1984',
                ifelse((year.numb=="2"),'1987',
                ifelse((year.numb=="3"),'1990',
                ifelse((year.numb=="4"),'1993',       
                ifelse((year.numb=="5"),'1996',
                ifelse((year.numb=="6"),'1999',
                ifelse((year.numb=="7"),'2001',
                ifelse((year.numb=="8"),'2003',
                ifelse((year.numb=="9"),'2005',
                ifelse((year.numb=="10"),'2007',
                ifelse((year.numb=="11"),'2009',
                ifelse((year.numb=="12"),'2011',
                ifelse((year.numb=="13"),'2013',
                ifelse((year.numb=="14"),'2015',        
                ""))))))))))))))
  )
head(SpByArea)

## Species richness by Area (11 areas) by Year (14 years) -- use presence/absence data here rather than mean densities.
## On hold until we figure out if we can actually calculate this with the data we have.
SpByArea<-SpByArea[SpByArea$area!="Total",]
year.numb<-SpByArea$year.numb
allSummary <- data.frame()


  for(i in c(1:14)){
    for(j in c(1:11)){
      spRichSum<-sum(SpByArea$Mean.avgPresence[SpByArea$year.numb==i&SpByArea$area==j])
      spRichOut<-as.data.frame(cbind(i,j,spRichSum))
      colnames(spRichOut)<-c("year","area","spRich")                                    
      allSummary <- rbind(allSummary, spRichOut)                             
      #write.csv(allSummary,"GoA_spRichByArea.csv")           
    }
  }

allSummary$area<-as.factor(allSummary$area)
SR <- ggplot(allSummary, aes(year, spRich, color=area)) + geom_line() + #facet_wrap(~area,ncol=2)
      scale_color_hue(breaks=c("1","2","3","4","5","6","7","8","9","10","11"))
SR


## Eric's code for calculating diversity metrics using occurrence data 
# Calculate average occurrence, by year in one core area. 
# Use MCMC samples to represent uncertainty
meanOccurrence = apply(projectedLatentGrid,c(2,3),mean)
occurrenceStats = (apply(meanOccurrence,2,quantile,c(0.025,0.5,0.975)))
occAll[,ii] = occurrenceStats[2,]
plot(seq(1981,2011), occurrenceStats[2,],type="b",lwd=2,ylab="Median occurrence (black=model,red=raw)",ylim=range(occurrenceStats),xlab="",main=this.spp)
lines(seq(1981,2011),occurrenceStats[1,],col="grey30",lty=3)
lines(seq(1981,2011),occurrenceStats[3,],col="grey30",lty=3)
points(rawOcc$Group.1+1980, rawOcc$x, col="red",lwd=2)
  
# Derive some stats related to richness and diversity
pdf("Expected species richness (core area).pdf")
plot(1981:2011, apply(occAll,1,sum), xlab="", ylab="Expected species", lwd=2, type="b")
dev.off()

pdf("Shannon diversity (core area).pdf")
plot(1981:2011, -apply(occAll * log(occAll), 1, sum), xlab="", ylab="Expected species", lwd=2, type="b")
dev.off()

pdf("Simpson diversity (core area).pdf")
plot(1981:2011, 1/apply(occAll^2, 1, sum), xlab="", ylab="Expected species", lwd=2, type="b")
dev.off()

