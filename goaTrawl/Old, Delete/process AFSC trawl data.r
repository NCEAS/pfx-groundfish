rm(list=ls())
library(INLA)


#Define directories for the data and for the plots
plot.dir	<-	"/Users/ole.shelton/GitHub/exxonValdez_nceas/goaTrawl/Output plots/"
data.dir	<-	"/Users/ole.shelton/GitHub/exxonValdez_nceas/goaTrawl/"

#### GO GET THE OBSERVED TRAWL DATA

# setwd(data.dir)
setwd("/Users/ole.shelton/Desktop/_TEMP/goaTrawl/")
dat = read.csv("goa_data.csv")

##	go get the species of interest list 
# setwd(paste(data.dir,"March-17-2015",sep=""))

setwd("/Users/ole.shelton/Desktop/_TEMP/goaTrawl/March-17-2015/")
dat.names	<-	read.csv("final species list.csv")
NAMES.sci	<-	dat.names$Scientific
NAMES.com	<-	dat.names$Common
#########################################################################################	
	
	
YEARS	<-	sort(unique(dat$YEAR))

all.years	<-	NULL

for(i in 1:length(YEARS)){
	dat.temp	<-	 dat[dat$YEAR == YEARS[i],]

	all.temp 	<- aggregate(dat.temp$WTCPUE,by=list(	
									Lat =	dat.temp$LATITUDE,
									Lon = 		dat.temp$LONGITUDE,
									Station = 	dat.temp$STATION,
									Year =	dat.temp$YEAR,
									DateTime =	dat.temp$DATETIME,
									Stratum = 	dat.temp$STRATUM,
								#	Common = 	dat.temp$COMMON,
								#	Scientific = dat.temp$SCIENTIFIC,
									BottomDepth = dat.temp$BOT_DEPTH,
									BottomTemp = dat.temp$BOT_TEMP,
									SurfTemp =	dat.temp$SURF_TEMP,
									Cruise	=	dat.temp$CRUISE,
									Haul	=	dat.temp$HAUL),
									sum)
	colnames(all.temp)[which(colnames(all.temp)=="x")]	<-	"sum.CPUE"
	all.temp	<-	all.temp[order(all.temp$Station),]

	for(j in 1:nrow(dat.names)){

			if(NAMES.sci[j] != ""){
				dat.sp	<-	dat.temp[dat.temp$SCIENTIFIC == as.character(NAMES.sci[j]),]
				dat.sp	<-	dat.sp[,c("LATITUDE","LONGITUDE","STATION","WTCPUE")]
				colnames(dat.sp)[4]	<-	as.character(NAMES.sci[j])
				
				all.temp	<- merge(all.temp,dat.sp, 
						by.y=c("LATITUDE","LONGITUDE","STATION"),by.x=c("Lat","Lon","Station"),all=T)
			}
			if(NAMES.sci[j] == ""){
				dat.sp	<-	dat.temp[dat.temp$COMMON == as.character(NAMES.com[j]),]
				dat.sp	<-	dat.sp[,c("LATITUDE","LONGITUDE","STATION","WTCPUE")]
				colnames(dat.sp)[4]	<-	as.character(NAMES.com[j])
				
				all.temp	<- merge(all.temp,dat.sp, 
						by.y=c("LATITUDE","LONGITUDE","STATION"),by.x=c("Lat","Lon","Station"),all=T)
			}
	}
		
	all.years	<-	rbind(all.years,all.temp)

}

#Replace NA with 0
temp	<-	all.years[,12:ncol(all.years)]
temp[is.na(temp==T)]	<-	0
all.years[,12:ncol(all.years)]	<- temp

write.csv(all.years,file="goa_trawl_top63.csv")




length(which(is.na(A[,as.character(NAMES.sci[j])])==F))



total.cpue	<-aggregate(dat$WTCPUE,by=list(	
									Common = 	dat$COMMON,
									Scientific = dat$SCIENTIFIC),
									sum)

total.occur	<-aggregate(dat$SCIENTIFIC,by=list(	
									Common = 	dat$COMMON,
									Scientific = dat$SCIENTIFIC),
									length)





	sort(unique(dat$COMMON[dat$SCIENTIFIC == ""]))
	






df	<-	df[order(df$Year,df$Lat),]
