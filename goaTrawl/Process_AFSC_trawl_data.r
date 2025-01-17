rm(list=ls())

# You need to install all the stuff associated with GDAL for this script to work.
library(rgdal)
library(mgcv)
library(INLA)
#library(dataone) 
library(httr)

#Define directories for the data and for the plots
plot.dir	<-	"/Users/ole.shelton/Documents/GitHub/exxonValdez_nceas/goaTrawl/Output plots/"
data.dir	<-	"/Users/ole.shelton/Documents/GitHub/exxonValdez_nceas/goaTrawl/Raw Trawl Data/"

#Define directories for the data and for the plots
#plot.dir	<-	"/Users/ole.shelton/GitHub/exxonValdez_nceas/goaTrawl/Output plots/"
#data.dir	<-	"/Users/ole.shelton/GitHub/exxonValdez_nceas/goaTrawl/"

#### GO GET THE OBSERVED TRAWL DATA from the GOA portal: goa.nceas.ucsb.edu 
mn_uri <- "https://goa.nceas.ucsb.edu/goa/d1/mn/v1"  ## define goa portal as DataONE member node
# mn <- MNode(mn_uri) ## needed for dataone client

pid <- "df35b.257.1"   # unique identifier for this data file
# obj <- get(mn, pid)
# dat <- read.csv(text=rawToChar(obj))  # read in file as text 
### temporary way to get the data using httr()
datGet <- GET(paste(mn_uri, "/object/", pid, sep=""))
dat <- content(datGet, as='parsed')
dat<-dat[!dat$NUMCPUE==-9999,] ## remove -9999 (NAs?)

##	Go get the species of interest list from portal: goa.nceas.ucsb.edu 
pidSIL <- "df35b.275.1"        # unique identifier for this data file
# objSIL <- get(mn, pidSIL)
# dat.names <- read.csv(text=rawToChar(objSIL))    # read in file as text 
### temporary way to get the data using httr()
dat.namesGet <- GET(paste(mn_uri, "/object/", pidSIL, sep=""))
dat.names <- content(dat.namesGet, as='parsed')



#### GO GET THE OBSERVED TRAWL DATA
 setwd(data.dir)
#setwd("/Users/ole.shelton/Desktop/_TEMP/goaTrawl/")
#dat 		<- read.csv("goa_data.csv")
#dat[is.na(dat)==T]	<-	-9999	# make sure missing data is -9999

#Summarize
#A		<- aggregate(dat$WTCPUE,by=list(cruise=dat$CRUISE,
#										lat = dat$LATITUDE,
#										long = dat$LONGITUDE,
#										vessel=dat$VESSEL,
#										haul=dat$HAUL),length)

# Go get the size data
dat.size		<-	read.csv("size_data_goa.csv")
dat.size.breaks	<-	read.csv("size_breaks_6_sp.csv")	
# Summarize
#A.size	<-	aggregate(dat.size$FREQUENCY,by=list(cruise=dat.size$CRUISE,
#										lat = dat.size$START_LATITUDE,
#										long = dat.size$START_LONGITUDE,
#										vessel=dat.size$VESSEL,
#										haul=dat.size$HAUL),length)

##	go get the species of interest list 
# setwd(paste(data.dir,"March-17-2015",sep=""))
# setwd("/Users/ole.shelton/Desktop/_TEMP/goaTrawl/March-17-2015/")

# dat.names	<-	read.csv("final species list.csv")
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
								#	DateTime =	dat.temp$DATETIME,
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

		print(dim(all.temp))
		print(length(unique(paste(dat.temp$LATITUDE,dat.temp$LONGITUDE))))
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

#########################################################################################
### Manually combine columns to combine species and species groups
#########################################################################################

all.mod	<- all.years

all.mod[,"Lepidopsettasp."]		 <-	all.mod[,"Lepidopsettasp."] + all.mod[,"Lepidopsettapolyxystra"] + all.mod[,"Lepidopsettabilineata"]
all.mod[,"DuskyandDarkRockfish"] <- all.mod[,"duskyanddarkrockfishesunid"] + all.mod[,"Sebastesvariabilis"] + all.mod[,"Sebastesciliatus"]
all.mod[,"RougheyeandBlackspottedRockfish"] <- all.mod[,"rougheyeandblackspottedrockfishunid."] + all.mod[,"Sebastesaleutianus"] + 
										all.mod[,"Sebastesmelanostictus"] + all.mod[,"Sebastes aleutianus"] + all.mod[,"Sebastes melanostictus"]

drops	<- c("Lepidopsettapolyxystra","Lepidopsettabilineata",
				"duskyanddarkrockfishesunid","Sebastesvariabilis","Sebastesciliatus",
				"rougheyeandblackspottedrockfishunid.", "Sebastesaleutianus", "Sebastesmelanostictus",
				"Sebastes aleutianus", "Sebastes melanostictus")

THESE	<- which(match(colnames(all.mod),drops,nomatch=0)>0)
all.mod	<-	all.mod[,-THESE]

########################################################################################
#### Split common species by size category using the haul-level data
########################################################################################

ID	<-	unique(dat.size$SPECIES_CODE)

for(i in 1:length(ID)){
	temp			<-	dat.size[dat.size$SPECIES_CODE == ID[i],]
	
	SPECIES			<-	as.character(dat.names$Scientific[is.na(match(dat.names$SpID,ID[i]))==F])
	temp$species	<-	SPECIES
	Break			<-	dat.size.breaks[dat.size.breaks$Species == SPECIES,"break.cm"]
	A				<-	dat.size.breaks[dat.size.breaks$Species == SPECIES,"fishbase.a"]
	B				<-	dat.size.breaks[dat.size.breaks$Species == SPECIES,"fishbase.b"]
	
	temp$per.capita.kg		<-	(A * (temp$LENGTH/10)^B) / 1000
	temp$mass.kg			<-	temp$per.capita.kg * temp$FREQUENCY
	
	temp.small		<-	temp[temp$LENGTH <= Break*10,]	
	temp.large		<-	temp[temp$LENGTH >  Break*10,]	

	small.agg		<-	aggregate(temp.small$mass.kg,by=list(Cruise=temp.small$CRUISE,
										Lat = temp.small$START_LATITUDE,
										Lon = temp.small$START_LONGITUDE,
										Haul = temp.small$HAUL),sum)
	colnames(small.agg)[5]	<-	"small.mass"
	large.agg		<-	aggregate(temp.large$mass.kg,by=list(Cruise=temp.large$CRUISE,
										Lat = temp.large$START_LATITUDE,
										Lon = temp.large$START_LONGITUDE,
										Haul = temp.large$HAUL),sum)
	colnames(large.agg)[5]	<-	"large.mass"

	both	<-	merge(small.agg,large.agg, all=T)
	both$small.mass[is.na(both$small.mass)==T]	<- 0
	both$large.mass[is.na(both$large.mass)==T]	<- 0
	both$frac.small	<- both$small.mass  / (both$small.mass + both$large.mass)
	both$frac.large	<- 1 - both$frac.small

	both$SpID <- ID[i]
	both$Scientific <- SPECIES

	#drop 2013 tows and 1985 tows
	both <- both[substr(both$Cruise,1,4) != "2013" & substr(both$Cruise,1,4) != "1985",]

	all.mod.temp	<-	data.frame(all.mod[,1:11],all.mod[,SPECIES])
	colnames(all.mod.temp)[12]	<-	SPECIES
	
	temp.2	<-	merge(all.mod.temp,both,by=c("Cruise","Lat","Lon","Haul"),all=T)

#	temp.2[is.na(temp.2[,SPECIES])==T,]

	# Get rid of hauls that appear only in the size database
	temp.3	<- temp.2[is.na(temp.2[,SPECIES])==F,] 
	temp.3$small.mass[temp.3[,SPECIES]==0]	<-	0
	temp.3$large.mass[temp.3[,SPECIES]==0]	<-	0

	# combine information for all species
	fin	<-	temp.3

	fin$small.cpue	<-	fin[,SPECIES] * fin$frac.small
	fin$large.cpue	<-	fin[,SPECIES] * fin$frac.large
	fin$small.cpue[fin[,SPECIES] == 0]	<-	0
	fin$large.cpue[fin[,SPECIES] == 0]	<-	0

	fin	<-	data.frame(fin[,1:12],fin[,c("small.cpue","large.cpue")])
	colnames(fin)[(ncol(fin)-1):ncol(fin)]	<- c(paste(SPECIES,colnames(fin)[(ncol(fin)-1):ncol(fin)],sep="."))
	if(i==1){all.sp	<-	fin}
	if(i>1){all.sp	<-	merge(all.sp,fin,all=T)}
}

### Replace missing values (NA) with -9999
all.sp[is.na(all.sp)==T]	<-	-9999


##########################################################################################
##### Add Albers Projection with standardized reference location.
##########################################################################################

aea.proj <- "+proj=aea +lat_1=51 +lat_2=62 +lon_0=-150 +x_0=0 +y_0=0 +datum=WGS84"

new				<- cbind(all.mod$Lon,all.mod$Lat)
new.sp			<- SpatialPoints(new,CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
new.sp.albers	<- spTransform(new.sp, CRS(aea.proj))
new.dat			<- as(new.sp.albers,"data.frame")
colnames(new.dat)	<- c("LonUTMAlbers","LatUTMAlbers")

all.albers	<-	data.frame(all.mod[,1:2],new.dat,all.mod[,3:ncol(all.mod)])

new				<- cbind(all.sp$Lon,all.sp$Lat)
new.sp			<- SpatialPoints(new,CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
new.sp.albers	<- spTransform(new.sp, CRS(aea.proj))
new.dat			<- as(new.sp.albers,"data.frame")
colnames(new.dat)	<- c("LonUTMAlbers","LatUTMAlbers")

all.size.albers	<-	data.frame(all.sp[,c("Lat","Lon")],new.dat,all.sp[,c(1,4:ncol(all.sp))])

##########################################################################################
### Add temperature data (spatial gam).
##########################################################################################

# Go get the projection points data frame.
setwd(data.dir)
project.dat	<-	read.csv(file="goa_projection_points.csv",header=T)

### GET RID OF NAs in the modified trawl data ("all.albers")
trawl.dat	<-	all.albers
trawl.dat2	<-	trawl.dat[is.na(trawl.dat$BottomDepth)==F,]
#trawl.dat2	<-	trawl.dat[trawl.dat$BottomDepth != -9999,]
trawl.trim	<-	trawl.dat[is.na(trawl.dat$BottomTemp)==F & 
							is.na(trawl.dat$BottomDepth)==F &
							trawl.dat$BottomTemp != -9999
							,1:12]
	
project.dat	<- project.dat[project.dat$NGDC24_M< -25,]
project.dat$log.BottomDepth	<-	log(-project.dat$NGDC24_M)
	
#### Fit a spatial GAM to the observed data for each year
YEAR	<- unique(trawl.dat2$Year)
all.dat	<- NULL

for(i in 1:length(YEAR)){
	dat.all.obs	<-	trawl.dat2[trawl.dat2$Year == YEAR[i],]
	dat.all.obs$log.BottomDepth	<- log(dat.all.obs$BottomDepth)

	dat	<-	trawl.trim[trawl.trim$Year == YEAR[i],]
	dat$log.BottomDepth	<- log(dat$BottomDepth)
	
 	out		<-	gam(BottomTemp ~ te(LonUTMAlbers,LatUTMAlbers,k=7)+s(BottomDepth),data=dat)
 	out.2	<-	gam(BottomTemp ~ te(LonUTMAlbers,LatUTMAlbers,k=7)+s(log.BottomDepth),data=dat)

	THESE	<-	c("LonUTMAlbers","LatUTMAlbers","log.BottomDepth")
	new.dat	<-	dat.all.obs[,THESE]
	new.points <- predict.gam(out.2,new.dat)
	
	dat.all.obs$Pred.bot.temp			<-	unlist(new.points)
	dat.all.obs$BottomTemp[dat.all.obs$BottomTemp == -9999]	<- NA
	dat.all.obs$Resid.bot.temp			<- dat.all.obs$BottomTemp - dat.all.obs$Pred.bot.temp

	all.dat	<- rbind(all.dat,dat.all.obs)
	
	### predict for all points in the projection set
	new.dat	<- project.dat[,c('LonUTMAlbers','LatUTMAlbers','log.BottomDepth')]
	project.dat[,paste("Bot.Temp.",YEAR[i],sep="")]		<- predict.gam(out.2,new.dat)
}

all.dat2		<-all.dat[,c(1:12,69:71,13:68)]
all.size.dat	<- merge(all.size.albers,all.dat2[,c(1:6,13:16)])

write.csv(all.dat2,file="goa_trawl_final_albers+temp.csv",row.names=F)
write.csv(all.size.dat,file="goa_trawl_final_size_albers+temp.csv",row.names=F)
write.csv(project.dat,file="goa_projection_points+temp.csv",row.names=F)
