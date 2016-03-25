# calculate groundfish temperature distributions
# Colette Ward, March 21 2016

library(dplyr)
library(mgcv)
library(ggplot2)

# lines 14-127 are data loading / prep from Ole Shelton's script "Process AFSC trawl data.r"


#### GO GET THE OBSERVED TRAWL DATA
setwd("~/Google Drive/GoA project/pfx-groundfish/goaTrawl/Raw Trawl Data")
# http://www.afsc.noaa.gov/RACE/groundfish/survey_data/data.htm
dat1 		<- read.csv("goa1984_1987.csv")
dat2 		<- read.csv("goa1990_1999.csv")
dat3 		<- read.csv("goa2001_2005.csv")
dat4 		<- read.csv("goa2007_2013.csv")
dat5 		<- read.csv("goa2015.csv")

dat     <- rbind(dat1,dat2)
dat     <- rbind(dat,dat3)
dat     <- rbind(dat,dat4)
dat     <- rbind(dat,dat5)

dat[is.na(dat)==T]	<-	-9999	# make sure missing data is -9999
# Get rid of some missing data in the haul 
dat     <- dat[dat$HAUL != -9999,]

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
  
  print(YEARS[i])
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

all.mod[,"Lepidopsetta sp."]		 <-	all.mod[,"Lepidopsetta sp."] + all.mod[,"Lepidopsetta polyxystra"] + all.mod[,"Lepidopsetta bilineata"]
all.mod[,"Dusky and Dark Rockfish"] <- all.mod[,"dusky and dark rockfishes unid."] + all.mod[,"Sebastes variabilis"] + all.mod[,"Sebastes ciliatus"]
all.mod[,"Rougheye and Blackspotted Rockfish"] <- all.mod[,"rougheye and blackspotted rockfish unid."] + all.mod[,"Sebastes aleutianus"] + 
  all.mod[,"Sebastes melanostictus"]

drops	<- c("Lepidopsetta polyxystra","Lepidopsetta bilineata",
           "dusky and dark rockfishes unid.","Sebastes variabilis","Sebastes ciliatus",
           "rougheye and blackspotted rockfish unid.","Sebastes aleutianus", "Sebastes melanostictus")

THESE	<- which(match(colnames(all.mod),drops,nomatch=0)>0)
all.mod	<-	all.mod[,-THESE]


#########################################################################################
#########################################################################################

# Fit GAMs to CPUE distributions across temperature

myTemps <- all.mod %>%
  dplyr::select(-Lat, -Lon, -Station, -Year, -Stratum, -BottomDepth, -SurfTemp, -Cruise, -Haul, -sum.CPUE) %>%
  filter(BottomTemp != -9999.0) %>%
  group_by(BottomTemp) %>%
  summarize_each(funs(sum)) %>% # sum CPUEs for each 0.1 degC temperature bin
  rename(AtheresthesStomias = `Atheresthes stomias`, DuskyDarkRockfish = `Dusky and Dark Rockfish`)
#View(myTemps)


# visualize distributions for a few species:
plot(log(myTemps$AtheresthesStomias + 1) ~ myTemps$BottomTemp, pch=16)
plot(myTemps$`Gadus chalcogrammus` ~ myTemps$BottomTemp, pch=16)
plot(myTemps$`Gadus macrocephalus` ~ myTemps$BottomTemp, pch=16)
plot(myTemps$`Hippoglossus stenolepis` ~ myTemps$BottomTemp, pch=16)
plot(myTemps$`Hippoglossoides elassodon` ~ myTemps$BottomTemp, pch=16)
plot(log(myTemps$`Raja binoculata`) ~ myTemps$BottomTemp, pch=16)
plot(log(myTemps$DuskyDarkRockfish) ~ myTemps$BottomTemp, pch=16)




# start with trial model for Dusky & Dark Rockfish
qqnorm(log(myTemps$DuskyDarkRockfish + 1)) # still not great
ggplot(myTemps, aes(log(DuskyDarkRockfish+1))) + geom_density()

# general model structure:
# species_mod <- gam(species ~ s(BottomTemp, k=4) + te(lat, lon), weights = 1/n ... , data = myTemps)

# conversation with B Williams:
# sampling is not all in the same area, so try accounting for lat/long
#lat & long are scaled differently; te (thin plate smoother or tensor product) 
#dont constrain te with k 
# consider weighting for inverse of number of samples at each temperature (1/n)
# temperature with max cpue is not informative because highest cpue will be a temperature sampled most (cpue is already adjusted for effort, but effort at temp is not necessarily accounted for ...)
# try a poisson distribution for CPUE instead of gaussian + log link


plot(myTemps$DuskyDarkRockfish ~ myTemps$BottomTemp, pch = 16, col = "red")
plot(log(myTemps$DuskyDarkRockfish) ~ myTemps$BottomTemp, pch = 16, col = "red")

dusky <- gam(myTemps$DuskyDarkRockfish +1 ~ s(myTemps$BottomTemp, k=3), data = myTemps, family=poisson)
plot(dusky$fitted.values ~ myTemps$BottomTemp)
summary(dusky)
plot(dusky, pages=1)
dusky$coefficients
gam.check(dusky)


