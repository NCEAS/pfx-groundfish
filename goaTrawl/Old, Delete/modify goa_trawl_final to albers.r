library(rgdal)

#### GO GET THE OBSERVED TRAWL DATA
setwd("/Users/ole.shelton/GitHub/exxonValdez_nceas/goaTrawl/")
df = read.csv("goa_trawl_final.csv")

aea.proj <- "+proj=aea +lat_1=51 +lat_2=62 +lon_0=-150 +x_0=0 +y_0=0 +datum=WGS84"

new	<-	 cbind(df$Lon,df$Lat)
new.sp	<- SpatialPoints(new,CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
new.sp.albers	<- spTransform(new.sp, CRS(aea.proj))

new.dat		<- as(new.sp.albers,"data.frame")

df$LonUTMAlbers	<-	new.dat$coords.x1
df$LatUTMAlbers	<-	new.dat$coords.x2


write.csv(df,"goa_trawl_top63_final_albers.csv",row.names=F)
dat<-read.csv("goa_500trawls_albers.csv")	





###################################################

dat.new	<- dat.centroid.2[,c("MASTER_ID","LONGITUDE","LATITUDE","x","y",  "WATER","SRTM_M","NGDC24_M")]
colnames(dat.new)[4:5]	<- c("LonUTMAlbers", "LatUTMAlbers")
write.csv(dat.new,file="goa_projection_points.csv",row.names=F)