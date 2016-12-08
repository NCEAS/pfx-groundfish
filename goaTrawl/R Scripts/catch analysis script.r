library(rgdal)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(viridis)

### Examine fish catches by stat6 area to include information 

setwd("/Users/ole.shelton/GitHub/pfx-groundfish/goaTrawl/_fishing areas gfish")
# this is the total catch of groundfish
# From Eric

# Yep -- I did it by permit, with only the "B" "C" and "M" permits. So it should be only longline and trawl, but would include a number of species as bycatch (including very small amounts of salmon, etc)
# 
# Codes here:
#   https://www.cfec.state.ak.us/misc/FshyDesC.htm
catch <- read.csv("pounds_by_stat6.csv")

##3 Read in mapping of stat areas to regions.
map.to.regions <- read.csv("/Users/ole.shelton/GitHub/pfx-commercial/GIS files/regions+OLE.csv")

catch     <- merge(catch,map.to.regions[,c("stat6","final.OLE")])
catch.sum <- catch %>% group_by(year,final.OLE) %>% summarize(tot.pound =sum(pounds))
catch.sum$met.ton <-  catch.sum$tot.pound * 0.000453592

# Region areas:
Area.est <- data.frame(
  matrix(c(
    "Alaska Peninsula",124367,
    "Cook Inlet",39443,
    "Kodiak",147333,
    "PWS",45136),
    4,2,byrow=T)
)
colnames(Area.est) <- c("final.OLE","km2")
Area.est$km2 <- as.numeric(as.character(Area.est$km2))

catch.sum <- merge(catch.sum,Area.est)
catch.sum$met.ton.km2 <- catch.sum$met.ton / catch.sum$km2
catch.sum$final.OLE <- as.character(catch.sum$final.OLE)
catch.sum$final.OLE[catch.sum$final.OLE == "PWS"] <- "Prince William Sound"
###########################################################################

COL <- viridis(4,begin=0,end=0.8)

biomass.plot <- ggplot(catch.sum) +
                  geom_point(aes(year,met.ton.km2,shape=final.OLE),size=2.5) +
                  geom_line(aes(year,met.ton.km2,group=final.OLE)) +
                  scale_shape(name="Region",solid=F) +
                  coord_cartesian(xlim=c(min(catch.sum$year)-0.75,max(catch.sum$year)+0.75),
                                  ylim=c(0,max(catch.sum$met.ton.km2)*1.05),
                                  expand=c(0)) +
               		labs(x="Year", y=expression("Catch (mt km"^-2*")")) +
                  theme_bw()


biomass.plot



quartz(file = "/Users/ole.shelton/GitHub/pfx-groundfish/goaTrawl/_fishing areas gfish/Plots/Catch by region.pdf",
       type="pdf",width=6,height=4)
  print(biomass.plot)
dev.off()


# catch[catch$stat6<240 & catch$stat6>230] <-230
# catch[catch$stat6<250 & catch$stat6>240] <-240
# catch[catch$stat6<260 & catch$stat6>250] <-250
# catch[catch$stat6<270 & catch$stat6>260] <-260
# catch[catch$stat6<280 & catch$stat6>270] <-270
# catch[catch$stat6<290 & catch$stat6>280] <-280
# catch[catch$stat6<300 & catch$stat6>290] <-290

### Read in shapefiles
setwd("/Users/ole.shelton/GitHub/pfx-commercial/GIS files/Groundfish Stat Areas")
shp.stat6	 <-	readOGR(dsn=".",layer="pvg_stat_2001")
dat.stat6   <- fortify(shp.stat6,"data.frame")

setwd("/Users/ole.shelton/GitHub/pfx-commercial/GIS files/Shapefile with region definitions")
shp.regions	 <-	readOGR(dsn=".",layer="NCEASDissolveNew")
dat.regions   <- fortify(shp.regions@data,"data.frame")


plot(shp.regions)

YEARS <- sort(unique(catch$year))
for(i in 1:length(YEARS)){
  nom <- paste(YEARS[i],"pounds",sep=".")
  
  match(shp.stat6$)
  
  
}



# summarise the catch by stat area and year.
catch %>% group_by(year)
