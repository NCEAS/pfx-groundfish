library(rgdal)
library(ggplot2)
library(RColorBrewer)

### Examine fish catches by stat6 area to include information 

setwd("/Users/ole.shelton/GitHub/pfx-groundfish/goaTrawl/_fishing areas gfish")
# this is the total catch of groundfish
# From Eric

# Yep -- I did it by permit, with only the "B" "C" and "M" permits. So it should be only longline and trawl, but would include a number of species as bycatch (including very small amounts of salmon, etc)
# 
# Codes here:
#   https://www.cfec.state.ak.us/misc/FshyDesC.htm
catch <- read.csv("pounds_by_stat6.csv")
catch[catch$stat6<240 & catch$stat6>230] <-230
catch[catch$stat6<250 & catch$stat6>240] <-240
catch[catch$stat6<260 & catch$stat6>250] <-250
catch[catch$stat6<270 & catch$stat6>260] <-260
catch[catch$stat6<280 & catch$stat6>270] <-270
catch[catch$stat6<290 & catch$stat6>280] <-280
catch[catch$stat6<300 & catch$stat6>290] <-290

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
