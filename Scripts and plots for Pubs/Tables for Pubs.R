# Script for making tables and calculating number in the manuscript.

# Calculate the area, location, of the regions used in my manuscript.

# go get the data
region.dat <- read.csv("/Users/ole.shelton/GitHub/pfx-groundfish/goaTrawl/Output Data/goa_discrete_areas_for_comparison(50_to_150m).csv")

region.area <- aggregate(region.dat$Area,by=list(Area=region.dat$Area),length)
region.area$km.2 = region.area$x * 4
colnames(region.area)[2] <- "n.proj.points"


write.csv(region.area,file="/Users/ole.shelton/GitHub/pfx-groundfish/Groundfish Writing/site characteristics.csv",row.names=F)