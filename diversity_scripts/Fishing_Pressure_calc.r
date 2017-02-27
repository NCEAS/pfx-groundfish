##############################################################################################
###  Script for Fishing Pressure Analysis                                                  ###
###  Originally written by Ole Shelton for PFX Groundfish paper #1                         ###
###       located pfx-groundfish/Scripts and plots for pubs/Catch and Temperature Script.r ###
###  Modified by Rachael E. Blake in January 2017 for                                      ###
###       groundfish paper #2 fishing pressure analysis                                    ###
##############################################################################################

# load libraries needed by this script
library(rgdal)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(viridis)

### Examine fish catches by stat6 area to include information 
catch <- read.csv("../goaTrawl/_fishing areas gfish/pounds_by_stat6.csv")
# this is the total catch of groundfish
# From Eric
# Yep -- I did it by permit, with only the "B" "C" and "M" permits. So it should be only longline and
# trawl, but would include a number of species as bycatch (including very small amounts of salmon, etc)
# Codes here:  #   https://www.cfec.state.ak.us/misc/FshyDesC.htm

### Read in mapping of stat areas to regions.
map_to_regions <- read.csv("../diversity-data/regions+OLE.csv")

catch1     <- merge(catch, map_to_regions[,c("stat6","final.OLE")])
catch_sum <- catch1 %>% 
             group_by(year,final.OLE) %>% 
             summarize(tot_pound =sum(pounds))
catch_sum$met_ton <-  catch_sum$tot_pound * 0.000453592

### Calculating the area (km2) of each fishery Region:
Area_est <- data.frame(
              matrix(c(
              "Alaska Peninsula",124367,
              "Cook Inlet",39443,
              "Kodiak",147333,
              "PWS",45136),
              4,2,byrow=T)
            )
colnames(Area_est) <- c("final.OLE","km2")
Area_est$km2 <- as.numeric(as.character(Area_est$km2))

catch_sum <- merge(catch_sum,Area_est)
catch_sum$met_ton_km2 <- catch_sum$met_ton / catch_sum$km2
catch_sum$final.OLE <- as.character(catch_sum$final.OLE)
catch_sum$final.OLE[catch_sum$final.OLE == "PWS"] <- "Prince William Sound"

###########################################################################
# Plotting time series of catch by region
COL <- viridis(4,begin=0,end=0.8)
biomass_plot <- ggplot(catch_sum) +
                  geom_point(aes(x=year, y=met_ton_km2, shape=final.OLE), size=2.5) +
                  geom_line(aes(x=year, y=met_ton_km2, group=final.OLE)) +
                  scale_shape(name="Region",solid=F) +
                  coord_cartesian(xlim=c(min(catch_sum$year)-0.75,max(catch_sum$year)+0.75),
                                  ylim=c(0,max(catch_sum$met_ton_km2)*1.05),
                                  expand=c(0)) +
               		labs(x="Year", y=expression("Catch (mt km"^-2*")")) +
                  theme_bw()
#biomass_plot



