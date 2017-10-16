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

###################################################
### MAKE MY OWN THEME TO SAVE LINES OF CODE
theme_boxplot <- function(base_size = 12){
  theme_bw(base_size) %+replace%
    theme(legend.key.size=unit(13,"points"),
          legend.text=element_text(size=I(11)),
          legend.key=element_blank(),
          legend.title=element_blank(),
          #legend.position="none",
          plot.margin=unit(c(0.75,1,0.75,1), "cm"), # respectively: top, right, bottom, left; refers to margin *outside* labels; default is c(1,1,0.5,0.5)
          panel.border=element_blank(),
          panel.spacing=unit(0,"lines"),
          axis.ticks.length=unit(1,"mm"),
          axis.text.x = element_text(margin=margin(5,0,0,0)),
          axis.text.y = element_text(margin=margin(0,5,0,0)),
          axis.text=element_text(size=11),
          axis.title.x=element_text(size=13, margin=margin(15,0,0,0)), 
          axis.title.y=element_text(size=13, angle=90, margin=margin(0,15,0,0)), 
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          strip.text.x=element_text(size=14),
          strip.background=element_rect(colour='black', fill='white'),
          #axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line = element_line(colour = 'black', size=0.5, linetype='solid')
          )
}
####################################################


### Examine fish catches by stat6 area to include information 
catch <- read.csv("../goaTrawl/_fishing areas gfish/pounds_by_stat6.csv")
# this is the total catch of groundfish
# From Eric
# Yep -- I did it by permit, with only the "B" "C" and "M" permits. So it should be only longline and
# trawl, but would include a number of species as bycatch (including very small amounts of salmon, etc)
# Codes here:  #   https://www.cfec.state.ak.us/misc/FshyDesC.htm

### Read in mapping of stat areas to regions.
map_to_regions <- read.csv("../diversity-data/regions+OLE.csv")

catch1 <- merge(catch, map_to_regions[,c("stat6","final.OLE")])

catch_sum <- catch1 %>% 
             group_by(year, final.OLE) %>% 
             summarize(tot_pound = sum(as.numeric(pounds)))

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

catch_sum <- merge(catch_sum, Area_est)
catch_sum$met_ton_km2 <- catch_sum$met_ton / catch_sum$km2
catch_sum$final.OLE <- as.character(catch_sum$final.OLE)
catch_sum$final.OLE[catch_sum$final.OLE == "PWS"] <- "Prince William Sound"

###########################################################################
# OLE's Figure for Groundfish Paper #1
# Plotting time series of catch by region
COL <- viridis(4,begin=0,end=0.8)
catch_sum$final.OLE1 <- factor(catch_sum$final.OLE, levels=c('Alaska Peninsula', 'Kodiak', 
                                                             'Cook Inlet', 'Prince William Sound')) # for ordering the plot

biomass_plot <- ggplot(catch_sum) +
                #geom_point(aes(x=year, y=met_ton_km2, shape=final.OLE1), size=2.5) +
                geom_line(aes(x=year, y=met_ton_km2, group=final.OLE1, color=final.OLE1), size=1.5) +
                scale_shape(solid=F) +
                coord_cartesian(xlim=c(min(catch_sum$year)-0.75,max(catch_sum$year)+0.75),
                                ylim=c(0,max(catch_sum$met_ton_km2)*1.05),
                                expand=c(0)) +
               	labs(x="Year", y=expression("Catch (mt km"^-2*")")) +
                theme_boxplot() + 
                scale_color_manual(values=c("#49006A","#DD3497","#FA9FB5","#E3C9C6")) + 
                theme(axis.text=element_text(size=15),
                      axis.title.x=element_text(size=15),
                      axis.title.y=element_text(size=15),
                      legend.text=element_text(size=15))
#biomass_plot
#############################################################################

# GLM Analysis

fp <- glm(met_ton_km2 ~ final.OLE, data=catch_sum, family=gaussian)
summary(fp)






