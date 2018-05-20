###############################################################
#####  Script for Primary Production Analysis             #####
#####  Created by Rachael E. Blake, March 2017 to test    #####
#####      the intermediate productivity hypothesis for   #####
#####      groundfish paper #2                            #####
###############################################################

# load libraries needed by this script
library(rgdal) ; library(plyr) ; library(ggplot2) ; library(reshape2)
library(scales) ; library(psych) ; library(rworldmap) ; library(rworldxtra) ;
library(sp) ; library(maptools) ; library(raster) ; library(rgeos) ;
library(maps) ; library(mapdata) ; library(mapproj) ; library(httr) ;
library(RColorBrewer) ; library(GISTools) ; library(stringi) ; library(stringr) ;
library(grDevices) ; library(dplyr)
#########################################

### NOTE: The files read in here are created by the Primary_Productivity_pre_calc.R script.  

# read in the SEAWIFS data subset to the 10 gfish study areas
SEAW_chl_gfish <- read.csv("../diversity-data/SEAW_chl_gfish_areas.csv", header=TRUE)


# MAP the clipped SEAWIFS data
# pinks <- c("#FFE6DA", "#E3C9C6", "#FCC5C0", "#FA9FB5", "#F768A1",  # colors for study areas
#            "#E7298A", "#DD3497", "#AE017E", "#7A0177",  "#49006A")   
# 
# SWFS_map <- akmap2 + 
#             geom_point(data=SEAW_chl_gfish, aes(x=longitude_degrees_east, y=latitude_degrees_north)) + 
#             geom_point(data=discrete_areas, 
#                        aes(LONGITUDE, LATITUDE, color=area1,
#                            fill=area1)) + 
#           # scale_color_manual(values=pinks, name="Study Area", guide=FALSE) +
#           # scale_fill_manual(values=pinks, name="Study Area", guide=FALSE) +
#             annotate("text", x=-146, y=59.9, label="1", size=4) + 
#             annotate("text", x=-147.6, y=59.6, label="2", size=4) +
#             annotate("text", x=-148.55, y=59.4, label="3", size=4) +
#             annotate("text", x=-149.5, y=59.4, label="4", size=4) +
#             annotate("text", x=-152.1, y=58.8, label="5", size=4) +
#             annotate("text", x=-150.3, y=58.4, label="6", size=4) +
#             annotate("text", x=-151.1, y=57.8, label="7", size=4, color="white") +
#             annotate("text", x=-152.85, y=56.8, label="8", size=4, color="white") +
#             annotate("text", x=-155.5, y=56.45, label="9", size=4, color="white") +
#             annotate("text", x=-157.25, y=55.75, label="10", size=4, color="white")
#
#SWFS_map


# read in the MODIS data subset to the 10 gfish study areas
MODIS_chl_gfish <- read.csv("../diversity-data/MODIS_chl_gfish_areas.csv", header=TRUE)


# # MAP the clipped MODIS data
# MDS_map <- akmap2 + 
#            geom_point(data=MODIS_chl_gfish, aes(x=longitude_degrees_east, y=latitude_degrees_north)) +
#            geom_point(data=discrete_areas, 
#                        aes(LONGITUDE, LATITUDE, color=area1,
#                            fill=area1)) + 
#           # scale_color_manual(values=pinks, name="Study Area", guide=FALSE) +
#           # scale_fill_manual(values=pinks, name="Study Area", guide=FALSE) +
#             annotate("text", x=-146, y=59.9, label="1", size=4) + 
#             annotate("text", x=-147.6, y=59.6, label="2", size=4) +
#             annotate("text", x=-148.55, y=59.4, label="3", size=4) +
#             annotate("text", x=-149.5, y=59.4, label="4", size=4) +
#             annotate("text", x=-152.1, y=58.8, label="5", size=4) +
#             annotate("text", x=-150.3, y=58.4, label="6", size=4) +
#             annotate("text", x=-151.1, y=57.8, label="7", size=4, color="white") +
#             annotate("text", x=-152.85, y=56.8, label="8", size=4, color="white") +
#             annotate("text", x=-155.5, y=56.45, label="9", size=4, color="white") +
#             annotate("text", x=-157.25, y=55.75, label="10", size=4, color="white")
# 
# MDS_map

############################################################
#  Combine the two datasets  
CHL_both <- SEAW_chl_gfish %>%
            full_join(MODIS_chl_gfish) %>%
            dplyr::select(-time_UTC, -optional) %>%
            filter(Year != 2016)
  
  

  
  
  
  
  





