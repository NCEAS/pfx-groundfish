###############################################################
#####  Script for Primary Production Analysis             #####
#####  Created by Rachael E. Blake, March 2017 to test    #####
#####      the intermediate productivity hypothesis for   #####
#####      groundfish paper #2                            #####
###############################################################

# load libraries needed by this script
library(rgdal) ; library(dplyr) ; library(ggplot2) ; library(reshape2)
library(scales) ; library(psych) ; library(rworldmap) ; library(rworldxtra)
library(sp) ; library(maptools) ; library(raster) ; library(rgeos) ; library(ggplot2) ; 
library(maps) ; library(mapdata) ; library(mapproj) ; library(httr) ; library(scales) ;
library(RColorBrewer) ; library(GISTools) ; library(stringi) ; library(stringr)

# read in some data
SEAWIFS_chl <- read.csv("./diversity-data/SEAWIFS_erdSW1chlamday_b482_a79d_4c9e.csv", header=TRUE)
head(SEAWIFS_chl)

# filter data to within our study area
SEAW_chl1 <- SEAWIFS_chl %>%
             mutate(Year = str_sub(time_UTC, 1,4),
                    Month = str_sub(time_UTC, 6,7),
                    Day = str_sub(time_UTC, 9,10)) %>%
             filter(latitude_degrees_north < 62 & latitude_degrees_north > 54) %>%
             filter(longitude_degrees_east < -142 & longitude_degrees_east > -165) %>%
             filter(!chlorophyll_mg_m.3 == "NaN") #%>%
             #filter(!time_UTC == "1997-10-16T00:00:00Z",
             #        !time_UTC == "1997-09-16T00:00:00Z",
             #        !time_UTC == "1998-05-16T00:00:00Z",
             #        !time_UTC == "1998-04-16T00:00:00Z",
             #        !time_UTC == "1998-03-16T00:00:00Z",
             #        !time_UTC == "1998-08-16T00:00:00Z",
             #        !time_UTC == "1998-06-16T00:00:00Z",
            #         !time_UTC == "1998-07-16T00:00:00Z")  # removing really high chla values...seems spurious

# read in some data
MODIS_chl <- read.csv("./diversity-data/MODIS_erdMH1chlamday_632a_18ce_f111.csv", header=TRUE)
head(MODIS_chl)




# map where the data are...
ak <- map_data('worldHires','USA:Alaska')

akmap2 <- ggplot() + 
          geom_polygon(data=ak, aes(long,lat,group=group), fill=8, 
                       color="black") +
          xlab(expression(paste(Longitude^o,~'W'))) +
          ylab(expression(paste(Latitude^o,~'N'))) +
          coord_map(xlim=c(-165, -142), ylim=c(54, 62)) +
          theme(panel.background=element_rect(fill='aliceblue')) 

SWFS_map <- akmap2 + 
            geom_point(data=SEAW_chl1, aes(x=longitude_degrees_east, y=latitude_degrees_north,
                                           color=chlorophyll_mg_m.3))

SWFS_map

MDS_map <- akmap2 + 
           geom_point(data= )


