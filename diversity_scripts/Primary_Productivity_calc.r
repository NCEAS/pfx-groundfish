###############################################################
#####  Script for Primary Production Analysis             #####
#####  Created by Rachael E. Blake, March 2017 to test    #####
#####      the intermediate productivity hypothesis for   #####
#####      groundfish paper #2                            #####
###############################################################

# load libraries needed by this script
library(rgdal) ; library(dplyr) ; library(ggplot2) ; library(reshape2)
library(scales) ; library(psych) ; library(rworldmap) ; library(rworldxtra) ;
library(sp) ; library(maptools) ; library(raster) ; library(rgeos) ;
library(maps) ; library(mapdata) ; library(mapproj) ; library(httr) ;
library(RColorBrewer) ; library(GISTools) ; library(tidyverse) ; library(stringr)

# read in the SEAWIFS data
SEAWIFS_chl <- read.csv("./diversity-data/SEAWIFS_erdSW1chlamday_b482_a79d_4c9e.csv", header=TRUE)
head(SEAWIFS_chl)

# filter data to within our study area
SEAW_chl1 <- SEAWIFS_chl %>%
             mutate(Year = str_sub(time_UTC, 1,4),   # This data appears to only be 1997 and 1998
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

# Study Areas
discrete_areas1 <- read.csv("./goaTrawl/Output Data/goa_discrete_areas_for_comparison(50_to_150m).csv")
discrete_areas1$area <- discrete_areas1$Area
discrete_areas1$Area <- factor(discrete_areas1$area)

# re-naming discrete areas for Diversity paper
discrete_areas <- discrete_areas1 %>%
                  mutate(area = ifelse(area == "8", "NULL", 
                                ifelse(area == "9", "8", 
                                ifelse(area == "10", "9", 
                                ifelse(area == "11", "10", area))))) %>%
                  filter(area != "NULL") %>%
                  mutate(area = factor(area))

discrete_areas$area1 <- factor(discrete_areas$area,
                               levels=c('1','2','3','4','5','6','7','8','9','10'))

pinks <- c("#FFE6DA", "#E3C9C6", "#FCC5C0", "#FA9FB5", "#F768A1",
           "#E7298A", "#DD3497", "#AE017E", "#7A0177",  "#49006A")   




# read in the MODIS data
MODIS_chl <- read.csv("./diversity-data/MODIS_erdMH1chlamday_632a_18ce_f111.csv", header=TRUE)
head(MODIS_chl)

MODIS_chl1 <- MODIS_chl %>%
              dplyr::rename(time_UTC = UTC, 
                     latitude_degrees_north = degrees_north,
                     longitude_degrees_east = degrees_east,
                     chlorophyll_mg_m3 = mg.m..3) %>%
              mutate(Year = str_sub(time_UTC, 1,4),
                     Month = str_sub(time_UTC, 6,7),
                     Day = str_sub(time_UTC, 9,10)) %>%
              filter(latitude_degrees_north < 62 & latitude_degrees_north > 54) %>%
              filter(longitude_degrees_east < -142 & longitude_degrees_east > -165) %>%
              filter(!chlorophyll_mg_m3 == "NaN")




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
                                           color="blue")) + 
            geom_point(data=discrete_areas, 
                       aes(LONGITUDE, LATITUDE, color=area1,
                           fill=area1)) + 
          # scale_color_manual(values=pinks, name="Study Area", guide=FALSE) +
          # scale_fill_manual(values=pinks, name="Study Area", guide=FALSE) +
            annotate("text", x=-146, y=59.9, label="1", size=4) + 
            annotate("text", x=-147.6, y=59.6, label="2", size=4) +
            annotate("text", x=-148.55, y=59.4, label="3", size=4) +
            annotate("text", x=-149.5, y=59.4, label="4", size=4) +
            annotate("text", x=-152.1, y=58.8, label="5", size=4) +
            annotate("text", x=-150.3, y=58.4, label="6", size=4) +
            annotate("text", x=-151.1, y=57.8, label="7", size=4, color="white") +
            annotate("text", x=-152.85, y=56.8, label="8", size=4, color="white") +
            annotate("text", x=-155.5, y=56.45, label="9", size=4, color="white") +
            annotate("text", x=-157.25, y=55.75, label="10", size=4, color="white")

SWFS_map


MDS_map <- akmap2 + 
           geom_point(data=MODIS_chl1, aes(x=longitude_degrees_east, y=latitude_degrees_north,
                                           color=chlorophyll_mg_m3))


