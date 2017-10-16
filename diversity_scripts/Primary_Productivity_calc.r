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

# set working dir
#setwd("~/NCEAS/GoA Portfolio Effects WG/pfx-groundfish")


# # Study Areas
# discrete_areas1 <- read.csv("./goaTrawl/Output Data/goa_discrete_areas_for_comparison(50_to_150m).csv")
# discrete_areas1$area <- discrete_areas1$Area
# discrete_areas1$Area <- factor(discrete_areas1$area)
# 
# # re-naming discrete areas for Diversity paper
# discrete_areas <- discrete_areas1 %>%
#                   mutate(area = ifelse(area == "8", "NULL", 
#                                 ifelse(area == "9", "8", 
#                                 ifelse(area == "10", "9", 
#                                 ifelse(area == "11", "10", area))))) %>%
#                   filter(area != "NULL") %>%
#                   mutate(area = factor(area))
# 
# discrete_areas$area1 <- factor(discrete_areas$area,
#                                levels=c('1','2','3','4','5','6','7','8','9','10'))
# 
# # map where the data are...
# # library(tidyverse) also has a function called map() which causes this to produce an error
# #ak <- map_data('worldHires','USA:Alaska')  
# ak <- fortify(maps::map('worldHires','USA:Alaska', plot=FALSE, fill=TRUE))
#   
# akmap2 <- ggplot() + 
#           geom_polygon(data=ak, aes(long,lat,group=group), fill=8, 
#                        color="black") +
#           xlab(expression(paste(Longitude^o,~'W'))) +
#           ylab(expression(paste(Latitude^o,~'N'))) +
#           coord_map(xlim=c(-165, -142), ylim=c(54, 62)) +
#           theme(panel.background=element_rect(fill='aliceblue')) 

##########################################
# Function for making convex hulls around study areas
#' Title
#'
#' @param area_number 
#' @param df 
#'
#' @return
#' @export
#'
#' @examples
make_cnvx_hull <- function(area_number, df){
  
                  # make the convex hull polygon
                  area_hull <- df %>% filter(area1 == area_number) %>% 
                               dplyr::select(LONGITUDE, LATITUDE)
                  hull_pos <- chull(area_hull)
                  hull_coords <- area_hull[c(hull_pos, hull_pos[1]), ] 
 
                  # check that it looks ok
                  #plot(area_hull, pch=19)
                  #lines(hull_coords, col="red")

                  # make the convex hull a SpatialPolygonsDataFrame, which can be saved as a shapefile 
                  # with rgdal::writeOGR()
                  sp_poly <- SpatialPolygons(list(Polygons(list(Polygon(hull_coords)), ID=1)), 
                                             proj4string=CRS("+proj=longlat +datum=WGS84"))
                  sp_poly_df <- SpatialPolygonsDataFrame(sp_poly, data=data.frame(ID=1))
                  rgdal::writeOGR(sp_poly_df,  "./diversity-data", layer=area_number, driver="ESRI Shapefile")
}

#########################################
#make a shapefile for each study area using convex hulls

# # use sapply to run this function over all values of discrete_areas$area1
# sapply(unique(discrete_areas$area1), function(x) make_cnvx_hull(area_number=x, df=discrete_areas))
# 
# # make individual shapefiles into one shapefile
# area_list <- list()
# area_list[[1]] <- readOGR(dsn="./diversity-data", layer="1")
# area_list[[2]] <- readOGR(dsn="./diversity-data", layer="2")
# area_list[[3]] <- readOGR(dsn="./diversity-data", layer="3")
# area_list[[4]] <- readOGR(dsn="./diversity-data", layer="4")
# area_list[[5]] <- readOGR(dsn="./diversity-data", layer="5")
# area_list[[6]] <- readOGR(dsn="./diversity-data", layer="6")
# area_list[[7]] <- readOGR(dsn="./diversity-data", layer="7")
# area_list[[8]] <- readOGR(dsn="./diversity-data", layer="8")
# area_list[[9]] <- readOGR(dsn="./diversity-data", layer="9")
# area_list[[10]] <- readOGR(dsn="./diversity-data", layer="10")
# 
# areas_shps1 <- raster::union(area_list[[1]],area_list[[2]])
# areas_shps2 <- raster::union(area_list[[3]],area_list[[4]])
# areas_shps3 <- raster::union(area_list[[5]],area_list[[6]])
# areas_shps4 <- raster::union(area_list[[7]],area_list[[8]])
# areas_shps5 <- raster::union(area_list[[9]],area_list[[10]])
# 
# areas_shps6 <- raster::union(areas_shps1,areas_shps2)
# areas_shps7 <- raster::union(areas_shps3,areas_shps4)
# areas_shps8 <- raster::union(areas_shps5,areas_shps6)
# areas_shps10 <- raster::union(areas_shps7,areas_shps8)
# 
# areas_shps10 # This is the shapefile of all 10 study area polygons
# 
# # makes a column that assigns area number to each polygon in the combined shapefile
# areas_shps10@data$area <- c(5,6,7,8,9,10,1,2,3,4)
# 
# # write the combined shapefile to the repo
# #rgdal::writeOGR(areas_shps10,  "./diversity-data", layer="ten", driver="ESRI Shapefile")
# 
# 
# #########################################
# # read in the SEAWIFS data
# SEAWIFS_chl <- read.csv("./diversity-data/SEAWIFS_erdSW1chlamday_b482_a79d_4c9e.csv", header=TRUE)
# head(SEAWIFS_chl)
# 
# # filter data to within our study area
# SEAW_chl1 <- SEAWIFS_chl %>%
#              mutate(Year = str_sub(time_UTC, 1,4),   # This data appears to only be 1997 and 1998
#                     Month = str_sub(time_UTC, 6,7),
#                     Day = str_sub(time_UTC, 9,10)) %>%
#              filter(latitude_degrees_north < 62 & latitude_degrees_north > 54) %>%
#              filter(longitude_degrees_east < -142 & longitude_degrees_east > -165) %>%
#              filter(!chlorophyll_mg_m3 == "NaN") %>%
#              mutate_if(is.character, as.integer) #%>%
#              #filter(!time_UTC == "1997-10-16T00:00:00Z",
#              #        !time_UTC == "1997-09-16T00:00:00Z",
#              #        !time_UTC == "1998-05-16T00:00:00Z",
#              #        !time_UTC == "1998-04-16T00:00:00Z",
#              #        !time_UTC == "1998-03-16T00:00:00Z",
#              #        !time_UTC == "1998-08-16T00:00:00Z",
#              #        !time_UTC == "1998-06-16T00:00:00Z",
#              #        !time_UTC == "1998-07-16T00:00:00Z")  # removing really high chla values...seems spurious
# 
# # Clip the SEAWIFS data to just what overlaps with our study areas
# SEAW_chl2 <- SEAW_chl1
# 
# # copying from this website: https://www.nceas.ucsb.edu/scicomp/usecases/point-in-polygon
# 
# coordinates(SEAW_chl2) <- c("longitude_degrees_east", "latitude_degrees_north") # makes area grid points into spatial polygon file I think
# 
# areas_shps10 # This is the shapefile of all 10 study area polygons
# 
# proj4string(SEAW_chl2) <- proj4string(areas_shps10)  # tell R that coordinates are in the same lat/lon reference system
# 
# # combine is.na() with over() to do the containment test; note that we
# # need to "demote" areas to a SpatialPolygons object first
# S_areas_chla <- !is.na(over(SEAW_chl2, as(areas_shps10, "SpatialPolygons")))
# 
# # use 'over' again, this time with areas_shps10 as a SpatialPolygonsDataFrame
# # object, to determine which areas (if any) contains each chla point, and
# # store the areas name as an attribute of the SEAW_chl2 data
# SEAW_chl2$area <- over(SEAW_chl2, areas_shps10)$area
# 
# # subset file to just those with area != NA
# SEAW_chl2_areas <- filter(data.frame(SEAW_chl2), !is.na(area))
# 
# write.csv(SEAW_chl2_areas, "./diversity-data/SEAW_chl_gfish_areas.csv", row.names=F)

# read in the data subset to the 10 gfish study areas
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

########################################################
# read in the MODIS data  

## NOTE: This file is 4 GB !!!
# MODIS_chl <- read.csv("./diversity-data/MODIS_erdMH1chlamday_632a_18ce_f111.csv", 
#                       skip=1, header=TRUE)
# head(MODIS_chl)
# 
# MODIS_chl1 <- MODIS_chl %>%
#               dplyr::rename(time_UTC = UTC, 
#                      latitude_degrees_north = degrees_north,
#                      longitude_degrees_east = degrees_east,
#                      chlorophyll_mg_m3 = mg.m..3) %>%
#               mutate(Year = str_sub(time_UTC, 1,4),
#                      Month = str_sub(time_UTC, 6,7),
#                      Day = str_sub(time_UTC, 9,10)) %>%
#               filter(latitude_degrees_north < 62 & latitude_degrees_north > 54) %>%
#               filter(longitude_degrees_east < -142 & longitude_degrees_east > -165) %>%
#               filter(!chlorophyll_mg_m3 == "NaN")
# 
# write.csv(MODIS_chl1, "./diversity-data/MODIS_Big_Box_Clipped.csv", row.names=F)
# 
# MODIS_chl_box1 <- read.csv("./diversity-data/MODIS_Big_Box_Clipped.csv", header=TRUE)
# head(MODIS_chl_box1)
# 
# MODIS_chl_box <- MODIS_chl_box1
# 
# 
# # clip the MODIS data to each study area
# # copying from this website: https://www.nceas.ucsb.edu/scicomp/usecases/point-in-polygon
# 
# coordinates(MODIS_chl_box) <- c("longitude_degrees_east", "latitude_degrees_north") # makes area grid points into spatial polygon file I think
# 
# areas_shps10 # This is the shapefile of all 10 study area polygons
# 
# proj4string(MODIS_chl_box) <- proj4string(areas_shps10)  # tell R that coordinates are in the same lat/lon reference system
# 
# # combine is.na() with over() to do the containment test; note that we
# # need to "demote" parks to a SpatialPolygons object first
# areas_chla <- !is.na(over(MODIS_chl_box, as(areas_shps10, "SpatialPolygons")))
# 
# # use 'over' again, this time with areas_shps10 as a SpatialPolygonsDataFrame
# # object, to determine which areas (if any) contains each chla point, and
# # store the areas name as an attribute of the MODIS_chl_box data
# MODIS_chl_box$area <- over(MODIS_chl_box, areas_shps10)$area
# 
# # subset file to just those with area != NA
# MODIS_chl_areas <- filter(data.frame(MODIS_chl_box), !is.na(area))
# 
# write.csv(MODIS_chl_areas, "./diversity-data/MODIS_chl_gfish_areas.csv", row.names=F)

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
  
  

  
  
  
  
  





