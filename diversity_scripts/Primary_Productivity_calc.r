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
 
 # map where the data are...
 # library(tidyverse) also has a function called map() which causes this to produce an error
 #ak <- map_data('worldHires','USA:Alaska')  
 ak <- fortify(maps::map('worldHires','USA:Alaska', plot=FALSE, fill=TRUE))
   
 akmap2 <- ggplot() + 
           geom_polygon(data=ak, aes(long,lat,group=group), fill=8, 
                        color="black") +
           xlab(expression(paste(Longitude^o,~'W'))) +
           ylab(expression(paste(Latitude^o,~'N'))) +
           coord_map(xlim=c(-165, -142), ylim=c(54, 62)) +
           theme(panel.background=element_rect(fill='aliceblue')) 

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


############################################################
############################################################

# read in the data subset to the 10 gfish study areas
SEAW_chl_gfish <- read.csv("../diversity-data/SEAW_chl_gfish_areas.csv", header=TRUE)



# read in the MODIS data subset to the 10 gfish study areas
MODIS_chl_gfish <- read.csv("../diversity-data/MODIS_chl_gfish_areas.csv", header=TRUE)



############################################################
#  Combine the two datasets  
CHL_both <- SEAW_chl_gfish %>%
            full_join(MODIS_chl_gfish) %>%
            dplyr::select(-time_UTC, -optional) %>%
            filter(Year != 2016)
  
  

  
  
  
  
  





