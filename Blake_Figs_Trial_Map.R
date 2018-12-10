
library(dplyr); library(tidyr) ; library(stringr) ; library(rgdal) ; 
library(sp) ; library(maptools) ; library(raster) ; library(ggplot2) ;  library(rgeos) ; 
library(maps) ; library(mapdata) ; library(mapproj) ; library(httr) ; library(extrafont) ; 
library(vegan) ; library(scales) ; library(gridExtra) ; library(RColorBrewer) ; 
library(GISTools) ; library(egg) 

library(marmap)


ak <- fortify(maps::map('worldHires','USA:Alaska', plot=FALSE, fill=TRUE))


ak_bathy <- getNOAA.bathy(-165, -142, 54, 62, res=1, keep=TRUE)


akmap4 <- autoplot(ak_bathy, geom=c("raster", "contour"), coast=FALSE) +
          scale_fill_gradient2(low="deepskyblue4", mid="lightskyblue4", high="white", guide=FALSE) +
          geom_polygon(data=ak, aes(long,lat,group=group), fill=8, color="black") +
          xlab(expression(paste(Longitude^o,~'W'))) +
          ylab(expression(paste(Latitude^o,~'N'))) +
          coord_cartesian(xlim=c(-165, -142), ylim=c(54, 62))
  
  
# Adding Study Areas to the map
discrete_areas1 <- read.csv("./goaTrawl/Output Data/goa_discrete_areas_for_comparison(50_to_150m).csv")

discrete_areas1$area <- discrete_areas1$Area
discrete_areas1$Area <- factor(discrete_areas1$area)

# Re-naming discrete areas for Diversity paper
discrete_areas <- discrete_areas1 %>%
                  mutate(area = ifelse(area == "8", "NULL", 
                                ifelse(area == "9", "8", 
                                ifelse(area == "10", "9", 
                                ifelse(area == "11", "10", area))))) %>%
                  filter(area != "NULL") %>%
                  mutate(area = factor(area))


discrete_areas$area1 <- factor(discrete_areas$area,
                               levels=c('1','2','3','4','5','6','7','8','9','10'))
# http://colorbrewer2.org/     brewer.pal(n=9, "PuRd")
pinks <- c("#FFE6DA", "#E3C9C6", "#FCC5C0", "#FA9FB5", "#F768A1",
           "#E7298A", "#DD3497", "#AE017E", "#7A0177",  "#49006A")     


ak_areas_all <- akmap4 +
                geom_point(data=discrete_areas, 
                           aes(LONGITUDE, LATITUDE, color=area1#,fill=area1
                               )) + 
                 scale_color_manual(values=pinks, name="Study Area", guide=FALSE) +
                # scale_fill_manual(values=pinks, name="Study Area", guide=FALSE) +
                 annotate("text", x=-146, y=59.9, label="1", size=4, fontface="bold") + 
                 annotate("text", x=-147.6, y=59.6, label="2", size=4, fontface="bold") +
                 annotate("text", x=-148.55, y=59.4, label="3", size=4, fontface="bold") +
                 annotate("text", x=-149.5, y=59.4, label="4", size=4, fontface="bold") +
                 annotate("text", x=-152.1, y=58.8, label="5", size=4, fontface="bold") +
                 annotate("text", x=-150.3, y=58.4, label="6", size=4, fontface="bold") +
                 annotate("text", x=-151.1, y=57.8, label="7", size=4, color="white", fontface="bold") +
                 annotate("text", x=-152.85, y=56.8, label="8", size=4, color="white", fontface="bold") +
                 annotate("text", x=-155.5, y=56.45, label="9", size=4, color="white", fontface="bold") +
                 annotate("text", x=-157.25, y=55.75, label="10", size=4, color="white", fontface="bold") +
                 annotate("segment", x=-156.25, y=54, xend=-156.25, yend=57.5, color="black", 
                          linetype="dashed", size=1) + 
                 annotate("segment", x=-149, y=58.8, xend=-149, yend=59.9, color="black", 
                          linetype="dashed", size=1) + 
                 annotate("segment", x=-149, y=58.95, xend=-142, yend=59.2, color="black", 
                          linetype="dashed", size=1) + 
                 annotate("segment", x=-149, y=58.8, xend=-151.25, yend=58.75, color="black", 
                          linetype="dashed", size=1) + 
                 annotate("segment", x=-151.25, y=58.75, xend=-151.25, yend=58.45, color="black", 
                          linetype="dashed", size=1) + 
                 annotate("segment", x=-151.25, y=58.45, xend=-153.35, yend=58.45, color="black", 
                          linetype="dashed", size=1) + 
                 annotate("segment", x=-153.35, y=58.45, xend=-153.35, yend=58.75, color="black", 
                          linetype="dashed", size=1) + 
                 annotate("text", x=-164, y=61.5, label="N", size=8, fontface="bold", color="black") + 
                 geom_segment(aes(x=-164,xend=-164,y=61.2,yend=61.8), size=1.5,
                              arrow=arrow(length = unit(0.5, "cm"))) +
                 annotate("text", x=-158, y=54.5, label="AP", fontface="bold", size=6) + 
                 annotate("text", x=-151, y=56.5, label="K", fontface="bold", size=6) + 
                 annotate("text", x=-152.5, y=59.5, label="CI", fontface="bold", size=6) + 
                 annotate("text", x=-144, y=59.5, label="PWS", fontface="bold", size=6) + 
                 theme(axis.text=element_text(size=15),
                       axis.title=element_text(size=15))
  