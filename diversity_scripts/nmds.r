##########################################################################################
###   title: "Groundfish multivariate analyses - Comparing Areas By Year"              ###
###   author: "Colette Ward"                                                           ###
###   date: "May 6, 2016"                                                              ###
###   NOTE:  This code was taken from Colette's nmds.Rmd file, but put                 ###
###   here so it could be sourced into the final Figures script.                       ###
##########################################################################################

# Load packages
library(vegan) ; library(mvnormtest) ; library(plyr) ; library(dplyr) ; library(tidyr)


# load & prep look-up table of common names
common <- read.csv("./diversity-data/trawl_species_control_file.csv", header = T, stringsAsFactors = FALSE)

common1 <- common %>%
           dplyr::select(database.name, common.name) %>%
           dplyr::rename(Species = database.name)
for (i in 1:nrow(common1)) { # add common names for Sebastes 1 & 2
  if(common1$Species[i] == "Dusky.and.Dark.Rockfish") {common1$common.name[i] <- "Sebastes 1"}
  if(common1$Species[i] == "Rougheye.and.Blackspotted.Rockfish") {common1$common.name[i] <- "Sebastes 2"}
}


# load mean annual CPUE data for Shallow Areas (these are Ole's 11 areas)
shallowCPUEArea <- read.csv("./diversity-data/All_sp_index_meanCPUEByArea.csv", header = T, stringsAsFactors = FALSE)

shallowCPUEArea2 <- left_join(shallowCPUEArea, common1, by = "Species") %>% # merge common names onto SPCPUEArea
                    dplyr::filter(area != "Total", area != "8") %>% # remove regional totals and Ole's area 8
                    mutate(area = revalue(area, c("9"="8", "10"="9", "11"="10"))) # renumber (old = new) shallow areas to account for splitting area 7 into 7, 8, 9 (but removing 8)


#which are the most dominant species?
s <- shallowCPUEArea2 %>%
     filter(area == "Total") %>%
     group_by(common.name) %>%
     summarize(sum(Mean.totalDensity)) %>%
     ungroup() %>%
     arrange(desc(`sum(Mean.totalDensity)`))
print(head(s))
#1 arrowtooth flounder                 900.2215
#2     Pacific halibut                 338.6928
#3         Pacific cod                 268.7385
#4     walleye pollock                 247.9965
#5 Pacific ocean perch                 165.1791
#6   northern rockfish                 109.9819


######################################################
######################################################
######################################################

# nMDS
# Transformation is Wisconsin-style double transformation (normalizes taxa to % abundance, then normalizes 
# abundances to the maximum for each species) of sqrt-transformed data.
# We do not use a Hellinger transformation because it reduces the influence of very rare species and 
# inflates the influence of very abundant species (renders the result driven by Arrowtooth).

# organize Shallow areas data for nMDS:
shallowCPUE_spread3 <- shallowCPUEArea2 %>%
                       dplyr::select(common.name, Mean.totalDensity, area, year) %>%
                       spread(common.name, Mean.totalDensity) %>% # make each species a column
                       mutate(Area = as.numeric(area)) %>% dplyr::select(-area) # convert area to numeric


### run the nMDS:
# Shallow areas:
ordShallowWiscSqrt <- metaMDS(shallowCPUE_spread3[2:54], distance='bray', trymax=1000)
ordShallowHellinger <- metaMDS(decostand(shallowCPUE_spread3[,2:54], "hell"), distance='bray', trymax=1000) # goodness of fit plots look better on data with Hellinger transformation
# also considered method = "log" (does not converge) and "max" (results are similar to wisc(sqrt( )))

par(mfrow=c(1,1))
stressplot(ordShallowWiscSqrt) # make sure stressplots show good fits

plot(ordShallowWiscSqrt, display = 'sites', type = 't', main = "Shallow Areas") # look at goodness-of-fit plots
points(ordShallowWiscSqrt, display = 'sites', cex = goodness (ordShallowWiscSqrt)*200)


### Ordination plots with grouping by areas:
# Create ordination plots for Wisconsin double standardization + Square root transformation:
# Plot shallow areas:
# ordiplot(ordShallowWiscSqrt, type="n", choices=c(1,2), main=NULL, xlab="NMDS 1", ylab="NMDS 2")#, ylim = c(-0.3, 0.2))
# #title(main="Shallow Areas, Wisconsin + Sqrt")
# #text(ordShallowWiscSqrt, display = "sites", cex = 0.5, col="black")
# text(ordShallowWiscSqrt, display = "spec", cex = 0.5, col="black")
# #ordiplot(ordShallowWiscSqrt, display = "sites") #, groups=shallowCPUE_spread3$Area)
# # colors from map and other plots ("#FFE6DA", "#E3C9C6", "#FCC5C0", "#FA9FB5", "#F768A1",
# #           "#E7298A", "#DD3497", "#AE017E", "#7A0177",  "#49006A")   
# ordihull(ordShallowWiscSqrt, groups=shallowCPUE_spread3$Area, show.groups=1, col="#FFE6DA", label=T, lwd=3, lty=1, cex=1.5)
# ordihull(ordShallowWiscSqrt, groups=shallowCPUE_spread3$Area, show.groups=2, col="#E3C9C6", label=T, lwd=3, lty=1, cex=1.5)
# ordihull(ordShallowWiscSqrt, groups=shallowCPUE_spread3$Area, show.groups=3, col="#FCC5C0", label=T, lwd=3, lty=1, cex=1.5)
# ordihull(ordShallowWiscSqrt, groups=shallowCPUE_spread3$Area, show.groups=4, col="#FA9FB5", label=T, lwd=3, lty=1, cex=1.5)
# ordihull(ordShallowWiscSqrt, groups=shallowCPUE_spread3$Area, show.groups=5, col="#F768A1", label=T, lwd=3, lty=1, cex=1.5)
# ordihull(ordShallowWiscSqrt, groups=shallowCPUE_spread3$Area, show.groups=6, col="#E7298A", label=T, lwd=3, lty=1, cex=1.5)
# ordihull(ordShallowWiscSqrt, groups=shallowCPUE_spread3$Area, show.groups=7, col="#DD3497", label=T, lwd=3, lty=1, cex=1.5)
# ordihull(ordShallowWiscSqrt, groups=shallowCPUE_spread3$Area, show.groups=8, col="#AE017E", label=T, lwd=3, lty=1, cex=1.5)
# ordihull(ordShallowWiscSqrt, groups=shallowCPUE_spread3$Area, show.groups=9, col="#7A0177", label=T, lwd=3, lty=1, cex=1.5)
# ordihull(ordShallowWiscSqrt, groups=shallowCPUE_spread3$Area, show.groups=10, col="#49006A", label=T, lwd=3, lty=1, cex=1.5)

# areas 8 & 9  are very different from the rest, and also from each other.
# there is also an east-west gradient in areas 1-7.
# and area 1 does not overlap with any others.
# is area 7 similar to other areas outside the Alaska Coastal Current?
# note that area 8 & 9 hulls look smaller than hulls for areas 1-4 & 7 ... suggests that temporal turnover in 8 & 9 was less than that in 1-4 & 7 ... use betadisper test (must use same transformation) and then anova(test) to test this. Only use betadisper on factors that permanova says are significant.


### Create ordination plots for Hellinger transformation:
# Plot shallow areas:
# ordiplot(ordShallowHellinger, type="n", choices=c(1,2), xlab="NMDS 1", ylab="NMDS 2")#, ylim = c(-0.3, 0.2)), main="Hellinger Transformation",
# #title(main="Shallow Areas, Hellinger")
# #text(ordShallowHellinger, display = "sites", cex = 0.5, col="black")
# text(ordShallowHellinger, display = "spec", cex = 0.5, col="black")
# #ordiplot(ordShallowHellinger, display = "sites") #, groups=shallowCPUE_spread3$Area)
# # colors from map and other plots ("#FFE6DA", "#E3C9C6", "#FCC5C0", "#FA9FB5", "#F768A1",
# #           "#E7298A", "#DD3497", "#AE017E", "#7A0177",  "#49006A")  
# ordihull(ordShallowHellinger, groups=shallowCPUE_spread3$Area, show.groups=1, col="#FFE6DA", label=T, lwd=3, lty=1, cex=1.5)
# ordihull(ordShallowHellinger, groups=shallowCPUE_spread3$Area, show.groups=2, col="#E3C9C6", label=T, lwd=3, lty=1, cex=1.5)
# ordihull(ordShallowHellinger, groups=shallowCPUE_spread3$Area, show.groups=3, col="#FCC5C0", label=T, lwd=3, lty=1, cex=1.5)
# ordihull(ordShallowHellinger, groups=shallowCPUE_spread3$Area, show.groups=4, col="#FA9FB5", label=T, lwd=3, lty=1, cex=1.5)
# ordihull(ordShallowHellinger, groups=shallowCPUE_spread3$Area, show.groups=5, col="#F768A1", label=T, lwd=3, lty=1, cex=1.5)
# ordihull(ordShallowHellinger, groups=shallowCPUE_spread3$Area, show.groups=6, col="#E7298A", label=T, lwd=3, lty=1, cex=1.5)
# ordihull(ordShallowHellinger, groups=shallowCPUE_spread3$Area, show.groups=7, col="#DD3497", label=T, lwd=3, lty=1, cex=1.5)
# ordihull(ordShallowHellinger, groups=shallowCPUE_spread3$Area, show.groups=8, col="#AE017E", label=T, lwd=3, lty=1, cex=1.5)
# ordihull(ordShallowHellinger, groups=shallowCPUE_spread3$Area, show.groups=9, col="#7A0177", label=T, lwd=3, lty=1, cex=1.5)
# ordihull(ordShallowHellinger, groups=shallowCPUE_spread3$Area, show.groups=10, col="#49006A", label=T, lwd=3, lty=1, cex=1.5)



### Ordination plots with grouping by years:
# Wisconsin transformation:
par(mfrow=c(1,2))
ordiplot(ordShallowWiscSqrt, type="n", choices=c(1,2), main=NULL, xlab="NMDS 1", ylab="NMDS 2")#, ylim = c(-0.3, 0.2))
title(main="Shallow Areas")
#text(ordShallowWiscSqrt, display = "sites", cex = 0.5, col="black")
text(ordShallowWiscSqrt, display = "spec", cex = 0.5, col="black")
#ordiplot(ordShallowWiscSqrt, display = "sites") #, groups=shallowCPUE_spread3$Area)
ordihull(ordShallowWiscSqrt, groups=shallowCPUE_spread3$year, show.groups=1984, col="red", label=T, lwd=3, lty=1)
ordihull(ordShallowWiscSqrt, groups=shallowCPUE_spread3$year, show.groups=1987, col="orange", label=T, lwd=3, lty=1)
ordihull(ordShallowWiscSqrt, groups=shallowCPUE_spread3$year, show.groups=1990, col="yellow", label=T, lwd=3, lty=1)
ordihull(ordShallowWiscSqrt, groups=shallowCPUE_spread3$year, show.groups=1993, col="green", label=T, lwd=3, lty=1)
ordihull(ordShallowWiscSqrt, groups=shallowCPUE_spread3$year, show.groups=1996, col="dark green", label=T, lwd=3, lty=1)
ordihull(ordShallowWiscSqrt, groups=shallowCPUE_spread3$year, show.groups=1999, col="turquoise", label=T, lwd=3, lty=1)
ordihull(ordShallowWiscSqrt, groups=shallowCPUE_spread3$year, show.groups=2001, col="blue", label=T, lwd=3, lty=1)
ordihull(ordShallowWiscSqrt, groups=shallowCPUE_spread3$year, show.groups=2003, col="dark blue", label=T, lwd=3, lty=1)
ordihull(ordShallowWiscSqrt, groups=shallowCPUE_spread3$year, show.groups=2005, col="lavender", label=T, lwd=3, lty=1)
ordihull(ordShallowWiscSqrt, groups=shallowCPUE_spread3$year, show.groups=2007, col="purple", label=T, lwd=3, lty=1)
ordihull(ordShallowWiscSqrt, groups=shallowCPUE_spread3$year, show.groups=2009, col="pink", label=T, lwd=3, lty=1)
ordihull(ordShallowWiscSqrt, groups=shallowCPUE_spread3$year, show.groups=2011, col="light grey", label=T, lwd=3, lty=1)
ordihull(ordShallowWiscSqrt, groups=shallowCPUE_spread3$year, show.groups=2013, col="dark grey", label=T, lwd=3, lty=1)
ordihull(ordShallowWiscSqrt, groups=shallowCPUE_spread3$year, show.groups=2015, col="black", label=T, lwd=3, lty=1)


# Hellinger Transformation:
par(mfrow=c(1,2))
ordiplot(ordShallowHellinger, type="n", choices=c(1,2), main=NULL, xlab="NMDS 1", ylab="NMDS 2")#, ylim = c(-0.3, 0.2))
title(main="Shallow Areas")
#text(ordShallowHellinger, display = "sites", cex = 0.5, col="black")
text(ordShallowHellinger, display = "spec", cex = 0.5, col="black")
#ordiplot(ordShallowHellinger, display = "sites") #, groups=shallowCPUE_spread3$Area)
ordihull(ordShallowHellinger, groups=shallowCPUE_spread3$year, show.groups=1984, col="red", label=T, lwd=3, lty=1)
ordihull(ordShallowHellinger, groups=shallowCPUE_spread3$year, show.groups=1987, col="orange", label=T, lwd=3, lty=1)
ordihull(ordShallowHellinger, groups=shallowCPUE_spread3$year, show.groups=1990, col="yellow", label=T, lwd=3, lty=1)
ordihull(ordShallowHellinger, groups=shallowCPUE_spread3$year, show.groups=1993, col="green", label=T, lwd=3, lty=1)
ordihull(ordShallowHellinger, groups=shallowCPUE_spread3$year, show.groups=1996, col="dark green", label=T, lwd=3, lty=1)
ordihull(ordShallowHellinger, groups=shallowCPUE_spread3$year, show.groups=1999, col="turquoise", label=T, lwd=3, lty=1)
ordihull(ordShallowHellinger, groups=shallowCPUE_spread3$year, show.groups=2001, col="blue", label=T, lwd=3, lty=1)
ordihull(ordShallowHellinger, groups=shallowCPUE_spread3$year, show.groups=2003, col="dark blue", label=T, lwd=3, lty=1)
ordihull(ordShallowHellinger, groups=shallowCPUE_spread3$year, show.groups=2005, col="lavender", label=T, lwd=3, lty=1)
ordihull(ordShallowHellinger, groups=shallowCPUE_spread3$year, show.groups=2007, col="purple", label=T, lwd=3, lty=1)
ordihull(ordShallowHellinger, groups=shallowCPUE_spread3$year, show.groups=2009, col="pink", label=T, lwd=3, lty=1)
ordihull(ordShallowHellinger, groups=shallowCPUE_spread3$year, show.groups=2011, col="light grey", label=T, lwd=3, lty=1)
ordihull(ordShallowHellinger, groups=shallowCPUE_spread3$year, show.groups=2013, col="dark grey", label=T, lwd=3, lty=1)
ordihull(ordShallowHellinger, groups=shallowCPUE_spread3$year, show.groups=2015, col="black", label=T, lwd=3, lty=1)


###############
# PERMANOVA

# Using Wisconsing transformation:
adonis(wisconsin(sqrt(shallowCPUE_spread3[,2:54])) ~ Area*year, data = shallowCPUE_spread3, permutations = 200, method = "bray")
#           Df SumsOfSqs MeanSqs F.Model      R2   Pr(>F)   
#Area        1    1.2116 1.21156  45.282 0.21624 0.004975 **
#year        1    0.6860 0.68597  25.638 0.12243 0.004975 **
#Area:year   1    0.0665 0.06655   2.487 0.01188 0.054726 . 
#Residuals 136    3.6388 0.02676         0.64945            
#Total     139    5.6029                 1.00000 

# Using Hellinger transformation:
adonis(decostand(shallowCPUE_spread3[2:54],"hell") ~ Area*year, data = shallowCPUE_spread3, permutations = 200, method = "bray")
# Note this significance is likely driven by differences between 8 & 9 vs all other sites, not sites 1-6
#           Df SumsOfSqs MeanSqs F.Model      R2   Pr(>F)   
#Area        1    0.8614 0.86138  39.894 0.20597 0.004975 **
#year        1    0.3261 0.32614  15.105 0.07799 0.004975 **
#Area:year   1    0.0579 0.05794   2.683 0.01385 0.039801 * 
#Residuals 136    2.9365 0.02159         0.70218            
#Total     139    4.1820                 1.00000   



###################
# Homogeneity of multivariate dispersion:

shallowDist <- vegdist(wisconsin(sqrt(shallowCPUE_spread3[,2:54])),"bray") # consider Hellinger transformation to deal with superabundant vs super rare species (see Legendre & Legendre)
#shallowDist = vegdist(wisconsin(sqrt(shallowCPUE_spread3[,2:54])),"jaccard") 


shallowTest_areas <- betadisper(shallowDist,shallowCPUE_spread3$Area) # this gives within-area variability (ie looks at temporal effects ... how different is each year from the mean across all years)
boxplot(shallowTest_areas) # indicates temporal turnover (beta diversity in time)
anova(shallowTest_areas) # p = 0.94   Shows that area of each hull in nmds is not significantly different ... ie no difference in temporal turnover between areas
#            Df   Sum Sq    Mean Sq F value Pr(>F)
# Groups      9 0.005072 0.00056356  0.3949 0.9357
# Residuals 130 0.185509 0.00142699 
# TukeyHSD(shallowTest_areas) # irrelevant


shallowTest_years <- betadisper(shallowDist,shallowCPUE_spread3$year) # within-year variability ... how different is each area from the mean across all areas?
boxplot(shallowTest_years)
anova(shallowTest_years) # p = 1 ... how does this add up with our signifant betwen-region result???
#            Df  Sum Sq   Mean Sq F value Pr(>F)
# Groups     13 0.00232 0.0001784  0.0378      1
# Residuals 126 0.59462 0.0047192
# TukeyHSD(shallowTest_years) # irrelevant











