# NMDS script
# Colette Ward, Apr 22 2016

# load packages
library(httr)
library(vegan)
library(plyr)
library(dplyr)
library(tidyr)



# *** Need to add new data for shallow areas (combining areas 7-9) and data for our new deep areas ***

# load CPUE data from our google drive (this is the old data for Shallow Areas 1-11)
URL_SPCPUEArea <- "https://drive.google.com/uc?export=download&id=0By1iaulIAI-udm1FT2trQUh5N1k"
SPCPUEArea_Get <- GET(URL_SPCPUEArea)
SPCPUEArea_1 <- content(SPCPUEArea_Get, as='text')
SPCPUEArea <- read.csv(file=textConnection(SPCPUEArea_1),stringsAsFactors=FALSE,head=TRUE)


######################################################


# prepare CPUE data for NMDS analysis:
cpue_spread <- SPCPUEArea %>%
  select(Species, Mean.totalDensity, area, year) %>%
  spread(Species, Mean.totalDensity) %>% # make each species a column
  mutate(area = revalue(area, c("Total" = "12")), # recode Total for looping later
         area = as.numeric(area)) # convert to numeric class


byArea_nmds_list <- split(cpue_spread, f = cpue_spread$area) # create a list of dataframes (one for each area; NB area 12 is Total)


byArea_nmds_list1 <- lapply(byArea_nmds_list, function(y){ row.names(y) <- y$year; y}) # create row names from year column
byArea_nmds_list2 <- lapply(byArea_nmds_list1, function(x) x[!(names(x) %in% c("area", "year"))]) # drop area & year


######################################################


# NMDS:
ord1 <- list()
for (i in seq_along(byArea_nmds_list2)) {
  ord1[[i]] <- metaMDS(byArea_nmds_list2[[i]], distance='bray', k=2, trymax=1000)
}


# examine Shepard's stressplots:
# shows the relationship between real distances between samples in resulting m dimensional ordination solution, 
# and their particular compositional dissimilarities expressed by the selected dissimilarity measure
par (mfrow = c(3,4), pty="m") 
stressplots <- list()
for (i in seq_along(ord1)) {
  stressplots[[i]] <- stressplot(ord1[[i]])
}



# examine goodness of fit plots:
# size represents goodness of fit (bigger = worse fit)
par (mfrow = c(3,4), pty="m") 
goodness <- list()
for (i in seq_along(ord1)) {
  goodness[[i]] <- plot (ord1[[i]], display = 'sites', type = 't', main = 'Goodness of fit')
  points (ord1[[i]], display = 'sites', cex = goodness (ord1[[i]])*200)
}



# ordination plots:
par (mfrow = c(3,4), pty="m") 
plots <- list()
for (i in seq_along(ord1)) {
  plots[[i]] <- plot(ord1[[i]], type = "n")
  points(ord1[[i]], display = "spec", cex = 0.8, pch=21, col="red", bg="yellow")
  text(ord1[[i]], display = "sites", cex=0.7, col="blue")
}