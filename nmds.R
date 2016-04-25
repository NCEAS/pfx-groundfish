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
  ord1[[i]] <- metaMDS(byArea_nmds_list2[[1]], distance='bray', k=2, trymax=1000, autotransform=F)
}


plots <- list()
for (j in seq_along(ord1)) {
  plots[[j]] <- plot(ord1[[j]], type = "n")
  points(ord1[[j]], display = "spec", cex = 0.8, pch=21, col="red", bg="yellow")
  text(ord1[[j]], display = "sites", cex=0.7, col="blue")
}

# show the plots:
plots[[1]]; plots[[2]]; plots[[3]]; plots[[4]]; plots[[5]]; plots[[6]]; plots[[7]]; plots[[8]]; plots[[9]]; plots[[10]]; plots[[11]]; plots[[12]]




