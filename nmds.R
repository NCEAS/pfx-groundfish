# NMDS script
# Colette Ward, Apr 22 2016

# load packages
library(httr)
library(vegan)
library(plyr)
library(dplyr)
library(tidyr)



# load & prepare Mean annual CPUE data for Shallow Areas:
URL_SPCPUEArea <- "https://drive.google.com/uc?export=download&id=0By1iaulIAI-uYzBOUFRtZklmX0U" # new data for shallow areas
SPCPUEArea_Get <- GET(URL_SPCPUEArea)
SPCPUEArea_1 <- content(SPCPUEArea_Get, as='text')
SPCPUEArea <- read.csv(file=textConnection(SPCPUEArea_1),stringsAsFactors=FALSE,head=TRUE)
#View(SPCPUEArea)


cpue_spread <- SPCPUEArea %>%
  select(Species, Mean.totalDensity, area, year) %>%
  spread(Species, Mean.totalDensity) %>% # make each species a column
  mutate(area = revalue(area, c("Total" = "10")), # recode Total for looping later
         area = as.numeric(area)) # convert to numeric class

shallowByArea_nmds_list <- split(cpue_spread, f = cpue_spread$area) # create a list of dataframes (one for each area; NB area 10 is Total)

shallowByArea_nmds_list1 <- lapply(shallowByArea_nmds_list, function(y){ row.names(y) <- y$year; y}) # create row names from year column
shallowByArea_nmds_list2 <- lapply(shallowByArea_nmds_list1, function(x) x[!(names(x) %in% c("area", "year"))]) # drop area & year



#########################


# load & prepare mean annual CPUE for Deep areas:
URL_deepCPUE <- "https://drive.google.com/uc?export=download&id=0By1iaulIAI-uVF9VWnNPX3Z3S3c"
deepCPUE_Get <- GET(URL_deepCPUE)
deepCPUE_1 <- content(deepCPUE_Get, as='text')
deepCPUE <- read.csv(file=textConnection(deepCPUE_1),stringsAsFactors=FALSE,head=TRUE)
#View(deepCPUE)

deepCPUE_spread <- deepCPUE %>%
  select(Species, Mean.totalDensity, area, year) %>%
  spread(Species, Mean.totalDensity) %>% # make each species a column
  mutate(area = revalue(area, c("Total" = "6")), # recode Total for looping later
         area = as.numeric(area)) # convert to numeric class

deepByArea_nmds_list <- split(deepCPUE_spread, f = deepCPUE_spread$area) # create a list of dataframes (one for each area; NB area 6 is Total)

deepByArea_nmds_list1 <- lapply(deepByArea_nmds_list, function(y){ row.names(y) <- y$year; y}) # create row names from year column
deepByArea_nmds_list2 <- lapply(deepByArea_nmds_list1, function(x) x[!(names(x) %in% c("area", "year"))]) # drop area & year



######################################################
######################################################
######################################################


# NMDS
# Methods info:
# transformation is Wisconsin-style double transformation (normalizes taxa to % abundance, 
# then normalizes abundances to the maximum for each species) of sqrt-transformed data

#########################

# 1. nMDS for Shallow areas:

ordShallow <- list()
for (i in seq_along(shallowByArea_nmds_list2)) {
  ordShallow[[i]] <- metaMDS(shallowByArea_nmds_list2[[i]], distance='bray', trymax=1000)
}


# examine Shepard's stressplots:
# goodness-of-fit test: plots ordination distances against community dissimilarities
par (mfrow = c(3,4), pty="m") 
stressplotsShallow <- list()
for (i in seq_along(ordShallow)) {
  stressplotsShallow[[i]] <- stressplot(ordShallow[[i]], main = paste("Area",i))
}


# examine goodness of fit plots:
# size represents goodness of fit (bigger = worse fit)
par (mfrow = c(3,4), pty="m") 
goodnessShallow <- list()
for (i in seq_along(ordShallow)) {
  goodnessShallow[[i]] <- plot (ordShallow[[i]], display = 'sites', type = 't', main = paste("Area",i))
  points (ordShallow[[i]], display = 'sites', cex = goodness (ordShallow[[i]])*200)
}


# Ordination plots:
# NB Area 10 is Total of our 9 areas
# >50 species requires small text for legibility; save as .pdf and enlarge to read
par(mfrow = c(2,2))
ordPlotsShallow1 <- list()
for (i in seq_along(ordShallow)) {
  ordPlotsShallow1[[i]] <- plot(ordShallow[[i]], type = "n", main = paste("Shallow Area",i))
  text(ordShallow[[i]], display = "spec", cex = 0.15, col="red")
  text(ordShallow[[i]], display = "sites", cex=0.7, col="blue")
}


#########################

# 2. nMDS for Deep areas:

ordDeep <- list()
for (i in seq_along(deepByArea_nmds_list2)) {
  ordDeep[[i]] <- metaMDS(deepByArea_nmds_list2[[i]], distance='bray', trymax=1000)
}


# examine Shepard's stressplots:
# goodness-of-fit test: plots ordination distances against community dissimilarities
par (mfrow = c(2,3), pty="m") 
stressplotsDeep <- list()
for (i in seq_along(ordDeep)) {
  stressplotsDeep[[i]] <- stressplot(ordDeep[[i]], main = paste("Area",i))
}


# examine goodness of fit plots:
# size represents goodness of fit (bigger = worse fit)
par (mfrow = c(2,3), pty="m") 
goodnessDeep <- list()
for (i in seq_along(ordDeep)) {
  goodnessDeep[[i]] <- plot (ordDeep[[i]], display = 'sites', type = 't', main = paste("Area",i))
  points (ordDeep[[i]], display = 'sites', cex = goodness (ordDeep[[i]])*200)
}


# Ordination plots:
# NB Area 6 is total of our 5 areas
# >50 species requires small text for legibility; save as .pdf and enlarge to read
par(mfrow = c(2,2))
ordPlotsDeep1 <- list()
for (i in seq_along(ordShallow)) {
  ordPlotsDeep1[[i]] <- plot(ordDeep[[i]], type = "n", main = paste("Deep Area",i))
  text(ordDeep[[i]], display = "spec", cex = 0.15, col="red")
  text(ordDeep[[i]], display = "sites", cex=0.7, col="blue")
}
