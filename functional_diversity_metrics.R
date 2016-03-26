# Functional Diversity Metrics Processing Script
# Colette Ward, Feb 29 2016

# load packages
library(httr)
library(plyr)
library(dplyr)
library(tidyr)
library(psych)
library(FD)

# load the functional trait data
URL_traits <- "https://drive.google.com/uc?export=download&id=0B1XbkXxdfD7uM2M1UnhtTzlGZGM"
traitsGet <- GET(URL_traits)
traits1 <- content(traitsGet, as='text')
traits_df <- read.csv(file=textConnection(traits1),stringsAsFactors=FALSE)

# or load it from local source:
traits_df <- read.csv("Groundfish-Functional-Diversity-Traits.csv", header=T, stringsAsFactors=FALSE)


######################################################
######################################################
######################################################


# minor dataframe cleaning:

traits_df1 <- traits_df %>%
  select(-reference, -database, -lengthType, -comments) # drop unnecessary columns


# assign location classes (GoA vs other)
for(i in 1:nrow(traits_df1)){
  if(traits_df1$region[i] %in% c("GoA", "GoA (Kodiak Is)", "GoA (Kodiak Island)", "GoA & Aleutian Islands", "Kodiak Is (GoA")) {
    traits_df1$location[i] <- "GoA"
  } else traits_df1$location[i] <- "other"
}



######################################################
######################################################
######################################################


# Process trait estimates:


horPos_df <- traits_df1[which(traits_df1$trait=='adultSlopeShelf'),]
# standardize adultSlopeShelf values
for(i in 1:nrow(horPos_df)) {
  if(horPos_df$estimate[i] %in% c("coastal", "inshore, intertidal, bays", "nearshore, coastal inlets & rivers")) {horPos_df$adultHorizPosition[i] <- "coastal"}
  else{if(horPos_df$estimate[i] %in% c("inshore, shallow, estuaries, rivers (summer), deep water (winter)",
                                                "ocean, coastal streams", "shelf deep edge (winter), shallow coastal water (summer)")) {horPos_df$adultHorizPosition[i] <- "oceanadromous"}
    else{if(horPos_df$estimate[i] %in% c("shallow shelf", "shelf nearshore")) {horPos_df$adultHorizPosition[i] <- "inner shelf"}
      else{if(horPos_df$estimate[i] %in% c("shelf", "shelf margin (winter), mid/outer shelf (summer)")) {horPos_df$adultHorizPosition[i] <- "shelf"}
        else{if(horPos_df$estimate[i] %in% c("outer shelf", "outer shelf, upper slope", "shelf, slope", "shelf, upper slope")) {horPos_df$adultHorizPosition[i]<- "outer shelf, upper slope"}
          else{if(horPos_df$estimate[i] %in% c("slope", "slope, shelf gullies, deep fjords")) {horPos_df$adultHorizPosition[i] <- "slope"}
            else{horPos_df$estimate[i] <- NA}
          }}}}}}
# for ease of binding with other dfs, remove original estimate column and rename adultHorizPosition to estimate:
horPos_df1 <- horPos_df %>%
  select(-estimate, -region) %>% rename(estimate = adultHorizPosition)




substrate_df <- traits_df1[which(traits_df1$trait=='adultSubstrate'),]
# standardize adultSubstrate values
for(i in 1:nrow(substrate_df)) {
  if(substrate_df$estimate[i] %in% c("mud", "mud, sand", "mud, sand, clay", "sand, silt", "soft")) {substrate_df$adultSubstrateCategory[i] <- "soft"}
  else{if(substrate_df$estimate[i] %in% c("rocky", "rocky?")) {substrate_df$adultSubstrateCategory[i] <- "rocky"}
    else{if(substrate_df$estimate[i] %in% c("hard, rocky, biogenic structure", "rocky, coral, sponges", "rocky, rough, corals")) {substrate_df$adultSubstrateCategory[i] <- "rocky and biogenic"}
      else{if(substrate_df$estimate[i] %in% c("muddy, hard, near rocks or gravel, biogenic structure")) {substrate_df$adultSubstrateCategory[i] <- "soft, hard, biogenic"}
        else{if(substrate_df$estimate[i] %in% c("rocky, kelp, sand")) {substrate_df$adultSubstrateCategory[i] <- "algal associate"}
          else{substrate_df$estimate[i] <- NA}
        }}}}}
# for ease of binding with other dfs, remove original estimate column and rename adultSubstrateCategory to estimate:
substrate_df1 <- substrate_df %>%
  select(-estimate, -region) %>% rename(estimate = adultSubstrateCategory) 




# other traits with a single character value:
trChr_df <- traits_df1[which(traits_df1$trait %in% c('adultWaterColumnPosition', 'trophicPosition', 'diet', 'guild', 'migratoryStatus')),]
trChr_df1 <- trChr_df %>% select(-region)




# calculate gender-specific means of von Bertalanffy K & L_infinity
#kl_df <- traits_df1[which(traits_df1$trait %in% c('K', 'Linfinity')),]
#kl_df1 <- kl_df %>%
#  mutate(estimate1 = as.numeric(estimate)) %>%
#  select(-estimate, -region) %>%
#  group_by(genus.species, common.name, trait, gender, location) %>%
#  summarise_each(funs(mean(., na.rm = TRUE))) %>%
#  ungroup %>%
#  rename(estimate = estimate1)




# retrieve maximum values of longevity, maximum observed size, and maximum depth regardless of gender
max_df <- traits_df1[which(traits_df1$trait %in% c('ageMaximum', 'lengthMaximum', 'depthMax')),]
max_df1 <- max_df %>%
  mutate(estimate1 = as.numeric(estimate)) %>%
  select(-estimate, -region) %>%
  mutate(gender = revalue(gender, c("f"="both", "m" = "both", "u" = "both"))) %>% # replace all values of gender with "both" (because we've calculated maximum values regardless of gender, and to facilitate binding with other dfs below)
  group_by(genus.species, common.name, trait, location) %>%
  summarise_each(funs(max(., na.rm = TRUE))) %>%
  ungroup %>%
  rename(estimate = estimate1)




# depth range: calculate mean breadth of depth range
depth_df <- traits_df1[which(traits_df1$trait %in% c('depthRange')),]
depth_df1 <- depth_df %>%
  mutate(estimate1=strsplit(estimate,split="-") %>% sapply(function(x) x[1])) %>% 
  mutate(estimate2=strsplit(estimate,split="-") %>% sapply(function(x) x[2])) %>%
  mutate(estimate1.num = as.numeric(estimate1), estimate2.num = as.numeric(estimate2)) %>% # change to numeric values
  rowwise() %>%
  mutate(estimate3 = estimate2.num-estimate1.num) %>% # calculate breadth of depth range
  ungroup() %>%
  select(-region, -estimate, -estimate1, -estimate2, -estimate1.num, -estimate2.num) %>%
  rename(estimate = estimate3) %>%
  group_by(genus.species, common.name, trait, gender, location) %>%
  summarise_each(funs(mean)) %>% # calculate means of depth range breadth for each taxa / gender / location combination
  ungroup 
  



# age- and size-at-maturity: calculate gender-specific means of given ranges
maturity_df <- traits_df1[which(traits_df1$trait %in% c('age50percentMaturity', 'firstMaturityAge', 'firstMaturityLength', 'length50percentMaturity')),]

maturity_df1 <- maturity_df %>%
  mutate(estimate = gsub("\\+", "", estimate)) %>% # remove "+"
  filter(estimate != "up to 12") %>% # remove this entry because it's not imformative for determining mean values
  
  # age- and size-at-maturity data here are often given as ranges
  # the following calculates the mean value of the min & max given for these ranges
  mutate(estimate1=strsplit(estimate,split="-") %>% sapply(function(x) x[1])) %>% 
  mutate(estimate2=strsplit(estimate,split="-") %>% sapply(function(x) x[2])) %>%
  # note that when there is a single estimate, rather than a range, this value is entered into estimate1 and estimate2 is NA
  mutate(estimate1.num = as.numeric(estimate1), estimate2.num = as.numeric(estimate2)) %>% # change these from character to numeric values
  rowwise() %>%
  mutate(estimate3 = mean(c(estimate1.num, estimate2.num), na.rm=T)) %>% # calculate mean of the 2 estimates
  ungroup() %>%
  select( -region, -estimate, -estimate1, -estimate2, -estimate1.num, -estimate2.num) %>%
  rename(estimate = estimate3) %>%
  group_by(genus.species, common.name, trait, gender, location) %>%
  summarise_each(funs(mean(., na.rm = TRUE))) %>% # calculate means of all estimates for each taxa / gender / location combination
  ungroup 
# warning re NAs is OK; it's just reporting that estimate2 is NA when the original cell is a single value rather than a range



# bind these dfs together (bind_rows). 
# This creates a dataframe of all the functional trait data we have (these are summarized values; eg means, maxima, etc) 
traits_df3 <- rbind(horPos_df1, substrate_df1, trChr_df1, max_df1, depth_df1, maturity_df1) #kl_df1,



######################################################
######################################################
######################################################

# merge in additional K & L_infinity data from Ben Williams:
KL_df <- read.csv("linf_k.csv", header=T, stringsAsFactors = F)
KL_df1 <- KL_df %>%
  select(-X) %>%
  rename(estimate = value) %>%
  mutate(genus.species = revalue(genus.species, c("Bathyraja.aleutica" = "Bathyraja aleutica", "Pleurogrammus.monopterygius" = "Pleurogrammus monopterygius", 
                                 "Raja.binoculata" = "Raja binoculata", "Sebastes.polyspinis" = "Sebastes polyspinis")))
for(i in 1:nrow(KL_df1)) { # add columns for gender & location
  KL_df1$gender[[i]] <- "goodEnough"
  KL_df1$location[[i]] <- "goodEnough"
}

traits_df4 <- rbind(traits_df3, KL_df1)



######################################################
######################################################
######################################################


# We want to use GoA data wherever it exists; if it doesn't, we'll use data from other locations.
GoA_df <- traits_df4 %>%
  filter(location %in% c("GoA", "goodEnough")) # create table of GoA data

other_df <- traits_df4 %>%
  filter(location == "other") # create table of data for location == "other"



combos <- unique(traits_df4[,c('genus.species','common.name','trait')]) # create a table of all species & trait combinations for which we have data
# specify which gender we'll use for each trait:
for(i in 1:nrow(combos)){
  if(combos$trait[i] %in% c("adultSlopeShelf", "adultSubstrate", "adultWaterColumnPosition", "depthMax", "depthRange", "ageMaximum", "diet", "guild",
                            "lengthMaximum", "migratoryStatus", "trophicPosition")) {combos$gender[i] <- "both"} # for these traits, we'll use data from males & females
  if(combos$trait[i] %in% c("age50percentMaturity", "firstMaturityAge", "firstMaturityLength", "length50percentMaturity" #"K", "Linfinity"
                            )) {combos$gender[i] <- "f"} # for life history traits, for now we'll use only female data
  if(combos$trait[i] %in% c("K", "Linfinity")) {combos$gender[i] <- "goodEnough"} 
  }
traitsGoA_df <- left_join(combos, GoA_df, by = c("genus.species", "common.name", "trait", "gender")) # merge GoA data onto combos




# now figure out where we have trait data for location = "other" but not for GoA:
GoA_df1 <- GoA_df %>% select(-location, -estimate) # remove location & estimate columns from GoA_df to facilitate set difference operation
other_df1 <- other_df %>% select(-location, -estimate) # remove location & estimate columns from other_df to facilitate set difference operation
diffs_df <- setdiff(other_df1, GoA_df1) # retrieve taxa-gender-trait combinations which exist for location == "other" but not GoA
diffsOther_df <- left_join(diffs_df, other_df, by = c("genus.species", "common.name", "trait", "gender")) # merge in other_df, only keeping rows that are in diffs_df



# now merge the table with non-GoA trait data onto the table with GoA data: 
traits_df5 <- left_join(traitsGoA_df, diffsOther_df, by = c("genus.species", "common.name", "trait", "gender")) %>% 
  transmute(genus.species, common.name, trait, gender, 
            location = ifelse(is.na(location.x), location.y, location.x),
            estimate = ifelse(is.na(estimate.x), estimate.y, estimate.x)) %>%
  filter(!is.na(estimate)) %>% # remove rows for which there is no data from either GoA or "other" for some of the desired trait-gender combinations
  # the next 6 lines make taxonomic names match style in abundance files
  mutate(genus.species = revalue(genus.species, c("Clupea pallasii" = "Clupea pallasi",
                                                  "Sebastes group 1" = "Dusky and Dark Rockfish", 
                                                  "Sebastes group 2" = "Rougheye and Blackspotted Rockfish", 
                                                  "Lepidopsetta spp." = "Lepidopsetta sp.", 
                                                  "Myctophidae spp." = "Myctophidae",
                                                  "Theragra chalcogramma" = "Gadus chalcogrammus"))) %>%
  mutate(genus.species = gsub(" ", ".", genus.species)) %>%
  rename(Species = genus.species)

for(i in 1:nrow(traits_df5)) {
  if(traits_df5$Species[i] == "Dusky.and.Dark.Rockfish") {traits_df5$common.name[i] <- "sebastes group 1"} 
}

View(traits_df5)

######################################################
######################################################
######################################################


# Organize functional trait data for analysis in FD package:

# Convert to wide format:
traits_df5$row <- 1:nrow(traits_df5) # create column of unique identifiers to facilitate data spread
traits_wide <- traits_df5 %>%
  spread(trait, estimate) %>%
  select(-gender, -location, -row) %>%
  group_by(Species, common.name) %>%
  summarize_each(funs(first(., order_by = is.na(.)))) %>%
  ungroup()
# fill in missing NAs?
cols = c(6:9, 11, 12, 14:17, 19); traits_wide[,cols] <- apply(traits_wide[,cols], 2, function(x) as.numeric(x)) # convert columns to numeric as needed
cols1 = c(3:5, 10, 13, 18); traits_wide[,cols1] <- lapply(traits_wide[,cols1] , factor) # convert columns to factor as needed
#View(traits_wide)
write.csv(traits_wide, file = "traits_wide.csv")



# we have the most data for:
# adultWaterColumnPosition (missing for 2 taxa)
#ageMax (missing for 9 taxa)
#diet (missing 0)
#sum(is.na(traits_wide$firstMaturityLength)) # missing 23
sum(is.na(traits_wide$K)) # missing 16
sum(is.na(traits_wide$Linfinity)) # missing 16
#guild (missing 0)
#lengthMaximum (missing 4)
sum(is.na(traits_wide$ageMaximum)) # missing 10
#trophicPosition (missing 0)
#depthRange (missing 6)
#depthMax (missing 6)


# which pairs of traits are correlated?
pairs.panels(traits_wide[,c(3:19)],smooth=F,density=T,ellipses=F,lm=T,digits=3,scale=T)
names(traits_wide)
# significant correlations:
# adultSlopeShelf, adultSubstrate
# K, firstMaturityAge, age50percentMaturity, ageMaximum
# Linfinity, trophicPosition, firstMaturityLength, length50percentMaturity, lengthMaximum
# depthRange & depthMax





ft_df <- traits_wide %>%
  select(Species, lengthMaximum, ageMaximum, depthMax, trophicPosition, adultWaterColumnPosition) %>% # select traits for which we have the most data
  filter(!(is.na(depthMax)), !(is.na(ageMaximum))) %>%
  filter(!(Species %in% c("Hydrolagus.colliei", "Merluccius.productus", "Sebastes.helvomaculatus"))) %>% # remove taxa for which there is no abundance data; there is also no abund data for Berryteuthis.magister
  arrange(Species)
rownames(ft_df) <- ft_df$Species # create row names from Species column
ft_df <- ft_df %>% select(-Species)
#View(ft_df)  # this is the dataframe we'll use for functional diversity analyses

#unique(sort(setdiff(ft_df$Species, SPCPUEArea$Species))) # sp in traits_df3 but not SPCPUEArea


######################################################
######################################################
######################################################


# load taxonomic relative presence data:
#URL_SpByArea <- "https://drive.google.com/uc?export=download&id=0By1iaulIAI-udlVNME9rQXEwZ1k"
#SpByArea_Get <- GET(URL_SpByArea)
#SpByArea_1 <- content(SpByArea_Get, as='text')
#SpByArea <- read.csv(file=textConnection(SpByArea_1),stringsAsFactors=FALSE,head=TRUE)
#View(SpByArea)

# load taxonomic abundance data:
URL_SPCPUEArea <- "https://drive.google.com/uc?export=download&id=0By1iaulIAI-udm1FT2trQUh5N1k"
SPCPUEArea_Get <- GET(URL_SPCPUEArea)
SPCPUEArea_1 <- content(SPCPUEArea_Get, as='text')
SPCPUEArea <- read.csv(file=textConnection(SPCPUEArea_1),stringsAsFactors=FALSE,head=TRUE)
#View(SPCPUEArea)

# NB  SPCPUEArea has only 53 taxa, not 57. Which ones are missing?
#spDiffs <- setdiff(traits_wide$Species, SPCPUEArea$Species); spDiffs
# "Berryteuthis.magister"   "Hydrolagus.colliei"      "Merluccius.productus"    "Sebastes.helvomaculatus"


# organize abundance data for analysis in FD package:


#unique(sort(setdiff(SPCPUEArea$Species, ft_df$Species))) # 

sp_df <- SPCPUEArea %>%
  select(area, year, Species, Mean.totalDensity) %>%
  filter(!(Species %in% c("Chionoecetes.bairdi", "Hemitripterus.bolini", "Hyas.lyratus", "Lycodes.brevipes", 
                          "Lycodes.palearis", "Lyopsetta.exilis", "Myctophidae", "Oncorhynchus.keta", 
                          "Oncorhynchus.tshawytscha"))) %>% # remove taxa for which we don't have all trait data
  mutate(area = revalue(area, c("Total" = "12")), # recode Total for looping later
         area = as.numeric(area)) # convert to numeric class
View(sp_df)



A <- sp_df %>% 
  select(Species, area, year, Mean.totalDensity) %>%
  arrange(Species) %>% # arrange in alphabetical order to match order in functional traits df (required by FD package)
  spread(Species, Mean.totalDensity) %>%
  select(-year)
byArea_list <- split(A, f = A$area) # create a list of dataframes (one for each area; NB area 12 is Total)

byArea_list
byArea_list[[1]]
byArea_list$`1`
byArea_list$`5`
byArea_list$`1`[2]

is.data.frame(byArea_list$`1`) # TRUE

byArea_list1 <- lapply(byArea_list, function(x) x[!(names(x) %in% c("area", "year"))]) # drop area & year

fd <- list()
for (i in seq_along(byArea_list1)) {
  fd[[i]] <- dbFD(ft_df, byArea_list1[[i]], calc.FRic = F, calc.CWM = F, calc.FDiv = F)
}

fd[[i]]$RaoQ

is.data.frame(fd[[1]]) # FALSE
is.data.frame(fd[[i]]$RaoQ) # FALSE
is.list(fd[[i]]$RaoQ) # FALSE

year <- unique(sort(SPCPUEArea$year))

Rao1 <- list()
for (i in seq_along(fd)) {
  Rao1[[i]] <- data.frame(fd[[i]]$RaoQ)
}
is.data.frame(Rao1) #FALSE

Rao <- data.frame(year, fd1$RaoQ, fd2$RaoQ, fd3$RaoQ, fd4$RaoQ, fd5$RaoQ, fd6$RaoQ, 
                  fd7$RaoQ, fd8$RaoQ, fd9$RaoQ, fd10$RaoQ, fd11$RaoQ)






# make a list of all dataframes
my.list <- list(d1, d2)
# or
mylist <- list()
mylist[[1]] <- mtcars
mylist[[2]] <- data.frame(a = rnorm(50), b = runif(50))


my_data <- list()
for (i in seq_along(my_files)) {
  my_data[[i]] <- read.csv(file = my_files[i])
}

my_data <- lapply(my_files, read.csv)
my_dfs <- lapply(byArea_list, createAreaDf)

names(my_data) <- gsub("\\.csv", "", my_files)
# or, if you prefer the consistent syntax of stringr
names(my_data) <- stringr::str_replace(my_files, pattern = ".csv", replacement = "")



Splitting a data frame into a list of data frames
This is super-easy, the base function split() does it for you. You can split by a column (or columns) of the data, or by anything else you want
mt_list = split(mtcars, f = mtcars$cyl)
# This gives a list of three data frames, one for each value of cyl




  


# create dataframes for each area
area1 <- sp_df %>% 
  filter(area == 1) %>%
  select(Species, year, Mean.totalDensity) %>%
  arrange(Species) %>% # arrange in alphabetical order to match order in functional traits df (required by FD package)
  spread(Species, Mean.totalDensity) %>%
  select(-year)

area2 <- sp_df %>% 
  filter(area == 2) %>%
  select(Species, year, Mean.totalDensity) %>%
  arrange(Species) %>% # arrange in alphabetical order to match order in functional traits df (required by FD package)
  spread(Species, Mean.totalDensity) %>%
  select(-year)

area3 <- sp_df %>% 
  filter(area == 3) %>%
  select(Species, year, Mean.totalDensity) %>%
  arrange(Species) %>% # arrange in alphabetical order to match order in functional traits df (required by FD package)
  spread(Species, Mean.totalDensity) %>%
  select(-year)

area4 <- sp_df %>% 
  filter(area == 4) %>%
  select(Species, year, Mean.totalDensity) %>%
  arrange(Species) %>% # arrange in alphabetical order to match order in functional traits df (required by FD package)
  spread(Species, Mean.totalDensity) %>%
  select(-year)

area5 <- sp_df %>% 
  filter(area == 5) %>%
  select(Species, year, Mean.totalDensity) %>%
  arrange(Species) %>% # arrange in alphabetical order to match order in functional traits df (required by FD package)
  spread(Species, Mean.totalDensity) %>%
  select(-year)

area6 <- sp_df %>% 
  filter(area == 6) %>%
  select(Species, year, Mean.totalDensity) %>%
  arrange(Species) %>% # arrange in alphabetical order to match order in functional traits df (required by FD package)
  spread(Species, Mean.totalDensity) %>%
  select(-year)

area7 <- sp_df %>% 
  filter(area == 7) %>%
  select(Species, year, Mean.totalDensity) %>%
  arrange(Species) %>% # arrange in alphabetical order to match order in functional traits df (required by FD package)
  spread(Species, Mean.totalDensity) %>%
  select(-year)

area8 <- sp_df %>% 
  filter(area == 8) %>%
  select(Species, year, Mean.totalDensity) %>%
  arrange(Species) %>% # arrange in alphabetical order to match order in functional traits df (required by FD package)
  spread(Species, Mean.totalDensity) %>%
  select(-year)

area9 <- sp_df %>% 
  filter(area == 9) %>%
  select(Species, year, Mean.totalDensity) %>%
  arrange(Species) %>% # arrange in alphabetical order to match order in functional traits df (required by FD package)
  spread(Species, Mean.totalDensity) %>%
  select(-year)

area10 <- sp_df %>% 
  filter(area == 10) %>%
  select(Species, year, Mean.totalDensity) %>%
  arrange(Species) %>% # arrange in alphabetical order to match order in functional traits df (required by FD package)
  spread(Species, Mean.totalDensity) %>%
  select(-year)

area11 <- sp_df %>% 
  filter(area == 11) %>%
  select(Species, year, Mean.totalDensity) %>%
  arrange(Species) %>% # arrange in alphabetical order to match order in functional traits df (required by FD package)
  spread(Species, Mean.totalDensity) %>%
  select(-year)



######################################################
######################################################
######################################################


# calculate & plot Rao's Q

fd1 <- dbFD(ft_df, area1, calc.FRic = F, calc.CWM = F, calc.FDiv = F)
fd2 <- dbFD(ft_df, area2, calc.FRic = F, calc.CWM = F, calc.FDiv = F)
fd3 <- dbFD(ft_df, area3, calc.FRic = F, calc.CWM = F, calc.FDiv = F)
fd4 <- dbFD(ft_df, area4, calc.FRic = F, calc.CWM = F, calc.FDiv = F)
fd5 <- dbFD(ft_df, area5, calc.FRic = F, calc.CWM = F, calc.FDiv = F)
fd6 <- dbFD(ft_df, area6, calc.FRic = F, calc.CWM = F, calc.FDiv = F)
fd7 <- dbFD(ft_df, area7, calc.FRic = F, calc.CWM = F, calc.FDiv = F)
fd8 <- dbFD(ft_df, area8, calc.FRic = F, calc.CWM = F, calc.FDiv = F)
fd9 <- dbFD(ft_df, area9, calc.FRic = F, calc.CWM = F, calc.FDiv = F)
fd10 <- dbFD(ft_df, area10, calc.FRic = F, calc.CWM = F, calc.FDiv = F)
fd11 <- dbFD(ft_df, area11, calc.FRic = F, calc.CWM = F, calc.FDiv = F)
# for each area:
# "Species x species distance matrix was not Euclidean. 'sqrt' correction was applied."


year <- unique(sort(SPCPUEArea$year))
Rao <- data.frame(year, fd1$RaoQ, fd2$RaoQ, fd3$RaoQ, fd4$RaoQ, fd5$RaoQ, fd6$RaoQ, 
                  fd7$RaoQ, fd8$RaoQ, fd9$RaoQ, fd10$RaoQ, fd11$RaoQ)



ggplot(data=Rao, aes(x=year, y = value)) + 
  geom_point(aes(y = fd1.RaoQ), size=2) +
  geom_point(aes(y = fd2.RaoQ), size=2) +
  geom_point(aes(y = fd3.RaoQ), size=2, col=2) +
  geom_point(aes(y = fd4.RaoQ), size=2, col=2) +
  geom_point(aes(y = fd5.RaoQ), size=2, col=2) +
  geom_point(aes(y = fd6.RaoQ), size=2) +
  geom_point(aes(y = fd7.RaoQ), size=2) +
  geom_point(aes(y = fd8.RaoQ), size=2) +
  geom_point(aes(y = fd9.RaoQ), size=2) +
  geom_point(aes(y = fd10.RaoQ), size=2) +
  geom_point(aes(y = fd11.RaoQ), size=2) +
  
  geom_line(aes(y = fd1.RaoQ), size=2) +
  geom_line(aes(y = fd2.RaoQ), size=2) +
  geom_line(aes(y = fd3.RaoQ), size=2, col=2) +
  geom_line(aes(y = fd4.RaoQ), size=2, col=2) +
  geom_line(aes(y = fd5.RaoQ), size=2, col=2) +
  geom_line(aes(y = fd6.RaoQ), size=2) +
  geom_line(aes(y = fd7.RaoQ), size=2) +
  geom_line(aes(y = fd8.RaoQ), size=2) +
  geom_line(aes(y = fd9.RaoQ), size=2) +
  geom_line(aes(y = fd10.RaoQ), size=2) +
  geom_line(aes(y = fd11.RaoQ), size=2) +
  
  theme(axis.line=element_line('black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  theme(axis.text.x = element_text(angle=90, size=18, colour = "black"))+
  theme(axis.text.y = element_text(size=22))+
  scale_x_continuous(breaks=c(year), labels=c(year)) +
  ylab("Rao's Q") +
  xlab("Year") 


# get Euclidean distance matrix from traits
#trait.dist <- dist(trait)
