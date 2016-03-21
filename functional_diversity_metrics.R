# Functional Diversity Metrics Processing Script
# Colette Ward, Feb 29 2016

# load packages
library(httr)
library(plyr)
library(dplyr)
library(tidyr)
library(FD)

# load the functional trait data
URL_traits <- "https://drive.google.com/uc?export=download&id=0B1XbkXxdfD7uM2M1UnhtTzlGZGM"
traitsGet <- GET(URL_traits)
traits1 <- content(traitsGet, as='text')
traits_df <- read.csv(file=textConnection(traits1),stringsAsFactors=FALSE)

# or load it from local source:
#traits_df <- read.csv("Groundfish-Functional-Diversity-Traits1.csv", header=T, stringsAsFactors=FALSE)


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
kl_df <- traits_df1[which(traits_df1$trait %in% c('K', 'Linfinity')),]
kl_df1 <- kl_df %>%
  mutate(estimate1 = as.numeric(estimate)) %>%
  select(-estimate, -region) %>%
  group_by(genus.species, common.name, trait, gender, location) %>%
  summarise_each(funs(mean(., na.rm = TRUE))) %>%
  ungroup %>%
  rename(estimate = estimate1)




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
traits_df3 <- rbind(horPos_df1, substrate_df1, trChr_df1, kl_df1, max_df1, depth_df1, maturity_df1)




######################################################
######################################################
######################################################


# We want to use GoA data wherever it exists; if it doesn't, we'll use data from other locations.
GoA_df <- traits_df3 %>%
  filter(location == "GoA") # create table of GoA data

other_df <- traits_df3 %>%
  filter(location == "other") # create table of data for location == "other"



combos <- unique(traits_df3[,c('genus.species','common.name','trait')]) # create a table of all species & trait combinations for which we have data
# specify which gender we'll use for each trait:
for(i in 1:nrow(combos)){
  if(combos$trait[i] %in% c("adultSlopeShelf", "adultSubstrate", "adultWaterColumnPosition", "depthMax", "depthRange", "ageMaximum", "diet", "guild",
                            "lengthMaximum", "migratoryStatus", "trophicPosition")) {combos$gender[i] <- "both"} # for these traits, we'll use data from males & females
  if(combos$trait[i] %in% c("age50percentMaturity", "firstMaturityAge", "firstMaturityLength", "length50percentMaturity",
                            "K", "Linfinity")) {combos$gender[i] <- "f"} # for life history traits, for now we'll use only female data
  }
traitsGoA_df <- left_join(combos, GoA_df, by = c("genus.species", "common.name", "trait", "gender")) # merge GoA data onto combos




# now figure out where we have trait data for location = "other" but not for GoA:
GoA_df1 <- GoA_df %>% select(-location, -estimate) # remove location & estimate columns from GoA_df to facilitate set difference operation
other_df1 <- other_df %>% select(-location, -estimate) # remove location & estimate columns from other_df to facilitate set difference operation
diffs_df <- setdiff(other_df1, GoA_df1) # retrieve taxa-gender-trait combinations which exist for location == "other" but not GoA
diffsOther_df <- left_join(diffs_df, other_df, by = c("genus.species", "common.name", "trait", "gender")) # merge in other_df, only keeping rows that are in diffs_df



# now merge the table with non-GoA trait data onto the table with GoA data: 
traits_df4 <- left_join(traitsGoA_df, diffsOther_df, by = c("genus.species", "common.name", "trait", "gender")) %>% 
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



######################################################
######################################################
######################################################


# Organize functional trait data for analysis in FD package:

# Convert to wide format:
traits_df4$row <- 1:nrow(traits_df4) # create column of unique identifiers to facilitate data spread
traits_wide <- traits_df4 %>%
  spread(trait, estimate) %>%
  select(-gender, -location, -row) %>%
  group_by(Species, common.name) %>%
  summarize_each(funs(first(., order_by = is.na(.)))) %>%
  ungroup()
# fill in missing NAs?
cols = c(6:9, 11, 12, 14:17, 19); traits_wide[,cols] <- apply(traits_wide[,cols], 2, function(x) as.numeric(x)) # convert columns to numeric as needed
cols1 = c(3:5, 10, 13, 18); traits_wide[,cols1] <- lapply(traits_wide[,cols1] , factor) # convert columns to factor as needed
#write.csv(traits_wide, file = "traits_wide.csv")


ft_df <- traits_wide %>%
  select(Species, lengthMaximum, trophicPosition, diet, guild, adultWaterColumnPosition) %>% # select traits for which we have the most data
  filter(!(Species %in% c("Berryteuthis.magister", "Chionoecetes.bairdi", "Hyas.lyratus", "Myctophidae", "Lepidopsetta.sp.", "Dusky.and.Dark.Rockfish"))) %>% # remove taxa for which some trait data is missing
  filter(!(Species %in% c("Hydrolagus.colliei", "Merluccius.productus", "Sebastes.helvomaculatus"))) %>% # remove taxa for which there is no abundance data; there is also no abund data for Berryteuthis.magister
  arrange(Species)
rownames(ft_df) <- ft_df$Species # create row names from Species column
ft_df <- ft_df %>% select(-Species)
# View(ft_df)  # this is the dataframe we'll use for functional diversity analyses


# we have the most data for:
# adultWaterColumnPosition (missing for 2 taxa)
#ageMax (missing for 9 taxa)
#diet (missing 0)
#sum(is.na(traits_wide$firstMaturityLength)) # missing 23
#sum(is.na(traits_wide$K)) # missing 27
#guild (missing 0)
#lengthMaximum (missing 4)
#trophicPosition (missing 0)
#depthRange (missing 6)
#depthMax (missing 6)



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

sp_df <- SPCPUEArea %>%
  select(area, year, Species, Mean.totalDensity) %>%
  filter(!(Species %in% c("Berryteuthis.magister", "Chionoecetes.bairdi", "Hyas.lyratus", 
                          "Myctophidae", "Lepidopsetta.sp.", "Dusky.and.Dark.Rockfish"))) %>% # remove taxa for which we don't have all trait data
  mutate(area = revalue(area, c("Total" = "12")), # recode Total for looping later
         area = as.numeric(area)) # convert to numeric class
#View(sp_df)




ar <- seq(1:12) 
createAreaDf <- function(sp_df){
  A <- sp_df %>% 
  filter(area == ar) %>%
  select(Species, year, Mean.totalDensity) %>%
  arrange(Species) %>% # arrange in alphabetical order to match order in functional traits df (required by FD package)
  spread(Species, Mean.totalDensity) %>%
  select(-year)
  return(A)
}

createAreaDf(sp_df)

for(i in 1:12) {
#for(i in 1:length(unique(sort(sp_df$area)))) {
  B = 
    sp_df %>% 
    filter(sp_df$area == i) %>%
    select(Species, year, Mean.totalDensity) %>%
    arrange(Species) %>% # arrange in alphabetical order to match order in functional traits df (required by FD package)
    spread(Species, Mean.totalDensity) %>%
    select(-year)
}


area11 <- sp_df %>% 
  filter(area == "11") %>%
  select(Species, year, Mean.totalDensity) %>%
  arrange(Species) %>% # arrange in alphabetical order to match order in functional traits df (required by FD package)
  spread(Species, Mean.totalDensity) %>%
  select(-year)
View(area11)



######################################################
######################################################
######################################################


# calculate Functional Diversity metric (Rao's Q):

# use function dbFD
# need dataframe of functional traits (x); species are rows
# need matrix of abundances of species in x; rows are sites, species are columns

a11 <- dbFD(ft_df, area11, calc.FRic = F, calc.CWM = F, calc.FDiv = F)
a11$RaoQ


