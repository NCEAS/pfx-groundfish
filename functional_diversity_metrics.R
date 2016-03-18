# Functional Diversity Metrics Processing Script
# Colette Ward, Feb 29 2016

# load packages
library(httr)
library(plyr)
library(dplyr)
library(tidyr)
#library(FD)

# load the data
URL_traits <- "https://drive.google.com/uc?export=download&id=0B1XbkXxdfD7uZjdLaHVRc3pveE0"
traitsGet <- GET(URL_traits)
traits1 <- content(traitsGet, as='text')
traits_df <- read.csv(file=textConnection(traits1),stringsAsFactors=FALSE)
#View(traits_df)

traits_df1 <- traits_df %>%
  select(-reference, -database, -lengthType, -comments) # drop unnecessary columns


# assign location classes (GoA vs other)
for(i in 1:nrow(traits_df1)){
  if(traits_df1$region[i] %in% c("GoA", "GoA (Kodiak Is)", "GoA (Kodiak Island)", "GoA & Aleutian Islands", "Kodiak Is (GoA")) {
    traits_df1$location[i] <- "GoA"
  } else traits_df1$location[i] <- "other"
}





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
            #else{horPos_df$estimate[i] <- NA}
          }}}}}#}
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
          #else{substrate_df$estimate[i] <- NA}
        }}}}#}
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




# retrieve maximum values of longevity and maximum observed size, regardless of gender
max_df <- traits_df1[which(traits_df1$trait %in% c('ageMaximum', 'lengthMaximum')),]
max_df1 <- max_df %>%
  mutate(estimate1 = as.numeric(estimate)) %>%
  select(-estimate, -region) %>%
  mutate(gender = revalue(gender, c("f"="both", "m" = "both", "u" = "both"))) %>% # replace all values of gender with "both" (because we've calculated maximum values regardless of gender, and to facilitate binding with other dfs below)
  group_by(genus.species, common.name, trait, location) %>%
  summarise_each(funs(max(., na.rm = TRUE))) %>%
  ungroup %>%
  rename(estimate = estimate1)
# warning re NAs occurs because of 80-100 value for Squalus acanthias 




# age- and size-at-maturity: calculate gender-specific means of given ranges
maturity_df <- traits_df1[which(traits_df1$trait %in% c('age50percentMaturity', 'firstMaturityAge', 'firstMaturityLength', 'length50percentMaturity')),]

maturty_df1 <- maturity_df %>%
  mutate(estimate = gsub("\\+", "", estimate)) %>% # remove "+"
  filter(estimate != "up to 12") %>% # remove this entry because it's not imformative for determining mean values
  
  # age- and size-at-maturity data here are often given as ranges
  # the following calculates the mean value of the min & max given for these ranges
  mutate(estimate1=strsplit(estimate,split="-") %>% sapply(function(x) x[1])) %>% 
  mutate(estimate2=strsplit(estimate,split="-") %>% sapply(function(x) x[2])) %>%
  # note that when there is a single estimate, rather than a range, this value is entered into estimate1 and estimate2 is NA
  mutate(estimate1.num = as.numeric(estimate1), estimate2.num = as.numeric(estimate2)) %>%
  rowwise() %>%
  mutate(estimate3 = mean(c(estimate1.num, estimate2.num), na.rm=T)) %>%
  ungroup() %>%
  select( -region, -estimate, -estimate1, -estimate2, -estimate1.num, -estimate2.num) %>%
  rename(estimate = estimate3)
# warning re NAs is OK; it's just reporting that estimate2 is NA when the original cell is a single value rather than a range




# bind these dfs together (bind_rows)
traits_df3 <- rbind(horPos_df1, substrate_df1, trChr_df1, kl_df1, max_df1, maturty_df1)
#View(traits_df3)


# we want to use GoA data wherever it exists; if it doesn't, we'll use data from other locations:
GoA_df <- traits_df3 %>%
  filter(location == "GoA") # create table of GoA data

other_df <- traits_df3 %>%
  filter(location == "other") #, subsetVector != T) # create table of data for location == "other"



combos <- unique(traits_df3[,c('genus.species','trait')]) # create df of all species & trait combinations for which we have data
# use a loop to add vector of desired gender for all traits
for(i in 1:nrow(combos)){
  if(combos$trait[i] %in% c("adultSlopeShelf", "adultSubstrate", "adultWaterColumnPosition", "ageMaximum", "diet", "guild",
                            "lengthMaximum", "migratoryStatus", "trophicPosition")) {combos$gender[i] <- "both"} # for these traits, we'll use data from males & females
  if(combos$trait[i] %in% c("age50percentMaturity", "firstMaturityAge", "firstMaturityLength", "length50percentMaturity",
                            "K", "Linfinity")) {combos$gender[i] <- "f"} # for life history traits, for now we'll use only female data
  }


View(combos)
# now merge in GoA values
traits_df4 <- left_join(combos, GoA_df)
View(traits_df4) # why is traits_df missing most values for common.name?

# then merge in location == other values, preferentially keep GoA value if it's already in there


diffs <- setdiff(other_df, GoA_df) # retrieve data which exist for location == "other" but not GoA
View(diffs)
traits_df5 <- left_join(traits_df4, diffs, by = c("genus.species", "common.name", "trait", "gender")) # does not work; contains multiple columns for location & estimate
View(traits_df5)
# add vector of T to every row in diffs







traits_df3$row <- 1:nrow(traits_df3) # create column of unique identifiers to facilitate data spread
traits_wide <- traits_df3 %>%
  #group_by(genus.species, common.name, gender, location) %>%
  spread(trait, estimate)
View(traits_wide)
# the data are transposed, but df is otherwise full of NAs! (rows did not collapse into single rows for each set of   genus.species, common.name, gender, location)



# also need to add:
# for traits %in% c(x, ..., n), use female data. if there is no female data, use "u" or "both"
# for traits %in% c(y, ..., n), use GoA data. if there is no GoA data, use "other" 




# problems following spread:
# merluccius productus trophic position
# squalus acanthias ageMaximum, both, GoA (NA) 
# limanda aspera length50Mat, f, other NaN


# later on, merge traits_wide1 with wide dataframes created below



traits_df1$row <- 1:nrow(traits_df1) # create column of unique identifiers to facilitate data spread
traits_wide2 <- spread(traits_df1, trait, estimate) # convert data to columns


# groups of traits that need the same treatment:
# 1. calculate means of female values of life history traits (usually have for either female or unknown gender or both)
trNum1 <- traits_wide2 %>% 
  select(genus.species, common.name, gender, location, K, Linfinity, 
         depthRangeShallow, depthRangeDeep, temperatureLower, temperatureUpper) %>%
  mutate_each(funs(as.numeric), K:temperatureUpper) %>%
  group_by(genus.species, common.name, gender, location) %>%
  summarise_each(funs(mean(., na.rm = TRUE))) %>%
  ungroup
View(trNum1)


# do (depthRangeShallow, depthRangeDeep, temperatureLower, temperatureUpper) separately because I will ignore gender
# also, calculate breadth of depth & temp ranges


# 2. calculate maximum value of ageMaximum and lengthMaximum, regardless of gender
trNum2 <- traits_wide2 %>% 
  select(genus.species, common.name, location, ageMaximum, lengthMaximum) %>%
  mutate_each(funs(as.numeric), ageMaximum:lengthMaximum) %>%
  group_by(genus.species, common.name, location) %>%
  summarise_each(funs(max(., na.rm = TRUE))) %>%
  ungroup
View(trNum2)


unique(sort(traits_wide$firstMaturityLength))

# numeric entries which are ranges separated by hyphens
trNum3 <- traits_wide2 %>% 
  select(genus.species, common.name, gender, location,
         age50percentMaturity, firstMaturityAge, firstMaturityLength, length50percentMaturity) %>%
  
  # remove "+" and "up to ..." values
  mutate(age50percentMaturity = gsub("\\+", "", age50percentMaturity)) %>% # remove "+"
  mutate(length50percentMaturity = gsub("\\+", "", length50percentMaturity)) %>%
  mutate(firstMaturityAge = gsub("\\+", "", firstMaturityAge)) %>%
  filter(firstMaturityAge != "up to 12") %>% # remove this entry because it's not imformative
  mutate(firstMaturityLength = gsub("\\+", "", firstMaturityLength)) %>%
  
  # for the function:
  # strip split
  # mutate to as.numeric
  # rowwise()
  # take means across rows
  # ungroup()
  
  
  
  # clean up age at 50% Maturity
  mutate(age50percentMaturity = gsub("\\+", "", age50percentMaturity)) %>% # remove "+"
  mutate(age50Mat1a=strsplit(age50percentMaturity,split="-") %>%
           sapply(function(x) x[1])) %>%
  mutate(age50Mat2a=strsplit(age50percentMaturity,split="-") %>%
           sapply(function(x) x[2])) %>%
  mutate(age50percentMaturity.a = as.numeric(age50percentMaturity), age50Mat1 = as.numeric(age50Mat1a), age50Mat2 = as.numeric(age50Mat2a)) %>%
  rowwise() %>%
  mutate(age50Mat = mean(c(age50percentMaturity.a, age50Mat1, age50Mat2), na.rm=T)) %>%
  

# clean up age at 50% Maturity
  mutate(length50percentMaturity = gsub("\\+", "", length50percentMaturity)) %>% # remove "+"
  mutate(length50Mat1a=strsplit(length50percentMaturity,split="-") %>%
           sapply(function(x) x[1])) %>%
  mutate(length50Mat2a=strsplit(length50percentMaturity,split="-") %>%
           sapply(function(x) x[2])) %>%
  mutate(length50percentMaturity.a = as.numeric(length50percentMaturity), length50Mat1 = as.numeric(length50Mat1a), length50Mat2 = as.numeric(length50Mat2a)) %>%
  rowwise() %>%
  mutate(length50Mat = mean(c(length50percentMaturity.a, length50Mat1, length50Mat2), na.rm=T)) %>%

  
  
  # clean up Age at First Maturity
  mutate(firstMaturityAge = gsub("\\+", "", firstMaturityAge)) %>% # remove "+"
  filter(firstMaturityAge != "up to 12") %>% # remove this entry because it's not imformative
  mutate(fMatAge1a=strsplit(firstMaturityAge,split="-") %>%
           sapply(function(x) x[1])) %>%
  mutate(fMatAge2a=strsplit(firstMaturityAge,split="-") %>%
           sapply(function(x) x[2])) %>%
  mutate(firstMaturityAge.a = as.numeric(firstMaturityAge), fMatAge1 = as.numeric(fMatAge1a), fMatAge2 = as.numeric(fMatAge2a)) %>%
  rowwise() %>%
  mutate(firstMatAge = mean(c(firstMaturityAge.a, fMatAge1, fMatAge2), na.rm=T)) 


# clean up Length at First Maturity
  mutate(firstMaturityLength = gsub("\\+", "", firstMaturityLength)) %>% # remove "+"
  mutate(fMatLength1a=strsplit(firstMaturityLength,split="-") %>%
           sapply(function(x) x[1])) %>%
  mutate(fMatLength2a=strsplit(firstMaturityLength,split="-") %>%
           sapply(function(x) x[2])) %>%
  mutate(firstMaturityLength.a = as.numeric(firstMaturityLength), fMatLength1 = as.numeric(fMatLength1a), fMatAge2 = as.numeric(fMatAge2a)) %>%
  rowwise() %>%
  mutate(firstMatAge = mean(c(firstMaturityLength.a, fMatAge1, fMatAge2), na.rm=T)) 




View(trNum3)
rm(trNum3)

    


# select which values to put into final dataset:
# for K and Linfinity, first and mean age at maturity, use gender == f
# for lengthMaximum and ageMaximum, use max value regardless of gender value

# if there is a value for GoA, use that, otherwise use value for location = other





df1 <- traits_wide %>%
  select(-row) %>%
  mutate(ageMaximum = as.numeric(ageMaximum))
ageMaximum, depthRangeDeep, depthRangeShallow))
View(df1)
rm(df1)

# for firstMaturityLength and firstMaturityAge, take the first value when there is a range
# for adultSubstrate, convert mud, sand, clay to soft


# select another df for which there is no GoA data
# if there is no GoA value (if region != "GoA"), then take mean of values in column



# calculate Rao's Q



# use function dbFD
# need dataframe of functional traits (x); species are rows
# need matrix of abundances of species in x; rows are sites, species are columns

ex1 <- dbFD(dummy$trait, dummy$abun)
ex1


#########################################################
# example from https://stat.ethz.ch/pipermail/r-sig-ecology/2010-November/001645.html
# get traits and abundances
trait <- runif(30) # trait with uniform values
abun <- rlnorm(30) # lognormal abundance distribution

# get names for both vectors
#names(trait) <- paste("sp", 1:30, sep = "")
#names(abun) <- paste("sp", 1:30, sep = "")

# get Euclidean distance matrix from traits
#trait.dist <- dist(trait)

# get Rao's Q and FDis using the dbFD function in the FD library
test1 <- dbFD(trait.dist, abun, calc.FRic = F, calc.CWM = F, calc.FDiv = F)
rao <- test1$RaoQ
fdis <- test1$FDis