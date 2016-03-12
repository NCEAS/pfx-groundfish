# Functional Diversity Metrics Processing Script
# Colette Ward, Feb 29 2016

# load packages
library(httr)
library(dplyr)
library(tidyr)
library(FD)

# load the data
URL_traits <- "https://drive.google.com/uc?export=download&id=0B1XbkXxdfD7uaDhWeV9tazRGLXM"
traitsGet <- GET(URL_traits)
traits1 <- content(traitsGet, as='text')
traits_df <- read.csv(file=textConnection(traits1),stringsAsFactors=FALSE)
View(traits_df)

# or call the data assembly script
#source("Groundfish-Functional-Diversity-Traits.csv")


traits_df1 <- traits_df %>%
  select(-reference, -database, -lengthType, -comments) # drop unnecessary columns


# assign location classes (GoA vs other)
for(i in 1:nrow(traits_df1)){
  if(traits_df1$region[i] %in% c("GoA", "GoA (Kodiak Is)", "GoA (Kodiak Island)", "GoA & Aleutian Islands", "Kodiak Is (GoA")) {
    traits_df1$location[i] <- "GoA"
  } else traits_df1$location[i] <- "other"
}


traits_df1$row <- 1:nrow(traits_df1) # create column of unique identifiers to facilitate data spread
traits_wide <- spread(traits_df1, trait, estimate) # convert data to columns


# standardize adultSubstrate values
for(i in 1:nrow(traits_wide)) {
  if(traits_wide$adultSubstrate[i] %in% c("mud", "mud, sand", "mud, sand, clay", "sand, silt", "soft")) {traits_wide$adultSubstrateCategory[i] <- "soft"}
  
  else{if(traits_wide$adultSubstrate[i] %in% c("rocky", "rocky?")) {traits_wide$adultSubstrateCategory[i] <- "hard"}
    
    else{if(traits_wide$adultSubstrate[i] %in% c("hard, rocky, biogenic structure", "rocky, coral, sponges", "rocky, rough, corals")) {traits_wide$adultSubstrateCategory[i] <- "hard and biogenic"}
      
      else{if(traits_wide$adultSubstrate[i] %in% c("muddy, hard, near rocks or gravel, biogenic structure", "rocky, kelp, sand")) {traits_wide$adultSubstrateCategory[i] <- "soft, hard, biogenic"}
        
        else{traits_wide$adultSubstrateCategory[i] <- NA}
      
    }}}}


# standardize adultSlopeShelf values
for(i in 1:nrow(traits_wide)) {
  if(traits_wide$adultSlopeShelf[i] %in% c("coastal", "inshore, intertidal, bays", "nearshore, coastal inlets & rivers")) {traits_wide$adultHorizPosition[i] <- "coastal"}
  
  else{if(traits_wide$adultSlopeShelf[i] %in% c("inshore, shallow, estuaries, rivers (summer), deep water (winter)",
                                                "ocean, coastal streams", "shelf deep edge (winter), shallow coastal water (summer)")) {traits_wide$adultHorizPosition[i] <- "oceanadromous"}
    
    else{if(traits_wide$adultSlopeShelf[i] %in% c("shallow shelf", "shelf nearshore")) {traits_wide$adultHorizPosition[i] <- "inner shelf"}
      
      else{if(traits_wide$adultSlopeShelf[i] %in% c("shelf", "shelf margin (winter), mid/outer shelf (summer)")) {traits_wide$adultHorizPosition[i] <- "shelf"}
        
        else{if(traits_wide$adultSlopeShelf[i] %in% c("outer shelf", "outer shelf, upper slope", "shelf, slope", "shelf, upper slope")) {traits_wide$adultHorizPosition[i] <- "outer shelf, upper slope"}
          
          else{if(traits_wide$adultSlopeShelf[i] %in% c("slope", "slope, shelf gullies, deep fjords")) {traits_wide$adultHorizPosition[i] <- "slope"}
            
            else{traits_wide$adultHorizPosition[i] <- NA}
        
      }}}}}}




# groups of traits that need the same treatment:

# character entries with only 1 value per species:
adultHorizPosition
adultSubstrateCategory
adultWaterColumnPosition
diet
guild
migratoryStatus


# take means of numeric entries
trNum1 <- traits_wide %>% 
  select(genus.species, common.name, gender, location, trophicPosition, K, Linfinity, 
         depthRangeShallow, depthRangeDeep, lengthConversionFLtoTL, temperatureLower, temperatureUpper,
         lengthConversionFLtoTL, lengthConversionSLtoTL) %>%
  mutate_each(funs(as.numeric),trophicPosition:lengthConversionSLtoTL) %>%
  group_by(genus.species, common.name, gender, location) %>%
  summarise_each(funs(mean(., na.rm = TRUE))) %>%
  ungroup

View(trNum1)
str(trNum1)

# take max value of ageMaximum and lengthMaximum
trNum2 <- traits_wide %>% 
  select(genus.species, common.name, gender, location, ageMaximum, lengthMaximum) %>%
  mutate_each(funs(as.numeric), ageMaximum:lengthMaximum) %>%
  group_by(genus.species, common.name, gender, location) %>%
  summarise_each(funs(max(., na.rm = TRUE))) %>%
  ungroup
View(trNum2)


unique(sort(traits_wide$firstMaturityLength))

# numeric entries which are ranges separated by hyphens
trNum3 <- traits_wide %>% 
  select(genus.species, common.name, gender, location,
         age50percentMaturity, firstMaturityAge, firstMaturityLength, length50percentMaturity) %>%
  
  # clean up age at 50% Maturity
  mutate(age50percentMaturity = gsub("\\+", "", age50percentMaturity)) %>% # remove "+"
  mutate(age50Mat1a=strsplit(age50percentMaturity,split="-") %>%
           sapply(function(x) x[1])) %>%
  mutate(age50Mat2a=strsplit(age50percentMaturity,split="-") %>%
           sapply(function(x) x[2])) %>%
  mutate(age50percentMaturity.a = as.numeric(age50percentMaturity), age50Mat1 = as.numeric(age50Mat1a), age50Mat2 = as.numeric(age50Mat2a)) %>%
  rowwise() %>%
  mutate(age50Mat = mean(c(age50percentMaturity.a, age50Mat1, age50Mat2), na.rm=T)) %>%
  ungroup() %>%
  

# clean up age at 50% Maturity
  mutate(length50percentMaturity = gsub("\\+", "", length50percentMaturity)) %>% # remove "+"
  mutate(length50Mat1a=strsplit(length50percentMaturity,split="-") %>%
           sapply(function(x) x[1])) %>%
  mutate(length50Mat2a=strsplit(length50percentMaturity,split="-") %>%
           sapply(function(x) x[2])) %>%
  mutate(length50percentMaturity.a = as.numeric(length50percentMaturity), length50Mat1 = as.numeric(length50Mat1a), length50Mat2 = as.numeric(length50Mat2a)) %>%
  rowwise() %>%
  mutate(length50Mat = mean(c(length50percentMaturity.a, length50Mat1, length50Mat2), na.rm=T)) %>%
  ungroup() %>%
  
  
  # clean up Age at First Maturity
  mutate(firstMaturityAge = gsub("\\+", "", firstMaturityAge)) %>% # remove "+"
  filter(firstMaturityAge != "up to 12") %>% # remove this entry because it's not imformative
  mutate(fMatAge1a=strsplit(firstMaturityAge,split="-") %>%
           sapply(function(x) x[1])) %>%
  mutate(fMatAge2a=strsplit(firstMaturityAge,split="-") %>%
           sapply(function(x) x[2])) %>%
  mutate(firstMaturityAge.a = as.numeric(firstMaturityAge), fMatAge1 = as.numeric(fMatAge1a), fMatAge2 = as.numeric(fMatAge2a)) %>%
  rowwise() %>%
  mutate(firstMatAge = mean(c(firstMaturityAge.a, fMatAge1, fMatAge2), na.rm=T)) %>%
  ungroup() %>%


# clean up Length at First Maturity
  mutate(firstMaturityLength = gsub("\\+", "", firstMaturityLength)) %>% # remove "+"
  mutate(fMatLength1a=strsplit(firstMaturityLength,split="-") %>%
           sapply(function(x) x[1])) %>%
  mutate(fMatLength2a=strsplit(firstMaturityLength,split="-") %>%
           sapply(function(x) x[2])) %>%
  mutate(firstMaturityLength.a = as.numeric(firstMaturityLength), fMatLength1 = as.numeric(fMatLength1a), fMatLength2 = as.numeric(fMatLength2a)) %>%
  rowwise() %>%
  mutate(firstMatLength = mean(c(firstMaturityLength.a, fMatLength1, fMatLength2), na.rm=T)) %>%
  ungroup() %>%
  
  select(genus.species, common.name, gender, location, age50Mat, length50Mat, firstMatAge, firstMatLength) %>%

# now group by ... and summarize_each(mean)  
#group_by(genus.species, common.name, gender, location) %>%

# when run all together, has all expected columns but all are short (len = 79) and only firstMatAge contains values
View(trNum3)
rm(trNum3)



    

  # sometimes the value is a range separated by hyphens, sometimes not
  # if value does not contain "-" then enter it into new column
  # if value does contain "-" then split it at "-" and enter first value into new column




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


# select just the rows I want and go from there ...
select(diet, guild, trophicPosition, lengthMaximum, ageMaximum)




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