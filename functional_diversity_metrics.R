# Functional Diversity Metrics Processing Script
# Colette Ward, Feb 29 2016

# load packages
library(httr)
library(plyr)
library(dplyr)
library(tidyr)
library(psych)
library(FD)
library(ggplot2)
library(gridExtra)


# load the functional trait data from local source in our repository:
traits_df <- read.csv("Groundfish-Functional-Diversity-Traits.csv", header=T, stringsAsFactors=FALSE)

# or load it from our google drive:
#URL_traits <- "https://drive.google.com/uc?export=download&id=0B1XbkXxdfD7uV0h5SG1UeC1lbjg"
#traitsGet <- GET(URL_traits)
#traits1 <- content(traitsGet, as='text')
#traits_df <- read.csv(file=textConnection(traits1),stringsAsFactors=FALSE)


######################################################
######################################################
######################################################


# minor dataframe cleaning:

traits_df1 <- traits_df %>%
  select(-reference, -database, -lengthType, -comments, -common.name) %>% # drop unnecessary columns
  mutate(genus.species = revalue(genus.species, c("Clupea pallasii" = "Clupea pallasi", # make Species names match those in CPUE file
                                                "Sebastes group 1" = "Dusky and Dark Rockfish", 
                                                "Sebastes group 2" = "Rougheye and Blackspotted Rockfish", 
                                                "Lepidopsetta spp." = "Lepidopsetta sp.", 
                                                "Myctophidae spp." = "Myctophidae",
                                                "Theragra chalcogramma" = "Gadus chalcogrammus"))) %>%
  mutate(genus.species = gsub(" ", ".", genus.species)) %>% # make Species names match those in CPUE file
  rename(Species = genus.species)


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
  select(-estimate, -region) %>% 
  rename(estimate = adultSubstrateCategory) 



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
  group_by(Species, trait, location) %>%
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
  group_by(Species, trait, gender, location) %>%
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
  group_by(Species, trait, gender, location) %>%
  summarise_each(funs(mean(., na.rm = TRUE))) %>% # calculate means of all estimates for each taxa / gender / location combination
  ungroup 
# warning re NAs is OK; it's just reporting that estimate2 is NA when the original cell is a single value rather than a range



###############

# pull in additional K & L_infinity values from Ben Williams:
KL_df <- read.csv("linf_k.csv", header=T, stringsAsFactors = F)
KL_df1 <- KL_df %>%
  select(-X, -common.name) %>%
  rename(Species = genus.species, estimate = value) %>%
  mutate(Species = revalue(Species, c("Clupea pallasii" = "Clupea pallasi", 
                                      "Lepidopsetta spp." = "Lepidopsetta sp.", 
                                      "Theragra chalcogramma" = "Gadus chalcogrammus",
                                      "Sebastes group 1" = "Dusky and Dark Rockfish"))) %>%
  mutate(Species = gsub(" ", ".", Species))

for(i in 1:nrow(KL_df1)) { # add columns for gender & location
  KL_df1$gender[[i]] <- "goodEnough"
  KL_df1$location[[i]] <- "goodEnough"
}


###############

# load depth coefficients from Ole ("All_sp_spatial_coef.csv")

URL_dCoef <- "https://drive.google.com/uc?export=download&id=0B1XbkXxdfD7ub1gtb09uLXFwQ2M"
dCoefGet <- GET(URL_dCoef)
dCoef1 <- content(dCoefGet, as='text')
dCoef_df <- read.csv(file=textConnection(dCoef1),stringsAsFactors=FALSE)

dCoef1 <- dCoef_df %>% 
  filter(Model == "pos") %>% # use positive model (vs binomial presence/absence)
  rename(estimate = Mean) %>%
  select(Species, estimate)

for(i in 1:nrow(dCoef1)) { # add columns for trait, gender, location
  dCoef1$trait[[i]] <- "depthCoefPos"
  dCoef1$gender[[i]] <- "both"
  dCoef1$location[[i]] <- "GoA"
}


###############

# bind these dfs together 
# This creates a dataframe of all the functional trait data we have. These are summarized values - ie means, maxima, etc
traits_df3 <- rbind(horPos_df1, substrate_df1, trChr_df1, max_df1, depth_df1, maturity_df1, KL_df1, dCoef1) #kl_df1,




######################################################
######################################################
######################################################


# This section creates a dataframe with GoA data wherever it exists; 
# where we don't have GoA data, we use data from other locations.
# we also select only female data for age / size at maturity

GoA_df <- traits_df3 %>%
  filter(location %in% c("GoA", "goodEnough")) # create table of GoA data

other_df <- traits_df3 %>%
  filter(location == "other") # create table of data for location == "other"



combos <- unique(traits_df3[,c('Species','trait')]) # create a table of all species & trait combinations for which we have data
# specify which gender we'll use for each trait:
for(i in 1:nrow(combos)){
  if(combos$trait[i] %in% c("adultSlopeShelf", "adultSubstrate", "adultWaterColumnPosition", "depthMax", "depthRange", "ageMaximum", "diet", "guild",
                            "lengthMaximum", "migratoryStatus", "trophicPosition", "depthCoefPos")) {combos$gender[i] <- "both"} # for these traits, we'll use data from males & females
  if(combos$trait[i] %in% c("age50percentMaturity", "firstMaturityAge", "firstMaturityLength", "length50percentMaturity" #"K", "Linfinity"
                            )) {combos$gender[i] <- "f"} # for life history traits, for now we'll use only female data
  if(combos$trait[i] %in% c("K", "Linfinity")) {combos$gender[i] <- "goodEnough"} 
  }
traitsGoA_df <- left_join(combos, GoA_df, by = c("Species", "trait", "gender")) # merge GoA data onto combos




# now figure out where we have trait data for location = "other" but not for GoA:
GoA_df1 <- GoA_df %>% select(-location, -estimate) # remove location & estimate columns from GoA_df to facilitate set difference operation
other_df1 <- other_df %>% select(-location, -estimate) # remove location & estimate columns from other_df to facilitate set difference operation
diffs_df <- setdiff(other_df1, GoA_df1) # retrieve taxa-gender-trait combinations which exist for location == "other" but not GoA
diffsOther_df <- left_join(diffs_df, other_df, by = c("Species", "trait", "gender")) # merge in other_df, only keeping rows that are in diffs_df



# now merge the table with non-GoA trait data onto the table with GoA data: 
traits_df4 <- left_join(traitsGoA_df, diffsOther_df, by = c("Species", "trait", "gender")) %>% 
  transmute(Species, trait, gender, 
            location = ifelse(is.na(location.x), location.y, location.x),
            estimate = ifelse(is.na(estimate.x), estimate.y, estimate.x)) %>%
  filter(!is.na(estimate)) # remove rows for which there is no data from either GoA or "other" for some of the desired trait-gender combinations


######################################################
######################################################
######################################################


# Organize functional trait data for analysis in FD package:

# Convert to wide format:
traits_df4$row <- 1:nrow(traits_df4) # create column of unique identifiers to facilitate data spread
traits_wide <- traits_df4 %>%
  spread(trait, estimate) %>%
  select(-gender, -location, -row) %>%
  group_by(Species) %>%
  summarize_each(funs(first(., order_by = is.na(.)))) %>%
  ungroup()


quant = c("age50percentMaturity", "ageMaximum", "depthCoefPos", "depthMax", "depthRange", 
         "firstMaturityAge", "firstMaturityLength", "K", "length50percentMaturity", 
         "lengthMaximum", "Linfinity", "trophicPosition")
traits_wide[,quant] <- apply(traits_wide[,quant], 2, function(x) as.numeric(x)) # convert columns to numeric as needed
traits_wide[,quant] <- apply(traits_wide[,quant], 2, function(y) log(y)) # log-transform the quantitative traits (as per recommendation in Botta-Dukat 2005)
head(traits_wide)

categ = c("adultSlopeShelf", "adultSubstrate", "adultWaterColumnPosition", "diet", "guild", "migratoryStatus") # convert columns to factors as needed
traits_wide[,categ] <- lapply(traits_wide[,categ] , factor)
#write.csv(traits_wide, file = "traits_wide.csv")



# which pairs of (log-transformed) quantitative traits are correlated?
pairs.panels(traits_wide[,quant],smooth=F,density=T,ellipses=F,lm=T,digits=3,scale=T)
names(traits_wide[,quant])
# significant correlations:
# K, firstMaturityAge, age50percentMaturity, ageMaximum
# Linfinity, trophicPosition, firstMaturityLength, length50percentMaturity, lengthMaximum
# depthRange & depthMax; note that depthCoefficient is not correlated with either



# we have the most data for:
# adultWaterColumnPosition (missing for 2 taxa)
#ageMax (missing for 9 taxa)
#diet (missing 0)
#sum(is.na(traits_wide$firstMaturityLength)) # missing 23
#sum(is.na(traits_wide$K)) # missing 16
#sum(is.na(traits_wide$Linfinity)) # missing 16
#guild (missing 0)
#lengthMaximum (missing 4)
#sum(is.na(traits_wide$ageMaximum)) # missing 10
#trophicPosition (missing 0)
#depthRange (missing 6)
#depthMax (missing 6)
#depthCoefPos (missing 0)


# Prep functional traits df to load into functional diversity analysis
ft_df <- traits_wide %>%
  select(Species, lengthMaximum, ageMaximum, depthMax, depthCoefPos) %>% # select log-transformed quantitative traits
  #select(Species, lengthMaximum, ageMaximum, depthMax, depthCoefPos, diet) %>% # select log-transformed quantitative traits & categorical diet trait
  filter(!(is.na(ageMaximum)), !(is.na(depthMax)), !is.na(depthCoefPos)) %>% # remove taxa for which trait data are missing
  arrange(Species)
rownames(ft_df) <- ft_df$Species # create row names from Species column
ft_df <- ft_df %>% select(-Species)
# View(ft_df)  # this is the dataframe we'll use for functional diversity analyses

#unique(sort(setdiff(ft_df$Species, SPCPUEArea$Species))) # sp in traits_df3 but not SPCPUEArea


######################################################
######################################################
######################################################


# load Mean annual CPUE data for Shallow Areas:
URL_SPCPUEArea <- "https://drive.google.com/uc?export=download&id=0By1iaulIAI-uYzBOUFRtZklmX0U" # new data for shallow areas
SPCPUEArea_Get <- GET(URL_SPCPUEArea)
SPCPUEArea_1 <- content(SPCPUEArea_Get, as='text')
SPCPUEArea <- read.csv(file=textConnection(SPCPUEArea_1),stringsAsFactors=FALSE,head=TRUE)
#View(SPCPUEArea)


# NB  SPCPUEArea and deepCPUE both have only 53 taxa, not 57. Which ones are missing?
# spDiffs <- setdiff(traits_wide$Species, deepCPUE$Species); spDiffs
# "Berryteuthis.magister"   "Hydrolagus.colliei"      "Merluccius.productus"    "Sebastes.helvomaculatus"


# organize Shallow Areas CPUE data for analysis in FD package:
# unique(sort(setdiff(SPCPUEArea$Species, ft_df$Species)))
spShallow_df <- SPCPUEArea %>%
  select(area, year, Species, Mean.totalDensity) %>%
  filter(!(Species %in% c("Chionoecetes.bairdi", "Hemitripterus.bolini", "Hyas.lyratus", "Lycodes.brevipes", 
                          "Lycodes.palearis", "Lyopsetta.exilis", "Myctophidae", "Oncorhynchus.keta", 
                          "Oncorhynchus.tshawytscha"))) %>% # remove taxa for which we don't have all trait data
  mutate(area = revalue(area, c("Total" = "10")), # recode Total for looping later
         area = as.numeric(area)) # convert to numeric class

Sh <- spShallow_df %>% 
  select(Species, area, year, Mean.totalDensity) %>%
  arrange(Species) %>% # arrange in alphabetical order to match order in functional traits df
  spread(Species, Mean.totalDensity) %>%
  select(-year)

shallowByArea_list <- split(Sh, f = Sh$area) # create a list of dataframes (one for each area; NB area 10 is Total)

shallowByArea_list1 <- lapply(shallowByArea_list, function(x) x[!(names(x) %in% c("area", "year"))]) # drop area & year



#########################


# load mean annual CPUE for Deep areas:
URL_deepCPUE <- "https://drive.google.com/uc?export=download&id=0By1iaulIAI-uVF9VWnNPX3Z3S3c"
deepCPUE_Get <- GET(URL_deepCPUE)
deepCPUE_1 <- content(deepCPUE_Get, as='text')
deepCPUE <- read.csv(file=textConnection(deepCPUE_1),stringsAsFactors=FALSE,head=TRUE)
#View(deepCPUE)


# organize Deep Areas CPUE data for analysis in FD package:
# unique(sort(setdiff(deepCPUE$Species, ft_df$Species)))
spDeep_df <- deepCPUE %>%
  select(area, year, Species, Mean.totalDensity) %>%
  filter(!(Species %in% c("Chionoecetes.bairdi", "Hemitripterus.bolini", "Hyas.lyratus", "Lycodes.brevipes", 
                          "Lycodes.palearis", "Lyopsetta.exilis", "Myctophidae", "Oncorhynchus.keta", 
                          "Oncorhynchus.tshawytscha"))) %>% # remove taxa for which we don't have all trait data
  mutate(area = revalue(area, c("Total" = "6")), # recode Total for looping later
         area = as.numeric(area)) # convert to numeric class

Dp <- spDeep_df %>% 
  select(Species, area, year, Mean.totalDensity) %>%
  arrange(Species) %>% # arrange in alphabetical order to match order in functional traits df
  spread(Species, Mean.totalDensity) %>%
  select(-year)

deepByArea_list <- split(Dp, f = Dp$area) # create a list of dataframes (one for each area; NB area 6 is Total)

deepByArea_list1 <- lapply(deepByArea_list, function(x) x[!(names(x) %in% c("area", "year"))]) # drop area & year


######################################################
######################################################
######################################################


# Calculate Functional Diversity (Rao's Q)


# Methods info from FD package documentation:
# "If not all traits are numeric, Gower’s (1971) standardization by the
#  range is automatically used; see gowdis for more details."

# "If x (trait df) is a matrix or a data frame that contains only continuous traits, no NAs, and no weights
# are specified (i.e. w is missing), a species-species Euclidean distance matrix is computed via dist.
# Otherwise, a Gower dissimilarity matrix is computed via gowdis."

# "Rao’s quadratic entropy (Q) is computed from the uncorrected species-species distance matrix via divc."

# gowdis computes the Gower (1971) similarity coefficient as described by Podani (1999),
# then converts it to a dissimilarity coefficient by using D = 1-S



# Rao's Q for Shallow Areas:
fdShallow <- list()
for (i in seq_along(shallowByArea_list1)) {
  fdShallow[[i]] <- dbFD(ft_df, shallowByArea_list1[[i]], calc.FRic = F, calc.CWM = F, calc.FDiv = F) # nb adding stand.x = T when using only quantitative traits makes no difference becasue data were already log-transformed
}
# when quantitative & categorical traits are used:
# "Species x species distance matrix was not Euclidean. 'sqrt' correction was applied."



# Rao's Q for Deep Areas:
fdDeep <- list()
for (i in seq_along(deepByArea_list1)) {
  fdDeep[[i]] <- dbFD(ft_df, deepByArea_list1[[i]], calc.FRic = F, calc.CWM = F, calc.FDiv = F)
}
# when quantitative & categorical traits are used:
# "Species x species distance matrix was not Euclidean. 'sqrt' correction was applied."


#############################

# Create dataframes of Rao's Q for plotting purposes:

# Shallow areas:
colsShallow <- paste("Area", 1:10, sep="")
shallowRao1 <- data.frame(matrix(NA_real_, nrow = 14, ncol = 10)); colnames(shallowRao1) <- colsShallow

for (i in seq_along(fdShallow)) {
  shallowRao1[,i] <- data.frame(as.data.frame(fdShallow[[i]]$RaoQ))
}

year <- as.data.frame(unique(sort(SPCPUEArea$year))); colnames(year) <- "year"
shallowRaoQ <- shallowRao1 %>% bind_cols(year) %>% 
  select(-Area10) # remove column for total areas combined




# Deep areas:
colsDeep <- paste("Area", 1:6, sep="")
deepRao1 <- data.frame(matrix(NA_real_, nrow = 14, ncol = 6)); colnames(deepRao1) <- colsDeep

for (i in seq_along(fdDeep)) {
  deepRao1[,i] <- data.frame(as.data.frame(fdDeep[[i]]$RaoQ))
}

year2 <- as.data.frame(unique(sort(deepCPUE$year))); colnames(year) <- "year"
deepRaoQ <- bind_cols(year, deepRao1) %>% 
  select(-Area6) # remove column for total areas combined

######################################################
######################################################
######################################################


# Temporal plots of Rao's Q

# 1. Shallow Areas:
year3 <- unique(sort(SPCPUEArea$year))

ggplot(data=shallowRaoQ, aes(x=year3, y = value)) + 
  geom_point(aes(y = Area1), size=2) +         geom_line(aes(y = Area1), size=2) +
  geom_point(aes(y = Area2), size=2) +         geom_line(aes(y = Area2), size=2) +
  geom_point(aes(y = Area3), size=2, col=2) +  geom_line(aes(y = Area3), size=2, col=2) +
  geom_point(aes(y = Area4), size=2, col=2) +  geom_line(aes(y = Area4), size=2, col=2) +
  geom_point(aes(y = Area5), size=2, col=2) +  geom_line(aes(y = Area5), size=2, col=2) +
  geom_point(aes(y = Area6), size=2) +         geom_line(aes(y = Area6), size=2) +
  geom_point(aes(y = Area7), size=2) +         geom_line(aes(y = Area7), size=2) +
  geom_point(aes(y = Area8), size=2) +         geom_line(aes(y = Area8), size=2) +
  geom_point(aes(y = Area9), size=2) +         geom_line(aes(y = Area9), size=2) +
  
  theme(axis.line=element_line('black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  theme(axis.text.x = element_text(angle=90, size=18, colour = "black"))+
  theme(axis.text.y = element_text(size=22))+
  scale_x_continuous(breaks=c(year3), labels=c(year3)) +
  labs(title = "Rao's Q for Shallow Areas", 
       x = "Year", y = "Rao's Q")



#########################

# 2. Deep Areas:
year4 <- unique(sort(deepCPUE$year))

ggplot(data=deepRaoQ, aes(x=year4, y = value)) + 
  geom_point(aes(y = Area1), size=2, col=2) +  geom_line(aes(y = Area1), size=2, col=2) +
  geom_point(aes(y = Area2), size=2, col=2) +  geom_line(aes(y = Area2), size=2, col=2) +
  geom_point(aes(y = Area3), size=2) +         geom_line(aes(y = Area3), size=2) +
  geom_point(aes(y = Area4), size=2) +         geom_line(aes(y = Area4), size=2) +
  geom_point(aes(y = Area5), size=2) +         geom_line(aes(y = Area5), size=2) +
  
  theme(axis.line=element_line('black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  theme(axis.text.x = element_text(angle=90, size=18, colour = "black"))+
  theme(axis.text.y = element_text(size=22))+
  scale_x_continuous(breaks=c(year4), labels=c(year4)) +
  labs(title = "Rao's Q for Deep Areas", 
       x = "Year", y = "Rao's Q")



######################################################
######################################################
######################################################


# Spatial boxplots of Rao's Q


# Create long-form dataframes of Rao's Q for Spatial Plots:
year <- as.data.frame(unique(sort(SPCPUEArea$year)))

shallowRao2 <- shallowRaoQ %>% 
  gather(key = area, value = RaosQ, Area1:Area9) %>%
  mutate(area = gsub("Area", "", area), 
         area = as.factor(area))

deepRao2 <- deepRaoQ %>% 
  gather(key = area, value = RaosQ, Area1:Area5) %>%
  mutate(area = gsub("Area", "", area), 
         area = as.factor(area))

#########################


# Load Mary & Rachael's boxplot theme:
theme_boxplot <- function(base_size = 12){
  theme_bw(base_size)%+replace%
    theme(legend.key.size=unit(12,"points"),
          legend.text=element_text(size=11),
          legend.key=element_blank(),
          legend.title=element_blank(),
          legend.background=element_rect(colour="white", fill="transparent"),
          plot.margin=unit(c(0.5,1,0.5,1), "lines"),  
          panel.border=element_blank(),
          panel.margin=unit(0,"lines"),
          panel.background=element_rect(fill=NA, colour=NA),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.ticks.length=unit(1,"mm"),
          axis.text.x = element_text(margin=margin(5,0,0,0), size=12), 
          axis.text.y = element_text(margin=margin(0,5,0,0), size=12),
          axis.title.x=element_text(size=12, margin=margin(15,0,0,0)),
          axis.title.y=element_text(size=12, angle=90, margin=margin(0,15,0,0)),
          strip.text.x=element_text(size=12),
          strip.background=element_rect(colour="black", fill='white'))
}


#########################


shallow_FD <- ggplot(data=shallowRao2, aes(x = area, y = RaosQ)) + 
  geom_boxplot() + theme_boxplot() +
  xlab("Area (West <-> East)") + ylab("Rao's Q") +
  xlim("9", "8", "7", "6", "5", "4", "3", "2", "1") + ylim(0.95, 2.6) +
  #theme_classic() +
  #theme(axis.text.x = element_text(size=15))
  theme(plot.background=element_blank(),
        axis.text.x = element_text(size=15))


deep_FD <- ggplot(data=deepRao2, aes(x = area, y = RaosQ)) + 
  geom_boxplot() + theme_boxplot() +
  xlab("Area (West <-> East)") + ylab("Rao's Q") +
  xlim("5", "4", "3", "2", "1") + ylim(0.95, 2.6) +
  #theme_classic() +
  #theme(axis.text.x = element_text(size=15))
  theme(plot.background=element_blank(),
        axis.text.x = element_text(size=15))


grid.arrange(shallow_FD, deep_FD, ncol=2)
