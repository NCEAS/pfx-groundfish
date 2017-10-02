##########################################################################
###  Testing statistical diff. in diversity metrics between study sites
###   Rachael E. Blake
###  August 9, 2017
##########################################################################

# load packages
library(tidyverse)

# load the files
shallowDat <- read.csv("./diversity-data/DiversityMetrics.Shallow.Bootstrapped.final.csv")
shallowDat$AREA_fac <- as.factor(shallowDat$AREA)


source("./diversity_scripts/functional_diversity_metrics.R") # shallowRao2 is the df I want.
shallowRao2$AREA_fac <- as.factor(shallowRao2$AREA)



# test differences between study sites in species richness, alpha diversity, saturation
sprich <- aov(Sp_rich~AREA_fac, data=shallowDat)
anova(sprich)
TukeyHSD(sprich, ordered=TRUE)

alpha <- aov(Eff_Num_Sp~AREA_fac, data=shallowDat)
anova(alpha)
TukeyHSD(alpha, ordered=TRUE)

beta <- aov(Exp_B_Div~AREA_fac, data=shallowDat)
anova(beta)
TukeyHSD(beta, ordered=TRUE)

RaoQ <- aov(RaosQ~AREA_fac, data=shallowRao2)
anova(RaoQ)
TukeyHSD(RaoQ, ordered=TRUE)


#####
# add column to shallowDat with one dummy value for areas 1-7, another dummy value for areas 8-10 
shallowDat2 <- shallowDat %>%
               mutate(AREA_group = ifelse(AREA %in% c(1:7), "low", "high"))

shallowRao22 <- shallowRao2 %>%
                mutate(AREA_group = ifelse(AREA %in% c(1:7), "low", "high"))

# test differences again, but using dummy area variable
sprich2 <- aov(Sp_rich~AREA_group, data=shallowDat2)
anova(sprich2)
TukeyHSD(sprich2, ordered=TRUE)

alpha2 <- aov(Eff_Num_Sp~AREA_group, data=shallowDat2)
anova(alpha2)
TukeyHSD(alpha2, ordered=TRUE)

beta2 <- aov(Exp_B_Div~AREA_group, data=shallowDat2)
anova(beta2)
TukeyHSD(beta2, ordered=TRUE)

RaoQ2 <- aov(RaosQ~AREA_group, data=shallowRao22)
anova(RaoQ2)
TukeyHSD(RaoQ2, ordered=TRUE)








