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



#####
# average the diversity metrics to get mean per area per year so that there's the same number of samples
# as for Rao's Q data

shallowDat3 <- shallowDat2 %>% 
               select(-Gamma, -SW, -Simp, -Exp_G_DivBoot, -AREA_group, -AREA_fac) %>% 
               group_by(AREA, YEAR) %>% 
               summarize_each(funs = "mean")  %>%
               mutate(AREA_group = ifelse(AREA %in% c(1:7), "low", "high"),
                      AREA_fac = as.factor(AREA)) %>% 
               ungroup() %>% 
               group_by(AREA_group) %>% 
               mutate(sp_rich_group_mn = mean(Sp_rich)) %>% 
               ungroup()


# test differences again, but using dummy area variable
sprich3 <- aov(Sp_rich~AREA_group, data=shallowDat3)
anova(sprich3)
TukeyHSD(sprich3, ordered=TRUE)

alpha3 <- aov(Eff_Num_Sp~AREA_group, data=shallowDat3)
anova(alpha3)
TukeyHSD(alpha3, ordered=TRUE)

beta3 <- aov(Exp_B_Div~AREA_group, data=shallowDat3)
anova(beta3)
TukeyHSD(beta3, ordered=TRUE)


# test differences again, but using dummy area variable
sprich4 <- aov(Sp_rich~AREA_fac, data=shallowDat3)
anova(sprich4)
TukeyHSD(sprich4, ordered=TRUE)

alpha4 <- aov(Eff_Num_Sp~AREA_fac, data=shallowDat3)
anova(alpha4)
TukeyHSD(alpha4, ordered=TRUE)

beta4 <- aov(Exp_B_Div~AREA_fac, data=shallowDat3)
anova(beta4)
TukeyHSD(beta4, ordered=TRUE)



