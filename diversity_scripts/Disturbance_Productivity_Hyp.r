##################################################################################
#####  Disturbance & Productivity script for PFX Groundfish Diversity paper  #####
#####  Script by Rachael E. Blake, created on Sept 12, 2017 from the         #####
#####  RMarkdown file of the same name.  This was created do the code can    #####
#####  be called in another RMarkdown file.                                  #####
##################################################################################

# Load packages
library(vegan) ; library(mvnormtest) ; library(plyr)
library(tidyverse) ; library(forcats) ; library(gridExtra)
library(ggplot2) ; library(psych) ; library(mgcv)
library(AICcmodavg) ; library(visreg) ; library(egg)


# source the diversity metrics csv file
div_met <- read.csv("../diversity-data/DiversityMetrics.Shallow.Bootstrapped.Final.csv")
#head(div_met)


# source the fishing pressure calculation script
source("../diversity_scripts/Fishing_Pressure_calc.r")

# want the catch_sum dataframe from here
#head(catch_sum)


# source the chlorophyll a calculation script
source("../diversity_scripts/Primary_Productivity_calc.r")

# want the CHL_both dataframe from here
#head(CHL_both)


# source the functional diversity calculation script
source("../diversity_scripts/functional_diversity_metrics.R")

# want the shallowRao2 dataframe from here
#head(shallowRao2)


#######################################
### MAKE MY OWN THEME TO SAVE LINES OF CODE
theme_boxplot <- function(base_size = 12){
  theme_bw(base_size) %+replace%
    theme(legend.key.size=unit(13,"points"),
          legend.text=element_text(size=I(11)),
          legend.key=element_blank(),
          #legend.title=element_blank(),
          #legend.position="none",
          plot.margin=unit(c(0.75,1,0.75,1), "cm"), # respectively: top, right, bottom, left; refers to margin *outside* labels; default is c(1,1,0.5,0.5)
          panel.border=element_blank(),
          panel.spacing=unit(0,"lines"),
          axis.ticks.length=unit(1,"mm"),
          axis.text.x = element_text(margin=margin(5,0,0,0)),
          axis.text.y = element_text(margin=margin(0,5,0,0)),
          axis.text=element_text(size=11),
          axis.title.x=element_text(size=13, margin=margin(15,0,0,0)), 
          axis.title.y=element_text(size=13, angle=90, margin=margin(0,15,0,0)), 
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          strip.text.x=element_text(size=14),
          strip.background=element_rect(colour='black', fill='white'),
          #axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line = element_line(colour = 'black', size=0.5, linetype='solid')
    )
}



# Fishing Pressure effects on Diversity Metrics

# clean the data
catch_sum2 <- catch_sum %>%
              dplyr::rename(Catch_area = final.OLE,
                            Year = year) %>%
              dplyr::select(Catch_area, Year, met_ton_km2) %>%
              arrange(Catch_area, Year)

div_met2 <- div_met %>%
            dplyr::rename(Year = YEAR) %>%
            mutate(Catch_area = plyr::mapvalues(AREA, c("1","2","3","4","5","6","7","8","9","10"), 
                                                c("Prince William Sound","Prince William Sound",
                                                  "Prince William Sound","Cook Inlet","Cook Inlet",
                                                  "Kodiak","Kodiak","Kodiak","Kodiak","Alaska Peninsula")),
            AREA = as.factor(AREA)) %>%
            dplyr::arrange(AREA, Year)

# make one large dataframe
dist_df <- merge(div_met2, catch_sum2, all.x=T)

# filter for complete cases
dist_df <- dist_df %>% full_join(shallowRao2, by=c("Year", "AREA")) %>% filter(complete.cases(.))  

# average to catch area for each year...
dist_df_cat_ar <- dist_df %>%
                  dplyr::select(-AREA, -Simp, -Invsimp) %>%
                  group_by(Catch_area, Year) %>%
                  summarize(Eff_Num_Sp_ca = mean(Eff_Num_Sp),
                            Sp_rich_ca = mean(Sp_rich),
                            Exp_B_Div_ca = mean(Exp_B_Div),
                            met_ton_km2_ca = mean(met_ton_km2),
                            RaosQ = mean(RaosQ)) %>%
                  ungroup() %>%
                  arrange(Catch_area, Year)




# color for plotting
colors_plot <- c("#49006A","#DD3497","#FA9FB5","#E3C9C6")



###with Linear fits

# Species Richness
#pairs.panels(dist_df[,-c(1:3)],smooth=F,density=T,ellipses=F,lm=T,digits=3,scale=T)

lmfit1 <- lm(Sp_rich_ca~met_ton_km2_ca, data=dist_df_cat_ar)  # lm() fit 

div_plot1 <- ggplot(dist_df_cat_ar, aes(x=met_ton_km2_ca, y=Sp_rich_ca, color=Catch_area)) + 
             geom_point(size=6) + theme_boxplot() + 
             stat_smooth(method="lm", col="blue") + ylab("Species Richness") + 
             xlab(expression(paste("Fishing Pressure (metric tonnes ", ~km^-2,")"))) + 
             theme(panel.background = element_blank(), legend.position = "none") +
             scale_color_manual(values=colors_plot) +
             labs(title=paste("Adj R^2 =",signif(summary(lmfit1)$adj.r.squared, 5),
                              #"Intercept =",signif(lmfit1$coef[[1]],5 ),
                              #" Slope =",signif(lmfit1$coef[[2]], 5),
                              " p =",signif(summary(lmfit1)$coef[2,4], 5)))

anova(lmfit1) # anova table of above lm() fit


# Alpha Diversity

lmfit2 <- lm(Eff_Num_Sp_ca~met_ton_km2_ca, data=dist_df_cat_ar)  # lm() fit 

div_plot2 <- ggplot(dist_df_cat_ar, aes(x=met_ton_km2_ca, y=Eff_Num_Sp_ca, color=Catch_area)) +
             geom_point(size=6) + theme_boxplot() +  
             stat_smooth(method="lm", col="blue") + ylab("Alpha Diversity") +
             xlab(expression(paste("Fishing Pressure (metric tonnes ", ~km^-2,")"))) + 
             theme(panel.background = element_blank(), legend.position = "none") +  
             scale_color_manual(values=colors_plot) +
             labs(title=paste("Adj R^2 =",signif(summary(lmfit2)$adj.r.squared, 5),
                              #"Intercept =",signif(lmfit2$coef[[1]],5 ),
                              #" Slope =",signif(lmfit2$coef[[2]], 5),
                              " p =",signif(summary(lmfit2)$coef[2,4], 5)))

anova(lmfit2) # anova table of above lm() fit


# Beta Diversity

lmfit3 <- lm(Exp_B_Div_ca~met_ton_km2_ca, data=dist_df_cat_ar)  # lm() fit 


div_plot3 <- ggplot(dist_df_cat_ar, aes(x=met_ton_km2_ca, y=Exp_B_Div_ca, color=Catch_area)) + 
             geom_point(size=6) + theme_boxplot() + scale_color_manual(values=colors_plot) +
             stat_smooth(method="lm", col="blue", size=2) + ylab("Beta Diversity ") + 
             xlab(expression(paste("Fishing Pressure (metric tonnes ", ~km^-2,")"))) +  
             theme(#legend.position = c(1.1,1.2), legend.justification = c(1,1),
                   legend.position = "none",
                   panel.background = element_blank()) + 
             labs(title=paste("Adj R^2 =", signif(summary(lmfit3)$adj.r.squared, 5),
                              #"Intercept =",signif(lmfit3$coef[[1]],5 ),
                              #" Slope =",signif(lmfit3$coef[[2]], 5),
                              " p =", signif(summary(lmfit3)$coef[2,4], 5)))

anova(lmfit3) # anova table of above lm() fit


# Functional Diversity     #shallowRao2

lmfit3a <- lm(RaosQ~met_ton_km2_ca, data=dist_df_cat_ar)  # lm() fit 

div_plot3a <- ggplot(dist_df_cat_ar, aes(x=met_ton_km2_ca, y=RaosQ, color=Catch_area)) + 
              geom_point(size=6) + theme_boxplot() + scale_color_manual(values=colors_plot) +
              stat_smooth(method="lm", col="blue", size=2) + ylab("Functional Diversity ") + 
              xlab(expression(paste("Fishing Pressure (metric tonnes ", ~km^-2,")"))) +  
              theme(legend.position = "none", panel.background = element_blank()) + 
              labs(title=paste("Adj R^2 =", signif(summary(lmfit3a)$adj.r.squared, 5),
                               #"Intercept =",signif(lmfit3a$coef[[1]],5 ),
                               #" Slope =",signif(lmfit3a$coef[[2]], 5),
                               " p =", signif(summary(lmfit3a)$coef[2,4], 5))) 

anova(lmfit3a) # anova table of above lm() fit


# arranging above plots
ggarrange(div_plot1, div_plot2, div_plot3, div_plot3a, ncol=2, nrow=2)



###With GAM fits

# Species Richness 
#assumes normally distributed data 
gam1 <- gam(Sp_rich~met_ton_km2, data=dist_df, family=gaussian) # essentially equal to glm()

#gam1a <- gam(Sp_rich~s(met_ton_km2, k=1), data=dist_df, family=gaussian) 

#gam1b <- gam(Sp_rich~s(met_ton_km2, k=2), data=dist_df, family=gaussian) 

gam1c <- mgcv::gam(Sp_rich~s(met_ton_km2, k=4), data=dist_df, family=gaussian) 

summary(gam1c) #; anova(gam1c)

#plot(gam1c, pages=1, residuals=TRUE, shade=T)

#gam.check(gam1c)

#AIC(gam1, gam1a, gam1b, gam1c)

# ggplot
gam1_plot <- ggplot(data=dist_df_cat_ar, aes(x=met_ton_km2_ca, y=Sp_rich_ca, color=Catch_area)) + 
             geom_point(size=6) + ylab("Species Richness") + theme_boxplot() +
             stat_smooth(method="gam", formula=y ~ s(x, k=4), col="red", size=2) + 
             xlab(expression(paste("Fishing Pressure (metric tonnes ", ~km^-2,")"))) + 
             scale_color_manual(values=colors_plot) + 
             theme(legend.position = "none", panel.background = element_blank()) +
             labs(title=paste("Adj R2 =", signif(summary(gam1c)$r.sq, 4),
                              "p =", signif(summary(gam1c)$s.table[1,4], 3)))
#gam1_plot


# Alpha Diversity
gam2 <- gam(Eff_Num_Sp~met_ton_km2, data=dist_df, family=gaussian)  # essentially equal to glm()

#gam2a <- gam(Eff_Num_Sp~s(met_ton_km2, k=1), data=dist_df, family=gaussian)

#gam2b <- gam(Eff_Num_Sp~s(met_ton_km2, k=2), data=dist_df, family=gaussian)

gam2c <- gam(Eff_Num_Sp~s(met_ton_km2, k=4), data=dist_df, family=gaussian)

#plot(gam2b, pages=1, residuals=TRUE, shade=T)

summary(gam2c) #; anova(gam2c) 
#AIC(gam2, gam2a, gam2b, gam2c)

# ggplot
# use plot = FALSE to get plot data from visreg without plotting
# plotdata2 <- visreg(gam2b, type="contrast", plot=FALSE)
# 
# # The output from visreg is a list of the same length as the number of 'x' variables,
# #   so we use ldply to pick the objects we want from the each list part and make a dataframe: 
# smooths2 <- data.frame(Variable = plotdata2$meta$x, 
#                        x=plotdata2$fit[[plotdata2$meta$x]], 
#                        smooth=plotdata2$fit$visregFit, 
#                        lower=plotdata2$fit$visregLwr, 
#                        upper=plotdata2$fit$visregUpr)
# 
# resid2 <- data.frame(Variable=plotdata2$meta$x,
#                      x=plotdata2$res[[plotdata2$meta$x]],
#                      y=plotdata2$res$visregRes)
# 
# # The ggplot:
# gam2_plota <- ggplot(smooths2, aes(x, smooth)) + 
#               geom_line() + theme_boxplot() +
#               geom_line(aes(y=lower), linetype="dashed") + 
#               geom_line(aes(y=upper), linetype="dashed") #+
#               #geom_point(data=resid2, aes(x, y)) +

#gam2_plota


# ggplot
gam2_plot <- ggplot(data=dist_df_cat_ar, aes(x=met_ton_km2_ca, y=Eff_Num_Sp_ca, color=Catch_area)) + 
             geom_point(size=6) +  ylab("Alpha Diversity") + 
             xlab(expression(paste("Fishing Pressure (metric tonnes ", ~km^-2,")"))) + 
             stat_smooth(method="gam", formula=y ~ s(x, k=4), col="red", size=2) + theme_boxplot() +
             scale_color_manual(values=colors_plot) + 
             theme(legend.position = c(0.87,0.95), panel.background = element_blank(),
                   legend.background = element_rect(fill="transparent")) +
             labs(title=paste("Adj R2 =", signif(summary(gam2c)$r.sq, 4),
                              "p =", signif(summary(gam2c)$s.table[1,4], 3)))

#gam2_plot



# Beta Diversity
gam3 <- gam(Exp_B_Div~met_ton_km2, data=dist_df, family=gaussian)  # essentially equal to glm()

#gam3a <- gam(Exp_B_Div~s(met_ton_km2, k=1), data=dist_df, family=gaussian)

#gam3b <- gam(Exp_B_Div~s(met_ton_km2, k=2), data=dist_df, family=gaussian)

gam3c <- gam(Exp_B_Div~s(met_ton_km2, k=4), data=dist_df, family=gaussian)

#plot(gam3b, pages=1, shade=T)

summary(gam3c) #; anova(gam3c)
#AIC(gam3, gam3a, gam3b, gam3c, gam3d)

# ggplot
gam3_plot <- ggplot(data=dist_df_cat_ar, aes(x=met_ton_km2_ca, y=Exp_B_Div_ca, color=Catch_area)) + 
             geom_point(size=6) +  ylab("Beta Diversity") + 
             xlab(expression(paste("Fishing Pressure (metric tonnes ", ~km^-2,")"))) + 
             stat_smooth(method="gam", formula=y ~ s(x, k=4), col="red", size=2) + theme_boxplot() +
             scale_color_manual(values=colors_plot) + 
             theme(legend.position = "none", panel.background = element_blank()) +
             labs(title=paste("Adj R2 =", signif(summary(gam3c)$r.sq, 4),
                              "p =", signif(summary(gam3c)$s.table[1,4], 3)))

#gam3_plot



# Functional Diversity
gam3a <- gam(RaosQ~met_ton_km2, data=dist_df, family=gaussian)  # essentially equal to glm()

#gam3aa <- gam(RaosQ~s(met_ton_km2, k=1), data=dist_df, family=gaussian)

#gam3bb <- gam(RaosQ~s(met_ton_km2, k=2), data=dist_df, family=gaussian)

gam3cc <- gam(RaosQ~s(met_ton_km2, k=4), data=dist_df, family=gaussian)

#plot(gam3bb, pages=1, shade=T)

summary(gam3cc) #; anova(gam3cc)
#AIC(gam3a, gam3aa, gam3bb, gam3cc)

# ggplot
gam3a_plot <- ggplot(data=dist_df_cat_ar, aes(x=met_ton_km2_ca, y=RaosQ, color=Catch_area)) + 
              geom_point(size=6) +  ylab("Functional Diversity") + 
              xlab(expression(paste("Fishing Pressure (metric tonnes ", ~km^-2,")"))) +  
              stat_smooth(method="gam", formula=y ~ s(x, k=4), col="red", size=2) + theme_boxplot() + 
              scale_color_manual(values=colors_plot) + 
              theme(legend.position = "none", panel.background = element_blank()) +
              labs(title=paste("Adj R2 =", signif(summary(gam3cc)$r.sq, 4),
                               "p =", signif(summary(gam3cc)$s.table[1,4], 5)))

#gam3a_plot



# arranging above plots
#grid.arrange(gam1_plot, gam2_plot, gam3_plot, gam3a_plot, ncol=2, nrow=2)

ggarrange(gam1_plot, gam2_plot, gam3_plot, gam3a_plot, ncol=2, nrow=2)




# arranging plots
#grid.arrange(gam1_plot, gam2_plot, div_plot3, gam3a_plot, ncol=2, nrow=2)




# Primary Productivity effects 

# # clean and prepare the data
# head(CHL_both)  # size is 19.6 Mb
# 
# # size is 7.8 Mb
# CHL_both_s <- CHL_both %>% dplyr::select(-Day, -latitude_degrees_north, -longitude_degrees_east) 
# 
# div_met3 <- div_met %>%    # size is 5.3 Mb
#             rename(Year = YEAR,
#                    area = AREA) %>%
#             dplyr::select(-Simp, -Invsimp, -SW, -Gamma) 
# 
# # make one large dataframe
# prod_df <- full_join(div_met3, CHL_both_s, by=c("Year","area"))   
# 
# # filter for complete cases
# prod_df <- prod_df %>% filter(complete.cases(.))   # size is 13.4 Gb
# 
# # average to chla for each year...
# prod_df_chlyr <- prod_df %>%                     # size is 4 Kb
#                  group_by(area, Year) %>%
#                  summarize(Eff_Num_Sp_mn = mean(Eff_Num_Sp),
#                            Sp_rich_mn = mean(Sp_rich),
#                            Exp_B_Div_mn = mean(Exp_B_Div),
#                            chl_mg_m3_mn = mean(chlorophyll_mg_m3)) %>%
#                  ungroup() %>%
#                  arrange(Year, area)

# average to chla for each month...
# prod_df_chlmth <- prod_df %>%
#                   group_by(area, Year, Month) %>%
#                   summarize(Eff_Num_Sp_mn = mean(Eff_Num_Sp),
#                             Sp_rich_mn = mean(Sp_rich),
#                             Exp_B_Div_mn = mean(Exp_B_Div),
#                             chl_mg_m3_mn = mean(chlorophyll_mg_m3)) %>%
#                   ungroup() %>%
#                   arrange(Year, Month, area)


# write.csv(prod_df_chlyr, "./diversity-data/Productivity_gfish_areas_yr.csv", row.names=F)
# write.csv(prod_df_chlmth, "./diversity-data/Productivity_gfish_areas_mth.csv", row.names=F)

prod_chl_yr <- read.csv("../diversity-data/Productivity_gfish_areas_yr.csv")
prod_chl_mth <- read.csv("../diversity-data/Productivity_gfish_areas_mth.csv")

# add in the functional diversity data
prod_chl_yr <- prod_chl_yr %>%
               dplyr::rename(AREA=area) %>%
               dplyr::mutate(AREA=as.factor(AREA)) %>%
               full_join(shallowRao2, by=c("AREA", "Year")) %>%
               filter(!is.na(Eff_Num_Sp_mn))

# Define order of the boxes (use this as "fill=" below)
prod_chl_yr$AREA1 <- factor(prod_chl_yr$AREA, levels=c("1","2","3","4","5","6","7","8","9","10"))

# colors for plots
pinks <- c("#FFE6DA", "#E3C9C6", "#FCC5C0", "#FA9FB5", "#F768A1",
           "#E7298A", "#DD3497", "#AE017E", "#7A0177",  "#49006A") 



### with Linear Fits
# Linear fits on Species Richness

lmfit4 <- lm(Sp_rich_mn~chl_mg_m3_mn, data=prod_chl_yr)  # lm() fit 

div_plot2c <- ggplot(prod_chl_yr, aes(x=chl_mg_m3_mn, y=Sp_rich_mn, color=as.factor(AREA1))) + 
              geom_point(size=4) + scale_color_manual(values=pinks) + theme_boxplot() +
              stat_smooth(method="lm", col="blue", size=2) + ylab("Species Richness") +
              xlab(expression(paste("Primary Productivity (chla mg ", ~m^3,")")))  + 
              theme(legend.position = "none", panel.background = element_blank()) +
              labs(title=paste("Adj R^2 =",signif(summary(lmfit4)$adj.r.squared, 5),
                               #"Intercept =",signif(lmfit2$coef[[1]],5 ),
                               #" Slope =",signif(lmfit2$coef[[2]], 5),
                               " p =",signif(summary(lmfit4)$coef[2,4], 5)))

anova(lmfit4) # anova table of above lm() fit



# Linear fits on alpha diversity

lmfit5 <- lm(Eff_Num_Sp_mn~chl_mg_m3_mn, data=prod_chl_yr)  # lm() fit 

div_plot3c <- ggplot(prod_chl_yr, aes(x=chl_mg_m3_mn, y=Eff_Num_Sp_mn, color=as.factor(AREA1))) + 
              geom_point(size=4) + scale_color_manual(values=pinks) + theme_boxplot() +
              stat_smooth(method="lm", col="blue", size=2) + ylab("Alpha Diversity") +
              xlab(expression(paste("Primary Productivity (chla mg ", ~m^3,")")))  + 
              theme(legend.position = "right", panel.background = element_blank()) +
              labs(title=paste("Adj R^2 =",signif(summary(lmfit5)$adj.r.squared, 5),
                               #"Intercept =",signif(lmfit2$coef[[1]],5 ),
                               #" Slope =",signif(lmfit2$coef[[2]], 5),
                               " p =",signif(summary(lmfit5)$coef[2,4], 5)))

anova(lmfit5) # anova table of above lm() fit



# Linear fits on beta diversity

lmfit6 <- lm(Exp_B_Div_mn~chl_mg_m3_mn, data=prod_chl_yr)  # lm() fit 

div_plot4c <- ggplot(prod_chl_yr, aes(x=chl_mg_m3_mn, y=Exp_B_Div_mn, color=as.factor(AREA1))) + 
              geom_point(size=4) + theme_boxplot() + scale_color_manual(values=pinks, name="Local\nArea") + 
              stat_smooth(method="lm", col="blue", size=2) + ylab("Beta Diversity") +
              xlab(expression(paste("Primary Productivity (chla mg ", ~m^3,")"))) + 
              theme(legend.position = c(1,0), legend.justification = c(0,0),
                    legend.background = element_rect(fill="transparent"), 
                    panel.background = element_blank()) +
              labs(title=paste("Adj R^2 =",signif(summary(lmfit6)$adj.r.squared, 5),
                               #"Intercept =",signif(lmfit2$coef[[1]],5 ),
                               #" Slope =",signif(lmfit2$coef[[2]], 5),
                               " p =",signif(summary(lmfit6)$coef[2,4], 5)))

anova(lmfit6) # anova table of above lm() fit


# Linear fits on functional diversity

lmfit6b <- lm(RaosQ~chl_mg_m3_mn, data=prod_chl_yr)  # lm() fit 

div_plot4d <- ggplot(prod_chl_yr, aes(x=chl_mg_m3_mn, y=RaosQ, color=as.factor(AREA1))) + 
              geom_point(size=4) + scale_color_manual(values=pinks) + theme_boxplot() +
              stat_smooth(method="lm", col="blue", size=2) + ylab("Functional Diversity") +
              xlab(expression(paste("Primary Productivity (chla mg ", ~m^3,")"))) + 
              theme(legend.position = "none", panel.background = element_blank()) +
              labs(title=paste("Adj R^2 =",signif(summary(lmfit6b)$adj.r.squared, 5),
                               #"Intercept =",signif(lmfit6b$coef[[1]],5 ),
                               #" Slope =",signif(lmfit6b$coef[[2]], 5),
                               " p =",signif(summary(lmfit6b)$coef[2,4], 5)))

anova(lmfit6b) # anova table of above lm() fit



# arranging above plots
#grid.arrange(div_plot2c, div_plot3c, div_plot4c, div_plot4d, ncol=2, nrow=2)

ggarrange(div_plot2c, div_plot3c, div_plot4c, div_plot4d, ncol=2, nrow=2)




### With GAM Fits

# Species Richness 
#assumes normally distributed data 
gam4 <- gam(Sp_rich_mn~chl_mg_m3_mn, data=prod_chl_yr, family=gaussian) # essentially equal to glm()

gam4a <- gam(Sp_rich_mn~s(chl_mg_m3_mn, k=4), data=prod_chl_yr, family=gaussian) 

summary(gam4a) #; anova(gam4a)

# ggplot
gam4_plot <- ggplot(data=prod_chl_yr, aes(x=chl_mg_m3_mn, y=Sp_rich_mn, color=as.factor(AREA1))) +
             geom_point(size=4) + scale_color_manual(values=pinks) + theme_boxplot() +
             stat_smooth(method="gam", formula=y ~ s(x, k=4), col="red", size=2) + ylab("Species Richness") +
             xlab(expression(paste("Primary Productivity (chla mg ", ~m^3,")")))  + 
             theme(legend.position = "none", panel.background = element_blank()) +
             labs(title=paste("Adj R2 =", signif(summary(gam4a)$r.sq, 4),
                              "p =", signif(summary(gam4a)$s.table[1,4], 3)))
#gam4_plot



# Alpha Diversity
#assumes normally distributed data 
gam5 <- gam(Eff_Num_Sp_mn~chl_mg_m3_mn, data=prod_chl_yr, family=gaussian) # essentially equal to glm()

gam5a <- gam(Eff_Num_Sp_mn~s(chl_mg_m3_mn, k=4), data=prod_chl_yr, family=gaussian) 

summary(gam5a) #; anova(gam5a)

# ggplot
gam5_plot <- ggplot(data=prod_chl_yr, aes(x=chl_mg_m3_mn, y=Eff_Num_Sp_mn, color=as.factor(AREA1))) + 
             geom_point(size=4) + scale_color_manual(values=pinks) + theme_boxplot() +
             stat_smooth(method="gam", formula=y ~ s(x, k=4), col="red", size=2) + ylab("Alpha Diversity") +
             xlab(expression(paste("Primary Productivity (chla mg ", ~m^3,")")))  + 
             theme(legend.position = "none", panel.background = element_blank()) +
             labs(title=paste("Adj R2 =", signif(summary(gam5a)$r.sq, 4),
                              "p =", signif(summary(gam5a)$s.table[1,4], 3)))
#gam5_plot



# Beta Diversity
#assumes normally distributed data 
gam6 <- gam(Exp_B_Div_mn~chl_mg_m3_mn, data=prod_chl_yr, family=gaussian) # essentially equal to glm()

gam6a <- gam(Exp_B_Div_mn~s(chl_mg_m3_mn, k=4), data=prod_chl_yr, family=gaussian) 

summary(gam6a) #; anova(gam6a)

# ggplot
gam6_plot <- ggplot(data=prod_chl_yr, aes(x=chl_mg_m3_mn, y=Exp_B_Div_mn, color=as.factor(AREA1))) + 
             geom_point(size=4) + scale_color_manual(values=pinks) + theme_boxplot() +
             stat_smooth(method="gam", formula=y ~ s(x, k=4), col="red", size=2) + ylab("Beta Diversity") +
             xlab(expression(paste("Primary Productivity (chla mg ", ~m^3,")")))  + 
             theme(legend.position = "none", panel.background = element_blank()) +
             labs(title=paste("Adj R2 =", signif(summary(gam6a)$r.sq, 4),
                              "p =", signif(summary(gam6a)$s.table[1,4], 3)))
#gam6_plot



# Functional Diversity
#assumes normally distributed data 
gam7 <- gam(RaosQ~chl_mg_m3_mn, data=prod_chl_yr, family=gaussian) # essentially equal to glm()

gam7a <- gam(RaosQ~s(chl_mg_m3_mn, k=4), data=prod_chl_yr, family=gaussian) 

summary(gam7a) #; anova(gam7a)

# ggplot
gam7_plot <- ggplot(data=prod_chl_yr, aes(x=chl_mg_m3_mn, y=RaosQ, color=as.factor(AREA1))) + 
             geom_point(size=4) + scale_color_manual(values=pinks) + theme_boxplot() +
             stat_smooth(method="gam", formula=y ~ s(x, k=4), col="red", size=2) + ylab("Functional Diversity") +
             xlab(expression(paste("Primary Productivity (chla mg ", ~m^3,")")))  + 
             theme(legend.position = "none", panel.background = element_blank()) +
             labs(title=paste("Adj R2 =", signif(summary(gam7a)$r.sq, 4),
                              "p =", signif(summary(gam7a)$s.table[1,4], 3)))
#gam7_plot



# arranging above plots
ggarrange(gam4_plot, gam5_plot, gam6_plot, gam7_plot, ncol=2, nrow=2)


# arranging above plots
ggarrange(gam4_plot, gam5_plot, div_plot4c, div_plot4d, ncol=2, nrow=2)


# correlation between primary productivity and fishing pressure
prod_over <- prod_chl_yr %>% 
             dplyr::filter(Year %in% c(2003:2013)) %>%
             dplyr::select(AREA, Year, chl_mg_m3_mn) %>% 
             dplyr::group_by(Year) %>%
             dplyr::summarize(ann_mn_chl = mean(chl_mg_m3_mn)) %>%
             dplyr::ungroup()

fish_over <- dist_df_cat_ar %>%
             dplyr::filter(Year %in% c(2003:2013)) %>%
             dplyr::select(Catch_area, Year, met_ton_km2_ca) %>%
             dplyr::group_by(Year) %>%
             dplyr::summarize(ann_mn_fish = mean(met_ton_km2_ca)) %>%
             dplyr::ungroup() 

prod_fish_over <- fish_over %>%
                  full_join(prod_over, by="Year")

ppfp <- ggplot(data=prod_fish_over, aes(x=ann_mn_fish, y=ann_mn_chl)) + 
        geom_point(size=3) + geom_text(aes(label=Year), hjust=-0.25, vjust=0) + 
        theme_boxplot() + xlab("Annual mean fishing pressure") + ylab("Annual mean Chla")
ppfp










