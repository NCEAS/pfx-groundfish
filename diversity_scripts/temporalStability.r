# -----------------------------------------------------
# title: "temporalStability"
# author: "Colette Ward"
# date: "10/16/2017"
# output: pdf_document
# -----------------------------------------------------
  
  
# Load packages
library(plyr); library(dplyr); library(tidyr); library(ggplot2)



# load & prep look-up table of common names
common <- read.csv("./diversity-data/trawl_species_control_file.csv", header = T, stringsAsFactors = FALSE)

common1 <- common %>%
  select(database.name, common.name) %>%
  rename(Species = database.name)
for (i in 1:nrow(common1)) { # add common names for Sebastes 1 & 2
  if(common1$Species[i] == "Dusky.and.Dark.Rockfish") {common1$common.name[i] <- "Sebastes 1"}
  if(common1$Species[i] == "Rougheye.and.Blackspotted.Rockfish") {common1$common.name[i] <- "Sebastes 2"}
}



# load & prep mean annual CPUE data

# load file for local communities:
shallowCPUEArea <- read.csv("./diversity-data/All_sp_index_meanCPUEByArea.Shallow.MH.final.csv", header = T, stringsAsFactors = FALSE)
shallowCPUEArea1 <- shallowCPUEArea %>% 
  select(-area_alph, -oldArea) %>%
  mutate(area = as.character(area))

# load file for Regional metacommunity:
shallowCPUETot <- read.csv("./diversity-data/All_sp_index_meanCPUEByArea.All.csv", header = T, stringsAsFactors = FALSE) 


# merge and clean up the files:
shallowCPUEArea2 <- bind_rows(shallowCPUEArea1, shallowCPUETot) %>% # merge local and metacommunity files
  left_join(common1, by = "Species") %>% # merge common names onto SPCPUEArea
  select(area, Mean.totalDensity, SD.totalDensity, year, Species, common.name) %>%
  mutate(area = revalue(area, c("Total" = "Region")))




# To address local vs regional community stability we can compare CVs of total CPUE between local communities and the regional metacommunity.
# The regional community does indeed show greater stability than local communities. Note that y-axis scales differ - deep areas are generally less stable than shallow areas.  
# **The magnitude of the portfolio effect can be assessed from the ratio of mean local CV / regional CV (ie the increase in stability that arises from spatial (beta) diversity; values >1 indicate a stabilizing effect): the ratio is ~2.3 for shallow areas.**


# Calculate CV
CV_shallow <- shallowCPUEArea2 %>% 
  group_by(area, year) %>%
  summarise(total = sum(Mean.totalDensity)) %>%
  ungroup() %>%
  
  group_by(area) %>%
  summarise(CV = sd(total) / mean(total)) %>%
  ungroup()


# mean CV of local areas
CV_shallow %>%
  filter(area != "Region") %>%
  summarise(mean = mean(CV), sd = sd(CV))
0.1744/0.0845 # ratio of mean local CV / regional CV = 2.06


# plot on a log scale
CV_shallow_plot <- ggplot(data=CV_shallow, aes(x=area, y = CV)) + 
  geom_point(aes(y = CV), size=5) +
  
  theme(panel.background = element_blank(),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        axis.text.x = element_text(margin=margin(5,0,0,0), size=12), 
        axis.text.y = element_text(margin=margin(0,5,0,0), size=12),
        axis.title.x=element_text(size=12, margin=margin(15,0,0,0)),
        axis.title.y=element_text(size=12, angle=90, margin=margin(0,15,0,0)),
        plot.title = element_blank()) +
  
  scale_x_discrete(limits=c(1:10, "Region")) +
  scale_y_log10(breaks = c(0.07, 0.08, 0.09, 0.1, 0.15, 0.2)) +
  labs(x = "Local Community", y = "CV (Total CPUE)", title = "Shallow Areas")

#CV_shallow_plot




