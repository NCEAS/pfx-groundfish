library(ggplot2)
theme_set(theme_bw(base_size=12)+ 
            theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank()))
library(dplyr)
library(tidyr)
library(reshape2)
library(ggtern)


#Pull in original data
meanCPUE <- read.csv("All_sp_index_meanCPUEByArea.csv") # load 
meanCPUE <- subset(meanCPUE,area!='Total') #remove total field
meanCPUE <- droplevels(meanCPUE) 
meanCPUE$area <- as.numeric(levels(meanCPUE$area))[meanCPUE$area]#use numeric/factor conventions
meanCPUE$Area <- factor(meanCPUE$area)
meanCPUE$vari <- meanCPUE$SD.totalDensity^2 #calculate variance
#meanCPUE$Species <- gsub("[.]","",meanCPUE$Species) #quit changing how names are spelled damnit!

meanOcc <- meanOcc %>%
  mutate(year = ifelse((year.numb=="1"),'1984',
                       ifelse((year.numb=="2"),'1987',
                      ifelse((year.numb=="3"),'1990',
                       ifelse((year.numb=="4"),'1993',       
                        ifelse((year.numb=="5"),'1996',
                        ifelse((year.numb=="6"),'1999',
                       ifelse((year.numb=="7"),'2001',
                       ifelse((year.numb=="8"),'2003',
                       ifelse((year.numb=="9"),'2005',
                        ifelse((year.numb=="10"),'2007',
                       ifelse((year.numb=="11"),'2009',
                       ifelse((year.numb=="12"),'2011',
                       ifelse((year.numb=="13"),'2013',
                       ifelse((year.numb=="14"),'2015',        
                     ""))))))))))))))
  )

#Input trawl species files
trawl_species <- read.csv("./goaTrawl/Output Data/trawl_species_control_file.csv")
s <- trawl_species$diet
s <- as.data.frame(s) %>% separate(s, into = paste("diet", 1:2, sep = ""))
trawl_species <- cbind(trawl_species,s)

temp <- trawl_species[,c("database.name","fish.invert","pelagic.benthic","total.biomass.fish","guild","diet1","diet2")]
colnames(temp)[1] <- "Species"

meanCPUE <- merge(meanCPUE,temp)

#Pull data for plotting
meanCPUE %>% 
  select(Species,Median.totalDensity, year, Area, diet1, guild) %>% 
  mutate(Spill = ifelse(year<1990, 'Pre', 'Post'),
         Year = factor(year))-> dat

#Diet data
dat %>% 
  group_by(Year, diet1, Area, Spill) %>% 
  summarize(density = sum(Median.totalDensity)) -> dat1

dat1 <- spread(dat1, diet1, density)

#Guild data
dat %>% 
  group_by(Year, guild, Area, Spill) %>% 
  summarize(density = sum(Median.totalDensity)) -> dat2

dat2 <- spread(dat2, guild, density)

#Diet by area

ggtern(data = dat1, aes(x = F, y = G, z = I)) + 
  geom_point(aes(fill = Area),
             size = 4, 
             shape = 21, 
             color = "black",
             alpha=.35) + 
  facet_wrap(~Year) +
  theme_hidegrid() + theme_hidelabels() + 
  theme(legend.justification=c(1,0), legend.position=c(1,0)) + 
  guides(fill=guide_legend(ncol=3))

#Diet by year

ggtern(data = dat1, aes(x = F, y = G, z = I)) + 
  geom_point(aes(fill = Year),
             size = 4, 
             shape = 21, 
             color = "black",
             alpha=.35) +facet_wrap(~Area)+theme_hidegrid()+theme_hidelabels()+
  theme(legend.justification=c(1,0), legend.position=c(1,0)) + 
  guides(fill=guide_legend(ncol=3))

#Diet pre/post spill

ggtern(data = dat1, aes(x = F, y = G, z = I)) + 
  geom_point(aes(fill = Spill),
             size = 4, 
             shape = 21, 
             color = "black",
             alpha=.35) +facet_wrap(~Area)+
  theme_hidegrid()+theme_hidelabels() +
  theme(legend.justification=c(1,0), legend.position=c(1,0))

#Guild by area
ggtern(data = dat2, aes(x = A, y = B, z = P)) + 
  geom_point(aes(fill = Area),
             size = 4, 
             shape = 21, 
             color = "black",
             alpha=.5) +facet_wrap(~Year)+
  theme_hidegrid()+theme_hidelabels()+ 
  theme(legend.justification=c(1,0), legend.position=c(1,0)) + 
  guides(fill=guide_legend(ncol=3))

#Guild by year
ggtern(data = dat2, aes(x = A, y = B, z = P)) + 
  geom_point(aes(fill = Year),
             size = 4, 
             shape = 21, 
             color = "black",
             alpha=.5) +facet_wrap(~Area)+
  theme_hidegrid()+theme_hidelabels()+ 
  theme(legend.justification=c(1,0), legend.position=c(1,0)) + 
  guides(fill=guide_legend(ncol=3))
#Guild by pre/post
ggtern(data = dat2, aes(x = A, y = B, z = P)) + 
  geom_point(aes(fill = Spill),
             size = 4, 
             shape = 21, 
             color = "black",
             alpha=.5) +facet_wrap(~Area) +theme_hidegrid()+theme_hidelabels()
