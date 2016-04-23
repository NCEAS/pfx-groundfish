##############################################################
### PFX Groundfish
### Script to get means per species per year per area
### Rachael Blake, April 2016
##############################################################

Sh_CPUE <- read.csv("C:/Users/rblake/Documents/NCEAS/GoA Portfolio Effects WG/pfx-groundfish/goaDat.uncond.Shallow.AllData.Combined.7.8.9.csv",
                    stringsAsFactors=FALSE, header=TRUE)

Dp_CPUE <- read.csv("C:/Users/rblake/Documents/NCEAS/GoA Portfolio Effects WG/pfx-groundfish/goaDat.uncond.Deep.AllData.csv",
                    stringsAsFactors=FALSE, header=TRUE)

Sh_OCCR <- read.csv("C:/Users/rblake/Documents/NCEAS/GoA Portfolio Effects WG/pfx-groundfish/goaDat.occ.Shallow.AllData.Combined.7.8.9.csv",
                    stringsAsFactors=FALSE, header=TRUE)

Dp_OCCR <- read.csv("C:/Users/rblake/Documents/NCEAS/GoA Portfolio Effects WG/pfx-groundfish/goaDat.occ.Deep.AllData.csv",
                    stringsAsFactors=FALSE, header=TRUE)







