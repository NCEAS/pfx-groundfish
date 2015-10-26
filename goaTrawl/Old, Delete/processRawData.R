library(INLA)

#setwd("/users/eric.ward/documents/exxonValdez_nceas/goaTrawl")
setwd("/users/ole.shelton/GitHub/exxonValdez_nceas/goaTrawl")



d = read.csv("goa_data.csv")

# Get rid of records with no bottom temp
d = d[which(d$BOT_TEMP != -9999 & d$SCIENTIFIC != ""),]
d$STATIONYEAR = paste(d$STATION,d$YEAR)

# Convert the raw data to more of the format used by NWFSC -- each spp in a column
species = unique(d$SCIENTIFIC)

# Only grab species with more than 500 observations = 78
spp500 = names(table(d$SCIENTIFIC)[which(table(d$SCIENTIFIC)>500)])

# convert data into data frame of tow / lat / lon / YEAR / 
tows = unique(as.character(d$STATION))

df = data.frame("StationYear"=unique(d$STATIONYEAR), "Lat"= d$LATITUDE[match(unique(d$STATIONYEAR),as.character(d$STATIONYEAR))], 
"Lon" = d$LONGITUDE[match(unique(d$STATIONYEAR),as.character(d$STATIONYEAR))],"Depth"= d$BOT_DEPTH[match(unique(d$STATIONYEAR),as.character(d$STATIONYEAR))],
"Temp"=d$BOT_TEMP[match(unique(d$STATIONYEAR),as.character(d$STATIONYEAR))])

# Add each of the species columns to the data.frame
df[spp500] = NA

for(i in 1:dim(df)[1]) {
  for(j in 7:dim(df)[2]) {
    indx = which(as.character(d$STATIONYEAR)==df$StationYear[i] & d$SCIENTIFIC == names(df)[j])
    if(length(indx)>0) {
      # Enter data
      df[i,j] = mean(d$WTCPUE[indx])
    }
  }
}

df$PID = 1
df$POS = seq(1,dim(df)[1])
df["X"]=df$Lon
df["Y"]=df$Lat

attr(df, "zone") <- 7
attr(df, "projection") ="LL"
df = convUL(xydata = df)
coords = cbind(utm$X,utm$Y)


L = unlist(strsplit(as.character(df[,1])," "))

df$Station = L[seq(1,length(L),2)]
df$Year = as.numeric(L[seq(2,length(L),2)])



