# calculate estimated optimal depth for each species
# occurrence and positive.
Name <- paste(getwd(),"/goaTrawl/_Estimated Coefs/All_sp_depth_coef.csv",sep="")
dat <- read.csv(Name)

anti.logit <- function(A){
  return(1/(1+exp(-A)))
}

quad.func <- function(dat,log.D,mod){
  if(mod=="binomial"){
    if(nrow(dat)>=3){
      Int   <- dat$mean[dat$Coef.name  == "Y1984"]
      Lin   <- dat$mean[dat$Coef.name  == "cent.log.depth"]
      Quad  <- dat$mean[dat$Coef.name == "cent.log.depth2"]
    }
    if(nrow(dat)>3){
      Int     <- dat$mean[dat$Coef.name  == "Y1984"]
        temp <- dat[grep("Y",temp$Coef.name),]
        temp <- temp[temp$Coef.name != "Y1984",]
        Int  <- mean(c(Int,temp$mean+Int))
  
      Lin   <- dat$mean[dat$Coef.name  == "cent.log.depth"]
      Quad  <- dat$mean[dat$Coef.name == "cent.log.depth2"]
    }
    Out   <- Int + Lin*log.D + Quad*log.D^2
  }
  if(mod=="pos"){
    if(nrow(dat)>=3){
      Int   <- dat$mean[dat$Coef.name  == "Y1984"]
      Lin   <- dat$mean[dat$Coef.name  == "cent.log.depth"]
      Quad  <- dat$mean[dat$Coef.name == "cent.log.depth2"]
    }
    if(nrow(dat)>3){
      Int     <- dat$mean[dat$Coef.name  == "Y1984"]
      temp <- dat[grep("Y",temp$Coef.name),]
      temp <- temp[temp$Coef.name != "Y1984",]
      Int  <- mean(c(Int,temp$mean+Int))
      
      Lin   <- dat$mean[dat$Coef.name  == "cent.log.depth"]
      Quad  <- dat$mean[dat$Coef.name == "cent.log.depth2"]
    }
    Out   <- Int + Lin*log.D + Quad*log.D^2
  }
  if(mod=="binomial"){Out <- anti.logit(Out)}
  if(mod=="pos"){Out <- exp(Out)}
  
  return(Out)
}


##### 
Bin <- NULL
Pos <- NULL 

species <- sort(unique(dat$Species))
mod     <- c("binomial","pos")
for(i in 1: length(species)){
  for(j in 1:length(mod)){
      temp  <- dat[dat$Species == species[i] & dat$Model == mod[j],] 
      D     <- seq(50,600,length.out=10000)
      log.D <- log(D) - mean(dat$log.center)
        if(nrow(temp)>3){
          X <- quad.func(temp,log.D,mod[j]) 
        }
        if(nrow(temp)<=3){
          X <- quad.func(temp,log.D,mod[j]) 
        }

      if(mod[j]=="binomial"){
        Bin <- cbind(Bin,X)
      }
      if(mod[j]=="pos"){
        Pos <- cbind(Pos,X)
      }
  }
}

colnames(Bin) <- species
colnames(Pos) <- species
Bin <- data.frame(depth.m=D,Bin)
Pos <- data.frame(depth.m=D,Pos)
Both <- list(Bin,Pos)

save(file=paste(getwd(),"/goaTrawl/_Estimated Coefs/Depth predictions mean.R",sep=""),Both)
 
