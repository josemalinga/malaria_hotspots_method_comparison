##################################################################
#Author: Josephine Malinga
#Date: Original (2019)
#
#Updated September 30th 2021
#
##################################################################
rm(list=ls())
getwd()

library(plyr); library(dplyr); library(data.table); library(maptools); library(gdata);library(ggplot2)
library(dbscan); library(fpc); library(tripack); library(spdep); library(reshape);library(spData)
library("Sim.DiffProc"); library(scales)

homesteads<-read.table("homesteads.txt", header=T, row.names = NULL)
maxNumPplehh<-8
numPple<-8
numSNPsPerInf<-1

#sample PEOPLE in the population who are infected
populate <- homesteads %>% slice(rep(1:n(), each = numPple))
names(populate)[1]<-"location.id"

setDT(populate)[, pID := seq_len(.N), by = location.id]
populate<-data.frame(populate)

#########################################################################################3
##random initial infections at the beginning of the simulation
# randomly select houses PEOPLE??? for the initial infections
# assuming the prevalence of 5% in the infection data

#CHECK HOW TO DEAL WITH MEASUREMENT ERROR SINCE WE ARE ASSUMING PERFECT DETECTION W PCR
#WE ASSUME THE INITIAL INFECTIONS - ONE CASE PER PERSON
prevLM <- 0.10
#moi<-1.2
#logoddsprevLM<-log(prevLM/(1-prevLM))
#logoddsprevpcr<-0.954+0.868*logoddsprevLM
#prevpcr<-exp(logoddsprevpcr)/(exp(logoddsprevpcr)+1)
#prevpcrmoi<-prevpcr*moi

numPpleWithInitInf<-trunc(nrow(populate)*prevLM)
maxNumPplehh<-8
cases<-sample_n(populate,numPpleWithInitInf,replace = F)
cases<-data.frame(cases)
write.table(cases,"initialInfs.txt", row.names = F)

#SINGLE HOTSPOT
#create the hotspots with different radii for the initial infections
homes<-homesteads[sample(nrow(homesteads), 1), ]
x <- homes$locx
y <- homes$locy
plot(x, y, ylim=c(0,10), xlim=c(0,10),pch=20)
#on whether its a single or multiple clusters - change the radii accordingly
sizeOfRadius=2.0

symbols(x,y,circles=sizeOfRadius,ylim=c(0,10), xlim=c(0,10), fg=c("red","blue","tan","green"), lwd=5, main = "hotspot1")
points(homes[,2:3],pch=20, cex=0.8,main="cases", col="red")

hotHme<-array(-9, dim=c(nrow(homesteads),2))
for (i in 1:nrow(homesteads)){
  dist<-sqrt((x-homesteads$locx[i])^2+(y-homesteads$locy[i])^2)
  if (dist<=sizeOfRadius){
    hotHme[i,1]<-homesteads$locx[i]
    hotHme[i,2]<-homesteads$locy[i]
  }
}

hotHme<-data.frame(hotHme[hotHme[,1]>0,])
colnames(hotHme)<-c("locx","locy")

hotHme$hotHme<-1
hotHme1<-merge(homesteads, hotHme, by=c("locx","locy"), all=T)
hotHme1$hotHme[is.na(hotHme1$hotHme)]<-0
names(hotHme1)[3:4]<-c("location.id","hotspot")

hotHme2<-data.frame(hotHme1[hotHme1$hotspot==1,])

write.table(hotHme1,"hotspots.txt", row.names = F)
write.table(hotHme2,"hotHomes.txt", row.names = F)


#MULTIPLE HOTSPOTS
homesteads<-read.table("homesteads.txt", header=T, row.names = NULL)
#create the hotspots with different radii for the initial infections
homes<-homesteads[sample(nrow(homesteads), 2), ]
x <- homes$locx
y <- homes$locy
plot(x, y, ylim=c(0,10), xlim=c(0,10),pch=20)
sizeOfRadius= c(1.35,1.35)
symbols(x,y,circles=sizeOfRadius,ylim=c(0,10), xlim=c(0,10), fg=c("red","blue","tan","green"), lwd=5, main = "hotspot1")
points(homes[,2:3],pch=20, cex=0.8,main="cases", col="red")

hot1<-array(-9, dim=c(nrow(homesteads),2))
for (j in 1:length(sizeOfRadius)){
  hotHme<-array(-9, dim=c(nrow(homesteads),2))
  for (i in 1:nrow(homesteads)){
    dist<-sqrt((x[j]-homesteads$locx[i])^2+(y[j]-homesteads$locy[i])^2)
    if (dist<=sizeOfRadius[j]){
      hotHme[i,1]<-homesteads$locx[i]
      hotHme[i,2]<-homesteads$locy[i]
    }
  }
  hot1<-rbind(hot1,hotHme)
}

hotHme0<-data.frame(hot1[hot1[,1]>0,])
hotHme<-hotHme0[!duplicated(hotHme0[,c(1,2)]),]
colnames(hotHme)<-c("locx","locy")

hotHme$hotHme<-1
hotHme1<-merge(homesteads, hotHme, by=c("locx","locy"), all=T)
hotHme1$hotHme[is.na(hotHme1$hotHme)]<-0
names(hotHme1)[3:4]<-c("location.id","hotspot")
hotHme2<-data.frame(hotHme1[hotHme1$hotspot==1,])

write.table(hotHme1,"hotspots_M.txt", row.names = F)
write.table(hotHme2,"hotHomes_M.txt", row.names = F)


#DECAY GRADIENT SINGLE HOTSPOT
homesteads<-read.table("homesteads.txt", header=T, row.names = NULL)
homes<-homesteads[sample(nrow(homesteads), 1), ]
x <- homes$locx
y <- homes$locy

sizeOfRadius=2.0

##############################################33
#ASSUMING A DECAY AT THE AGE OF THE GRID (DO WE RESTRICT DISTANCE?)
hotspotFactor=2
hotHme<-array(-9, dim=c(nrow(homesteads),3))
for (i in 1:nrow(homesteads)){
  dist<-sqrt((x-homesteads$locx[i])^2+(y-homesteads$locy[i])^2)
  if (dist<=sizeOfRadius){
    hotHme[i,1]<-homesteads$locx[i]
    hotHme[i,2]<-homesteads$locy[i]
    hotHme[i,3]<-1
  }
  
  if (dist>sizeOfRadius & dist<(sizeOfRadius*1.2)){
    hotHme[i,1]<-homesteads$locx[i]
    hotHme[i,2]<-homesteads$locy[i]
    #hotHme[i,3]<-exp(-(dist^2))
    hotHme[i,3]<-exp(-(dist^2)/(2*(hotspotFactor^2)))
  }
}

#remove arrays with -9
hotHme<-data.frame(hotHme[hotHme[,1]>0,])
colnames(hotHme)<-c("locx","locy","hotFactor")
#hotHme$hotFactor<-rescale(hotHme$hotFactor, to=c(0.5,1))

#hotHme$hotHme<-1
hotHme1<-merge(homesteads, hotHme, by=c("locx","locy"), all=T)
#names(hotHme)[3:5]<-c("location.id","hotFactor","hotspot")
hotHme1$hotFactor[is.na(hotHme1$hotFactor)]<-0.00001
hotHme1$hotFactor[hotHme1$hotFactor>0.001 & hotHme1$hotFactor<1]<-rescale(hotHme1$hotFactor[hotHme1$hotFactor>0.001 & hotHme1$hotFactor<1], to=c(0,1))

hotHme1$hotspot[hotHme1$hotFactor>0.001]<-1
hotHme1$hotspot[is.na(hotHme1$hotspot)]<-0

hotHme2<-data.frame(hotHme1[hotHme1$hotspot==1,])

#table(hotHme1$hotFactor[hotHme1$hotFactor==1])
#table(hotHme1$hotspot)
plot(hotHme1$locx, hotHme1$locy, col="blue")
points(hotHme2$locx, hotHme2$locy, col=hotHme2$hotFactor)

write.table(hotHme1,"hotspots_D.txt", row.names = F)
write.table(hotHme2,"hotHomes_D.txt", row.names = F)


#IRREGULAR HOTSPOTS - MAYBE ALONG A RIVER AT FIRST
homesteads<-read.table("homesteads.txt", header=T, row.names = NULL)
#simulate a TS object and mimic this pattern
npint=50

x1<-BM(N = npint, t0 = 0, T = 0.5, C = 2)
peri<-time(x1)
x2<-data.frame(x1)
#ran1<-sample(0:10,1)
#ran2<-sample(0:10,1)
x<-rescale(as.numeric(peri), to=c(0,10))
y<-rescale(x2$X1, to=c(0,10))
x2$x<-x;x2$y<-y
plot(x2$x, x2$y,type="o",xlim=c(0,10),ylim=c(0,10))

sizeOfRadius= rep(0.3, npint)   #runif(2)
#symbols(x,y,circles=sizeOfRadius,ylim=c(0,10), xlim=c(0,10), lwd=5, main = "hotspot1")
#points(homes[,2:3],pch=20, cex=0.8,main="cases", col="red")

hot1<-array(-9, dim=c(nrow(homesteads),2))
for (j in 1:length(sizeOfRadius)){
  hotHme<-array(-9, dim=c(nrow(homesteads),2))
  for (i in 1:nrow(homesteads)){
    dist<-sqrt((x[j]-homesteads$locx[i])^2+(y[j]-homesteads$locy[i])^2)
    if (dist<=sizeOfRadius[j]){
      hotHme[i,1]<-homesteads$locx[i]
      hotHme[i,2]<-homesteads$locy[i]
    }
  }
  hot1<-rbind(hot1,hotHme)
}

hotHme0<-data.frame(hot1[hot1[,1]>0,])
hotHme0<-hotHme0[,1:2]

hotHme<-hotHme0[!duplicated(hotHme0[,c(1,2)]),]
colnames(hotHme)<-c("locx","locy")

hotHme$hotHme<-1

hotHme1<-merge(homesteads, hotHme, by=c("locx","locy"),all=T)
hotHme1$hotHme[is.na(hotHme1$hotHme)]<-0
names(hotHme1)[3:4]<-c("location.id","hotspot")

hotHme2<-data.frame(hotHme1[hotHme1$hotspot==1,])

table(hotHme1$hotspot)

write.table(hotHme1,"hotspots_R.txt", row.names = F)
write.table(hotHme2,"hotHomes_R.txt", row.names = F)


