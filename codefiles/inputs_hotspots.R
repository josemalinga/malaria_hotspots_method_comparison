##################################################################
#Author: Josephine Malinga
#Date: Original (2019)
#
#Updated September 30th 2021
#
##################################################################

rm=lis=ls()
setwd(getwd())
library(plyr); library(dplyr); library(data.table); library(maptools); library(gdata);library(ggplot2)
library(dbscan); library(fpc); library(tripack); library(spdep); library(reshape);library(spData)
library("Sim.DiffProc"); library(scales); library(qpcR)

# Homesteads estimated from the Kilifi HDSS
# simulate for now - 2000 homesteads on a 10by10 grid with 8 people per house 
num.homesteads<-2000 
numPple<-8
homesteads<-array(-9, dim=c(num.homesteads,3))
colnames(homesteads)<-c("location.id","locx","locy")
homesteads[,1]<-seq(1:num.homesteads)
homesteads[,2]<-runif(num.homesteads,0.00,10.00) 
homesteads[,3]<-runif(num.homesteads,0.00,10.00) 
homesteads<-data.frame(homesteads)
plot(homesteads$locx, homesteads$locy, col="red", pch=1, main="cases")

write.table(homesteads,"homesteads.txt", row.names = F)

# Compute distances between pairs of homesteads
# Run the function to calculate distance between pairs in kms
distance_km<-function(i,j,data_km){
  dist<-sqrt((data_km$x[i]-data_km$x[j])^2+(data_km$y[i]-data_km$y[j])^2)
  return(dist)
}
dist_km<-Vectorize(distance_km,vectorize.args=list("i","j"))

# data on grid map to match the 10km*10km simulated area
data_km<-data.frame(x=homesteads$locx,y=homesteads$locy)
time_km<-system.time(res_km<-outer(1:num.homesteads,1:num.homesteads,dist_km,data_km=data_km))
#maximum distance is  km
hist(res_km)

id<-1:nrow(res_km); res_km<-cbind(res_km, id) 

#calculate the pair distance between all houses
#to get the distance difference for i<j
#remove half the triangle
res_km1<-res_km
res_km1[upper.tri(res_km1,diag=FALSE)] <- NA

#reshaping the data into long format
res_km1<-as.data.frame(res_km1)
res_km1$id<-1:nrow(res_km1)

diffDist<-melt(res_km1,id=c("id"))

dist_diff<-diffDist[complete.cases(diffDist), ]
dist_diff$id_x<-as.numeric(dist_diff$variable); dist_diff$variable<-NULL

houses<-homesteads; names(houses)[1]<-"id"
dist.Diff<-merge(houses, dist_diff, by = "id")

colnames(houses)<-c("id_x", "locx_x", "locy_x")
dist.Diff<-merge(houses, dist.Diff, by = "id_x"); names(dist.Diff)[7]<-"dist.diff"
dist.Diff <- subset(dist.Diff, select=c(id_x,id,dist.diff,locx,locy,locx_x,locy_x))

write.table(dist.Diff,"dist_diff.txt", row.names = F)

#to identify neigh of points by eucledian distance 
#to identify neigh of points by great circle distance in km (longlat=T)
#d1=min;d2=max
xycoords <- cbind(homesteads$locx, homesteads$locy)
nears <- dnearneigh(xycoords, 0, 1.0)
nears

neig<-array(NA, dim=c(1,length(nears)))
for (i in 1:length(nears)){
  sam<-nears[[i]]
  neig<-qpcR:::rbind.na(neig, sam)
  
}

neigh1<-neig[-1,]
neigh<-data.frame(seq(1:num.homesteads),neigh1)

# na="" when saving replaces the "NA" values with blanks
write.table(neigh,"near_matrix.txt", row.names = F,col.names = F, na = "")






