##################################################################
#Author: Josephine Malinga
#Date: Original (2019)
#
#Updated September 30th 2021
#
##################################################################

############################
rm(list=ls())
#setwd("/Users/dukeghi/Documents/SCICORE Files/Malaria Hotspots/simulations_sep/rep64")

getwd() #23,28,37,3,15, 42

# CIRCULAR HOTSPOT
#single hotspot
library(dplyr)
homesteads<-read.table("homesteads.txt", header = T)
# after obtaining the results from the saTScan model
data<-read.table("res0.gis.txt", header = T, sep="")
data<-data[data$P_VALUE<0.05,]
data1<-read.table("res0.col.txt", header = T, sep="")
data1<-data1[data1$P_VALUE<0.05,]
data2<-merge(data, data1, by=c("LOC_ID", "CLUSTER"), all = T)
colnames(homesteads)[1]<-c("LOC_ID")
data3<-merge(data2, homesteads, by=c("LOC_ID"), all=T)
data3$REL_RISK[is.na(data3$REL_RISK)]<-1.0
data3<-data3[!duplicated(data3$LOC_ID),]
data3$high[data3$REL_RISK>1.0]<-1
data3$high[data3$REL_RISK<1.0]<-0
data3$high[is.na(data3$high)]<-0
data3$CLUSTER[is.na(data3$CLUSTER)]<-0
table(data3$CLUSTER)
data3$hotspot<-0

# the significant clusters
data3$hotspot[data3$CLUSTER>0]<-1
#names(data3)[1]<-"LOC_ID"
data4<-data3[,c("LOC_ID","hotspot")]; 
data5<-merge(data4, homesteads, by="LOC_ID")
names(data5)[1:2]<-c("location.id", "hotspots_satO")
# from observed
hotspots<-read.table("hotspots.txt", header = T)
hotnew<-merge(hotspots, data5, by=c("location.id", "locx", "locy") )
par(mfrow=c(2,2))
write.table(hotnew,"hotspots_single.txt", row.names = F)

############################
#multiple hotspot
homesteads<-read.table("homesteads.txt", header = T)
# after obtaining the results from the saTScan model
data<-read.table("res0_M.gis.txt", header = T, sep="")
data<-data[data$P_VALUE<0.05,]
data1<-read.table("res0_M.col.txt", header = T, sep="")
data1<-data1[data1$P_VALUE<0.05,]
data2<-merge(data, data1, by=c("LOC_ID", "CLUSTER"), all = T)
colnames(homesteads)[1]<-c("LOC_ID")
data3<-merge(data2, homesteads, by=c("LOC_ID"), all=T)
data3$REL_RISK[is.na(data3$REL_RISK)]<-1.0
data3<-data3[!duplicated(data3$LOC_ID),]
data3$high[data3$REL_RISK>1.0]<-1
data3$high[data3$REL_RISK<1.0]<-0
data3$high[is.na(data3$high)]<-0
data3$CLUSTER[is.na(data3$CLUSTER)]<-0
table(data3$CLUSTER)
data3$hotspot<-0
# the significant clusters
data3$hotspot[data3$CLUSTER>0]<-1
#names(data3)[1]<-"LOC_ID"
data4<-data3[,c("LOC_ID","hotspot")]; 
data5<-merge(data4, homesteads, by="LOC_ID")
names(data5)[1:2]<-c("location.id", "hotspots_satO")
# from observed
hotspots<-read.table("hotspots_M.txt", header = T)
hotnew<-merge(hotspots,data5, by=c("location.id","locx","locy"))
write.table(hotnew,"hotspots_multiple.txt", row.names = F)

############################
#decay hotspot
homesteads<-read.table("homesteads.txt", header = T)
# after obtaining the results from the saTScan model
data<-read.table("res0_D.gis.txt", header = T, sep="")
data<-data[data$P_VALUE<0.05,]
data1<-read.table("res0_D.col.txt", header = T, sep="")
data1<-data1[data1$P_VALUE<0.05,]
data2<-merge(data, data1, by=c("LOC_ID", "CLUSTER"), all = T)
colnames(homesteads)[1]<-c("LOC_ID")
data3<-merge(data2, homesteads, by=c("LOC_ID"), all=T)
data3$REL_RISK[is.na(data3$REL_RISK)]<-1.0
data3<-data3[!duplicated(data3$LOC_ID),]
data3$high[data3$REL_RISK>1.0]<-1
data3$high[data3$REL_RISK<1.0]<-0
data3$high[is.na(data3$high)]<-0
data3$CLUSTER[is.na(data3$CLUSTER)]<-0
table(data3$CLUSTER)
data3$hotspot<-0
# the significant clusters
data3$hotspot[data3$CLUSTER>0]<-1
#names(data3)[1]<-"LOC_ID"
data4<-data3[,c("LOC_ID","hotspot")]; 
data5<-merge(data4, homesteads, by="LOC_ID")
names(data5)[1:2]<-c("location.id", "hotspots_satO")
# from observed
hotspots<-read.table("hotspots_D.txt", header = T)
hotnew<-merge(hotspots, data5, by=c("location.id", "locx", "locy") )
write.table(hotnew,"hotspots_decay.txt", row.names = F)

############################
#river hotspot
homesteads<-read.table("homesteads.txt", header = T)
# after obtaining the results from the saTScan model
data<-read.table("res0_R.gis.txt", header = T, sep="")
data<-data[data$P_VALUE<0.05,]
data1<-read.table("res0_R.col.txt", header = T, sep="")
data1<-data1[data1$P_VALUE<0.05,]
data2<-merge(data, data1, by=c("LOC_ID", "CLUSTER"), all = T)
colnames(homesteads)[1]<-c("LOC_ID")
data3<-merge(data2, homesteads, by=c("LOC_ID"), all=T)
data3$REL_RISK[is.na(data3$REL_RISK)]<-1.0
data3<-data3[!duplicated(data3$LOC_ID),]
data3$high[data3$REL_RISK>1.0]<-1
data3$high[data3$REL_RISK<1.0]<-0
data3$high[is.na(data3$high)]<-0
data3$CLUSTER[is.na(data3$CLUSTER)]<-0
table(data3$CLUSTER)
data3$hotspot<-0
# the significant clusters
data3$hotspot[data3$CLUSTER>0]<-1
#names(data3)[1]<-"LOC_ID"
data4<-data3[,c("LOC_ID","hotspot")]; 
data5<-merge(data4, homesteads, by="LOC_ID")
names(data5)[1:2]<-c("location.id", "hotspots_satO")
# from observed
hotspots<-read.table("hotspots_R.txt", header = T)
hotnew<-merge(hotspots, data5, by=c("location.id", "locx", "locy") )
write.table(hotnew,"hotspots_river.txt", row.names = F)

############################
#short hotspot
homesteads<-read.table("homesteads.txt", header = T)
# after obtaining the results from the saTScan model
data<-read.table("res0_sh.gis.txt", header = T, sep="")
data<-data[data$P_VALUE<0.05,]
data1<-read.table("res0_sh.col.txt", header = T, sep="")
data1<-data1[data1$P_VALUE<0.05,]
data2<-merge(data, data1, by=c("LOC_ID", "CLUSTER"), all = T)
colnames(homesteads)[1]<-c("LOC_ID")
data3<-merge(data2, homesteads, by=c("LOC_ID"), all=T)
data3$REL_RISK[is.na(data3$REL_RISK)]<-1.0
data3<-data3[!duplicated(data3$LOC_ID),]
data3$high[data3$REL_RISK>1.0]<-1
data3$high[data3$REL_RISK<1.0]<-0
data3$high[is.na(data3$high)]<-0
data3$CLUSTER[is.na(data3$CLUSTER)]<-0
table(data3$CLUSTER)
data3$hotspot<-0
# the significant clusters
data3$hotspot[data3$CLUSTER>0]<-1
#names(data3)[1]<-"LOC_ID"
data4<-data3[,c("LOC_ID","hotspot")]; 
data5<-merge(data4, homesteads, by="LOC_ID")
names(data5)[1:2]<-c("location.id", "hotspots_satO")
# from observed
hotspots<-read.table("hotspots.txt", header = T)
hotnew<-merge(hotspots, data5, by=c("location.id", "locx", "locy") )
write.table(hotnew,"hotspots_short.txt", row.names = F)

############################
#long hotspot
homesteads<-read.table("homesteads.txt", header = T)
# after obtaining the results from the saTScan model
data<-read.table("res0_lo.gis.txt", header = T, sep="")
data<-data[data$P_VALUE<0.05,]
data1<-read.table("res0_lo.col.txt", header = T, sep="")
data1<-data1[data1$P_VALUE<0.05,]
data2<-merge(data, data1, by=c("LOC_ID", "CLUSTER"), all = T)
colnames(homesteads)[1]<-c("LOC_ID")
data3<-merge(data2, homesteads, by=c("LOC_ID"), all=T)
data3$REL_RISK[is.na(data3$REL_RISK)]<-1.0
data3<-data3[!duplicated(data3$LOC_ID),]
data3$high[data3$REL_RISK>1.0]<-1
data3$high[data3$REL_RISK<1.0]<-0
data3$high[is.na(data3$high)]<-0
data3$CLUSTER[is.na(data3$CLUSTER)]<-0
table(data3$CLUSTER)
data3$hotspot<-0
# the significant clusters
data3$hotspot[data3$CLUSTER>0]<-1
#names(data3)[1]<-"LOC_ID"
data4<-data3[,c("LOC_ID","hotspot")]; 
data5<-merge(data4, homesteads, by="LOC_ID")
names(data5)[1:2]<-c("location.id", "hotspots_satO")
# from observed
hotspots<-read.table("hotspots.txt", header = T)
hotnew<-merge(hotspots, data5, by=c("location.id", "locx", "locy") )
write.table(hotnew,"hotspots_long.txt", row.names = F)

############################
#seasonal hotspot
homesteads<-read.table("homesteads.txt", header = T)
# after obtaining the results from the saTScan model
data<-read.table("res0_ss.gis.txt", header = T, sep="")
data<-data[data$P_VALUE<0.05,]
data1<-read.table("res0_ss.col.txt", header = T, sep="")
data1<-data1[data1$P_VALUE<0.05,]
data2<-merge(data, data1, by=c("LOC_ID", "CLUSTER"), all = T)
colnames(homesteads)[1]<-c("LOC_ID")
data3<-merge(data2, homesteads, by=c("LOC_ID"), all=T)
data3$REL_RISK[is.na(data3$REL_RISK)]<-1.0
data3<-data3[!duplicated(data3$LOC_ID),]
data3$high[data3$REL_RISK>1.0]<-1
data3$high[data3$REL_RISK<1.0]<-0
data3$high[is.na(data3$high)]<-0
data3$CLUSTER[is.na(data3$CLUSTER)]<-0
table(data3$CLUSTER)
data3$hotspot<-0
# the significant clusters
data3$hotspot[data3$CLUSTER>0]<-1
#names(data3)[1]<-"LOC_ID"
data4<-data3[,c("LOC_ID","hotspot")]; 
data5<-merge(data4, homesteads, by="LOC_ID")
names(data5)[1:2]<-c("location.id", "hotspots_satO")
# from observed
hotspots<-read.table("hotspots.txt", header = T)
hotnew<-merge(hotspots, data5, by=c("location.id", "locx", "locy") )
write.table(hotnew,"hotspots_ss.txt", row.names = F)


# ELLIPTIC HOTSPOT
###########################################3
#single hotspot
library(dplyr)
homesteads<-read.table("homesteads.txt", header = T)
# after obtaining the results from the saTScan model
data<-read.table("resE0.gis.txt", header = T, sep="")
data<-data[data$P_VALUE<0.05,]
data1<-read.table("resE0.col.txt", header = T, sep="")
data1<-data1[data1$P_VALUE<0.05,]
data2<-merge(data, data1, by=c("LOC_ID", "CLUSTER"), all = T)
colnames(homesteads)[1]<-c("LOC_ID")
data3<-merge(data2, homesteads, by=c("LOC_ID"), all=T)
data3$REL_RISK[is.na(data3$REL_RISK)]<-1.0
data3<-data3[!duplicated(data3$LOC_ID),]
data3$high[data3$REL_RISK>1.0]<-1
data3$high[data3$REL_RISK<1.0]<-0
data3$high[is.na(data3$high)]<-0
data3$CLUSTER[is.na(data3$CLUSTER)]<-0
table(data3$CLUSTER)
data3$hotspot<-0
# the significant clusters
data3$hotspot[data3$CLUSTER>0]<-1
#names(data3)[1]<-"LOC_ID"
data4<-data3[,c("LOC_ID","hotspot")]; 
data5<-merge(data4, homesteads, by="LOC_ID")
names(data5)[1:2]<-c("location.id", "hotspots_satE")
# from observed
hotspots<-read.table("hotspots_single.txt", header = T)
hotnew<-merge(hotspots, data5, by=c("location.id", "locx", "locy") )
write.table(hotnew,"hotspots_single.txt", row.names = F)

############################
#multiple hotspot
homesteads<-read.table("homesteads.txt", header = T)
# after obtaining the results from the saTScan model
data<-read.table("resE0_M.gis.txt", header = T, sep="")
data<-data[data$P_VALUE<0.05,]
data1<-read.table("resE0_M.col.txt", header = T, sep="")
data1<-data1[data1$P_VALUE<0.05,]
data2<-merge(data, data1, by=c("LOC_ID", "CLUSTER"), all = T)
colnames(homesteads)[1]<-c("LOC_ID")
data3<-merge(data2, homesteads, by=c("LOC_ID"), all=T)
data3$REL_RISK[is.na(data3$REL_RISK)]<-1.0
data3<-data3[!duplicated(data3$LOC_ID),]
data3$high[data3$REL_RISK>1.0]<-1
data3$high[data3$REL_RISK<1.0]<-0
data3$high[is.na(data3$high)]<-0
data3$CLUSTER[is.na(data3$CLUSTER)]<-0
table(data3$CLUSTER)
data3$hotspot<-0
# the significant clusters
data3$hotspot[data3$CLUSTER>0]<-1
#names(data3)[1]<-"LOC_ID"
data4<-data3[,c("LOC_ID","hotspot")]; 
data5<-merge(data4, homesteads, by="LOC_ID")
names(data5)[1:2]<-c("location.id", "hotspots_satE")
# from observed
hotspots<-read.table("hotspots_multiple.txt", header = T)
hotnew<-merge(hotspots, data5, by=c("location.id","locx","locy") )
write.table(hotnew,"hotspots_multiple.txt", row.names = F)

############################
#decay hotspot
homesteads<-read.table("homesteads.txt", header = T)
# after obtaining the results from the saTScan model
data<-read.table("resE0_D.gis.txt", header = T, sep="")
data<-data[data$P_VALUE<0.05,]
data1<-read.table("resE0_D.col.txt", header = T, sep="")
data1<-data1[data1$P_VALUE<0.05,]
data2<-merge(data, data1, by=c("LOC_ID", "CLUSTER"), all = T)
colnames(homesteads)[1]<-c("LOC_ID")
data3<-merge(data2, homesteads, by=c("LOC_ID"), all=T)
data3$REL_RISK[is.na(data3$REL_RISK)]<-1.0
data3<-data3[!duplicated(data3$LOC_ID),]
data3$high[data3$REL_RISK>1.0]<-1
data3$high[data3$REL_RISK<1.0]<-0
data3$high[is.na(data3$high)]<-0
data3$CLUSTER[is.na(data3$CLUSTER)]<-0
table(data3$CLUSTER)
data3$hotspot<-0
# the significant clusters
data3$hotspot[data3$CLUSTER>0]<-1
#names(data3)[1]<-"LOC_ID"
data4<-data3[,c("LOC_ID","hotspot")]; 
data5<-merge(data4, homesteads, by="LOC_ID")
names(data5)[1:2]<-c("location.id", "hotspots_satE")
# from observed
hotspots<-read.table("hotspots_decay.txt", header = T)
hotnew<-merge(hotspots, data5, by=c("location.id", "locx", "locy") )
write.table(hotnew,"hotspots_decay.txt", row.names = F)

############################
#river hotspot
homesteads<-read.table("homesteads.txt", header = T)
# after obtaining the results from the saTScan model
data<-read.table("resE0_R.gis.txt", header = T, sep="")
data<-data[data$P_VALUE<0.05,]
data1<-read.table("resE0_R.col.txt", header = T, sep="")
data1<-data1[data1$P_VALUE<0.05,]
data2<-merge(data, data1, by=c("LOC_ID", "CLUSTER"), all = T)
colnames(homesteads)[1]<-c("LOC_ID")
data3<-merge(data2, homesteads, by=c("LOC_ID"), all=T)
data3$REL_RISK[is.na(data3$REL_RISK)]<-1.0
data3<-data3[!duplicated(data3$LOC_ID),]
data3$high[data3$REL_RISK>1.0]<-1
data3$high[data3$REL_RISK<1.0]<-0
data3$high[is.na(data3$high)]<-0
data3$CLUSTER[is.na(data3$CLUSTER)]<-0
table(data3$CLUSTER)
data3$hotspot<-0
# the significant clusters
data3$hotspot[data3$CLUSTER>0]<-1
#names(data3)[1]<-"LOC_ID"
data4<-data3[,c("LOC_ID","hotspot")]; 
data5<-merge(data4, homesteads, by="LOC_ID")
names(data5)[1:2]<-c("location.id", "hotspots_satE")
# from observed
hotspots<-read.table("hotspots_river.txt", header = T)
hotnew<-merge(hotspots, data5, by=c("location.id", "locx", "locy") )
write.table(hotnew,"hotspots_river.txt", row.names = F)

############################
#short hotspot
homesteads<-read.table("homesteads.txt", header = T)
# after obtaining the results from the saTScan model
data<-read.table("resE0_sh.gis.txt", header = T, sep="")
data<-data[data$P_VALUE<0.05,]
data1<-read.table("resE0_sh.col.txt", header = T, sep="")
data1<-data1[data1$P_VALUE<0.05,]
data2<-merge(data, data1, by=c("LOC_ID", "CLUSTER"), all = T)
colnames(homesteads)[1]<-c("LOC_ID")
data3<-merge(data2, homesteads, by=c("LOC_ID"), all=T)
data3$REL_RISK[is.na(data3$REL_RISK)]<-1.0
data3<-data3[!duplicated(data3$LOC_ID),]
data3$high[data3$REL_RISK>1.0]<-1
data3$high[data3$REL_RISK<1.0]<-0
data3$high[is.na(data3$high)]<-0
data3$CLUSTER[is.na(data3$CLUSTER)]<-0
table(data3$CLUSTER)
data3$hotspot<-0
# the significant clusters
data3$hotspot[data3$CLUSTER>0]<-1
#names(data3)[1]<-"LOC_ID"
data4<-data3[,c("LOC_ID","hotspot")]; 
data5<-merge(data4, homesteads, by="LOC_ID")
names(data5)[1:2]<-c("location.id", "hotspots_satE")
# from observed
hotspots<-read.table("hotspots_short.txt", header = T)
hotnew<-merge(hotspots, data5, by=c("location.id", "locx", "locy") )
write.table(hotnew,"hotspots_short.txt", row.names = F)

############################
#long hotspot
homesteads<-read.table("homesteads.txt", header = T)
# after obtaining the results from the saTScan model
data<-read.table("resE0_lo.gis.txt", header = T, sep="")
data<-data[data$P_VALUE<0.05,]
data1<-read.table("resE0_lo.col.txt", header = T, sep="")
data1<-data1[data1$P_VALUE<0.05,]
data2<-merge(data, data1, by=c("LOC_ID", "CLUSTER"), all = T)
colnames(homesteads)[1]<-c("LOC_ID")
data3<-merge(data2, homesteads, by=c("LOC_ID"), all=T)
data3$REL_RISK[is.na(data3$REL_RISK)]<-1.0
data3<-data3[!duplicated(data3$LOC_ID),]
data3$high[data3$REL_RISK>1.0]<-1
data3$high[data3$REL_RISK<1.0]<-0
data3$high[is.na(data3$high)]<-0
data3$CLUSTER[is.na(data3$CLUSTER)]<-0
table(data3$CLUSTER)
data3$hotspot<-0
# the significant clusters
data3$hotspot[data3$CLUSTER>0]<-1
#names(data3)[1]<-"LOC_ID"
data4<-data3[,c("LOC_ID","hotspot")]; 
data5<-merge(data4, homesteads, by="LOC_ID")
names(data5)[1:2]<-c("location.id", "hotspots_satE")
# from observed
hotspots<-read.table("hotspots_long.txt", header = T)
hotnew<-merge(hotspots, data5, by=c("location.id", "locx", "locy") )
write.table(hotnew,"hotspots_long.txt", row.names = F)

############################
#long hotspot
homesteads<-read.table("homesteads.txt", header = T)
# after obtaining the results from the saTScan model
data<-read.table("resE0_ss.gis.txt", header = T, sep="")
data<-data[data$P_VALUE<0.05,]
data1<-read.table("resE0_ss.col.txt", header = T, sep="")
data1<-data1[data1$P_VALUE<0.05,]
data2<-merge(data, data1, by=c("LOC_ID", "CLUSTER"), all = T)
colnames(homesteads)[1]<-c("LOC_ID")
data3<-merge(data2, homesteads, by=c("LOC_ID"), all=T)
data3$REL_RISK[is.na(data3$REL_RISK)]<-1.0
data3<-data3[!duplicated(data3$LOC_ID),]
data3$high[data3$REL_RISK>1.0]<-1
data3$high[data3$REL_RISK<1.0]<-0
data3$high[is.na(data3$high)]<-0
data3$CLUSTER[is.na(data3$CLUSTER)]<-0
table(data3$CLUSTER)
data3$hotspot<-0
# the significant clusters
data3$hotspot[data3$CLUSTER>0]<-1
#names(data3)[1]<-"LOC_ID"
data4<-data3[,c("LOC_ID","hotspot")]; 
data5<-merge(data4, homesteads, by="LOC_ID")
names(data5)[1:2]<-c("location.id", "hotspots_satE")
# from observed
hotspots<-read.table("hotspots_ss.txt", header = T)
hotnew<-merge(hotspots, data5, by=c("location.id", "locx", "locy") )
write.table(hotnew,"hotspots_ss.txt", row.names = F)


#Draw the graphs for original and simulated
par(mfrow=c(3,5))

hotspots<-read.table("hotspots_single.txt", header = T)
plot(hotspots$locx, xlab="",ylab="", hotspots$locy, col=as.factor(hotspots$hotspot), main="single hotspot", pch=20,cex=1.5)
plot(hotspots$locx, xlab="",ylab="", hotspots$locy, col=as.factor(hotspots$hotspots_satO), pch=20,cex=1.5)
plot(hotspots$locx, hotspots$locy, col=as.factor(hotspots$hotspots_satE), pch=20,cex=1.5)
plot(hotspots$locx, hotspots$locy, col=as.factor(hotspots$hotspots_flexO), pch=20,cex=1.5)
plot(hotspots$locx, hotspots$locy, col=as.factor(hotspots$hotspots_mc), pch=20,cex=1.5)

hotspots<-read.table("hotspots_multiple.txt", header = T)
plot(hotspots$locx, hotspots$locy, col=as.factor(hotspots$hotspot), main="multiple hotspots", pch=20,cex=1.5)
plot(hotspots$locx, xlab="",ylab="", hotspots$locy, col=as.factor(hotspots$hotspots_satO), pch=20,cex=1.5)
plot(hotspots$locx, hotspots$locy, col=as.factor(hotspots$hotspots_satE), pch=20,cex=1.5)
plot(hotspots$locx, hotspots$locy, col=as.factor(hotspots$hotspots_flexO), pch=20,cex=1.5)
plot(hotspots$locx, hotspots$locy, col=as.factor(hotspots$hotspots_mc), pch=20,cex=1.0)

hotspots<-read.table("hotspots_decay.txt", header = T)
plot(hotspots$locx, hotspots$locy, col=as.factor(hotspots$hotspot), main="decay gradient", pch=20,cex=1.5)
plot(hotspots$locx, hotspots$locy, col=as.factor(hotspots$hotspots_satO), pch=20,cex=1.5)
plot(hotspots$locx, hotspots$locy, col=as.factor(hotspots$hotspots_satE), pch=20,cex=1.5)
plot(hotspots$locx, hotspots$locy, col=as.factor(hotspots$hotspots_flexO), pch=20,cex=1.5)
plot(hotspots$locx, hotspots$locy, col=as.factor(hotspots$hotspots_mc), pch=20,cex=1.0)

hotspots<-read.table("hotspots_river.txt", header = T)
plot(hotspots$locx, hotspots$locy, col=as.factor(hotspots$hotspot), main="irregular hotspot", pch=20,cex=1.5)
plot(hotspots$locx, hotspots$locy, col=as.factor(hotspots$hotspots_satO), pch=20,cex=1.5)
plot(hotspots$locx, hotspots$locy, col=as.factor(hotspots$hotspots_satE), pch=20,cex=1.5)
plot(hotspots$locx, hotspots$locy, col=as.factor(hotspots$hotspots_flexO), pch=20,cex=1.5)
plot(hotspots$locx, hotspots$locy, col=as.factor(hotspots$hotspots_mc), pch=20,cex=1.0)

hotspots<-read.table("hotspots_short.txt", header = T)
plot(hotspots$locx, hotspots$locy, col=as.factor(hotspots$hotspot), main="short distance", pch=20,cex=1.5)
plot(hotspots$locx, hotspots$locy, col=as.factor(hotspots$hotspots_satO), pch=20,cex=1.5)
plot(hotspots$locx, hotspots$locy, col=as.factor(hotspots$hotspots_satE), pch=20,cex=1.5)
plot(hotspots$locx, hotspots$locy, col=as.factor(hotspots$hotspots_flexO), pch=20,cex=1.5)
plot(hotspots$locx, hotspots$locy, col=as.factor(hotspots$hotspots_mc), pch=20,cex=1.0)

hotspots<-read.table("hotspots_long.txt", header = T)
plot(hotspots$locx, hotspots$locy, col=as.factor(hotspots$hotspot), main="long distance", pch=20,cex=1.5)
plot(hotspots$locx, hotspots$locy, col=as.factor(hotspots$hotspots_satO), pch=20,cex=1.5)
plot(hotspots$locx, hotspots$locy, col=as.factor(hotspots$hotspots_satE), pch=20,cex=1.5)
plot(hotspots$locx, hotspots$locy, col=as.factor(hotspots$hotspots_flexO), pch=20,cex=1.5)
plot(hotspots$locx, hotspots$locy, col=as.factor(hotspots$hotspots_mc), pch=20,cex=1.0)
