#data preparation
#set seed in the code files R folder - for different replicates


#________________________________________________________________________________________________________________
##############  Function to check  new infections against current infections ##########################################################
#_______________and delete if necessary_______________________________________________________________________________

checkNewInfections<-function(newInf) { 
  
  # newInf<-data.frame(newInf)
  numberNew<-length(newInf[,1])
  if (numberNew>0) {  
    # check total number of current infections per person per house
    for (n3 in 1:length(newInf[,1])) {
      currentMatch<-current_genotypes[((current_genotypes[,"location.id"]==newInf[n3,"location.id"]) & (current_genotypes[,"pID"]==newInf[n3,"pID"])),]
      if ( length(currentMatch[,1])>numCurrentInfAllowed )  {
        newInf[n3,"clear"]<-1
      }
    }
    
  } # closes if numberNew>0
  
  # remove marked new infections 
  if (is.array(newInf)) {
    newInf<-newInf[newInf[,"clear"]==0,]
  }
  
  return(newInf)
}

# ________________________________________________________________________________________________________________
##############  Function to get new infections per timestep  ###############################################################
#________________________________________________________________________________________________________________
#********************
getNewInfections<-function(t) {
  currentNumInf<-length(current_genotypes[,1])
  # dummy row for newInf, omit later after loop
  newInf<-current_genotypes[1,]  
  
  for (i in 1:currentNumInf){
    # number of new infections for current infection
    meanNew<-meanNumNewInfPerTimestep
    
    #5% for the infections at the edge of the grid
    if ((current_genotypes[i, "locx"]<0.5) | (current_genotypes[i, "locx"]>9.5) | (current_genotypes[i, "locy"]<0.5) | (current_genotypes[i, "locy"]>9.5)){
      meanNew<-(meanNumNewInfPerTimestep/2)
    }
    
    currenttime<-current_genotypes[i,"time"]
    ageInf<- t - currenttime
    
    currenthome<-current_genotypes[i,"location.id"]
    currentperson<-current_genotypes[i,"pID"]	
    
    for (q in 1:nrow(infHousesHot)){    
      if (currenthome == infHousesHot[q,"location.id"]) {
        meanNew<-meanNew * hotspotFactor
      }
    }
    
    #infectivity over time - decays with an exponential function (an individual is fully infectious for "infsTime" days)
    infsTime<-75
    rho<-1/infsTime
    
    meanNew[currenthome]<-meanNew*exp(-rho*ageInf)
    #meanNew[currenthome]<-meanNew[currenthome]*(dist_kernel[currenthome,1])/max(dist_kernel$sum_kernel)
    
    numNewInfPerTimestep<-rpois(1,meanNew[currenthome])
    
    # assign new infections
    if (numNewInfPerTimestep>0){
      for (j in 1:numNewInfPerTimestep) {
        newInf<-rbind(newInf,current_genotypes[i,])
        
        # get new locations
        newInf[length(newInf[,1]), "location.id"]<-getNewLoc(currenthome)
        
        #condition for person ID - if it is in the same house at the same time then a different person is infected
        maxNumpplehh<-8
        numPple<-1:maxNumpplehh
        
        if (newInf[length(newInf[,1]), "location.id"]==currenthome & newInf[length(newInf[,1]), "time"]==currenttime){
          numPple<-subset(numPple, numPple!=currentperson)}
        
        newInf[length(newInf[,1]), "locx"]<-homesteads$locx[homesteads$id==newInf[length(newInf[,1]), "location.id"]]
        newInf[length(newInf[,1]), "locy"]<-homesteads$locy[homesteads$id==newInf[length(newInf[,1]), "location.id"]]
        
        if (newInf[length(newInf[,1]), "location.id"]!=currenthome & newInf[length(newInf[,1]), "time"]==currenttime){
          numPple<-1:maxNumpplehh}
        
        newInf[length(newInf[,1]), "locx"]<-homesteads$locx[homesteads$id==newInf[length(newInf[,1]), "location.id"]]
        newInf[length(newInf[,1]), "locy"]<-homesteads$locy[homesteads$id==newInf[length(newInf[,1]), "location.id"]]
        
        newInf[length(newInf[,1]), "pID"]<-sample(numPple,1)
        
        newInf[,"status"]<-1
        
        #				set clearance to 0 for new infections
        newInf[length(newInf[,1]),"clear"]<-0 
        # get timestep
        newInf[,"time"]<-rep(t,length(newInf[,1]))
        # col1=clear, col2=time, col3=locx, col4=locy, col5-104=snps 
      }
    }
  } # closes i loop
  
  # tidy new infections up
  # remove first dummy row of newInf
  # added 'if' to prevent errors: todo - see what happened
  if (is.array(newInf)) {
    newInf<-newInf[-1,]
  }  
  
  return(newInf)
} 


#________________________________________________________________________________________________________________
##############  Function to effect select new household  ###############################################################
#________________________________________________________________________________________________________________
getNewLoc<-function(currenthome) {
  specificHomes <- dist_diff[ which(dist_diff$id==currenthome | dist_diff$id_x ==currenthome), ]
  specificHomes[,"scaledProp"]<-specificHomes[,"prop"]/sum(specificHomes[,"prop"])
  specificHomes <- specificHomes[order(-specificHomes[,"scaledProp"]),]
  specificHomes[,"cumProp"]<-cumsum(specificHomes[,"scaledProp"])
  rand_new<-runif(1)
  
  specificHomes<-specificHomes[c(1,1:nrow(specificHomes)),]
  specificHomes[1,"cumProp"]<-0
  
  smallHomes<-specificHomes[specificHomes[,"cumProp"]<rand_new,]
  homePair<-smallHomes[nrow(smallHomes),1:2]
  
  newloc<-homePair[1]
  if (newloc==currenthome) newloc<-homePair[2]
  
  return(newloc)
}


# __________________________________________________________________________
#
# getPotentialImportedInfections
# function to create a set of imported infections which can potetnially be added in to the simulated infections
#*************

getPotentialImportedInfections<-function() {
  
  nImp<-2
  numSNPsPerInf<-1
  
  potAddInf <- array(-9, dim=c(nImp, numSNPsPerInf))

 # for (i in 1:nImp) {
   # for (j in 1:numSNPsPerInf) {
    #  potAddInf[i,j]<-getBaselineSNPs(domAlleleFreqList[j]) 
    #}
  #}
  # add also random locations 1-200 in homesteads.txt
  # to do: make this generalisable
  
  potAdd<-sample(homesteads[,"id"], nImp, replace=TRUE)
  potAdd<-homesteads[potAdd,]
  potAdd$pID<-sample(1:maxNumpplehh, nImp, replace=T)
  colnames(potAdd)<-c("potAddH", "locx", "locy","pID")
  
  # add a list of random times (0-2500) 
  potAddt <- sample(1:400, nImp, replace=TRUE) #time changed t reflect new timesteps in obsInf
  potAddt <- sort(potAddt)
  # make last time very long to avoid missing value problem
  potAddt[length(potAddt)]<-999999
  clear<- rep(0,nImp)
  numCommInf <- rep(0,nImp)
  
  potAddInf <- cbind(potAdd$potAddH, potAdd$pID, clear, potAddt, numCommInf, potAdd$locx, potAdd$locy, potAddInf)
  potAddInf<-data.frame(potAddInf)
  names(potAddInf)[1]<-"potAddH";names(potAddInf)[2]<-"pID";names(potAddInf)[6]<-"locx";names(potAddInf)[7]<-"locy"
  return(potAddInf)
}

# _______________________________________________________________________
# setUpInitialCurrentGenotypes
# function to set up initial current_genotypes from read-in initial infections
#*************
numPreSNPcol<-7
setUpInitialCurrentGenotypes<-function() {
  num.init<-dim(infHousesInitInf)[1]
  clear<-rep(0,num.init)
  time<-rep(0,num.init)
  numCommInf<-rep(num.init, num.init)
  current_genotypes <- array(-9, dim=c(length(infHousesInitInf[,1]), numSNPsPerInf+7))
  current_genotypes[,1] <- infHousesInitInf[,"location.id"]
  current_genotypes[,2] <- infHousesInitInf[,"pID"]
  current_genotypes[,3] <- clear
  current_genotypes[,4] <- time
  current_genotypes[,5] <- numCommInf
  current_genotypes[,6] <- infHousesInitInf[,"locx"]
  current_genotypes[,7] <- infHousesInitInf[,"locy"]
  #for (i in 1:numSNPsPerInf) {current_genotypes[,(i+7)] <- infHousesInitInf[,(i+4)]}
  current_genotypes<-data.frame(current_genotypes)
  current_genotypes[,numSNPsPerInf+7]<-1
  colnames(current_genotypes)[1:8]<-c("location.id","pID","clear","time","numCommInf", "locx", "locy","status")
  #firstSNPcol<-numPreSNPcol+1
  #current_genotypes<-addGenoID(current_genotypes)
  return(current_genotypes)
}


#________________________________________________________________________________________________________________
##############   updateCurrentGenotypes for one time step     ##############################################
#________________________________________________________________________________________________________________

updateCurrentGenotypes<-function(t){
  
  # mark for clearance
  for (md in 1:length(current_genotypes[,1])) {    
    rand<-runif(1)
    if (rand<probClearancePerTimestep){current_genotypes[md,"clear"]<-1}
    if (t-current_genotypes[md,"time"]>maxDuration) {current_genotypes[md,"clear"]<-1}    
  }
  
  # clear doomed infections 
  current_genotypes<-current_genotypes[current_genotypes[,"clear"]==0,]
  
  # make new "current" dataset
  current_genotypes<-rbind(current_genotypes,newInf)
  
  # to speedup the iteration set ceiling
  if (length(current_genotypes[,1])>25000){
   nt<-seq(1:length(current_genotypes[,1]))
   nt <- sample(nt, 25000,replace=FALSE)
    current_genotypes<-current_genotypes[nt,]
  cat("ceiling reached at time",t,"\n")
  }
  
  return(current_genotypes)
}

# End of functionsHotspots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Load required libraries
library(SpatialEpi);library(data.table); require(utils)
library(plyr);library(dplyr); library(dtplyr)
library(foreign);library("scales")

############################################################################################\
rate<-0.314

############
setwd(getwd())
# HOUSE-PAIRS AND DISTANCES 
# allpairs is a dataset with simulated house locations and obsInf locations and distances between every possible pair of house-obsInf
# to find closest house with an infection at time t 
# allpairs<-read.table("allpairs_dist.txt", header=TRUE, row.names=NULL)

# HOMESTEADS
# simulate homesteads and calculate distances for infection 'movement', if not previously done and no homestead data
num.homesteads<-2000

# maximum number of people in a household (average)
maxNumpplehh<-8

# the factor for being in a hotpot

hotspotFactor<-1.2

#read in simulated homesteads
homesteads<-read.table("homesteads.txt", header=T, row.names = NULL)
colnames(homesteads)[1:3]<-c("id","locx","locy")

# source("../tempSimHomesteads.r")
dist_diff<-read.table("dist_diff.txt", header=TRUE, row.names=NULL)

#rate<-4.0
mu<-1/(2*(rate^2))
dist_diff[,"prop"]<-exp(-mu*(dist_diff[,"dist.diff"]^2)) #the normal kernel

#seonality/EIR data
#eir_data<-read.table("eir_data.txt", header=T, row.names = NULL)

# number of SNPs
numSNPsPerInf<-1

# number of total infections allowed per house
numCurrentInfAllowed<-8

# max duration of infections (otherwise exponential) todo ev: assign duration in getNewInfections 
maxDuration<-4000

# initial infections: version with structured initial infections

infHousesInitInf<-read.table("initialInfs.txt",header=TRUE)
colnames(infHousesInitInf)[1:4]<-c("location.id","locx","locy", "pID")

# read files with hotspot homes

infHousesHot<-read.table("hotHomes.txt",header=TRUE)
colnames(infHousesHot)[1:4]<-c("locx","locy","location.id", "hotHme")

current_genotypes<-setUpInitialCurrentGenotypes() 

# set up array of potential imported infections to be added in to the simulation
potAddInf<-getPotentialImportedInfections()

# UPDATE TIMESTEPS 
currentNumInf<-length(current_genotypes[,1])
meanNumNewInfPerTimestep<-0.036
probClearancePerTimestep<-0.025
pt<-1

# wamr-up period to get a substantial number of infections
for (t in 1:146) {  
  #meanNumNewInfPerTimestep<-meanNumNewInfPerTimestep*eir_data$EIRscaled[eir_data$time==t]
  newInf<-getNewInfections(t) 
  #newInf<-checkNewInfections(newInf)
  current_genotypes<-updateCurrentGenotypes(t)
  if (t==potAddInf$potAddt[pt]) {
    potImport<-potAddInf[pt,]; names(potImport)<-names(current_genotypes)
    current_genotypes<-rbind(current_genotypes,potImport)
    #currently  included to test 
    pt<-pt+1
  }
  
}

meanNumNewInfPerTimestep<-0.0354
probClearancePerTimestep<-0.025

# running the model

for (t in 147:219) {  
  #meanNumNewInfPerTimestep<-meanNumNewInfPerTimestep*eir_data$EIRscaled[eir_data$time==t]
  newInf<-getNewInfections(t) 
  #newInf<-checkNewInfections(newInf)
  current_genotypes<-updateCurrentGenotypes(t)
  if (t==potAddInf$potAddt[pt]) {
    potImport<-potAddInf[pt,]; names(potImport)<-names(current_genotypes)
    current_genotypes<-rbind(current_genotypes,potImport)
   #currently  included to test 
  pt<-pt+1
  }
  if (t==165) {
    # write to text file for investigations
    infs0<-current_genotypes
  }
  
  if (t==218) {
    # write to text file for investigations
    infs1<-current_genotypes
  }

}

#############################
#output files
write.table(infs0,"crossT0.txt", row.names = F)
write.table(infs1,"crossT1.txt", row.names = F)

meanNumNewInfPerTimestep<-0.0352
probClearancePerTimestep<-0.025

for (t in 220:366) {  
  #meanNumNewInfPerTimestep<-meanNumNewInfPerTimestep*eir_data$EIRscaled[eir_data$time==t]
  newInf<-getNewInfections(t) 
  #newInf<-checkNewInfections(newInf)
  current_genotypes<-updateCurrentGenotypes(t)
  if (t==potAddInf$potAddt[pt]) {
    potImport<-potAddInf[pt,]; names(potImport)<-names(current_genotypes)
    current_genotypes<-rbind(current_genotypes,potImport)
   #currently  included to test 
  pt<-pt+1
  }
  
  if (t==293) {
    # write to text file for investigations
    infs2<-current_genotypes
  }
  
  if (t==365) {
    # write to text file for investigations
    infs3<-current_genotypes
  }
}

write.table(infs2,"crossT2.txt", row.names = F)
write.table(infs3,"crossT3.txt", row.names = F)


#############:::::::::::::::::::::::::::::::::::::::::::::::::::::::#####################
#############:::::::::::::::::::::::::::::::::::::::::::::::::::::::#####################

# RUN IN THE SIMULATION AND OBTAIN CROSSECTIONAL FILES 
library(ggplot2)
library(RColorBrewer)
homesteads<-read.table("homesteads.txt", header=T, row.names = NULL)
names(homesteads)[1]<-"location.id"
numPple<-8

# create case files from both the crossectional surveys
#3 months
cases2<-infs0
# remove duplicates
cases3<-cases2[!duplicated(cases2[,c("location.id", "pID","locx","locy")]),]
#plot(cases3$locx, cases3$locy, col="darkblue", pch=20, cex=2, main="cases")
#save case file
cases3$status<-1
cases4<-subset(cases3, select = c(location.id, pID, locx, locy, status))
write.table(cases4,"casesT03.txt", row.names = F)

# get population estimates and create control files
populate <- homesteads %>% slice(rep(1:n(), each = numPple))
names(populate)[1]<-"location.id"
setDT(populate)[, pID := seq_len(.N), by = location.id]
populate<-data.frame(populate)
#cases plus controls
names(homesteads)[1]<-"location.id"
cases6<-anti_join(homesteads, cases4, by = c("location.id"))
cases5<-merge(cases4, populate, by=c("location.id","locx","locy", "pID"), all=T)
cases5$status[is.na(cases5$status)]<-0
write.table(cases5,"case_ctrlT03.txt", row.names = F, col.names = F)

#controls
controls2<-anti_join(populate, cases4, by = c("location.id", "pID"))
controls2$status<-1
write.table(controls2,"ctrlsT03.txt", row.names = F)

#population in each household
populate$num<-1
pop_file<-ddply(populate, .(location.id), summarize,  pop=sum(num))
write.table(pop_file,"pop_file.txt", row.names = F)

# case file for flexscan
cases_total2<-ddply(cases4, .(location.id), summarize,  case=sum(status))
cases_total2<-merge(cases_total2, pop_file, by="location.id", all=T)
cases_total2$case[is.na(cases_total2$case)]<-0
write.table(cases_total2,"cases_popT03.txt", row.names = F, col.names = F)

# create case files from both the crossectional surveys
#1 year
cases2<-infs1
# remove duplicates
cases3<-cases2[!duplicated(cases2[,c("location.id", "pID","locx","locy")]),]
#plot(cases3$locx, cases3$locy, col="darkblue", pch=20, cex=2, main="cases")
#save case file
cases3$status<-1
cases4<-subset(cases3, select = c(location.id, pID, locx, locy, status))
write.table(cases4,"casesT.txt", row.names = F)
# get population estimates and create control files
populate <- homesteads %>% slice(rep(1:n(), each = numPple))
names(populate)[1]<-"location.id"
setDT(populate)[, pID := seq_len(.N), by = location.id]
populate<-data.frame(populate)
#cases plus controls
names(homesteads)[1]<-"location.id"
cases6<-anti_join(homesteads, cases4, by = c("location.id"))
cases5<-merge(cases4, populate, by=c("location.id","locx","locy", "pID"), all=T)
cases5$status[is.na(cases5$status)]<-0
write.table(cases5,"case_ctrlT.txt", row.names = F, col.names = F)
#controls
controls2<-anti_join(populate, cases4, by = c("location.id", "pID"))
controls2$status<-1
write.table(controls2,"ctrlsT.txt", row.names = F)
#populcation in each household
populate$num<-1
pop_file<-ddply(populate, .(location.id), summarize,  pop=sum(num))
write.table(pop_file,"pop_file.txt", row.names = F)
# case file for flexscan
cases_total2<-ddply(cases4, .(location.id), summarize,  case=sum(status))
cases_total2<-merge(cases_total2, pop_file, by="location.id", all=T)
cases_total2$case[is.na(cases_total2$case)]<-0
write.table(cases_total2,"cases_popT.txt", row.names = F, col.names = F)

# create case files from both the crossectional surveys
#3 years
cases2<-infs3
# remove duplicates
cases3<-cases2[!duplicated(cases2[,c("location.id", "pID","locx","locy")]),]
#plot(cases3$locx, cases3$locy, col="darkblue", pch=20, cex=2, main="cases")
#save case file
cases3$status<-1
cases4<-subset(cases3, select = c(location.id, pID, locx, locy, status))
write.table(cases4,"casesT3y.txt", row.names = F)
# get population estimates and create control files
populate <- homesteads %>% slice(rep(1:n(), each = numPple))
names(populate)[1]<-"location.id"
setDT(populate)[, pID := seq_len(.N), by = location.id]
populate<-data.frame(populate)
#cases plus controls
names(homesteads)[1]<-"location.id"
cases6<-anti_join(homesteads, cases4, by = c("location.id"))
cases5<-merge(cases4, populate, by=c("location.id","locx","locy", "pID"), all=T)
cases5$status[is.na(cases5$status)]<-0
write.table(cases5,"case_ctrlT3y.txt", row.names = F, col.names = F)
#controls
controls2<-anti_join(populate, cases4, by = c("location.id", "pID"))
controls2$status<-1
write.table(controls2,"ctrlsT3y.txt", row.names = F)
#populcation in each household
populate$num<-1
pop_file<-ddply(populate, .(location.id), summarize,  pop=sum(num))
write.table(pop_file,"pop_file.txt", row.names = F)
# case file for flexscan
cases_total2<-ddply(cases4, .(location.id), summarize,  case=sum(status))
cases_total2<-merge(cases_total2, pop_file, by="location.id", all=T)
cases_total2$case[is.na(cases_total2$case)]<-0
write.table(cases_total2,"cases_popT3y.txt", row.names = F, col.names = F)

