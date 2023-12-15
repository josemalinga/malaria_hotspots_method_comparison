
rm(list=ls())
library(SpatialEpi)
library(dplyr)
library(tidyr)
library("FRK")
#detach(package:plyr)

set.seed(42)

# Create file names for all replicate folders

# file_names = list.files(path = paste0(getwd()), pattern = "rep*", full.names = FALSE)
# file_names=file_names[file_names!="measure_hotspots.R" & file_names!="rep0" & file_names!="repold"
#                       & file_names!="replicates.txt" & file_names!="measure_hotspots2.R"]
# file_names=file_names[file_names!="rep33" & file_names!="rep63" & file_names!="rep93"
#                       & file_names!="rep62" & file_names!="rep92"]
# length(file_names)
# file_names
# write.table(file_names, paste0(getwd(), "/file_names.txt"), row.names=FALSE, col.names=FALSE)


#for (i in 1:28){
# User 
  user = strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
  
  # Working directory
   setwd(paste0("/scicore/home/penny/",user,"/malaria_hotspots/simulations_nov"))
  
   file_names <- read.table(paste0(getwd(),"/file_names.txt"), header = F)
   for (i in length(file_names)){
     file_names$V1[i] <- as.character(file_names$V1[i] )
   }
   file_names
   
   
   task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
   task_id <- as.numeric(task_id)
  
   experiment <- file_names$V1[task_id]
   experiment
  # experiment <- as.character(experiment)
  exp <- levels(experiment)[task_id]
  exp
  
  # experiment <- paste0("rep",i); experiment
  setwd(paste0(getwd(),"/", exp))
  getwd()
  
  
# SINGLE
exp_all = c("T03","T","T3y")

ncenter = 15
hh <- 2000

exp = exp_all[1]

hot <- read.table("hotspots.txt", header = T)

# load data with cases and population
case_pop <- read.table(paste0("cases_pop", exp, ".txt"))
names(case_pop) <- c("location.id","cases","population")

## Process geographical information and convert to grid
geo <- read.table("homesteads.txt", header = T)
geo1 <- geo[,-1]
names(geo1) <- c("x", "y")

# create data clustered by distance
ncluster <- kmeans(cbind(geo1$x, geo1$y), centers = ncenter)
#plot(geo1$x, geo1$y, col = ncluster$cluster, pch = 20)

clust.mean <- data.frame(ncluster$centers)
names(clust.mean) <- c("x", "y")

clust.all = clust.mean[FALSE,]
clust.all$numclust <- clust.all[,1]

for (i in 1:ncenter){
  geof <- geo1[ncluster$cluster==i,]
  geof <- grid2latlong(geof)
  
  geof$numclust <- i
  
  sb <- geof[1,]
  geof <- rbind(geof,sb)
  
  clust.all <- rbind(clust.all, geof)
}

clust.all$numclust <- as.factor(clust.all$numclust)

#extract cluster allocations
nclust <- data.frame(ncluster$cluster)
nclust$location.id <- seq(1:2000)
nclust2 <- inner_join(nclust, geo, by="location.id")

# merge with cases and population
nclust3 <- inner_join(nclust2, case_pop, by="location.id")

# aggregate cases and population by cluster
nclust4 <- nclust3 %>% group_by(ncluster.cluster) %>%
  summarize(cases = sum(cases), population = sum(population))

## Get aggregated counts of population and cases for each county
population <- tapply(nclust3$population,nclust3$ncluster.cluster,sum)
cases <- tapply(nclust3$cases,nclust3$ncluster.cluster,sum)

## Set Parameters
pop.upper.bound <- 0.15
n.simulations <- 9999
alpha.level <- 0.05
plot <- FALSE

## Kulldorff using Binomial likelihoods
binomial <- kulldorff(geo, cases, population, NULL, pop.upper.bound, n.simulations,
                      alpha.level, plot)
cluster <- binomial$most.likely.cluster$location.IDs.included
cluster

cl_data <- geo[geo$location.id %in% cluster, ]
cl_data$cluster <- 1
names(cl_data)[1] <- "ncluster.cluster"

nclust5 <- left_join(nclust2, cl_data[,-c(2:3)], by="ncluster.cluster")
nclust5$cluster[is.na(nclust5$cluster)] <- 0
nclust5$cluster <- as.factor(nclust5$cluster)

plot(nclust5$locx, nclust5$locy, col=nclust5$cluster, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")
cluster
table(nclust5$cluster)

plot(geo$locx, geo$locy, col="black", pch=1, main="simulated", xlab="",ylab="",xaxt = "n",yaxt = "n")
points(hot$locx, hot$locy, col=hot$hotspot, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")

names(nclust5)[5] <- paste0("kull",exp)
hotnew <- inner_join(hot, nclust5[,c(2,5)], by="location.id")
write.table(hotnew, "hotspots_single.txt")

#at 1 year
exp = exp_all[2]

hot <- read.table("hotspots_single.txt", header = T)

# load data with cases and population
case_pop <- read.table(paste0("cases_pop", exp, ".txt"))
names(case_pop) <- c("location.id","cases","population")

# merge with cases and population
nclust3 <- inner_join(nclust2, case_pop, by="location.id")

# aggregate cases and population by cluster
nclust4 <- nclust3 %>% group_by(ncluster.cluster) %>%
  summarize(cases = sum(cases), population = sum(population))

## Get aggregated counts of population and cases for each county
population <- tapply(nclust3$population,nclust3$ncluster.cluster,sum)
cases <- tapply(nclust3$cases,nclust3$ncluster.cluster,sum)

## Set Parameters
pop.upper.bound <- 0.15
n.simulations <- 9999
alpha.level <- 0.05
plot <- FALSE

## Kulldorff using Binomial likelihoods
binomial <- kulldorff(geo, cases, population, NULL, pop.upper.bound, n.simulations,
                      alpha.level, plot)
cluster <- binomial$most.likely.cluster$location.IDs.included
cluster

cl_data <- geo[geo$location.id %in% cluster, ]
cl_data$cluster <- 1
names(cl_data)[1] <- "ncluster.cluster"

nclust5 <- left_join(nclust2, cl_data[,-c(2:3)], by="ncluster.cluster")
nclust5$cluster[is.na(nclust5$cluster)] <- 0
nclust5$cluster <- as.factor(nclust5$cluster)

#par(mfrow=c(2,2))
plot(nclust5$locx, nclust5$locy, col=nclust5$cluster, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")
cluster
table(nclust5$cluster)

plot(geo$locx, geo$locy, col="black", pch=1, main="simulated", xlab="",ylab="",xaxt = "n",yaxt = "n")
points(hot$locx, hot$locy, col=hot$hotspot, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")

names(nclust5)[5] <- paste0("kull",exp)
hotnew <- inner_join(hot, nclust5[,c(2,5)], by="location.id")
write.table(hotnew, "hotspots_single.txt")

#at 3 year
exp = exp_all[3]

hot <- read.table("hotspots_single.txt", header = T)

# load data with cases and population
case_pop <- read.table(paste0("cases_pop", exp, ".txt"))
names(case_pop) <- c("location.id","cases","population")

# merge with cases and population
nclust3 <- inner_join(nclust2, case_pop, by="location.id")

# aggregate cases and population by cluster
nclust4 <- nclust3 %>% group_by(ncluster.cluster) %>%
  summarize(cases = sum(cases), population = sum(population))

## Get aggregated counts of population and cases for each county
population <- tapply(nclust3$population,nclust3$ncluster.cluster,sum)
cases <- tapply(nclust3$cases,nclust3$ncluster.cluster,sum)

## Set Parameters
pop.upper.bound <- 0.15
n.simulations <- 9999
alpha.level <- 0.05
plot <- FALSE

## Kulldorff using Binomial likelihoods
binomial <- kulldorff(geo, cases, population, NULL, pop.upper.bound, n.simulations,
                      alpha.level, plot)
cluster <- binomial$most.likely.cluster$location.IDs.included
cluster

cl_data <- geo[geo$location.id %in% cluster, ]
cl_data$cluster <- 1
names(cl_data)[1] <- "ncluster.cluster"

nclust5 <- left_join(nclust2, cl_data[,-c(2:3)], by="ncluster.cluster")
nclust5$cluster[is.na(nclust5$cluster)] <- 0
nclust5$cluster <- as.factor(nclust5$cluster)

#par(mfrow=c(2,2))
plot(nclust5$locx, nclust5$locy, col=nclust5$cluster, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")
cluster
table(nclust5$cluster)

plot(geo$locx, geo$locy, col="black", pch=1, main="simulated", xlab="",ylab="",xaxt = "n",yaxt = "n")
points(hot$locx, hot$locy, col=hot$hotspot, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")

names(nclust5)[5] <- paste0("kull",exp)
hotnew <- inner_join(hot, nclust5[,c(2,5)], by="location.id")
write.table(hotnew, "hotspots_single.txt")


# MULTIPLE
exp_all = c("T03","T","T3y")

ncenter = 15

hh <- 2000

exp = exp_all[1]

hot <- read.table("hotspots_M.txt", header = T)

# load data with cases and population
case_pop <- read.table(paste0("cases_pop", exp, ".txt"))
names(case_pop) <- c("location.id","cases","population")

## Process geographical information and convert to grid
geo <- read.table("homesteads.txt", header = T)
geo1 <- geo[,-1]
names(geo1) <- c("x", "y")

# create data clustered by distance
ncluster <- kmeans(cbind(geo1$x, geo1$y), centers = ncenter)
#plot(geo1$x, geo1$y, col = ncluster$cluster, pch = 20)

clust.mean <- data.frame(ncluster$centers)
names(clust.mean) <- c("x", "y")

clust.all = clust.mean[FALSE,]
clust.all$numclust <- clust.all[,1]

for (i in 1:ncenter){
  geof <- geo1[ncluster$cluster==i,]
  geof <- grid2latlong(geof)
  
  geof$numclust <- i
  
  sb <- geof[1,]
  geof <- rbind(geof,sb)
  
  clust.all <- rbind(clust.all, geof)
}

clust.all$numclust <- as.factor(clust.all$numclust)

#extract cluster allocations
nclust <- data.frame(ncluster$cluster)
nclust$location.id <- seq(1:2000)
nclust2 <- inner_join(nclust, geo, by="location.id")

# merge with cases and population
nclust3 <- inner_join(nclust2, case_pop, by="location.id")

# aggregate cases and population by cluster
nclust4 <- nclust3 %>% group_by(ncluster.cluster) %>%
  summarize(cases = sum(cases), population = sum(population))

## Get aggregated counts of population and cases for each county
population <- tapply(nclust3$population,nclust3$ncluster.cluster,sum)
cases <- tapply(nclust3$cases,nclust3$ncluster.cluster,sum)

## Set Parameters
pop.upper.bound <- 0.15
n.simulations <- 9999
alpha.level <- 0.05
plot <- FALSE

## Kulldorff using Binomial likelihoods
binomial <- kulldorff(geo, cases, population, NULL, pop.upper.bound, n.simulations,
                      alpha.level, plot)
cluster <- binomial$most.likely.cluster$location.IDs.included
cluster

cl_data <- geo[geo$location.id %in% cluster, ]
cl_data$cluster <- 1
names(cl_data)[1] <- "ncluster.cluster"

nclust5 <- left_join(nclust2, cl_data[,-c(2:3)], by="ncluster.cluster")
nclust5$cluster[is.na(nclust5$cluster)] <- 0
nclust5$cluster <- as.factor(nclust5$cluster)

# plot(nclust5$locx, nclust5$locy, col=nclust5$cluster, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")
# cluster
# table(nclust5$cluster)
# 
# plot(geo$locx, geo$locy, col="black", pch=1, main="simulated", xlab="",ylab="",xaxt = "n",yaxt = "n")
# points(hot$locx, hot$locy, col=hot$hotspot, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")

names(nclust5)[5] <- paste0("kull",exp)
hotnew <- inner_join(hot, nclust5[,c(2,5)], by="location.id")
write.table(hotnew, "hotspots_multiple.txt")

#at 1 year
exp = exp_all[2]

hot <- read.table("hotspots_multiple.txt", header = T)

# load data with cases and population
case_pop <- read.table(paste0("cases_pop", exp, ".txt"))
names(case_pop) <- c("location.id","cases","population")

# merge with cases and population
nclust3 <- inner_join(nclust2, case_pop, by="location.id")

# aggregate cases and population by cluster
nclust4 <- nclust3 %>% group_by(ncluster.cluster) %>%
  summarize(cases = sum(cases), population = sum(population))

## Get aggregated counts of population and cases for each county
population <- tapply(nclust3$population,nclust3$ncluster.cluster,sum)
cases <- tapply(nclust3$cases,nclust3$ncluster.cluster,sum)

## Set Parameters
pop.upper.bound <- 0.15
n.simulations <- 9999
alpha.level <- 0.05
plot <- FALSE

## Kulldorff using Binomial likelihoods
binomial <- kulldorff(geo, cases, population, NULL, pop.upper.bound, n.simulations,
                      alpha.level, plot)
cluster <- binomial$most.likely.cluster$location.IDs.included
cluster

cl_data <- geo[geo$location.id %in% cluster, ]
cl_data$cluster <- 1
names(cl_data)[1] <- "ncluster.cluster"

nclust5 <- left_join(nclust2, cl_data[,-c(2:3)], by="ncluster.cluster")
nclust5$cluster[is.na(nclust5$cluster)] <- 0
nclust5$cluster <- as.factor(nclust5$cluster)

#par(mfrow=c(2,2))
# plot(nclust5$locx, nclust5$locy, col=nclust5$cluster, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")
# cluster
# table(nclust5$cluster)
# 
# plot(geo$locx, geo$locy, col="black", pch=1, main="simulated", xlab="",ylab="",xaxt = "n",yaxt = "n")
# points(hot$locx, hot$locy, col=hot$hotspot, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")

names(nclust5)[5] <- paste0("kull",exp)
hotnew <- inner_join(hot, nclust5[,c(2,5)], by="location.id")
write.table(hotnew, "hotspots_multiple.txt")

#at 3 year
exp = exp_all[3]

hot <- read.table("hotspots_multiple.txt", header = T)

# load data with cases and population
case_pop <- read.table(paste0("cases_pop", exp, ".txt"))
names(case_pop) <- c("location.id","cases","population")

# merge with cases and population
nclust3 <- inner_join(nclust2, case_pop, by="location.id")

# aggregate cases and population by cluster
nclust4 <- nclust3 %>% group_by(ncluster.cluster) %>%
  summarize(cases = sum(cases), population = sum(population))

## Get aggregated counts of population and cases for each county
population <- tapply(nclust3$population,nclust3$ncluster.cluster,sum)
cases <- tapply(nclust3$cases,nclust3$ncluster.cluster,sum)

## Set Parameters
pop.upper.bound <- 0.15
n.simulations <- 9999
alpha.level <- 0.05
plot <- FALSE

## Kulldorff using Binomial likelihoods
binomial <- kulldorff(geo, cases, population, NULL, pop.upper.bound, n.simulations,
                      alpha.level, plot)
cluster <- binomial$most.likely.cluster$location.IDs.included
cluster

cl_data <- geo[geo$location.id %in% cluster, ]
cl_data$cluster <- 1
names(cl_data)[1] <- "ncluster.cluster"

nclust5 <- left_join(nclust2, cl_data[,-c(2:3)], by="ncluster.cluster")
nclust5$cluster[is.na(nclust5$cluster)] <- 0
nclust5$cluster <- as.factor(nclust5$cluster)

#par(mfrow=c(2,2))
# plot(nclust5$locx, nclust5$locy, col=nclust5$cluster, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")
# cluster
# table(nclust5$cluster)
# 
# plot(geo$locx, geo$locy, col="black", pch=1, main="simulated", xlab="",ylab="",xaxt = "n",yaxt = "n")
# points(hot$locx, hot$locy, col=hot$hotspot, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")

names(nclust5)[5] <- paste0("kull",exp)
hotnew <- inner_join(hot, nclust5[,c(2,5)], by="location.id")
write.table(hotnew, "hotspots_multiple.txt")


# DECAY
exp_all = c("T03","T","T3y")

ncenter = 15

hh <- 2000

exp = exp_all[1]

hot <- read.table("hotspots_D.txt", header = T)

# load data with cases and population
case_pop <- read.table(paste0("cases_pop", exp, ".txt"))
names(case_pop) <- c("location.id","cases","population")

## Process geographical information and convert to grid
geo <- read.table("homesteads.txt", header = T)
geo1 <- geo[,-1]
names(geo1) <- c("x", "y")

# create data clustered by distance
ncluster <- kmeans(cbind(geo1$x, geo1$y), centers = ncenter)
#plot(geo1$x, geo1$y, col = ncluster$cluster, pch = 20)

clust.mean <- data.frame(ncluster$centers)
names(clust.mean) <- c("x", "y")

clust.all = clust.mean[FALSE,]
clust.all$numclust <- clust.all[,1]

for (i in 1:ncenter){
  geof <- geo1[ncluster$cluster==i,]
  geof <- grid2latlong(geof)
  
  geof$numclust <- i
  
  sb <- geof[1,]
  geof <- rbind(geof,sb)
  
  clust.all <- rbind(clust.all, geof)
}

clust.all$numclust <- as.factor(clust.all$numclust)

#extract cluster allocations
nclust <- data.frame(ncluster$cluster)
nclust$location.id <- seq(1:2000)
nclust2 <- inner_join(nclust, geo, by="location.id")

# merge with cases and population
nclust3 <- inner_join(nclust2, case_pop, by="location.id")

# aggregate cases and population by cluster
nclust4 <- nclust3 %>% group_by(ncluster.cluster) %>%
  summarize(cases = sum(cases), population = sum(population))

## Get aggregated counts of population and cases for each county
population <- tapply(nclust3$population,nclust3$ncluster.cluster,sum)
cases <- tapply(nclust3$cases,nclust3$ncluster.cluster,sum)

## Set Parameters
pop.upper.bound <- 0.15
n.simulations <- 9999
alpha.level <- 0.05
plot <- FALSE

## Kulldorff using Binomial likelihoods
binomial <- kulldorff(geo, cases, population, NULL, pop.upper.bound, n.simulations,
                      alpha.level, plot)
cluster <- binomial$most.likely.cluster$location.IDs.included
cluster

cl_data <- geo[geo$location.id %in% cluster, ]
cl_data$cluster <- 1
names(cl_data)[1] <- "ncluster.cluster"

nclust5 <- left_join(nclust2, cl_data[,-c(2:3)], by="ncluster.cluster")
nclust5$cluster[is.na(nclust5$cluster)] <- 0
nclust5$cluster <- as.factor(nclust5$cluster)

# plot(nclust5$locx, nclust5$locy, col=nclust5$cluster, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")
# cluster
# table(nclust5$cluster)
# 
# plot(geo$locx, geo$locy, col="black", pch=1, main="simulated", xlab="",ylab="",xaxt = "n",yaxt = "n")
# points(hot$locx, hot$locy, col=hot$hotspot, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")

names(nclust5)[5] <- paste0("kull",exp)
hotnew <- inner_join(hot, nclust5[,c(2,5)], by="location.id")
write.table(hotnew, "hotspots_decay.txt")

#at 1 year
exp = exp_all[2]

hot <- read.table("hotspots_decay.txt", header = T)

# load data with cases and population
case_pop <- read.table(paste0("cases_pop", exp, ".txt"))
names(case_pop) <- c("location.id","cases","population")

# merge with cases and population
nclust3 <- inner_join(nclust2, case_pop, by="location.id")

# aggregate cases and population by cluster
nclust4 <- nclust3 %>% group_by(ncluster.cluster) %>%
  summarize(cases = sum(cases), population = sum(population))

## Get aggregated counts of population and cases for each county
population <- tapply(nclust3$population,nclust3$ncluster.cluster,sum)
cases <- tapply(nclust3$cases,nclust3$ncluster.cluster,sum)

## Set Parameters
pop.upper.bound <- 0.15
n.simulations <- 9999
alpha.level <- 0.05
plot <- FALSE

## Kulldorff using Binomial likelihoods
binomial <- kulldorff(geo, cases, population, NULL, pop.upper.bound, n.simulations,
                      alpha.level, plot)
cluster <- binomial$most.likely.cluster$location.IDs.included
cluster

cl_data <- geo[geo$location.id %in% cluster, ]
cl_data$cluster <- 1
names(cl_data)[1] <- "ncluster.cluster"

nclust5 <- left_join(nclust2, cl_data[,-c(2:3)], by="ncluster.cluster")
nclust5$cluster[is.na(nclust5$cluster)] <- 0
nclust5$cluster <- as.factor(nclust5$cluster)

#par(mfrow=c(2,2))
# plot(nclust5$locx, nclust5$locy, col=nclust5$cluster, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")
# cluster
# table(nclust5$cluster)
# 
# plot(geo$locx, geo$locy, col="black", pch=1, main="simulated", xlab="",ylab="",xaxt = "n",yaxt = "n")
# points(hot$locx, hot$locy, col=hot$hotspot, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")

names(nclust5)[5] <- paste0("kull",exp)
hotnew <- inner_join(hot, nclust5[,c(2,5)], by="location.id")
write.table(hotnew, "hotspots_decay.txt")

#at 3 year
exp = exp_all[3]

hot <- read.table("hotspots_decay.txt", header = T)

# load data with cases and population
case_pop <- read.table(paste0("cases_pop", exp, ".txt"))
names(case_pop) <- c("location.id","cases","population")

# merge with cases and population
nclust3 <- inner_join(nclust2, case_pop, by="location.id")

# aggregate cases and population by cluster
nclust4 <- nclust3 %>% group_by(ncluster.cluster) %>%
  summarize(cases = sum(cases), population = sum(population))

## Get aggregated counts of population and cases for each county
population <- tapply(nclust3$population,nclust3$ncluster.cluster,sum)
cases <- tapply(nclust3$cases,nclust3$ncluster.cluster,sum)

## Set Parameters
pop.upper.bound <- 0.15
n.simulations <- 9999
alpha.level <- 0.05
plot <- FALSE

## Kulldorff using Binomial likelihoods
binomial <- kulldorff(geo, cases, population, NULL, pop.upper.bound, n.simulations,
                      alpha.level, plot)
cluster <- binomial$most.likely.cluster$location.IDs.included
cluster

cl_data <- geo[geo$location.id %in% cluster, ]
cl_data$cluster <- 1
names(cl_data)[1] <- "ncluster.cluster"

nclust5 <- left_join(nclust2, cl_data[,-c(2:3)], by="ncluster.cluster")
nclust5$cluster[is.na(nclust5$cluster)] <- 0
nclust5$cluster <- as.factor(nclust5$cluster)

#par(mfrow=c(2,2))
# plot(nclust5$locx, nclust5$locy, col=nclust5$cluster, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")
# cluster
# table(nclust5$cluster)
# 
# plot(geo$locx, geo$locy, col="black", pch=1, main="simulated", xlab="",ylab="",xaxt = "n",yaxt = "n")
# points(hot$locx, hot$locy, col=hot$hotspot, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")

names(nclust5)[5] <- paste0("kull",exp)
hotnew <- inner_join(hot, nclust5[,c(2,5)], by="location.id")
write.table(hotnew, "hotspots_decay.txt")


# RIVER
exp_all = c("T03","T","T3y")

ncenter = 15

hh <- 2000

exp = exp_all[1]

hot <- read.table("hotspots_R.txt", header = T)

# load data with cases and population
case_pop <- read.table(paste0("cases_pop", exp, ".txt"))
names(case_pop) <- c("location.id","cases","population")

## Process geographical information and convert to grid
geo <- read.table("homesteads.txt", header = T)
geo1 <- geo[,-1]
names(geo1) <- c("x", "y")

# create data clustered by distance
ncluster <- kmeans(cbind(geo1$x, geo1$y), centers = ncenter)
#plot(geo1$x, geo1$y, col = ncluster$cluster, pch = 20)

clust.mean <- data.frame(ncluster$centers)
names(clust.mean) <- c("x", "y")

clust.all = clust.mean[FALSE,]
clust.all$numclust <- clust.all[,1]

for (i in 1:ncenter){
  geof <- geo1[ncluster$cluster==i,]
  geof <- grid2latlong(geof)
  
  geof$numclust <- i
  
  sb <- geof[1,]
  geof <- rbind(geof,sb)
  
  clust.all <- rbind(clust.all, geof)
}

clust.all$numclust <- as.factor(clust.all$numclust)

#extract cluster allocations
nclust <- data.frame(ncluster$cluster)
nclust$location.id <- seq(1:2000)
nclust2 <- inner_join(nclust, geo, by="location.id")

# merge with cases and population
nclust3 <- inner_join(nclust2, case_pop, by="location.id")

# aggregate cases and population by cluster
nclust4 <- nclust3 %>% group_by(ncluster.cluster) %>%
  summarize(cases = sum(cases), population = sum(population))

## Get aggregated counts of population and cases for each county
population <- tapply(nclust3$population,nclust3$ncluster.cluster,sum)
cases <- tapply(nclust3$cases,nclust3$ncluster.cluster,sum)

## Set Parameters
pop.upper.bound <- 0.15
n.simulations <- 9999
alpha.level <- 0.05
plot <- FALSE

## Kulldorff using Binomial likelihoods
binomial <- kulldorff(geo, cases, population, NULL, pop.upper.bound, n.simulations,
                      alpha.level, plot)
cluster <- binomial$most.likely.cluster$location.IDs.included
cluster

cl_data <- geo[geo$location.id %in% cluster, ]
cl_data$cluster <- 1
names(cl_data)[1] <- "ncluster.cluster"

nclust5 <- left_join(nclust2, cl_data[,-c(2:3)], by="ncluster.cluster")
nclust5$cluster[is.na(nclust5$cluster)] <- 0
nclust5$cluster <- as.factor(nclust5$cluster)

# plot(nclust5$locx, nclust5$locy, col=nclust5$cluster, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")
# cluster
# table(nclust5$cluster)
# 
# plot(geo$locx, geo$locy, col="black", pch=1, main="simulated", xlab="",ylab="",xaxt = "n",yaxt = "n")
# points(hot$locx, hot$locy, col=hot$hotspot, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")

names(nclust5)[5] <- paste0("kull",exp)
hotnew <- inner_join(hot, nclust5[,c(2,5)], by="location.id")
write.table(hotnew, "hotspots_river.txt")

#at 1 year
exp = exp_all[2]

hot <- read.table("hotspots_river.txt", header = T)

# load data with cases and population
case_pop <- read.table(paste0("cases_pop", exp, ".txt"))
names(case_pop) <- c("location.id","cases","population")

# merge with cases and population
nclust3 <- inner_join(nclust2, case_pop, by="location.id")

# aggregate cases and population by cluster
nclust4 <- nclust3 %>% group_by(ncluster.cluster) %>%
  summarize(cases = sum(cases), population = sum(population))

## Get aggregated counts of population and cases for each county
population <- tapply(nclust3$population,nclust3$ncluster.cluster,sum)
cases <- tapply(nclust3$cases,nclust3$ncluster.cluster,sum)

## Set Parameters
pop.upper.bound <- 0.15
n.simulations <- 9999
alpha.level <- 0.05
plot <- FALSE

## Kulldorff using Binomial likelihoods
binomial <- kulldorff(geo, cases, population, NULL, pop.upper.bound, n.simulations,
                      alpha.level, plot)
cluster <- binomial$most.likely.cluster$location.IDs.included
cluster

cl_data <- geo[geo$location.id %in% cluster, ]
cl_data$cluster <- 1
names(cl_data)[1] <- "ncluster.cluster"

nclust5 <- left_join(nclust2, cl_data[,-c(2:3)], by="ncluster.cluster")
nclust5$cluster[is.na(nclust5$cluster)] <- 0
nclust5$cluster <- as.factor(nclust5$cluster)

#par(mfrow=c(2,2))
# plot(nclust5$locx, nclust5$locy, col=nclust5$cluster, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")
# cluster
# table(nclust5$cluster)
# 
# plot(geo$locx, geo$locy, col="black", pch=1, main="simulated", xlab="",ylab="",xaxt = "n",yaxt = "n")
# points(hot$locx, hot$locy, col=hot$hotspot, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")

names(nclust5)[5] <- paste0("kull",exp)
hotnew <- inner_join(hot, nclust5[,c(2,5)], by="location.id")
write.table(hotnew, "hotspots_river.txt")

#at 3 year
exp = exp_all[3]

hot <- read.table("hotspots_river.txt", header = T)

# load data with cases and population
case_pop <- read.table(paste0("cases_pop", exp, ".txt"))
names(case_pop) <- c("location.id","cases","population")

# merge with cases and population
nclust3 <- inner_join(nclust2, case_pop, by="location.id")

# aggregate cases and population by cluster
nclust4 <- nclust3 %>% group_by(ncluster.cluster) %>%
  summarize(cases = sum(cases), population = sum(population))

## Get aggregated counts of population and cases for each county
population <- tapply(nclust3$population,nclust3$ncluster.cluster,sum)
cases <- tapply(nclust3$cases,nclust3$ncluster.cluster,sum)

## Set Parameters
pop.upper.bound <- 0.15
n.simulations <- 9999
alpha.level <- 0.05
plot <- FALSE

## Kulldorff using Binomial likelihoods
binomial <- kulldorff(geo, cases, population, NULL, pop.upper.bound, n.simulations,
                      alpha.level, plot)
cluster <- binomial$most.likely.cluster$location.IDs.included
cluster

cl_data <- geo[geo$location.id %in% cluster, ]
cl_data$cluster <- 1
names(cl_data)[1] <- "ncluster.cluster"

nclust5 <- left_join(nclust2, cl_data[,-c(2:3)], by="ncluster.cluster")
nclust5$cluster[is.na(nclust5$cluster)] <- 0
nclust5$cluster <- as.factor(nclust5$cluster)

#par(mfrow=c(2,2))
# plot(nclust5$locx, nclust5$locy, col=nclust5$cluster, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")
# cluster
# table(nclust5$cluster)
# 
# plot(geo$locx, geo$locy, col="black", pch=1, main="simulated", xlab="",ylab="",xaxt = "n",yaxt = "n")
# points(hot$locx, hot$locy, col=hot$hotspot, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")

names(nclust5)[5] <- paste0("kull",exp)
hotnew <- inner_join(hot, nclust5[,c(2,5)], by="location.id")
write.table(hotnew, "hotspots_river.txt")


# LONG
exp_all = c("T03","T","T3y")

ncenter = 15

hh <- 2000

exp = exp_all[1]

hot <- read.table("hotspots.txt", header = T)

# load data with cases and population
case_pop <- read.table(paste0("cases_pop", exp, ".txt"))
names(case_pop) <- c("location.id","cases","population")

## Process geographical information and convert to grid
geo <- read.table("homesteads.txt", header = T)
geo1 <- geo[,-1]
names(geo1) <- c("x", "y")

# create data clustered by distance
ncluster <- kmeans(cbind(geo1$x, geo1$y), centers = ncenter)
#plot(geo1$x, geo1$y, col = ncluster$cluster, pch = 20)

clust.mean <- data.frame(ncluster$centers)
names(clust.mean) <- c("x", "y")

clust.all = clust.mean[FALSE,]
clust.all$numclust <- clust.all[,1]

for (i in 1:ncenter){
  geof <- geo1[ncluster$cluster==i,]
  geof <- grid2latlong(geof)
  
  geof$numclust <- i
  
  sb <- geof[1,]
  geof <- rbind(geof,sb)
  
  clust.all <- rbind(clust.all, geof)
}

clust.all$numclust <- as.factor(clust.all$numclust)

#extract cluster allocations
nclust <- data.frame(ncluster$cluster)
nclust$location.id <- seq(1:2000)
nclust2 <- inner_join(nclust, geo, by="location.id")

# merge with cases and population
nclust3 <- inner_join(nclust2, case_pop, by="location.id")

# aggregate cases and population by cluster
nclust4 <- nclust3 %>% group_by(ncluster.cluster) %>%
  summarize(cases = sum(cases), population = sum(population))

## Get aggregated counts of population and cases for each county
population <- tapply(nclust3$population,nclust3$ncluster.cluster,sum)
cases <- tapply(nclust3$cases,nclust3$ncluster.cluster,sum)

## Set Parameters
pop.upper.bound <- 0.15
n.simulations <- 9999
alpha.level <- 0.05
plot <- FALSE

## Kulldorff using Binomial likelihoods
binomial <- kulldorff(geo, cases, population, NULL, pop.upper.bound, n.simulations,
                      alpha.level, plot)
cluster <- binomial$most.likely.cluster$location.IDs.included
cluster

cl_data <- geo[geo$location.id %in% cluster, ]
cl_data$cluster <- 1
names(cl_data)[1] <- "ncluster.cluster"

nclust5 <- left_join(nclust2, cl_data[,-c(2:3)], by="ncluster.cluster")
nclust5$cluster[is.na(nclust5$cluster)] <- 0
nclust5$cluster <- as.factor(nclust5$cluster)

# plot(nclust5$locx, nclust5$locy, col=nclust5$cluster, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")
# cluster
# table(nclust5$cluster)
# 
# plot(geo$locx, geo$locy, col="black", pch=1, main="simulated", xlab="",ylab="",xaxt = "n",yaxt = "n")
# points(hot$locx, hot$locy, col=hot$hotspot, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")

names(nclust5)[5] <- paste0("kull",exp)
hotnew <- inner_join(hot, nclust5[,c(2,5)], by="location.id")
write.table(hotnew, "hotspots_long.txt")

#at 1 year
exp = exp_all[2]

hot <- read.table("hotspots_long.txt", header = T)

# load data with cases and population
case_pop <- read.table(paste0("cases_pop", exp, ".txt"))
names(case_pop) <- c("location.id","cases","population")

# merge with cases and population
nclust3 <- inner_join(nclust2, case_pop, by="location.id")

# aggregate cases and population by cluster
nclust4 <- nclust3 %>% group_by(ncluster.cluster) %>%
  summarize(cases = sum(cases), population = sum(population))

## Get aggregated counts of population and cases for each county
population <- tapply(nclust3$population,nclust3$ncluster.cluster,sum)
cases <- tapply(nclust3$cases,nclust3$ncluster.cluster,sum)

## Set Parameters
pop.upper.bound <- 0.15
n.simulations <- 9999
alpha.level <- 0.05
plot <- FALSE

## Kulldorff using Binomial likelihoods
binomial <- kulldorff(geo, cases, population, NULL, pop.upper.bound, n.simulations,
                      alpha.level, plot)
cluster <- binomial$most.likely.cluster$location.IDs.included
cluster

cl_data <- geo[geo$location.id %in% cluster, ]
cl_data$cluster <- 1
names(cl_data)[1] <- "ncluster.cluster"

nclust5 <- left_join(nclust2, cl_data[,-c(2:3)], by="ncluster.cluster")
nclust5$cluster[is.na(nclust5$cluster)] <- 0
nclust5$cluster <- as.factor(nclust5$cluster)

#par(mfrow=c(2,2))
# plot(nclust5$locx, nclust5$locy, col=nclust5$cluster, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")
# cluster
# table(nclust5$cluster)
# 
# plot(geo$locx, geo$locy, col="black", pch=1, main="simulated", xlab="",ylab="",xaxt = "n",yaxt = "n")
# points(hot$locx, hot$locy, col=hot$hotspot, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")

names(nclust5)[5] <- paste0("kull",exp)
hotnew <- inner_join(hot, nclust5[,c(2,5)], by="location.id")
write.table(hotnew, "hotspots_long.txt")

#at 3 year
exp = exp_all[3]

hot <- read.table("hotspots_long.txt", header = T)

# load data with cases and population
case_pop <- read.table(paste0("cases_pop", exp, ".txt"))
names(case_pop) <- c("location.id","cases","population")

# merge with cases and population
nclust3 <- inner_join(nclust2, case_pop, by="location.id")

# aggregate cases and population by cluster
nclust4 <- nclust3 %>% group_by(ncluster.cluster) %>%
  summarize(cases = sum(cases), population = sum(population))

## Get aggregated counts of population and cases for each county
population <- tapply(nclust3$population,nclust3$ncluster.cluster,sum)
cases <- tapply(nclust3$cases,nclust3$ncluster.cluster,sum)

## Set Parameters
pop.upper.bound <- 0.15
n.simulations <- 9999
alpha.level <- 0.05
plot <- FALSE

## Kulldorff using Binomial likelihoods
binomial <- kulldorff(geo, cases, population, NULL, pop.upper.bound, n.simulations,
                      alpha.level, plot)
cluster <- binomial$most.likely.cluster$location.IDs.included
cluster

cl_data <- geo[geo$location.id %in% cluster, ]
cl_data$cluster <- 1
names(cl_data)[1] <- "ncluster.cluster"

nclust5 <- left_join(nclust2, cl_data[,-c(2:3)], by="ncluster.cluster")
nclust5$cluster[is.na(nclust5$cluster)] <- 0
nclust5$cluster <- as.factor(nclust5$cluster)

#par(mfrow=c(2,2))
# plot(nclust5$locx, nclust5$locy, col=nclust5$cluster, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")
# cluster
# table(nclust5$cluster)
# 
# plot(geo$locx, geo$locy, col="black", pch=1, main="simulated", xlab="",ylab="",xaxt = "n",yaxt = "n")
# points(hot$locx, hot$locy, col=hot$hotspot, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")

names(nclust5)[5] <- paste0("kull",exp)
hotnew <- inner_join(hot, nclust5[,c(2,5)], by="location.id")
write.table(hotnew, "hotspots_long.txt")


# SHORT
exp_all = c("T03","T","T3y")

ncenter = 15

hh <- 2000

exp = exp_all[1]

hot <- read.table("hotspots.txt", header = T)

# load data with cases and population
case_pop <- read.table(paste0("cases_pop", exp, ".txt"))
names(case_pop) <- c("location.id","cases","population")

## Process geographical information and convert to grid
geo <- read.table("homesteads.txt", header = T)
geo1 <- geo[,-1]
names(geo1) <- c("x", "y")

# create data clustered by distance
ncluster <- kmeans(cbind(geo1$x, geo1$y), centers = ncenter)
#plot(geo1$x, geo1$y, col = ncluster$cluster, pch = 20)

clust.mean <- data.frame(ncluster$centers)
names(clust.mean) <- c("x", "y")

clust.all = clust.mean[FALSE,]
clust.all$numclust <- clust.all[,1]

for (i in 1:ncenter){
  geof <- geo1[ncluster$cluster==i,]
  geof <- grid2latlong(geof)
  
  geof$numclust <- i
  
  sb <- geof[1,]
  geof <- rbind(geof,sb)
  
  clust.all <- rbind(clust.all, geof)
}

clust.all$numclust <- as.factor(clust.all$numclust)

#extract cluster allocations
nclust <- data.frame(ncluster$cluster)
nclust$location.id <- seq(1:2000)
nclust2 <- inner_join(nclust, geo, by="location.id")

# merge with cases and population
nclust3 <- inner_join(nclust2, case_pop, by="location.id")

# aggregate cases and population by cluster
nclust4 <- nclust3 %>% group_by(ncluster.cluster) %>%
  summarize(cases = sum(cases), population = sum(population))

## Get aggregated counts of population and cases for each county
population <- tapply(nclust3$population,nclust3$ncluster.cluster,sum)
cases <- tapply(nclust3$cases,nclust3$ncluster.cluster,sum)

## Set Parameters
pop.upper.bound <- 0.15
n.simulations <- 9999
alpha.level <- 0.05
plot <- FALSE

## Kulldorff using Binomial likelihoods
binomial <- kulldorff(geo, cases, population, NULL, pop.upper.bound, n.simulations,
                      alpha.level, plot)
cluster <- binomial$most.likely.cluster$location.IDs.included
cluster

cl_data <- geo[geo$location.id %in% cluster, ]
cl_data$cluster <- 1
names(cl_data)[1] <- "ncluster.cluster"

nclust5 <- left_join(nclust2, cl_data[,-c(2:3)], by="ncluster.cluster")
nclust5$cluster[is.na(nclust5$cluster)] <- 0
nclust5$cluster <- as.factor(nclust5$cluster)

# plot(nclust5$locx, nclust5$locy, col=nclust5$cluster, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")
# cluster
# table(nclust5$cluster)
# 
# plot(geo$locx, geo$locy, col="black", pch=1, main="simulated", xlab="",ylab="",xaxt = "n",yaxt = "n")
# points(hot$locx, hot$locy, col=hot$hotspot, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")

names(nclust5)[5] <- paste0("kull",exp)
hotnew <- inner_join(hot, nclust5[,c(2,5)], by="location.id")
write.table(hotnew, "hotspots_short.txt")

#at 1 year
exp = exp_all[2]

hot <- read.table("hotspots_short.txt", header = T)

# load data with cases and population
case_pop <- read.table(paste0("cases_pop", exp, ".txt"))
names(case_pop) <- c("location.id","cases","population")

# merge with cases and population
nclust3 <- inner_join(nclust2, case_pop, by="location.id")

# aggregate cases and population by cluster
nclust4 <- nclust3 %>% group_by(ncluster.cluster) %>%
  summarize(cases = sum(cases), population = sum(population))

## Get aggregated counts of population and cases for each county
population <- tapply(nclust3$population,nclust3$ncluster.cluster,sum)
cases <- tapply(nclust3$cases,nclust3$ncluster.cluster,sum)

## Set Parameters
pop.upper.bound <- 0.15
n.simulations <- 9999
alpha.level <- 0.05
plot <- FALSE

## Kulldorff using Binomial likelihoods
binomial <- kulldorff(geo, cases, population, NULL, pop.upper.bound, n.simulations,
                      alpha.level, plot)
cluster <- binomial$most.likely.cluster$location.IDs.included
cluster

cl_data <- geo[geo$location.id %in% cluster, ]
cl_data$cluster <- 1
names(cl_data)[1] <- "ncluster.cluster"

nclust5 <- left_join(nclust2, cl_data[,-c(2:3)], by="ncluster.cluster")
nclust5$cluster[is.na(nclust5$cluster)] <- 0
nclust5$cluster <- as.factor(nclust5$cluster)

#par(mfrow=c(2,2))
# plot(nclust5$locx, nclust5$locy, col=nclust5$cluster, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")
# cluster
# table(nclust5$cluster)
# 
# plot(geo$locx, geo$locy, col="black", pch=1, main="simulated", xlab="",ylab="",xaxt = "n",yaxt = "n")
# points(hot$locx, hot$locy, col=hot$hotspot, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")

names(nclust5)[5] <- paste0("kull",exp)
hotnew <- inner_join(hot, nclust5[,c(2,5)], by="location.id")
write.table(hotnew, "hotspots_short.txt")

#at 3 year
exp = exp_all[3]

hot <- read.table("hotspots_short.txt", header = T)

# load data with cases and population
case_pop <- read.table(paste0("cases_pop", exp, ".txt"))
names(case_pop) <- c("location.id","cases","population")

# merge with cases and population
nclust3 <- inner_join(nclust2, case_pop, by="location.id")

# aggregate cases and population by cluster
nclust4 <- nclust3 %>% group_by(ncluster.cluster) %>%
  summarize(cases = sum(cases), population = sum(population))

## Get aggregated counts of population and cases for each county
population <- tapply(nclust3$population,nclust3$ncluster.cluster,sum)
cases <- tapply(nclust3$cases,nclust3$ncluster.cluster,sum)

## Set Parameters
pop.upper.bound <- 0.15
n.simulations <- 9999
alpha.level <- 0.05
plot <- FALSE

## Kulldorff using Binomial likelihoods
binomial <- kulldorff(geo, cases, population, NULL, pop.upper.bound, n.simulations,
                      alpha.level, plot)
cluster <- binomial$most.likely.cluster$location.IDs.included
cluster

cl_data <- geo[geo$location.id %in% cluster, ]
cl_data$cluster <- 1
names(cl_data)[1] <- "ncluster.cluster"

nclust5 <- left_join(nclust2, cl_data[,-c(2:3)], by="ncluster.cluster")
nclust5$cluster[is.na(nclust5$cluster)] <- 0
nclust5$cluster <- as.factor(nclust5$cluster)

#par(mfrow=c(2,2))
# plot(nclust5$locx, nclust5$locy, col=nclust5$cluster, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")
# cluster
# table(nclust5$cluster)
# 
# plot(geo$locx, geo$locy, col="black", pch=1, main="simulated", xlab="",ylab="",xaxt = "n",yaxt = "n")
# points(hot$locx, hot$locy, col=hot$hotspot, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")

names(nclust5)[5] <- paste0("kull",exp)
hotnew <- inner_join(hot, nclust5[,c(2,5)], by="location.id")
write.table(hotnew, "hotspots_short.txt")


# WET SEASON
exp_all = c("T03","T","T3y")

ncenter = 15

hh <- 2000

exp = exp_all[1]

hot <- read.table("hotspots.txt", header = T)

# load data with cases and population
case_pop <- read.table(paste0("cases_pop", exp, ".txt"))
names(case_pop) <- c("location.id","cases","population")

## Process geographical information and convert to grid
geo <- read.table("homesteads.txt", header = T)
geo1 <- geo[,-1]
names(geo1) <- c("x", "y")

# create data clustered by distance
ncluster <- kmeans(cbind(geo1$x, geo1$y), centers = ncenter)
#plot(geo1$x, geo1$y, col = ncluster$cluster, pch = 20)

clust.mean <- data.frame(ncluster$centers)
names(clust.mean) <- c("x", "y")

clust.all = clust.mean[FALSE,]
clust.all$numclust <- clust.all[,1]

for (i in 1:ncenter){
  geof <- geo1[ncluster$cluster==i,]
  geof <- grid2latlong(geof)
  
  geof$numclust <- i
  
  sb <- geof[1,]
  geof <- rbind(geof,sb)
  
  clust.all <- rbind(clust.all, geof)
}

clust.all$numclust <- as.factor(clust.all$numclust)

#extract cluster allocations
nclust <- data.frame(ncluster$cluster)
nclust$location.id <- seq(1:2000)
nclust2 <- inner_join(nclust, geo, by="location.id")

# merge with cases and population
nclust3 <- inner_join(nclust2, case_pop, by="location.id")

# aggregate cases and population by cluster
nclust4 <- nclust3 %>% group_by(ncluster.cluster) %>%
  summarize(cases = sum(cases), population = sum(population))

## Get aggregated counts of population and cases for each county
population <- tapply(nclust3$population,nclust3$ncluster.cluster,sum)
cases <- tapply(nclust3$cases,nclust3$ncluster.cluster,sum)

## Set Parameters
pop.upper.bound <- 0.15
n.simulations <- 9999
alpha.level <- 0.05
plot <- FALSE

## Kulldorff using Binomial likelihoods
binomial <- kulldorff(geo, cases, population, NULL, pop.upper.bound, n.simulations,
                      alpha.level, plot)
cluster <- binomial$most.likely.cluster$location.IDs.included
cluster

cl_data <- geo[geo$location.id %in% cluster, ]
cl_data$cluster <- 1
names(cl_data)[1] <- "ncluster.cluster"

nclust5 <- left_join(nclust2, cl_data[,-c(2:3)], by="ncluster.cluster")
nclust5$cluster[is.na(nclust5$cluster)] <- 0
nclust5$cluster <- as.factor(nclust5$cluster)

# plot(nclust5$locx, nclust5$locy, col=nclust5$cluster, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")
# cluster
# table(nclust5$cluster)
# 
# plot(geo$locx, geo$locy, col="black", pch=1, main="simulated", xlab="",ylab="",xaxt = "n",yaxt = "n")
# points(hot$locx, hot$locy, col=hot$hotspot, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")

names(nclust5)[5] <- paste0("kull",exp)
hotnew <- inner_join(hot, nclust5[,c(2,5)], by="location.id")
write.table(hotnew, "hotspots_season.txt")

#at 1 year
exp = exp_all[2]

hot <- read.table("hotspots_season.txt", header = T)

# load data with cases and population
case_pop <- read.table(paste0("cases_pop", exp, ".txt"))
names(case_pop) <- c("location.id","cases","population")

# merge with cases and population
nclust3 <- inner_join(nclust2, case_pop, by="location.id")

# aggregate cases and population by cluster
nclust4 <- nclust3 %>% group_by(ncluster.cluster) %>%
  summarize(cases = sum(cases), population = sum(population))

## Get aggregated counts of population and cases for each county
population <- tapply(nclust3$population,nclust3$ncluster.cluster,sum)
cases <- tapply(nclust3$cases,nclust3$ncluster.cluster,sum)

## Set Parameters
pop.upper.bound <- 0.15
n.simulations <- 9999
alpha.level <- 0.05
plot <- FALSE

## Kulldorff using Binomial likelihoods
binomial <- kulldorff(geo, cases, population, NULL, pop.upper.bound, n.simulations,
                      alpha.level, plot)
cluster <- binomial$most.likely.cluster$location.IDs.included
cluster

cl_data <- geo[geo$location.id %in% cluster, ]
cl_data$cluster <- 1
names(cl_data)[1] <- "ncluster.cluster"

nclust5 <- left_join(nclust2, cl_data[,-c(2:3)], by="ncluster.cluster")
nclust5$cluster[is.na(nclust5$cluster)] <- 0
nclust5$cluster <- as.factor(nclust5$cluster)

#par(mfrow=c(2,2))
# plot(nclust5$locx, nclust5$locy, col=nclust5$cluster, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")
# cluster
# table(nclust5$cluster)
# 
# plot(geo$locx, geo$locy, col="black", pch=1, main="simulated", xlab="",ylab="",xaxt = "n",yaxt = "n")
# points(hot$locx, hot$locy, col=hot$hotspot, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")

names(nclust5)[5] <- paste0("kull",exp)
hotnew <- inner_join(hot, nclust5[,c(2,5)], by="location.id")
write.table(hotnew, "hotspots_season.txt")

#at 3 year
exp = exp_all[3]

hot <- read.table("hotspots_season.txt", header = T)

# load data with cases and population
case_pop <- read.table(paste0("cases_pop", exp, ".txt"))
names(case_pop) <- c("location.id","cases","population")

# merge with cases and population
nclust3 <- inner_join(nclust2, case_pop, by="location.id")

# aggregate cases and population by cluster
nclust4 <- nclust3 %>% group_by(ncluster.cluster) %>%
  summarize(cases = sum(cases), population = sum(population))

## Get aggregated counts of population and cases for each county
population <- tapply(nclust3$population,nclust3$ncluster.cluster,sum)
cases <- tapply(nclust3$cases,nclust3$ncluster.cluster,sum)

## Set Parameters
pop.upper.bound <- 0.15
n.simulations <- 9999
alpha.level <- 0.05
plot <- FALSE

## Kulldorff using Binomial likelihoods
binomial <- kulldorff(geo, cases, population, NULL, pop.upper.bound, n.simulations,
                      alpha.level, plot)
cluster <- binomial$most.likely.cluster$location.IDs.included
cluster

cl_data <- geo[geo$location.id %in% cluster, ]
cl_data$cluster <- 1
names(cl_data)[1] <- "ncluster.cluster"

nclust5 <- left_join(nclust2, cl_data[,-c(2:3)], by="ncluster.cluster")
nclust5$cluster[is.na(nclust5$cluster)] <- 0
nclust5$cluster <- as.factor(nclust5$cluster)

#par(mfrow=c(2,2))
# plot(nclust5$locx, nclust5$locy, col=nclust5$cluster, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")
# cluster
# table(nclust5$cluster)
# 
# plot(geo$locx, geo$locy, col="black", pch=1, main="simulated", xlab="",ylab="",xaxt = "n",yaxt = "n")
# points(hot$locx, hot$locy, col=hot$hotspot, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")

names(nclust5)[5] <- paste0("kull",exp)
hotnew <- inner_join(hot, nclust5[,c(2,5)], by="location.id")
write.table(hotnew, "hotspots_season.txt")


# DRY SEASON
exp_all = c("T03","T","T3y")

ncenter = 15

hh <- 2000

exp = exp_all[1]

hot <- read.table("hotspots.txt", header = T)

# load data with cases and population
case_pop <- read.table(paste0("cases_pop", exp, ".txt"))
names(case_pop) <- c("location.id","cases","population")

## Process geographical information and convert to grid
geo <- read.table("homesteads.txt", header = T)
geo1 <- geo[,-1]
names(geo1) <- c("x", "y")

# create data clustered by distance
ncluster <- kmeans(cbind(geo1$x, geo1$y), centers = ncenter)
#plot(geo1$x, geo1$y, col = ncluster$cluster, pch = 20)

clust.mean <- data.frame(ncluster$centers)
names(clust.mean) <- c("x", "y")

clust.all = clust.mean[FALSE,]
clust.all$numclust <- clust.all[,1]

for (i in 1:ncenter){
  geof <- geo1[ncluster$cluster==i,]
  geof <- grid2latlong(geof)
  
  geof$numclust <- i
  
  sb <- geof[1,]
  geof <- rbind(geof,sb)
  
  clust.all <- rbind(clust.all, geof)
}

clust.all$numclust <- as.factor(clust.all$numclust)

#extract cluster allocations
nclust <- data.frame(ncluster$cluster)
nclust$location.id <- seq(1:2000)
nclust2 <- inner_join(nclust, geo, by="location.id")

# merge with cases and population
nclust3 <- inner_join(nclust2, case_pop, by="location.id")

# aggregate cases and population by cluster
nclust4 <- nclust3 %>% group_by(ncluster.cluster) %>%
  summarize(cases = sum(cases), population = sum(population))

## Get aggregated counts of population and cases for each county
population <- tapply(nclust3$population,nclust3$ncluster.cluster,sum)
cases <- tapply(nclust3$cases,nclust3$ncluster.cluster,sum)

## Set Parameters
pop.upper.bound <- 0.15
n.simulations <- 9999
alpha.level <- 0.05
plot <- FALSE

## Kulldorff using Binomial likelihoods
binomial <- kulldorff(geo, cases, population, NULL, pop.upper.bound, n.simulations,
                      alpha.level, plot)
cluster <- binomial$most.likely.cluster$location.IDs.included
cluster

cl_data <- geo[geo$location.id %in% cluster, ]
cl_data$cluster <- 1
names(cl_data)[1] <- "ncluster.cluster"

nclust5 <- left_join(nclust2, cl_data[,-c(2:3)], by="ncluster.cluster")
nclust5$cluster[is.na(nclust5$cluster)] <- 0
nclust5$cluster <- as.factor(nclust5$cluster)

# plot(nclust5$locx, nclust5$locy, col=nclust5$cluster, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")
# cluster
# table(nclust5$cluster)
# 
# plot(geo$locx, geo$locy, col="black", pch=1, main="simulated", xlab="",ylab="",xaxt = "n",yaxt = "n")
# points(hot$locx, hot$locy, col=hot$hotspot, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")

names(nclust5)[5] <- paste0("kull",exp)
hotnew <- inner_join(hot, nclust5[,c(2,5)], by="location.id")
write.table(hotnew, "hotspots_dry.txt")

#at 1 year
exp = exp_all[2]

hot <- read.table("hotspots_dry.txt", header = T)

# load data with cases and population
case_pop <- read.table(paste0("cases_pop", exp, ".txt"))
names(case_pop) <- c("location.id","cases","population")

# merge with cases and population
nclust3 <- inner_join(nclust2, case_pop, by="location.id")

# aggregate cases and population by cluster
nclust4 <- nclust3 %>% group_by(ncluster.cluster) %>%
  summarize(cases = sum(cases), population = sum(population))

## Get aggregated counts of population and cases for each county
population <- tapply(nclust3$population,nclust3$ncluster.cluster,sum)
cases <- tapply(nclust3$cases,nclust3$ncluster.cluster,sum)

## Set Parameters
pop.upper.bound <- 0.15
n.simulations <- 9999
alpha.level <- 0.05
plot <- FALSE

## Kulldorff using Binomial likelihoods
binomial <- kulldorff(geo, cases, population, NULL, pop.upper.bound, n.simulations,
                      alpha.level, plot)
cluster <- binomial$most.likely.cluster$location.IDs.included
cluster

cl_data <- geo[geo$location.id %in% cluster, ]
cl_data$cluster <- 1
names(cl_data)[1] <- "ncluster.cluster"

nclust5 <- left_join(nclust2, cl_data[,-c(2:3)], by="ncluster.cluster")
nclust5$cluster[is.na(nclust5$cluster)] <- 0
nclust5$cluster <- as.factor(nclust5$cluster)

#par(mfrow=c(2,2))
# plot(nclust5$locx, nclust5$locy, col=nclust5$cluster, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")
# cluster
# table(nclust5$cluster)
# 
# plot(geo$locx, geo$locy, col="black", pch=1, main="simulated", xlab="",ylab="",xaxt = "n",yaxt = "n")
# points(hot$locx, hot$locy, col=hot$hotspot, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")

names(nclust5)[5] <- paste0("kull",exp)
hotnew <- inner_join(hot, nclust5[,c(2,5)], by="location.id")
write.table(hotnew, "hotspots_dry.txt")

#at 3 year
exp = exp_all[3]

hot <- read.table("hotspots_dry.txt", header = T)

# load data with cases and population
case_pop <- read.table(paste0("cases_pop", exp, ".txt"))
names(case_pop) <- c("location.id","cases","population")

# merge with cases and population
nclust3 <- inner_join(nclust2, case_pop, by="location.id")

# aggregate cases and population by cluster
nclust4 <- nclust3 %>% group_by(ncluster.cluster) %>%
  summarize(cases = sum(cases), population = sum(population))

## Get aggregated counts of population and cases for each county
population <- tapply(nclust3$population,nclust3$ncluster.cluster,sum)
cases <- tapply(nclust3$cases,nclust3$ncluster.cluster,sum)

## Set Parameters
pop.upper.bound <- 0.15
n.simulations <- 9999
alpha.level <- 0.05
plot <- FALSE

## Kulldorff using Binomial likelihoods
binomial <- kulldorff(geo, cases, population, NULL, pop.upper.bound, n.simulations,
                      alpha.level, plot)
cluster <- binomial$most.likely.cluster$location.IDs.included
cluster

cl_data <- geo[geo$location.id %in% cluster, ]
cl_data$cluster <- 1
names(cl_data)[1] <- "ncluster.cluster"

nclust5 <- left_join(nclust2, cl_data[,-c(2:3)], by="ncluster.cluster")
nclust5$cluster[is.na(nclust5$cluster)] <- 0
nclust5$cluster <- as.factor(nclust5$cluster)

#par(mfrow=c(2,2))
# plot(nclust5$locx, nclust5$locy, col=nclust5$cluster, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")
# cluster
# table(nclust5$cluster)
# 
# plot(geo$locx, geo$locy, col="black", pch=1, main="simulated", xlab="",ylab="",xaxt = "n",yaxt = "n")
# points(hot$locx, hot$locy, col=hot$hotspot, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")

names(nclust5)[5] <- paste0("kull",exp)
hotnew <- inner_join(hot, nclust5[,c(2,5)], by="location.id")
write.table(hotnew, "hotspots_dry.txt")

#}
