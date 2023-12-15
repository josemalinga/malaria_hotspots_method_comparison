
rm(list=ls())
library(SpatialEpi)
library(dplyr)
library(smerc)
library(tidyr)
library("FRK")
#detach(package:plyr)

remove.packages(SpatialEpi)
install.packages("remotes")
library(remotes)
install_version("SpatialEpi", "1.2.1")

set.seed(42)

for (i in 1:30){
  # c(1,5,7,9,11,13,15,17,19,20)
  # User 
  user = strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
  
  # Working directory
  setwd(paste0("/scicore/home/penny/",user,"/malaria_hotspots/simulations_nov"))
  
  experiment <- paste0("rep",i);experiment
  sed <- i*42
  set.seed(sed)
  
  setwd(paste0(getwd(),"/",experiment))
  getwd()
  
  
# SINGLE
exp_all = c("T03","T","T3y")

ncenter = n.strata = 5
hh <- 2000

exp = exp_all[2]

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
plot(geo1$x, geo1$y, col = ncluster$cluster, pch = 20)

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

## Get aggregated counts of population and cases for each county - CHECK THE EXPECTED CASES AND SEE VALUES
population <- tapply(nclust3$population,nclust3$ncluster.cluster,sum)
cases <- tapply(nclust3$cases,nclust3$ncluster.cluster,sum)

expected.cases <- expected(nclust3$population, nclust3$cases, n.strata)

## Set Parameters
pop.upper.bound <- 0.1
n.simulations <- 9999
alpha.level <- 0.05
plot <- FALSE
k <- 10

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



## Kulldorff using Poisson likelihoods
#poisson <- kulldorff(geo, cases, population, expected.cases, pop.upper.bound, 
#                     n.simulations, alpha.level, plot)
#cluster <- poisson$most.likely.cluster$location.IDs.included

#bayesian analysis results
bayes_res <- besag.newell(geo, population, cases, expected.cases=NULL, k, alpha.level)





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
type <- "_M"
case_pop <- read.table(paste0("cases_pop", exp, type, ".txt"))
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
type <- "_M"
case_pop <- read.table(paste0("cases_pop", exp, type, ".txt"))
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
type <- "_M"
case_pop <- read.table(paste0("cases_pop", exp, type, ".txt"))
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
type <- "_D"
case_pop <- read.table(paste0("cases_pop", exp, type, ".txt"))
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
type <- "_D"
case_pop <- read.table(paste0("cases_pop", exp, type, ".txt"))
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
type <- "_D"
case_pop <- read.table(paste0("cases_pop", exp, type, ".txt"))
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
type <- "_R"
case_pop <- read.table(paste0("cases_pop", exp, type, ".txt"))
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
type <- "_R"
case_pop <- read.table(paste0("cases_pop", exp, type, ".txt"))
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
type <- "_R"
case_pop <- read.table(paste0("cases_pop", exp, type, ".txt"))
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
type <- "_long"
case_pop <- read.table(paste0("cases_pop", exp, type, ".txt"))
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
type <- "_long"
case_pop <- read.table(paste0("cases_pop", exp, type, ".txt"))
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
type <- "_long"
case_pop <- read.table(paste0("cases_pop", exp, type, ".txt"))
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
type <- "_short"
case_pop <- read.table(paste0("cases_pop", exp, type, ".txt"))
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
type <- "_short"
case_pop <- read.table(paste0("cases_pop", exp, type, ".txt"))
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
type <- "_short"
case_pop <- read.table(paste0("cases_pop", exp, type, ".txt"))
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


# # WET SEASON
# exp_all = c("T03","T","T3y")
# 
# ncenter = 15
# 
# hh <- 2000
# 
# exp = exp_all[1]
# 
# hot <- read.table("hotspots.txt", header = T)
# 
# # load data with cases and population
# type <- "_seas"
# case_pop <- read.table(paste0("cases_pop", exp, type, ".txt"))
# names(case_pop) <- c("location.id","cases","population")
# 
# ## Process geographical information and convert to grid
# geo <- read.table("homesteads.txt", header = T)
# geo1 <- geo[,-1]
# names(geo1) <- c("x", "y")
# 
# # create data clustered by distance
# ncluster <- kmeans(cbind(geo1$x, geo1$y), centers = ncenter)
# #plot(geo1$x, geo1$y, col = ncluster$cluster, pch = 20)
# 
# clust.mean <- data.frame(ncluster$centers)
# names(clust.mean) <- c("x", "y")
# 
# clust.all = clust.mean[FALSE,]
# clust.all$numclust <- clust.all[,1]
# 
# for (i in 1:ncenter){
#   geof <- geo1[ncluster$cluster==i,]
#   geof <- grid2latlong(geof)
#   
#   geof$numclust <- i
#   
#   sb <- geof[1,]
#   geof <- rbind(geof,sb)
#   
#   clust.all <- rbind(clust.all, geof)
# }
# 
# clust.all$numclust <- as.factor(clust.all$numclust)
# 
# #extract cluster allocations
# nclust <- data.frame(ncluster$cluster)
# nclust$location.id <- seq(1:2000)
# nclust2 <- inner_join(nclust, geo, by="location.id")
# 
# # merge with cases and population
# nclust3 <- inner_join(nclust2, case_pop, by="location.id")
# 
# # aggregate cases and population by cluster
# nclust4 <- nclust3 %>% group_by(ncluster.cluster) %>%
#   summarize(cases = sum(cases), population = sum(population))
# 
# ## Get aggregated counts of population and cases for each county
# population <- tapply(nclust3$population,nclust3$ncluster.cluster,sum)
# cases <- tapply(nclust3$cases,nclust3$ncluster.cluster,sum)
# 
# ## Set Parameters
# pop.upper.bound <- 0.15
# n.simulations <- 9999
# alpha.level <- 0.05
# plot <- FALSE
# 
# ## Kulldorff using Binomial likelihoods
# binomial <- kulldorff(geo, cases, population, NULL, pop.upper.bound, n.simulations,
#                       alpha.level, plot)
# cluster <- binomial$most.likely.cluster$location.IDs.included
# cluster
# 
# cl_data <- geo[geo$location.id %in% cluster, ]
# cl_data$cluster <- 1
# names(cl_data)[1] <- "ncluster.cluster"
# 
# nclust5 <- left_join(nclust2, cl_data[,-c(2:3)], by="ncluster.cluster")
# nclust5$cluster[is.na(nclust5$cluster)] <- 0
# nclust5$cluster <- as.factor(nclust5$cluster)
# 
# # plot(nclust5$locx, nclust5$locy, col=nclust5$cluster, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")
# # cluster
# # table(nclust5$cluster)
# # 
# # plot(geo$locx, geo$locy, col="black", pch=1, main="simulated", xlab="",ylab="",xaxt = "n",yaxt = "n")
# # points(hot$locx, hot$locy, col=hot$hotspot, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")
# 
# names(nclust5)[5] <- paste0("kull",exp)
# hotnew <- inner_join(hot, nclust5[,c(2,5)], by="location.id")
# write.table(hotnew, "hotspots_season.txt")
# 
# #at 1 year
# exp = exp_all[2]
# 
# hot <- read.table("hotspots_season.txt", header = T)
# 
# # load data with cases and population
# type <- "_seas"
# case_pop <- read.table(paste0("cases_pop", exp, type, ".txt"))
# names(case_pop) <- c("location.id","cases","population")
# 
# # merge with cases and population
# nclust3 <- inner_join(nclust2, case_pop, by="location.id")
# 
# # aggregate cases and population by cluster
# nclust4 <- nclust3 %>% group_by(ncluster.cluster) %>%
#   summarize(cases = sum(cases), population = sum(population))
# 
# ## Get aggregated counts of population and cases for each county
# population <- tapply(nclust3$population,nclust3$ncluster.cluster,sum)
# cases <- tapply(nclust3$cases,nclust3$ncluster.cluster,sum)
# 
# ## Set Parameters
# pop.upper.bound <- 0.15
# n.simulations <- 9999
# alpha.level <- 0.05
# plot <- FALSE
# 
# ## Kulldorff using Binomial likelihoods
# binomial <- kulldorff(geo, cases, population, NULL, pop.upper.bound, n.simulations,
#                       alpha.level, plot)
# cluster <- binomial$most.likely.cluster$location.IDs.included
# cluster
# 
# cl_data <- geo[geo$location.id %in% cluster, ]
# cl_data$cluster <- 1
# names(cl_data)[1] <- "ncluster.cluster"
# 
# nclust5 <- left_join(nclust2, cl_data[,-c(2:3)], by="ncluster.cluster")
# nclust5$cluster[is.na(nclust5$cluster)] <- 0
# nclust5$cluster <- as.factor(nclust5$cluster)
# 
# #par(mfrow=c(2,2))
# # plot(nclust5$locx, nclust5$locy, col=nclust5$cluster, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")
# # cluster
# # table(nclust5$cluster)
# # 
# # plot(geo$locx, geo$locy, col="black", pch=1, main="simulated", xlab="",ylab="",xaxt = "n",yaxt = "n")
# # points(hot$locx, hot$locy, col=hot$hotspot, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")
# 
# names(nclust5)[5] <- paste0("kull",exp)
# hotnew <- inner_join(hot, nclust5[,c(2,5)], by="location.id")
# write.table(hotnew, "hotspots_season.txt")
# 
# #at 3 year
# exp = exp_all[3]
# 
# hot <- read.table("hotspots_season.txt", header = T)
# 
# # load data with cases and population
# type <- "_seas"
# case_pop <- read.table(paste0("cases_pop", exp, type, ".txt"))
# names(case_pop) <- c("location.id","cases","population")
# 
# # merge with cases and population
# nclust3 <- inner_join(nclust2, case_pop, by="location.id")
# 
# # aggregate cases and population by cluster
# nclust4 <- nclust3 %>% group_by(ncluster.cluster) %>%
#   summarize(cases = sum(cases), population = sum(population))
# 
# ## Get aggregated counts of population and cases for each county
# population <- tapply(nclust3$population,nclust3$ncluster.cluster,sum)
# cases <- tapply(nclust3$cases,nclust3$ncluster.cluster,sum)
# 
# ## Set Parameters
# pop.upper.bound <- 0.15
# n.simulations <- 9999
# alpha.level <- 0.05
# plot <- FALSE
# 
# ## Kulldorff using Binomial likelihoods
# binomial <- kulldorff(geo, cases, population, NULL, pop.upper.bound, n.simulations,
#                       alpha.level, plot)
# cluster <- binomial$most.likely.cluster$location.IDs.included
# cluster
# 
# cl_data <- geo[geo$location.id %in% cluster, ]
# cl_data$cluster <- 1
# names(cl_data)[1] <- "ncluster.cluster"
# 
# nclust5 <- left_join(nclust2, cl_data[,-c(2:3)], by="ncluster.cluster")
# nclust5$cluster[is.na(nclust5$cluster)] <- 0
# nclust5$cluster <- as.factor(nclust5$cluster)
# 
# #par(mfrow=c(2,2))
# # plot(nclust5$locx, nclust5$locy, col=nclust5$cluster, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")
# # cluster
# # table(nclust5$cluster)
# # 
# # plot(geo$locx, geo$locy, col="black", pch=1, main="simulated", xlab="",ylab="",xaxt = "n",yaxt = "n")
# # points(hot$locx, hot$locy, col=hot$hotspot, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")
# 
# names(nclust5)[5] <- paste0("kull",exp)
# hotnew <- inner_join(hot, nclust5[,c(2,5)], by="location.id")
# write.table(hotnew, "hotspots_season.txt")
# 
# 
# # DRY SEASON
# exp_all = c("T03","T","T3y")
# 
# ncenter = 15
# 
# hh <- 2000
# 
# exp = exp_all[1]
# 
# hot <- read.table("hotspots.txt", header = T)
# 
# # load data with cases and population
# type <- "_dry"
# case_pop <- read.table(paste0("cases_pop", exp, type, ".txt"))
# names(case_pop) <- c("location.id","cases","population")
# 
# ## Process geographical information and convert to grid
# geo <- read.table("homesteads.txt", header = T)
# geo1 <- geo[,-1]
# names(geo1) <- c("x", "y")
# 
# # create data clustered by distance
# ncluster <- kmeans(cbind(geo1$x, geo1$y), centers = ncenter)
# #plot(geo1$x, geo1$y, col = ncluster$cluster, pch = 20)
# 
# clust.mean <- data.frame(ncluster$centers)
# names(clust.mean) <- c("x", "y")
# 
# clust.all = clust.mean[FALSE,]
# clust.all$numclust <- clust.all[,1]
# 
# for (i in 1:ncenter){
#   geof <- geo1[ncluster$cluster==i,]
#   geof <- grid2latlong(geof)
#   
#   geof$numclust <- i
#   
#   sb <- geof[1,]
#   geof <- rbind(geof,sb)
#   
#   clust.all <- rbind(clust.all, geof)
# }
# 
# clust.all$numclust <- as.factor(clust.all$numclust)
# 
# #extract cluster allocations
# nclust <- data.frame(ncluster$cluster)
# nclust$location.id <- seq(1:2000)
# nclust2 <- inner_join(nclust, geo, by="location.id")
# 
# # merge with cases and population
# nclust3 <- inner_join(nclust2, case_pop, by="location.id")
# 
# # aggregate cases and population by cluster
# nclust4 <- nclust3 %>% group_by(ncluster.cluster) %>%
#   summarize(cases = sum(cases), population = sum(population))
# 
# ## Get aggregated counts of population and cases for each county
# population <- tapply(nclust3$population,nclust3$ncluster.cluster,sum)
# cases <- tapply(nclust3$cases,nclust3$ncluster.cluster,sum)
# 
# ## Set Parameters
# pop.upper.bound <- 0.15
# n.simulations <- 9999
# alpha.level <- 0.05
# plot <- FALSE
# 
# ## Kulldorff using Binomial likelihoods
# binomial <- kulldorff(geo, cases, population, NULL, pop.upper.bound, n.simulations,
#                       alpha.level, plot)
# cluster <- binomial$most.likely.cluster$location.IDs.included
# cluster
# 
# cl_data <- geo[geo$location.id %in% cluster, ]
# cl_data$cluster <- 1
# names(cl_data)[1] <- "ncluster.cluster"
# 
# nclust5 <- left_join(nclust2, cl_data[,-c(2:3)], by="ncluster.cluster")
# nclust5$cluster[is.na(nclust5$cluster)] <- 0
# nclust5$cluster <- as.factor(nclust5$cluster)
# 
# # plot(nclust5$locx, nclust5$locy, col=nclust5$cluster, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")
# # cluster
# # table(nclust5$cluster)
# # 
# # plot(geo$locx, geo$locy, col="black", pch=1, main="simulated", xlab="",ylab="",xaxt = "n",yaxt = "n")
# # points(hot$locx, hot$locy, col=hot$hotspot, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")
# 
# names(nclust5)[5] <- paste0("kull",exp)
# hotnew <- inner_join(hot, nclust5[,c(2,5)], by="location.id")
# write.table(hotnew, "hotspots_dry.txt")
# 
# #at 1 year
# exp = exp_all[2]
# 
# hot <- read.table("hotspots_dry.txt", header = T)
# 
# # load data with cases and population
# type <- "_dry"
# case_pop <- read.table(paste0("cases_pop", exp, type, ".txt"))
# names(case_pop) <- c("location.id","cases","population")
# 
# # merge with cases and population
# nclust3 <- inner_join(nclust2, case_pop, by="location.id")
# 
# # aggregate cases and population by cluster
# nclust4 <- nclust3 %>% group_by(ncluster.cluster) %>%
#   summarize(cases = sum(cases), population = sum(population))
# 
# ## Get aggregated counts of population and cases for each county
# population <- tapply(nclust3$population,nclust3$ncluster.cluster,sum)
# cases <- tapply(nclust3$cases,nclust3$ncluster.cluster,sum)
# 
# ## Set Parameters
# pop.upper.bound <- 0.15
# n.simulations <- 9999
# alpha.level <- 0.05
# plot <- FALSE
# 
# ## Kulldorff using Binomial likelihoods
# binomial <- kulldorff(geo, cases, population, NULL, pop.upper.bound, n.simulations,
#                       alpha.level, plot)
# cluster <- binomial$most.likely.cluster$location.IDs.included
# cluster
# 
# cl_data <- geo[geo$location.id %in% cluster, ]
# cl_data$cluster <- 1
# names(cl_data)[1] <- "ncluster.cluster"
# 
# nclust5 <- left_join(nclust2, cl_data[,-c(2:3)], by="ncluster.cluster")
# nclust5$cluster[is.na(nclust5$cluster)] <- 0
# nclust5$cluster <- as.factor(nclust5$cluster)
# 
# #par(mfrow=c(2,2))
# # plot(nclust5$locx, nclust5$locy, col=nclust5$cluster, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")
# # cluster
# # table(nclust5$cluster)
# # 
# # plot(geo$locx, geo$locy, col="black", pch=1, main="simulated", xlab="",ylab="",xaxt = "n",yaxt = "n")
# # points(hot$locx, hot$locy, col=hot$hotspot, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")
# 
# names(nclust5)[5] <- paste0("kull",exp)
# hotnew <- inner_join(hot, nclust5[,c(2,5)], by="location.id")
# write.table(hotnew, "hotspots_dry.txt")
# 
# #at 3 year
# exp = exp_all[3]
# 
# hot <- read.table("hotspots_dry.txt", header = T)
# 
# # load data with cases and population
# type <- "_dry"
# case_pop <- read.table(paste0("cases_pop", exp, type, ".txt"))
# names(case_pop) <- c("location.id","cases","population")
# 
# # merge with cases and population
# nclust3 <- inner_join(nclust2, case_pop, by="location.id")
# 
# # aggregate cases and population by cluster
# nclust4 <- nclust3 %>% group_by(ncluster.cluster) %>%
#   summarize(cases = sum(cases), population = sum(population))
# 
# ## Get aggregated counts of population and cases for each county
# population <- tapply(nclust3$population,nclust3$ncluster.cluster,sum)
# cases <- tapply(nclust3$cases,nclust3$ncluster.cluster,sum)
# 
# ## Set Parameters
# pop.upper.bound <- 0.15
# n.simulations <- 9999
# alpha.level <- 0.05
# plot <- FALSE
# 
# ## Kulldorff using Binomial likelihoods
# binomial <- kulldorff(geo, cases, population, NULL, pop.upper.bound, n.simulations,
#                       alpha.level, plot)
# cluster <- binomial$most.likely.cluster$location.IDs.included
# cluster
# 
# cl_data <- geo[geo$location.id %in% cluster, ]
# cl_data$cluster <- 1
# names(cl_data)[1] <- "ncluster.cluster"
# 
# nclust5 <- left_join(nclust2, cl_data[,-c(2:3)], by="ncluster.cluster")
# nclust5$cluster[is.na(nclust5$cluster)] <- 0
# nclust5$cluster <- as.factor(nclust5$cluster)
# 
# #par(mfrow=c(2,2))
# # plot(nclust5$locx, nclust5$locy, col=nclust5$cluster, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")
# # cluster
# # table(nclust5$cluster)
# # 
# # plot(geo$locx, geo$locy, col="black", pch=1, main="simulated", xlab="",ylab="",xaxt = "n",yaxt = "n")
# # points(hot$locx, hot$locy, col=hot$hotspot, pch=1, main="cases", xlab="",ylab="",xaxt = "n",yaxt = "n")
# 
# names(nclust5)[5] <- paste0("kull",exp)
# hotnew <- inner_join(hot, nclust5[,c(2,5)], by="location.id")
# write.table(hotnew, "hotspots_dry.txt")

# MID HOTSPOT
exp_all = c("T03","T","T3y")

ncenter = 15

hh <- 2000

exp = exp_all[1]

hot <- read.table("hotspots.txt", header = T)

# load data with cases and population
type <- "_mid"
case_pop <- read.table(paste0("cases_pop", exp, type, ".txt"))
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
write.table(hotnew, "hotspots_mid.txt")

#at 1 year
exp = exp_all[2]

hot <- read.table("hotspots_mid.txt", header = T)

# load data with cases and population
type <- "_mid"
case_pop <- read.table(paste0("cases_pop", exp, type, ".txt"))
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
write.table(hotnew, "hotspots_mid.txt")

#at 3 year
exp = exp_all[3]

hot <- read.table("hotspots_mid.txt", header = T)

# load data with cases and population
type <- "_mid"
case_pop <- read.table(paste0("cases_pop", exp, type, ".txt"))
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
write.table(hotnew, "hotspots_mid.txt")

# LAG HOTSPOT
exp_all = c("T03","T","T3y")

ncenter = 15

hh <- 2000

exp = exp_all[1]

hot <- read.table("hotspots.txt", header = T)

# load data with cases and population
type <- "_lag"
case_pop <- read.table(paste0("cases_pop", exp, type, ".txt"))
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
write.table(hotnew, "hotspots_lag.txt")

#at 1 year
exp = exp_all[2]

hot <- read.table("hotspots_lag.txt", header = T)

# load data with cases and population
type <- "_lag"
case_pop <- read.table(paste0("cases_pop", exp, type, ".txt"))
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
write.table(hotnew, "hotspots_lag.txt")

#at 3 year
exp = exp_all[3]

hot <- read.table("hotspots_lag.txt", header = T)

# load data with cases and population
type <- "_lag"
case_pop <- read.table(paste0("cases_pop", exp, type, ".txt"))
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
write.table(hotnew, "hotspots_lag.txt")

}

