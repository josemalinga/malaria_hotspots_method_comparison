##################################################################
#Author: Josephine Malinga
#Date: Original (2019)
#
#Updated September 30th 2021
#
##################################################################

# RUN THIS THROUGH THE CLUSTER (sbatch batchShell1.sh) - check if this can be run in parallel or not

rm(list=ls())
#install.packages('geoR', dependencies = TRUE)

#create folders
#create replicates and their files
for (i in 33:35){
    #i=2
    # User 
    user = strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
    
    # Working directory
    setwd(paste0("/scicore/home/penny/",user,"/malaria_hotspots/simulations_nov"))
    
    exp <- paste0("rep",i);exp
    sed <- i*42
    set.seed(sed)
    
    #create the foldres to store the simulation codes + results
    source("codefiles/create_folders.R")
    
    setwd(paste0(getwd(),"/",exp))
    getwd()
    par(mfrow=c(2,2))
    
    #create the datasets we need for the simulations.
    #output includes a grid of homesteads dataset and the calculated distance matrices
    source("codefiles/inputs_hotspots.R")
    
    #define the household considered as hotspots - only run once unless changing everything
    #output includes datasets defining the different shapes and sizes of hotspots
    source("codefiles/create_hotspots.R")
    
    #Run the simulation model with the different hotspot types
    source("codefiles/single_hotspots.R")
    source("codefiles/multiple_hotspots.R")
    source("codefiles/decay_hotspots.R")
    source("codefiles/river_hotspots.R")
    source("codefiles/short_hotspots.R")
    source("codefiles/long_hotspots.R")
    # source("codefiles/season_hotspots.R")
    source("codefiles/mid_hotspots.R")
    source("codefiles/lag_hotspots.R")
    
    #Run the simulation model with the different hotspot types 
    #- to run Spatial Epi 
    
    # source("codefiles/prevmap_bayes.R")
    # source("codefiles/bayes_M.R")
    # source("codefiles/bayes_D.R")
    # source("codefiles/bayes_R.R")
    # source("codefiles/bayes_short.R")
    # source("codefiles/bayes_long.R")
    # source("codefiles/bayes_seas.R")
    # source("codefiles/bayes_dryseas.R")
    
    #Running satSCAN in R - is this possible
    
   
    #Running FlexSCAN in R - is this possible
    
}
    
