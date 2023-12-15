#!/bin/bash

#SBATCH --job-name=OM_results                   														
#SBATCH --cpus-per-task=1                  													
#SBATCH --mem-per-cpu=5G              														

#SBATCH --qos=1day          																

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/scicore/home/penny/malinga/malaria_hotspots/simulations_nov/error_out/outn.txt     							
#SBATCH --error=/scicore/home/penny/malinga/malaria_hotspots/simulations_nov/error_out/errn.txt

#load your required modules below
#################################


#export your required environment variables below
#################################################


#add your command lines below
#############################
ml R/3.6.0-foss-2018b

Rscript /scicore/home/penny/malinga/malaria_hotspots/simulations_nov/statscans_hotspots.R 
