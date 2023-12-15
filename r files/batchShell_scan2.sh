#!/bin/bash

#SBATCH --job-name=OM_results                   														
#SBATCH --cpus-per-task=1                  													
#SBATCH --mem-per-cpu=10G              														

#SBATCH --qos=6hours          																

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/scicore/home/penny/malinga/malaria_hotspots/simulations_nov/error_out/outn%A_%a.out   							
#SBATCH --error=/scicore/home/penny/malinga/malaria_hotspots/simulations_nov/error_out/errn%A_%a.err

#SBATCH --array=1-54%54

ml R/3.6.0-foss-2018b

SEEDFILE=$HOME/malaria_hotspots/simulations_nov/file_names.txt
SEED=$(sed -n ${SLURM_ARRAY_TASK_ID}p $SEEDFILE)

Rscript /scicore/home/penny/malinga/malaria_hotspots/simulations_nov/statscans_hotspots2.R $SLURM_ARRAY_TASK_ID

