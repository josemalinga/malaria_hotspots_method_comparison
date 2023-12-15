#!/bin/bash -l

#SBATCH --cpus-per-task=1                  													
#SBATCH --mem-per-cpu=50G              														

#SBATCH --time=00:25:00        																
#SBATCH --qos=30min   																

#SBATCH --output=/scicore/home/schindle/malinga/mosq_movement/simDatahse_2P_10/out1.txt     							
#SBATCH --error=/scicore/home/schindle/malinga/mosq_movement/simDatahse_2P_10/err1.txt     							

#SBATCH --array=1-6500%1

#SBATCH -n 1

name=$(sed -n "$SLURM_ARRAY_TASK_ID"p namelist1000)

echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID


cp $name aggTrapSimul.txt

g++ -c -g -l include mosqmove_nelder_new.cpp
g++ mosqmove_nelder_new.o asa047.o -lm

./a.out

cp results.txt findings/$name




