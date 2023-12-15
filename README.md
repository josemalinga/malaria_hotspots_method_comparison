# malaria_hotspots_method_comparison
Using a simulation model to compare different statistical methods to detect hotspots of malaria transmission

## Code for running the hotspots simulation models
Updated by 
Name: Josephine Malinga;
Date: 15.12.2023;
Email: josephine.malinga@unibas.ch

## Steps 0
For each replicate the code to RUN is the main code “r files/code_files.R” whose components are outlined below
 - Set seed for each run
 
## Steps 1
- Runs the “inputs_hotspots.R” code: this code creates the datasets we need for the simulations.
- The output includes a grid of homesteads dataset and the distance matrices
 
## Steps 2
- Runs the “create_hotspots.R” code: this code creates the the baseline areas of higher transmission that we need (defines the household considered as hotspots). 
- The output include datasets defining the different shapes and sizes.
 
## Step 3
- Runs the simulation models according to hotspot characteristic. There are 6 types of hotspots; single, multiple, decay, river, short, long, seasonal
- The output includes the cross sectional datasets at 3mo, 1y and 3y

## Step 4
- Runs the simulation models using the case files from Step 3 to obtain the bayesian results. 
- The output includes the cross sectional datasets for all the shapes/sizes with a new variables for hotspot status for each household

## Step 5
- Using SaTScan detect the households in the hotspots for each shape/size
- Run the code “extract_hotspots.R” to get new variables for hotspot status for each household as defined 
 
