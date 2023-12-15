# malaria_hotspots_method_comparison
Using a simulation model to compare different statistical methods to detect hotspots of malaria transmission

# Code for running the hotspots simulation models
J Malinga
15.02.2021

For each replicate the code to RUN is the main code “code_files.R” whose components are outlined below
 - Set seed for each run
 
# Steps 1
- Runs the “inputs_hotspots.R” code: this code creates the datasets we need for the simulations.
- The output includes a grid of homesteads dataset and the distance matrices
 
# Steps 2
- Runs the “create_hotspots.R” code: this code creates the the baseline areas of higher transmission that we need (defines the household considered as hotspots). 
- The output include datasets defining the different shapes and sizes.
 
# Step 3
- Runs the simulation models according to hotspot characteristic. There are 6 types of hotspots; single, multiple, decay, river, short, long, seasonal
- The output includes the cross sectional datasets at 3mo, 1y and 3y 
 
