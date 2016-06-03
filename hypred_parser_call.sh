#!/bin/bash

## Code to run a parsing function on the simulation results
# This code should output 6 "RData" files, 1 for each population:TP formation
## combination. There are 3 populations (MN, ND, MNxND) and 2 TP formations
## (random and cumulative)

# Load modules
module load R

# Set working directory
cd /panfs/roc/groups/6/smithkp/neyhartj/Genomic_Selection/Simulations/BarleySimGS-TPUpdate

# Set the experiment
exp=MAF_experiment
# Set the directory containing the files to parse
parsedir=Files/$exp

# Set window or cumulative to control flow
justwindow=true

# Control flow based on experiment
if [ "$exp" = "Base_experiment" ]; then

  # For each combination, declare a filename, then launch the R script
  
  # MN and random
  filename="simulation_results_q100_sel0.1_popmakeup-MN_tpformation-window_collective_300416.RData"
  files=$(find $parsedir -maxdepth 1 -name "simulation_results*_popmakeup-MN_*tpformation-window*")
  
  Rscript Code/hypred_simulation_parse.R $(echo $filename $files)
  
  # ND and random
  #filename="simulation_results_q100_sel0.1_popmakeup-ND_tpformation-window_collective_300416.RData"
  #files=$(find $parsedir -maxdepth 1 -name "simulation_results*_popmakeup-ND_*tpformation-window*")
  
  #Rscript Code/hypred_simulation_parse.R $(echo $filename $files)
  
  # MNxND and random
  #filename="simulation_results_q100_sel0.1_popmakeup-MNxND_tpformation-window_collective_300416.RData"
  #files=$(find $parsedir -maxdepth 1 -name "simulation_results*_popmakeup-MNxND_*tpformation-window*")
  
  #Rscript Code/hypred_simulation_parse.R $(echo $filename $files)

  # Exit if just the window sub-experiment is completed
  if [ "$justwindow" = true ]; then

	echo "Running the parse script on the window results. Exiting..." && exit
  else
  
  # MN and cumulative
  filename="simulation_results_q100_sel0.1_popmakeup-MN_tpformation-cumulative_collective_300416.RData"
  files=$(find $parsedir -maxdepth 1 -name "simulation_results*_popmakeup-MN_*tpformation-cumulative*")
  
  Rscript Code/hypred_simulation_parse.R $(echo $filename $files)
  
  # ND and cumulative
  filename="simulation_results_q100_sel0.1_popmakeup-ND_tpformation-cumulative_collective_300416.RData"
  files=$(find $parsedir -maxdepth 1 -name "simulation_results*_popmakeup-ND_*tpformation-cumulative*")
  
  Rscript Code/hypred_simulation_parse.R $(echo $filename $files)
  
  # MNxND and cumulative
  filename="simulation_results_q100_sel0.1_popmakeup-MNxND_tpformation-cumulative_collective_300416.RData"
  files=$(find $parsedir -maxdepth 1 -name "simulation_results*_popmakeup-MNxND_*tpformation-cumulative*")
  
  Rscript Code/hypred_simulation_parse.R $(echo $filename $files)

  fi
  
else 
  
  if [ "$exp" = "MAF_experiment" ]; then
  
    # MN and random
    filename="simulation_results_q100_sel0.1_popmakeup-MN_tpformation-window_maf_collective.RData"
    files=$(find $parsedir -maxdepth 1 -name "simulation_results*_popmakeup-MN_*tpformation-window_maf*")
    
    Rscript Code/hypred_simulation_parse.R $(echo $filename $files)
    
    # ND and random
    filename="simulation_results_q100_sel0.1_popmakeup-ND_tpformation-window_maf_collective.RData"
    files=$(find $parsedir -maxdepth 1 -name "simulation_results*_popmakeup-ND_*tpformation-window_maf*")
    
    Rscript Code/hypred_simulation_parse.R $(echo $filename $files)
    
    # MNxND and random
    filename="simulation_results_q100_sel0.1_popmakeup-MNxND_tpformation-window_maf_collective.RData"
    files=$(find $parsedir -maxdepth 1 -name "simulation_results*_popmakeup-MNxND_*tpformation-window_maf*")
    
    Rscript Code/hypred_simulation_parse.R $(echo $filename $files)
  
  # Close the else-if statement
  fi
  
fi

  
