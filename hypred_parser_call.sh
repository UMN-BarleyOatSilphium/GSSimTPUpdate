#!/bin/bash

#PBS -l walltime=24:00:00,mem=22gb,nodes=1:ppn=1
#PBS -N hypred_simulations_parse_data
#PBS -M neyha001@umn.edu
#PBS -m abe
#PBS -r n


## Code to run a parsing function on the simulation results
# This code should output 6 "RData" files, 1 for each population:TP formation
## combination. There are 3 populations (MN, ND, MNxND) and 2 TP formations
## (random and cumulative)

# Load modules
module load R

# Set working directory
cd /panfs/roc/groups/6/smithkp/neyhartj/Genomic_Selection/Simulations/BarleySimGS-TPUpdate

# Set the directory containing the files to parse
parsedir=Files/Base_experiment

# For each combination, declare a filename, then launch the R script

# MN and random
filename="simulation_results_q100_sel0.1_popmakeup-MN_tpformation-window_collective_300416.RData"
files=$(find $parsedir -maxdepth 1 -name "simulation_results*_popmakeup-MN_*tpformation-window*")

Rscript Code/hypred_simulation_parse.R $(echo $filename $files)

# ND and random
filename="simulation_results_q100_sel0.1_popmakeup-ND_tpformation-window_collective_300416.RData"
files=$(find $parsedir -maxdepth 1 -name "simulation_results*_popmakeup-ND_*tpformation-window*")

Rscript Code/hypred_simulation_parse.R $(echo $filename $files)

# MNxND and random
filename="simulation_results_q100_sel0.1_popmakeup-MNxND_tpformation-window_collective_300416.RData"
files=$(find $parsedir -maxdepth 1 -name "simulation_results*_popmakeup-MNxND_*tpformation-window*")

Rscript Code/hypred_simulation_parse.R $(echo $filename $files)

exit

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
