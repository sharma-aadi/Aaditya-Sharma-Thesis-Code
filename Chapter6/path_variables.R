## Set PATH to project folder in case of R PATH configuration error
setwd("./")

## Required packages
required_packages <- paste0("./required_packages.R")

## File locations for all the scripts to be sourced 
Simulation <- paste0("./Simulation/Model_v1_Simulate_Data.R")
In_Silico_Digestion <- paste0("./Simulation/Model_v1_In_Silico_Digest.R")
Diagnostics <- paste0("./Simulation/Model_v1_Diagnostics.R")
Inference <- paste0("./Inference/Model_v1_Inference.R") 
Visualization <- paste0("./Visualization/Model_v1_Visualization_Main.R") 

## This is where the proteomes for simulating protein digestion are
Proteome_Directory <- "./" 

## This is where the samplers for generating MCMC chains live
Sampler_Directory <- paste0("./Inference/Samplers/")
