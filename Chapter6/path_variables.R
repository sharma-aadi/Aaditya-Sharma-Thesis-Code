## Main is where all the non-archived code lives
Main <- "/Users/Aadi/Desktop/PhD\ Code\ Repository/Main"

## Directory where all the modelling code is
Modelling_Directory <- paste0(Main, "/Year_2/Combined_Repository/Proteomics/Modelling/")

## Project is the directory where the current project lives
Project <-paste0(Modelling_Directory,"Metropolis_Hastings_Implementations/Model_v1/Metropolis_Hastings_Implementation")

## Set PATH to project folder in case of R PATH configuration error
setwd(Project)

## Required packages
required_packages <- paste0(Project,"/required_packages.R")

## File locations for all the scripts to be sourced 
Simulation <- paste0(Project, "/Simulation/Model_v1_Simulate_Data.R")
In_Silico_Digestion <- paste0(Project, "/Simulation/Model_v1_In_Silico_Digest.R")
Diagnostics <- paste0(Project, "/Simulation/Model_v1_Diagnostics.R")
Inference <- paste0(Project, "/Inference/Model_v1_Inference.R") 
Visualization <- paste0(Project, "/Visualization/Model_v1_Visualization_Main.R") 

## This is where the proteomes for simulating protein digestion are
Proteome_Directory <- "/Users/Aadi/Desktop/PhD_Project/Data/Proteomes/" 

## Processed MaxQuant data files used for estimating model hyperparameters
Data <- "/Users/Aadi/Desktop/PhD_Project/Data/Proteomics_Sanquin_LFQ/PLT_2" # Unused

## This is where the samplers for generating MCMC chains live
Sampler_Directory <- paste0(Project, "/Inference/Samplers/")
