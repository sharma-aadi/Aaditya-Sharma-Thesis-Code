#######################################################################################################################
## set seed for simulation 
#######################################################################################################################

set.seed(1000)

#######################################################################################################################
## Clear Global Environment
#######################################################################################################################

# rm(list=ls()[! ls() %in% c("load_packages","save_counter", "todays_date")])


#######################################################################################################################
## Set Path variables
#######################################################################################################################

source(file = "path_variables.R")

#######################################################################################################################
## Load Required Packages
####################################################################################################################### 

source(required_packages)

#######################################################################################################################
## Simulation Parameters
#######################################################################################################################

## Proteome file name (this should be in the proteomes directory)
Simulation_Proteome <- "Shared_Peptide_Dataset_4_1"

## Set the number of Samples and Proteins in the simulated dataset
Number_of_Samples <- 4                   # Specify the number of samples in the experiment
Number_of_Proteins <- 4                 # Specify the max number of proteins in the experiment 

### Protein abundance parameters [88,5 are defaults]
shape_protein <- 88          
rate_protein <- 5

## Specify the minimum and maximum length of the peptides in the simulated data table
Min_Length <-  7
Max_Length <- 40

Peptide_Length_Range <- c(Min_Length, Max_Length)  

## Specify the maximum number of missed cleavages in the simulated tryptic peptides
Max_Missed_Cleavages <- 0                   

## Peptide-level sd (0.1 is default)
Peptide_SD <- 0.25

#######################################################################################################################
## Sampling Parameters
#######################################################################################################################

## Choose type of Metropolis-Hastings sampler 
Sampler_Type <- "Adaptive_Uniform_Real_Flat"
Parameter_Initialization <- "Random"    # Either Mean or Random

## Seed Parameters 
Simulation_Seed <- 1
Sampling_Seed <- 1

## MCMC Parameters
Number_of_Iterations <- 6000  # default is 6000
Burn_in_Iterations <- 1000       # default is 2
Thinning <- 1
Number_of_Chains <- 1

## Metropolis Hastings Parameters 
Proposal_Standard_Deviation <- 1     # Intialization for all proposal variances [0.1 is default]
Tuning_Speed <- 0.1

#######################################################################################################################
## Simulate Data (Generates data matrix Y & model matrices eta and M)
#######################################################################################################################

source(Simulation)

#######################################################################################################################
## Conduct Inference (Runs MCMC sampler for posterior inference)
#######################################################################################################################

source(Inference)

#######################################################################################################################
## Load all the visualization functions
#######################################################################################################################

source(Visualization)
