#######################################################################################################################
## Prepare Data [Note: The Y, M, and eta matrices are all transposes of what we have in the STAN implementation]
#######################################################################################################################

model.data <- list(P = length(unique(Peptide_Table$Protein)),
                   T = length(unique(Peptide_Table$Peptide)),
                   S = Number_of_Samples,
                   Y = t(Y), 
                   M = t(M))

#######################################################################################################################
## Load Metropolis-Hastings Sampler
#######################################################################################################################

source(paste0(Sampler_Directory, Sampler_Type, ".R"))

#######################################################################################################################
## Set Sampling Seed
#######################################################################################################################

## Set Seed and sample
if(Sampling_Seed != 0){
    cat("\n")
    
    print(paste0("Setting sampling seed: ", Sampling_Seed))
    
    cat("\n")
    set.seed(Sampling_Seed)
}

#######################################################################################################################
## Run MCMC
#######################################################################################################################

## model_fit individual runs

model_fit <- MH_Sampler(data = model.data,
                   init = Parameter_Initialization,
                   iter = Number_of_Iterations,
                   warmup = Burn_in_Iterations,
                   thin = Thinning,
                   chains = Number_of_Chains,
                   proposal_sd = Proposal_Standard_Deviation,
                   tuning_speed = Tuning_Speed)

#######################################################################################################################
## Print MCMC run information
#######################################################################################################################

## Print tuned proposal variance matrix
cat("\nMatrix of tuned variances [Sample X Protein]:\n")
print(model_fit$tuned_proposal_variance)

## Print acceptance rate matrix
cat("\nMatrix of average acceptance rates [Sample X Protein]:\n")
print(model_fit$acceptance_rate)

## Add space
cat("\n")
