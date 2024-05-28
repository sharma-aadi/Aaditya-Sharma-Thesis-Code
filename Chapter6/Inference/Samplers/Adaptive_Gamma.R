#######################################################################################################################
## Metropolis-Hastings Sampler
#######################################################################################################################

MH_Sampler <- function(data,
                       init,
                       iter = 2000,
                       warmup = 0,
                       thin = 1,
                       chains = 1,
                       proposal_sd = 1,
                       tuning_speed = 0.25){
    
    # Parameters [Type in later]
    ## data 
    ## init
    ## iter
    ## warmup
    ## thin
    ## chains
    ## proposal_sd

    Samples <- data$S
    Proteins <- data$P
    Peptides <- data$T
    
    # Initialize parameter array
    eta <- array(data = 0, dim = c(Samples, Proteins, iter))
    
    # Initialize matrix of proposal variances
    proposal_sd <- array(data = proposal_sd, dim = c(Samples, Proteins))
    
    
    if(init == "Mean"){
        cat("\n")
        print(paste0("Parameter initialization is set to ", tolower(Parameter_Initialization), "."))
        cat("\n")
        ## Initialize all the eta's at the mean of the simulating distribution
        eta[,,1] <- shape_protein*(1/rate_protein)
    }
    
    if(init == "Random"){
        cat("\n")
        print(paste0("Parameter initialization is set to ", tolower(Parameter_Initialization), "."))
        cat("\n")
        for(i in 1:Samples){
            for(j in 1:Proteins){
                eta[i,j,1] <- runif(n = 1, min = 15, max = 20)
            }
        }
        cat("\nInitial parameter matrix [Sample X Protein]:\n")
        print(eta[,,1])
        cat("\n")
    }
    
    # Create a variable for current proposal
    current_proposal <- numeric(length = 1)
    
    # Create a variable for MH ratio
    MH_Ratio <- numeric(length = 1)
    
    # Create a vector to keep track of the acceptance ratio [First iteration is the initialization so probability == 1]
    acceptance_probability <-  array(data = 0, dim = c(Samples, Proteins, iter))
    acceptance_probability[,,1] <- 1
    
    # Create other variables required for sampling
    
    log_term_old <- log_term_new <- numeric(length = 1)
    
    eta_vector_old <- eta_vector_new <- numeric(length = Proteins)
    
    multiplicity_vector <- integer(length = Proteins)
    
    acceptance_monitor <- array(data = 0, dim = c(Samples, Proteins))
    
    # Run MCMC
    for(i in 2:(iter)){
        
        if(i%%100 == 0){
            
            if(i%%1000 == 0){
                print(i)
            }
            
            acceptance_monitor <- apply(acceptance_probability[,,((i-100 + 1):(i))], c(1,2), mean)
            
            for(sam_check in 1:Samples){
                for (prot_check in 1:Proteins) {
                    if(acceptance_monitor[sam_check, prot_check] < 0.4){
                        proposal_sd[sam_check, prot_check] <- 
                            (proposal_sd[sam_check, prot_check])*(1-tuning_speed) #reduction ratio
                    } else if(acceptance_monitor[sam_check, prot_check] > 0.5){
                        proposal_sd[sam_check, prot_check] <- 
                            (proposal_sd[sam_check, prot_check])*(1+tuning_speed) #increment ratio
                    }
                }
            }
        }
        
        for(s in 1:(Samples)){
            
            for(p in 1:(Proteins)){
                
                ## Set terms to be computed to 0!
                denominator <- numerator <- numeric(length = 1) 
                protein_prior_term <- numeric(length = 1) 
                
                ## Propose a new value
                current_proposal <- rnorm(n = 1, mean = eta[s,p,i-1], sd = proposal_sd[s,p])
                
                ## Create current and proposed eta vectors
                eta_vector_old <- eta_vector_new <- eta[s,,i-1]
                eta_vector_new[p] <- current_proposal
                
                ## Compute the difference in the gamma prior term [does not depend on the]
                
                protein_prior_term <- -5*(eta_vector_new[p] - eta_vector_old[p]) + 
                    87*(log(eta_vector_new[p]) - log(eta_vector_old[p]))
                
                ## Compute the log terms for all the peptides
                for(t in 1:(Peptides)){
                    
                    multiplicity_vector <- data$M[,t]
                    
                    if(data$M[p,t] != 0){
                        
                        log_term_old <- log(
                            dot(
                                x = exp(eta_vector_old),
                                y = multiplicity_vector
                            )
                        )
                        
                        log_term_new <- log(
                            dot(
                                x = exp(eta_vector_new),
                                y = multiplicity_vector
                            )
                        )
                        
                        denominator <- denominator - 50*((log_term_old)^2) + 
                            100*(log(data$Y[s,t])*log_term_old) + 0.5*log_term_old
                        
                        
                        numerator <- numerator - 50*((log_term_new)^2) +
                            100*(log(data$Y[s,t])*log_term_new) + 0.5*log_term_new
                        
                    }
                    
                }
                
                MH_Ratio <- numerator - denominator + protein_prior_term
                
                if(log(runif(n = 1)) <= MH_Ratio){
                    eta[s,p,i] <- current_proposal
                    acceptance_probability[s,p,i] <- 1
                } else {
                    eta[s,p,i] <- eta[s,p,i-1]
                    acceptance_probability[s,p,i] <- 0
                }
                
            }
            
        }
    }
    
    return(list(samples = eta,
                acceptance_rate = apply(acceptance_probability[,,(warmup+1):iter], c(1,2), mean), 
                tuned_proposal_variance = proposal_sd)
           )
}


