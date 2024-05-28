#######################################################################################################################
## Run several chains
#######################################################################################################################

candidate_proposal_sd <- c(1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.003, 0.0005, 0.0001)

sample_chain_sd_repeat<-array(data = 0, dim = c(length(candidate_proposal_sd), Number_of_Chains, Number_of_Iterations))
ap_chain_sd_repeat <- array(data = 0, dim = c(length(candidate_proposal_sd), Number_of_Chains, Number_of_Iterations))

for(sd in 1:length(candidate_proposal_sd)){
    print(sd)
    for(c in 1:Number_of_Chains){
        print(candidate_proposal_sd[sd])
        print(c)
        
        samples <- MH_Sampler(data = model.data,
                              init = "Auto",
                              iter = Number_of_Iterations,
                              warmup = Burn_in_Iterations,
                              thin = 1,
                              chains = 1,
                              seed = c,
                              proposal_sd = candidate_proposal_sd[sd])
        
        sample_chain_sd_repeat[sd,c,] <- samples$samples[1,1,]
        ap_chain_sd_repeat[sd,c,] <- samples$am[1,1,]
        
        
    }
    
}

sample_tibble <- driver::gather_array(sample_chain_sd_repeat, Sampled_Log_Intensity, Proposal_SD, Chain, Iteration)

## Create proposal sd labels
Proposal_SD.labs <- paste0("sd_",candidate_proposal_sd)
names(Proposal_SD.labs) <- 1:length(candidate_proposal_sd)

## Create chain labels
Chain.labs <- paste0("chain_",1:Number_of_Chains)
names(Chain.labs) <- 1:Number_of_Chains


ggplot(data = (sample_tibble %>% 
                   filter(Proposal_SD %in% 1:9)),
       aes(x = Iteration, y = Sampled_Log_Intensity, color = as.factor(Chain))) +
    facet_grid(Proposal_SD~Chain, scales = "free",
               labeller = labeller(Proposal_SD = Proposal_SD.labs,
                                   Chain = Chain.labs)) +
    viridis::scale_color_viridis(name = 'Chain_ID', na.value="grey", option = "viridis", discrete=TRUE) +
    geom_hline(aes(yintercept = eta[1, 1]), size = 1, color = c("blue")) +
    labs(title = "MCMC Traces for several chains generated using different proposal standard deviation",
         subtitle = "Blue line denotes true value") +
    geom_line(size = 0.75) +
    theme(axis.text.x = element_text(angle = 90))

ggplot(data = (sample_tibble %>%
                   filter(Proposal_SD %in% 1:9)),
       aes(x = Sampled_Log_Intensity, color = as.factor(Chain), fill = as.factor(Chain))) +
    facet_grid(Proposal_SD~Chain, scales = "free",
               labeller = labeller(Proposal_SD = Proposal_SD.labs,
                                   Chain = Chain.labs)) +
    viridis::scale_color_viridis(name = 'Chain_ID', na.value="grey", option = "viridis", discrete=TRUE) +
    viridis::scale_fill_viridis(name = 'Chain_ID', na.value="grey", option = "viridis", discrete=TRUE) +
    geom_vline(aes(xintercept = eta[1, 1]), size = 1, color = c("blue")) +
    labs(title = "Posterior histograms for several chains generated using different proposal standard deviation",
         subtitle = "Blue line denotes true value") +
    geom_histogram(size = 0.75) +
    theme(axis.text.x = element_text(angle = 90))

acf(x = model_fit$samples[1,1,])

#######################################################################################################################
## Plotting peptide intensities for a given protein
#######################################################################################################################

test <- Y[!(rownames(Y) == "EIKPLYEQL"),]
test <- tibble(rownames_to_column(as.data.frame(test), "Sequence"))
test <- test %>% pivot_longer(-c(Sequence), names_to = "Sample", values_to = "Intensity") %>% 
            mutate(Intensity = log2(Intensity))
ggplot(data = test, aes(x = as_factor(Sample), y = Intensity))+
    geom_boxplot()+
        theme_linedraw()
