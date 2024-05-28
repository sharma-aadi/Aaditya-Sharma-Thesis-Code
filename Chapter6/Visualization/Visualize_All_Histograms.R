######################################################################################################################
## Function - 2: Visualize_All_Histograms()
#######################################################################################################################

Visualize_All_Histograms <- function(){
    
    plotting_eta <- eta
    
    rownames(plotting_eta) <- 1:Number_of_Proteins
    colnames(plotting_eta) <- 1:Number_of_Samples
    
    true_value_tibble <- driver::gather_array(plotting_eta, True_Value, Protein, Sample)
    
    ## Create Sample labels
    Sample.labs <- paste0("Sample ", 1:Number_of_Samples)
    names(Sample.labs) <- 1:(Number_of_Samples)
    
    ## Create Protein labels
    # Protein.labs <- paste0("Protein ",1:Number_of_Proteins) 
    Protein.labs <- rownames(eta)
    names(Protein.labs) <- 1:Number_of_Proteins
    
    plotting_tibble <- driver::gather_array(model_fit$samples, Sampled_Log_Intensity, Sample, Protein, Iteration)
    
    posterior_statistics <- plotting_tibble %>% filter(Iteration > Burn_in_Iterations) %>% 
        group_by(Protein, Sample) %>% 
        summarise(Mean_Sampled_Log_Intensity = mean(Sampled_Log_Intensity),
                  Median_Sampled_Log_Intensity = median(Sampled_Log_Intensity))

    All_Histograms <- ggplot(data = (plotting_tibble %>% filter(Iteration > Burn_in_Iterations)),
                             aes(x = Sampled_Log_Intensity,
                                 color = as.factor(Sample),
                                 fill =  as.factor(Sample))) +
        geom_histogram(binwidth = 0.01) +
        geom_vline(data = (true_value_tibble), aes(xintercept = True_Value), size  = 1) +
        geom_vline(data = (posterior_statistics),
                   aes(xintercept = Mean_Sampled_Log_Intensity), size = 1, linetype = 2, color = "red")+
        facet_grid(Sample~Protein, scales = "free", 
                   labeller = labeller(Sample = Sample.labs,
                                       Protein = Protein.labs)) +
        viridis::scale_color_viridis(name = 'Count', na.value="grey", option = "viridis", discrete=TRUE) +
        viridis::scale_fill_viridis(name = 'Sample', na.value="grey", option = "viridis", discrete=TRUE) +
        theme_linedraw() +
        theme(axis.text.x = element_text(angle = 90),
              legend.position = "none")+
        labs(title = "Post-warmup histogram for protein log-intensity.",
             subtitle = "True value is in black and posterior mean is in dashed red.",
             y = "Count", x = "Sampled log-intensity") +
        theme(
            plot.title = element_text(size=20, face = "bold"),
            plot.subtitle = element_text(size=18, face = "bold"),
            axis.title.y = element_text(color="grey20", size=16, face="bold"),
            axis.title.x = element_text(color="grey20", size=16, face="bold"),
            axis.text.y = element_text(color="grey20", size=14, face="bold"),
            axis.text.x = element_text(color="grey20", size=14, face="bold", angle = 90),
            strip.text.x = element_text(size = 12),
            strip.text.y = element_text(size = 12))
    
    return(All_Histograms)
    
}
