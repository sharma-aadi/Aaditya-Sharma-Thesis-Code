#######################################################################################################################
##  Visualize_Single_Protein()
#######################################################################################################################

Visualize_Single_Protein <- function(Protein_Number){

    Protein_Index <- Protein_Number
    
    plotting_eta <- eta
    
    rownames(plotting_eta) <- 1:Number_of_Proteins
    colnames(plotting_eta) <- 1:Number_of_Samples
    
    true_value_tibble <- driver::gather_array(plotting_eta, True_Value, Protein, Sample)
    
    ## Create Plotting Tibble from Sampling array
    
    ## Create Sample labels
    Sample.labs <- paste0("Sample ", 1:Number_of_Samples)
    names(Sample.labs) <- 1:(Number_of_Samples)
    
    ## Create Protein labels
    Protein.labs <- paste0("Protein ",1:Number_of_Proteins)
    names(Protein.labs) <- 1:Number_of_Proteins
    
    plotting_tibble <- driver::gather_array(model_fit$samples, Sampled_Log_Intensity, Sample, Protein, Iteration)
    
    ## Sample-Wise traces
    
    Individual_Trace_Plots <- ggplot(data = (plotting_tibble %>% filter(Protein == Protein_Index) %>% 
                                                 filter(Iteration > Burn_in_Iterations)),
                                     aes(y = Sampled_Log_Intensity, x = Iteration,
                                         color = as.factor(Sample))) +
        geom_line() +
        geom_hline(data = (true_value_tibble %>% filter(Protein == Protein_Index)),
                   aes(yintercept = True_Value))+
        facet_grid(Sample~Protein, scales = "free", 
                   labeller = labeller(Sample = Sample.labs,
                                       Protein = Protein.labs)) +
        viridis::scale_color_viridis(name = 'Sample', na.value="grey", option = "viridis", discrete=TRUE) +
        theme_linedraw() +
        labs(y = "Sampled log-intensity") +
        theme(axis.text.x = element_text(angle = 90),
              legend.position = "none")
    
    ## Sample-Wise histograms
    
    Individual_Histograms <- ggplot(data = (plotting_tibble %>% filter(Protein == Protein_Index) %>%
                                                filter(Iteration > Burn_in_Iterations)),
                                    aes(x = Sampled_Log_Intensity, 
                                        color = as.factor(Sample),
                                        fill =  as.factor(Sample))) +
        geom_histogram(binwidth = 0.01) +
        geom_vline(data = (true_value_tibble %>% filter(Protein == Protein_Index)), aes(xintercept = True_Value)) +
        facet_grid(Sample~Protein, scales = "free", 
                   labeller = labeller(Sample = Sample.labs,
                                       Protein = Protein.labs)) +
        viridis::scale_color_viridis(name = 'Count', na.value="grey", option = "viridis", discrete=TRUE) +
        viridis::scale_fill_viridis(name = 'Sample', na.value="grey", option = "viridis", discrete=TRUE) +
        theme_linedraw() +
        labs(x = "Sampled log-intensity") +
        theme(axis.text.x = element_text(angle = 90),
              legend.position = "none")
    
    # Return GGPlot 
    return(Individual_Trace_Plots + Individual_Histograms + plot_layout(widths = c(5, 1)))
}