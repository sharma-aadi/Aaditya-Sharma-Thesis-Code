#######################################################################################################################
## Function - 1: Visualize_Single_Protein()
#######################################################################################################################

Visualize_Single_Protein <- function(Protein_Number = 1){
    
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
    # Protein.labs <- paste0("Protein ",1:Number_of_Proteins) 
    Protein.labs <- rownames(eta)
    names(Protein.labs) <- 1:Number_of_Proteins
    
    plotting_tibble <- driver::gather_array(model_fit$samples, Sampled_Log_Intensity, Sample, Protein, Iteration)
    
    posterior_statistics <- plotting_tibble %>% filter(Iteration > Burn_in_Iterations) %>% 
        group_by(Protein, Sample) %>% 
        summarise(Mean_Sampled_Log_Intensity = mean(Sampled_Log_Intensity),
                  Median_Sampled_Log_Intensity = median(Sampled_Log_Intensity))
    
    ## Sample-Wise traces
    
    Individual_Trace_Plots <- ggplot(data = (plotting_tibble %>% filter(Protein == Protein_Index) %>% 
                                                 filter(Iteration > Burn_in_Iterations)),
                                     aes(y = Sampled_Log_Intensity, x = Iteration,
                                         color = as.factor(Sample))) +
        geom_line() +
        geom_hline(data = (true_value_tibble %>% filter(Protein == Protein_Index)),
                   aes(yintercept = True_Value))+
        geom_hline(data = (posterior_statistics%>% filter(Protein == Protein_Index)),
                   aes(yintercept = Mean_Sampled_Log_Intensity), size = 1, linetype = 2, color = "red")+
        facet_grid(Protein~Sample, scales = "free", 
                   labeller = labeller(Sample = Sample.labs,
                                       Protein = Protein.labs)) +
        viridis::scale_color_viridis(name = 'Sample', na.value="grey", option = "viridis", discrete=TRUE) +
        theme_linedraw() +
        labs(y = "Sampled log-intensity", title = paste0("Trace plots for protein ", Protein_Index ," from chain 1"),
             subtitle = "True value is in solid black and posterior mean is in dashed red.") +
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
        geom_vline(data = (posterior_statistics%>% filter(Protein == Protein_Index)),
                   aes(xintercept = Mean_Sampled_Log_Intensity), size = 1, linetype = 2, color = "red")+
        facet_grid(Protein~Sample, scales = "free", 
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

######################################################################################################################
## Function - 2: Visualize_All_Histograms()
#######################################################################################################################

Visualize_All_Histograms <- function(Show.Info = F, Scale = "Logarithmic"){
    
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
    
    if(Scale == "Logarithmic"){
        All_Histograms <- ggplot(data = (plotting_tibble %>% filter(Iteration > Burn_in_Iterations)),
                                 aes(x = Sampled_Log_Intensity,
                                     color = as.factor(Sample),
                                     fill =  as.factor(Sample))) +
            geom_histogram(binwidth = 0.1) +
            geom_vline(data = (true_value_tibble), aes(xintercept = True_Value), size  = 1) +
            geom_vline(data = (posterior_statistics),
                       aes(xintercept = Mean_Sampled_Log_Intensity), size = 1, linetype = 2, color = "red")+
            facet_grid(Protein~Sample, scales = "free", 
                       labeller = labeller(Sample = Sample.labs,
                                           Protein = Protein.labs)) +
            labs(x = "Sampled log-intensity")+
            # labs(x = "Sampled log-intensity",
            #      title = "Post-warmup histograms for protein sampled log-intensity.")+
            viridis::scale_color_viridis(name = 'Count', na.value="grey", option = "viridis", discrete=TRUE) +
            viridis::scale_fill_viridis(name = 'Sample', na.value="grey", option = "viridis", discrete=TRUE) +
            theme_linedraw() +
            theme(axis.text.x = element_text(angle = 90),
                  legend.position = "none")+
            theme(
                plot.title = element_text(size=20, face = "bold"),
                plot.subtitle = element_text(size=18, face = "bold"),
                axis.title.y = element_text(size=14),
                axis.title.x = element_text(size=14),
                axis.text.y = element_text(size=12),
                axis.text.x = element_text(size=12, angle = 0),
                strip.text = element_text(size = 14, color = "black"),
                strip.background = element_rect(fill = "white", color = "black"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank())
    } else if(Scale == "Real"){
        All_Histograms <- ggplot(data = (plotting_tibble %>% filter(Iteration > Burn_in_Iterations)),
                                 aes(x = exp(Sampled_Log_Intensity),
                                     color = as.factor(Sample),
                                     fill =  as.factor(Sample))) +
            geom_histogram(bins = 15) +
            geom_vline(data = (true_value_tibble), aes(xintercept = exp(True_Value)), size  = 1) +
            geom_vline(data = (posterior_statistics),
                       aes(xintercept = exp(Mean_Sampled_Log_Intensity)), size = 1, linetype = 2, color = "red")+
            facet_grid(Protein~Sample, scales = "free", 
                       labeller = labeller(Sample = Sample.labs,
                                           Protein = Protein.labs)) +
            labs(x = "Sampled log-intensity")+
            # labs(x = "Sampled intensity",
            #      title = "Post-warmup histograms for sampled protein intensity.")+
            viridis::scale_color_viridis(name = 'Count', na.value="grey", option = "viridis", discrete=TRUE) +
            viridis::scale_fill_viridis(name = 'Sample', na.value="grey", option = "viridis", discrete=TRUE) +
            theme_linedraw() +
            theme(axis.text.x = element_text(angle = 90),
                  legend.position = "none")+
            theme(
                plot.title = element_text(size=20, face = "bold"),
                plot.subtitle = element_text(size=18, face = "bold"),
                axis.title.y = element_text(size=14),
                axis.title.x = element_text(size=14),
                axis.text.y = element_text(size=12),
                axis.text.x = element_text(size=12, angle = 0),
                strip.text = element_text(size = 14, color = "black"),
                strip.background = element_rect(fill = "white", color = "black"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank())
    }
    
    if(Show.Info == F){
        return(All_Histograms+ labs(y = "Count") )
    } else{
        return(All_Histograms+ labs(subtitle = "True value is in solid black and posterior mean is in dashed red.",
                                    y = "Count",
                                    caption = paste0("Simulation seed is ", Simulation_Seed, " and sampling seed is ",
                                                     Sampling_Seed, 
                                                     ". Initialization is set to ", Parameter_Initialization,
                                                     ".") ))
    }
}

#######################################################################################################################
## Function - 3: Visualize_All_Traces()
#######################################################################################################################

Visualize_All_Traces <- function(Show.Info = F, Show.Warmup = F, Scale = "Logarithmic"){
    
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
    
    if(Scale == "Logarithmic"){
        if(Show.Warmup){
            posterior_statistics <- plotting_tibble %>% filter(Iteration > 0) %>% 
                group_by(Protein, Sample) %>% 
                summarise(Mean_Sampled_Log_Intensity = mean(Sampled_Log_Intensity),
                          Median_Sampled_Log_Intensity = median(Sampled_Log_Intensity))
            
            All_Traces <- ggplot(data = (plotting_tibble %>% filter(Iteration > 0)),
                                 aes(y = (Sampled_Log_Intensity), x = Iteration,
                                     color = as.factor(Sample))) +
                geom_line() +
                geom_hline(data = (true_value_tibble), aes(yintercept = (True_Value)), size  = 1) +
                geom_hline(data = (posterior_statistics),
                           aes(yintercept = (Mean_Sampled_Log_Intensity)), size = 1, linetype = 2, color = "red")+
                facet_grid(Protein~Sample, scales = "free", 
                           labeller = labeller(Sample = Sample.labs,
                                               Protein = Protein.labs)) +
                labs(title = "Post-warmup traceplots for protein log-intensity.",
                     y = "Sampled log-intensity") +
                viridis::scale_color_viridis(name = 'Count', na.value="grey", option = "viridis", discrete=TRUE) +
                viridis::scale_fill_viridis(name = 'Sample', na.value="grey", option = "viridis", discrete=TRUE) +
                theme_linedraw() +
                theme(axis.text.x = element_text(angle = 0),
                      legend.position = "none") +
                theme(
                    plot.title = element_text(size=18, face = "bold"),
                    plot.subtitle = element_text(size=16, face = "bold"),
                    axis.title.y = element_text(color="grey20", size=16, face="bold"),
                    axis.title.x = element_text(color="grey20", size=16, face="bold"),
                    axis.text.y = element_text(color="grey20", size=14, face="bold"),
                    axis.text.x = element_text(color="grey20", size=14, face="bold", angle = 90),
                    strip.text.x = element_text(size = 12),
                    strip.text.y = element_text(size = 12))
        } else{
            posterior_statistics <- plotting_tibble %>% filter(Iteration > Burn_in_Iterations) %>% 
                group_by(Protein, Sample) %>% 
                summarise(Mean_Sampled_Log_Intensity = mean(Sampled_Log_Intensity),
                          Median_Sampled_Log_Intensity = median(Sampled_Log_Intensity))
            
            All_Traces <- ggplot(data = (plotting_tibble %>% filter(Iteration > Burn_in_Iterations)),
                                 aes(y = Sampled_Log_Intensity, x = Iteration,
                                     color = as.factor(Sample))) +
                geom_line() +
                geom_hline(data = (true_value_tibble), aes(yintercept = True_Value), size  = 1) +
                geom_hline(data = (posterior_statistics),
                           aes(yintercept = Mean_Sampled_Log_Intensity), size = 1, linetype = 2, color = "red")+
                facet_grid(Protein~Sample, scales = "free", 
                           labeller = labeller(Sample = Sample.labs,
                                               Protein = Protein.labs)) +
                labs(title = "Post-warmup traceplots for protein log-intensity.",
                     y = "Sampled log-intensity") +
                viridis::scale_color_viridis(name = 'Count', na.value="grey", option = "viridis", discrete=TRUE) +
                viridis::scale_fill_viridis(name = 'Sample', na.value="grey", option = "viridis", discrete=TRUE) +
                theme_linedraw() +
                theme(axis.text.x = element_text(angle = 0),
                      legend.position = "none") +
                theme(
                    plot.title = element_text(size=18, face = "bold"),
                    plot.subtitle = element_text(size=16, face = "bold"),
                    axis.title.y = element_text(color="grey20", size=16, face="bold"),
                    axis.title.x = element_text(color="grey20", size=16, face="bold"),
                    axis.text.y = element_text(color="grey20", size=14, face="bold"),
                    axis.text.x = element_text(color="grey20", size=14, face="bold", angle = 90),
                    strip.text.x = element_text(size = 12),
                    strip.text.y = element_text(size = 12))
        }
    } else if(Scale == "Real"){
        if(Show.Warmup){
            posterior_statistics <- plotting_tibble %>% filter(Iteration > 0) %>% 
                group_by(Protein, Sample) %>% 
                summarise(Mean_Sampled_Log_Intensity = mean(Sampled_Log_Intensity),
                          Median_Sampled_Log_Intensity = median(Sampled_Log_Intensity))
            
            All_Traces <- ggplot(data = (plotting_tibble %>% filter(Iteration > 0)),
                                 aes(y = exp(Sampled_Log_Intensity), x = Iteration,
                                     color = as.factor(Sample))) +
                geom_line() +
                geom_hline(data = (true_value_tibble), aes(yintercept = exp(True_Value)), size  = 1) +
                geom_hline(data = (posterior_statistics),
                           aes(yintercept = exp(Mean_Sampled_Log_Intensity)), size = 1, linetype = 2, color = "red")+
                facet_grid(Protein~Sample, scales = "free", 
                           labeller = labeller(Sample = Sample.labs,
                                               Protein = Protein.labs)) +
                labs(title = "Post-warmup traceplots for protein intensity.",
                     y = "Sampled intensity") +
                viridis::scale_color_viridis(name = 'Count', na.value="grey", option = "viridis", discrete=TRUE) +
                viridis::scale_fill_viridis(name = 'Sample', na.value="grey", option = "viridis", discrete=TRUE) +
                theme_linedraw() +
                theme(axis.text.x = element_text(angle = 0),
                      legend.position = "none") +
                theme(
                    plot.title = element_text(size=18, face = "bold"),
                    plot.subtitle = element_text(size=16, face = "bold"),
                    axis.title.y = element_text(color="grey20", size=16, face="bold"),
                    axis.title.x = element_text(color="grey20", size=16, face="bold"),
                    axis.text.y = element_text(color="grey20", size=14, face="bold"),
                    axis.text.x = element_text(color="grey20", size=14, face="bold", angle = 90),
                    strip.text.x = element_text(size = 12),
                    strip.text.y = element_text(size = 12))
        } else{
            posterior_statistics <- plotting_tibble %>% filter(Iteration > Burn_in_Iterations) %>% 
                group_by(Protein, Sample) %>% 
                summarise(Mean_Sampled_Log_Intensity = mean(Sampled_Log_Intensity),
                          Median_Sampled_Log_Intensity = median(Sampled_Log_Intensity))
            
            All_Traces <- ggplot(data = (plotting_tibble %>% filter(Iteration > Burn_in_Iterations)),
                                 aes(y = exp(Sampled_Log_Intensity), x = Iteration,
                                     color = as.factor(Sample))) +
                geom_line() +
                geom_hline(data = (true_value_tibble), aes(yintercept = exp(True_Value)), size  = 1) +
                geom_hline(data = (posterior_statistics),
                           aes(yintercept = exp(Mean_Sampled_Log_Intensity)), size = 1, linetype = 2, color = "red")+
                facet_grid(Protein~Sample, scales = "free", 
                           labeller = labeller(Sample = Sample.labs,
                                               Protein = Protein.labs)) +
                labs(title = "Post-warmup traceplots for protein intensity.",
                     y = "Sampled intensity") +
                viridis::scale_color_viridis(name = 'Count', na.value="grey", option = "viridis", discrete=TRUE) +
                viridis::scale_fill_viridis(name = 'Sample', na.value="grey", option = "viridis", discrete=TRUE) +
                theme_linedraw() +
                theme(axis.text.x = element_text(angle = 0),
                      legend.position = "none") +
                theme(
                    plot.title = element_text(size=18, face = "bold"),
                    plot.subtitle = element_text(size=16, face = "bold"),
                    axis.title.y = element_text(color="grey20", size=16, face="bold"),
                    axis.title.x = element_text(color="grey20", size=16, face="bold"),
                    axis.text.y = element_text(color="grey20", size=14, face="bold"),
                    axis.text.x = element_text(color="grey20", size=14, face="bold", angle = 90),
                    strip.text.x = element_text(size = 12),
                    strip.text.y = element_text(size = 12))
        }
    }
    
    if(Show.Info == F){
        return(All_Traces+
                   labs(subtitle = "True value is in solid black and posterior mean is in dashed red.",
                        x = "Iteration") )
    } else{
        return(All_Traces +
                   labs(subtitle = "True value is in solid black and posterior mean is in dashed red.",
                        x = "Iteration",
                        caption = paste0("Simulation seed is ", Simulation_Seed, " and sampling seed is ",
                                         Sampling_Seed, 
                                         ". Initialization is set to ", Parameter_Initialization,
                                         "."))) 
    }
}

#######################################################################################################################
## Function - 4: Trace_Summation()
#######################################################################################################################

Trace_Summation <- function(Protein_Numbers = c(1)){
    
    Protein_Index <- Protein_Numbers
    
    plotting_eta <- eta
    
    rownames(plotting_eta) <- 1:Number_of_Proteins
    colnames(plotting_eta) <- 1:Number_of_Samples
    
    true_value_tibble <- driver::gather_array(plotting_eta, True_Value, Protein, Sample) %>%
                            filter(Protein %in% Protein_Numbers) %>% group_by(Sample) %>% 
                                summarise(Summed_Intensity = log(sum(exp(True_Value))))
    
    ## Create Plotting Tibble from Sampling array
    
    ## Create Sample labels
    Sample.labs <- paste0("Sample ", 1:Number_of_Samples)
    names(Sample.labs) <- 1:(Number_of_Samples)
    
    ## Create Protein labels
    # Protein.labs <- paste0("Protein ",1:Number_of_Proteins) 
    Protein.labs <- rownames(eta)
    names(Protein.labs) <- 1:Number_of_Proteins
    
    plotting_tibble <- driver::gather_array(model_fit$samples, Sampled_Log_Intensity, Sample, Protein, Iteration)
    
    # posterior_statistics <- plotting_tibble %>% filter(Iteration > Burn_in_Iterations) %>% 
    #     group_by(Protein, Sample) %>% 
    #     summarise(Mean_Sampled_Log_Intensity = mean(Sampled_Log_Intensity),
    #               Median_Sampled_Log_Intensity = median(Sampled_Log_Intensity))
    
    ## Summation posterior statistics
    summed_posterior <- (plotting_tibble %>% filter(Protein %in% Protein_Index) %>% 
                             filter(Iteration > Burn_in_Iterations) %>% 
                             group_by(Sample, Iteration) %>% 
                             summarise(Sampled_Log_Intensity = log(sum(exp(Sampled_Log_Intensity)))) %>% 
                             ungroup())
    
    summed_posterior <- summed_posterior %>% group_by(Sample) %>%
                         summarise(posterior_mean = mean(Sampled_Log_Intensity)) %>% ungroup()
    ## Sample-Wise traces
    
    Summed_Trace_Plots <- ggplot(data = (plotting_tibble %>% filter(Protein %in% Protein_Index) %>% 
                                                 filter(Iteration > Burn_in_Iterations) %>% 
                                             group_by(Sample, Iteration) %>% 
                    summarise(Sampled_Log_Intensity = log(sum(exp(Sampled_Log_Intensity)))) %>% 
                                             ungroup()),
                                     aes(y = Sampled_Log_Intensity, x = Iteration,
                                         color = as.factor(Sample))) +
        geom_line() +
        geom_hline(data = true_value_tibble, aes(yintercept = Summed_Intensity), size = 1) +
        geom_hline(data = summed_posterior, aes(yintercept = posterior_mean, group = Sample),
                   size = 1,linetype = 2, color = "red") +
        facet_wrap(Sample~., scales = "free", 
                   labeller = labeller(Sample = Sample.labs), nrow = Number_of_Samples) +
        viridis::scale_color_viridis(name = 'Sample', na.value="grey", option = "viridis", discrete=TRUE) +
        theme_linedraw() +
        labs(y = "Summed sampled log-intensity")+
        # labs(y = "Summed sampled log-intensity",
        # title = paste0("Traceplots and histograms of summed sampled abundance values of all proteins from chain 1"),
        #      subtitle = "True total abundance is in solid black and mean of summed posterior in dashed red.") +
        theme(legend.position = "none",
              axis.title.y = element_text(size=14),
              axis.title.x = element_text(size=14),
              axis.text.y = element_text(size=12),
              axis.text.x = element_text(size=12, angle = 0),
              strip.text = element_text(size = 14, color = "black"),
              strip.background = element_rect(fill = "white", color = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
    
    ## Sample-Wise histograms
    
    Summed_Histograms <- ggplot(data = (plotting_tibble %>% filter(Protein %in% Protein_Index) %>% 
                                                filter(Iteration > Burn_in_Iterations) %>% 
                                                group_by(Sample, Iteration) %>% 
                                        summarise(Sampled_Log_Intensity = log(sum(exp(Sampled_Log_Intensity)))) %>% 
                                                ungroup()),
                                    aes(x = Sampled_Log_Intensity,
                                        color = as.factor(Sample),
                                        fill = as.factor(Sample))) +
        geom_histogram(binwidth = 0.01) +
        geom_vline(data = true_value_tibble, aes(xintercept = Summed_Intensity), size = 1) +
        geom_vline(data = summed_posterior, aes(xintercept = posterior_mean, group = Sample),
                   size = 1,linetype = 2, color = "red") +
        facet_wrap(.~Sample, scales = "free", nrow = Number_of_Samples,
                   labeller = labeller(Sample = Sample.labs)) +
        viridis::scale_color_viridis(name = 'Count', na.value="grey", option = "viridis", discrete=TRUE) +
        viridis::scale_fill_viridis(name = 'Sample', na.value="grey", option = "viridis", discrete=TRUE) +
        theme_linedraw() +
        labs(x = "Summed log-intensity", y = "Count") +
        theme(legend.position = "none",
              axis.title.y = element_text(size=14),
              axis.title.x = element_text(size=14),
              axis.text.y = element_text(size=12),
              axis.text.x = element_text(size=12, angle = 0),
              strip.text = element_text(size = 14, color = "black"),
              strip.background = element_rect(fill = "white", color = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
    
    # Return GGPlot 
    return(Summed_Trace_Plots + Summed_Histograms + plot_layout(widths = c(5, 1)))
}


#######################################################################################################################
## Function - 5: Two_Protein_Scatterplot()
#######################################################################################################################

Two_Protein_Scatterplot <- function(Protein.1 = 1, Protein.2 = 2){
    
    Two_Protein_Indices <- c(Protein.1, Protein.2)
    
    plotting_eta <- eta
    
    rownames(plotting_eta) <- 1:Number_of_Proteins
    colnames(plotting_eta) <- 1:Number_of_Samples
    
    true_value_tibble <- driver::gather_array(plotting_eta, True_Value, Protein, Sample)
    
    ## Create Sample labels
    Sample.labs <- paste0("Sample ", 1:Number_of_Samples)
    names(Sample.labs) <- 1:(Number_of_Samples)
    
    ## Create Protein labels
    Protein.labs <- rownames(eta)
    names(Protein.labs) <- 1:Number_of_Proteins
    
    plotting_tibble <- driver::gather_array(model_fit$samples, Sampled_Log_Intensity, Sample, Protein, Iteration)
    
    posterior_statistics <- plotting_tibble %>% filter(Iteration > Burn_in_Iterations) %>% 
        group_by(Protein, Sample) %>% 
        summarise(Mean_Sampled_Log_Intensity = mean(Sampled_Log_Intensity),
                  Median_Sampled_Log_Intensity = median(Sampled_Log_Intensity))
    
    scatterplot_tibble <- plotting_tibble %>% filter(Protein %in% Two_Protein_Indices) %>% 
        filter(Iteration > Burn_in_Iterations) %>% mutate(Sampled_Log_Intensity = exp(Sampled_Log_Intensity)) %>% 
        mutate(Iteration = Iteration - Burn_in_Iterations) %>% 
        pivot_wider(names_from = Protein, values_from = Sampled_Log_Intensity, names_prefix = "Protein_")
    
    Two_Protein_Scatterplot <- ggplot(data = scatterplot_tibble, 
                                      aes(x = get(paste0("Protein_", Two_Protein_Indices[1])),
                                          y = get(paste0("Protein_", Two_Protein_Indices[2])),
                                          color = as.factor(Sample))) + 
        geom_point() +
        facet_wrap(~Sample, nrow = 2, scales = "free",
                   labeller = labeller(Sample = Sample.labs)) +
        labs(x = Protein.labs[Protein.1], y = Protein.labs[Protein.2]) +
        # labs(x = Protein.labs[Protein.1], y = Protein.labs[Protein.2],
        #      title = paste0("Scatterplot of sampled intensities on the real scale of ",
        #                     Protein.labs[Protein.1] ," against ", Protein.labs[Protein.2])) +
        viridis::scale_color_viridis(name = 'Sample', na.value="grey", option = "viridis", discrete=TRUE) +
        geom_smooth(method='lm', formula= y~x, color = "black") +
        # stat_cor(label.x = min(scatterplot_tibble[,3]),
        #          label.y = (max(scatterplot_tibble[,4]) - min(scatterplot_tibble[,4]))/2) + 
        stat_regline_equation(label.x =  0,
                              label.y =  0, size = 5, color = "black") +
        theme_linedraw()+
        theme(
            plot.title = element_text(size=20, face = "bold"),
            plot.subtitle = element_text(size=18, face = "bold"),
            axis.title.y = element_text(size=14),
            axis.title.x = element_text(size=14),
            axis.text.y = element_text(size=12),
            axis.text.x = element_text(size=12, angle = 0),
            strip.text = element_text(size = 14, color = "black"),
            strip.background = element_rect(fill = "white", color = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "none")
    
    return(Two_Protein_Scatterplot)
}

