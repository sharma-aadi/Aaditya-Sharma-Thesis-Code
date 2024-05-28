#######################################################################################################################
##  Load model
#######################################################################################################################

model_fit <- model_fit_Shared_2_4_v2

#######################################################################################################################
##  Visualize Model Fit - 1 : Create True Value Tibble for plotting
#######################################################################################################################
plotting_eta <- eta

rownames(plotting_eta) <- 1:Number_of_Proteins
colnames(plotting_eta) <- 1:Number_of_Samples

true_value_tibble <- driver::gather_array(plotting_eta, True_Value, Protein, Sample)

#######################################################################################################################
##  Visualize Model Fit - 2.1: Create Plotting Tibble from Sampling array
#######################################################################################################################

## Create Sample labels
Sample.labs <- paste0("Sample ", 1:Number_of_Samples)
names(Sample.labs) <- 1:(Number_of_Samples)

## Create Protein labels
Protein.labs <- paste0("Protein ",1:Number_of_Proteins)
names(Protein.labs) <- 1:Number_of_Proteins

plotting_tibble <- driver::gather_array(model_fit$samples, Sampled_Log_Intensity, Sample, Protein, Iteration)

#######################################################################################################################
##  Visualize Model Fit - 2.2: Create Plotting Tibble for posterior statistics
#######################################################################################################################

posterior_statistics <- plotting_tibble %>% filter(Iteration > Burn_in_Iterations) %>% 
    group_by(Protein, Sample) %>% 
    summarise(Mean_Sampled_Log_Intensity = mean(Sampled_Log_Intensity),
              Median_Sampled_Log_Intensity = median(Sampled_Log_Intensity))

#######################################################################################################################
##  Visualize Model Fit - 1 : Scatterplots
#######################################################################################################################
Two_Protein_Indices <- c(1,2)

scatterplot_tibble <- plotting_tibble %>% filter(Protein %in% Two_Protein_Indices) %>% 
                        filter(Iteration > Burn_in_Iterations) %>% 
                            mutate(Iteration = Iteration - Burn_in_Iterations) %>% 
                                pivot_wider(names_from = Protein, values_from = Sampled_Log_Intensity)


Two_Protein_Scatterplot <- ggplot(data = scatterplot_tibble, aes(x = `1`, y = `2`)) +
                            geom_point() +
                             facet_wrap(~Sample, nrow = 2, scales = "free",
                                        labeller = labeller(Sample = Sample.labs)) +
                                labs(x = "Protein 1", y = "Protein 2",
                                     title = "Scatterplot of sampled log-intensities of protein 1 against protein 2") +
    theme_linedraw()+
    theme( 
        plot.title = element_text(size=20, face = "bold"),
        plot.subtitle = element_text(size=18, face = "bold"),
        axis.title.y = element_text(color="grey20", size=16, face="bold"),
        axis.title.x = element_text(color="grey20", size=16, face="bold"),
        axis.text.y = element_text(color="grey20", size=14, face="bold"),
        axis.text.x = element_text(color="grey20", size=14, face="bold", angle = 90),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12))

Two_Protein_Scatterplot
