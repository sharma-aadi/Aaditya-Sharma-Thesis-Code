#######################################################################################################################
## Main Plot
#######################################################################################################################

## Plot - 1: All Proteins Histograms
mcmc_1 <- Visualize_All_Histograms(Show.Info = F, Scale = "Logarithmic")    
# Visualize_All_Histograms(Show.Info = T, Scale = "Real")
mcmc_1
## Plot - 2: All Proteins Traces
Visualize_All_Traces(Show.Info = F, Show.Warmup = T, Scale = "Logarithmic")
# Visualize_All_Traces(Show.Info = T, Show.Warmup = F, Scale = "Real")

## Plot - 3: Sum of all traces (How well can the model estimate total protein amount)
mcmc_2 <-Trace_Summation(Protein_Numbers = c(1:Number_of_Proteins)) 

## Plot - 4: Scatterplot of sampled abundances for two proteins [On real scale by default]
mcmc_3 <- Two_Protein_Scatterplot(Protein.1 = 1, Protein.2 = 4)

#######################################################################################################################
## Arrange plots for thesis
#######################################################################################################################

ggsave(plot = ggarrange(mcmc_1, mcmc_2, nrow = 2, labels = c("A", "B")),
       filename = paste0("mcmc_figure_thesis.pdf"),
       path = paste0("./"),
       width = 11.69, height = 12, units = "in",
       dpi = 320)

ggsave(plot = mcmc_1,
       filename = paste0("mcmc_figure_thesis_1.pdf"),
       path = paste0("./"),
       width = 11.69, height = 6, units = "in",
       dpi = 320)

ggsave(plot = mcmc_2,
       filename = paste0("mcmc_figure_thesis_2.pdf"),
       path = paste0("./"),
       width = 11.69, height = 6, units = "in",
       dpi = 320)

ggsave(plot = mcmc_3,
       filename = paste0("mcmc_figure_thesis_3.pdf"),
       path = paste0("./"),
       width = 11.69, height = 6, units = "in",
       dpi = 320)


#######################################################################################################################
## Other Plots
#######################################################################################################################

## 


## Other plots - 1: Individual protein (Traces + Histograms)
Visualize_Single_Protein(Protein_Number = 2)

## Other plots - 2: Sum of traces of specific proteins
Trace_Summation(Protein_Numbers = c(1,3,5))
Trace_Summation(Protein_Numbers = c(2,4,6))

#######################################################################################################################
## Extras
#######################################################################################################################

gamma_distribution <- function(x,
                               alpha = shape_protein, beta = rate_protein){
    return(((beta^alpha)/(gamma(alpha)))*(x^(alpha-1))*(exp(-beta*x)))
}

curve(gamma_distribution(x), xlim = c(0,30), main = "Prior probability density function", ylab = "", xlab = "")

