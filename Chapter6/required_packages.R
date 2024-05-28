if(!exists("load_packages")){
    # Data handling and manipulation
    library(tidyverse, quietly = TRUE)       # Tidyverse for handling data frames
    library(reshape2, quietly = TRUE)
    library(driver, quietly = TRUE)          # gather_array function
    
    # In-silico digestion
    library(seqinr, quietly = TRUE)          # For reading in fasta files
    library(stringr, quietly = TRUE)         # str_sub and str_detect function implementations
    
    # Data Simulation
    library(fitdistrplus, quietly = TRUE)    # Distribution Fitting
    library(MASS, quietly = TRUE)            # Maximum-Likelihood based curve fitting
    library(pracma, quietly = TRUE)          # Dot product
    
    # Sampling
    library(rstan, quietly = TRUE)           # RStan
    
    # Plotting
    
    library(bayesplot, quietly = TRUE)
    library(ggpubr, quietly = TRUE)
    library(patchwork, quietly = TRUE)       # compose R plots
    
    load_packages <- TRUE
}
