#######################################################################################################################
## Notes
#######################################################################################################################

# 1. number_of_cells is the total number of cells obtained
# 2. number_of_proteins is the total number of proteins present in the sample
# 3. up_regulated_proteins and down_regulated_proteins is the number of proteins with positive and negative 
#    fold-changes, respectively. Note that up_regulated_proteins + down_regulated_proteins should be less than
#    number_of_proteins

#######################################################################################################################
## initial parameters
#######################################################################################################################

simulation_seed <- 1

#######################################################################################################################
## sample parameters
#######################################################################################################################

number_of_cells <- 10000
number_of_conditions <- 2            # keep these fixed for now
number_of_technical_replicates <- 3  # keep these fixed for now
number_of_biological_replicates <- 1 # keep these fixed for now

#######################################################################################################################
## protein parameters
#######################################################################################################################


number_of_proteins <- 100     ## Total number of proteins in the samples
number_of_de_proteins <- 20   ## Number of differentially expressed proteins in condition 2

# protein_copy_number_mean <- 4
# protein_copy_number_sd <- 0.75

protein_shape <- 18.24*0.75
protein_rate <-  6*0.75

#######################################################################################################################
## peptide parameters
#######################################################################################################################

peptide_detectability <- 0.75 # proportion of peptides in the tryptic digest that are detectable

#######################################################################################################################
## biological parameters
#######################################################################################################################

biological_error <- 0.2   # sd of the mean-zero normal distribution modelling biological error (0.15)

#######################################################################################################################
## digestion parameters
#######################################################################################################################

number_of_digestion_simulations <- 1000
missed_cleavage_probability <- 0.1       # 0.35 makes 33% of the count data missed cleavage == 1 (0.1)

#######################################################################################################################
## technological variables
#######################################################################################################################

technological_error <- 0.1 # fixed sd of technical variation (0.1)

#######################################################################################################################
## fold_change vector
#######################################################################################################################
    
fold_change <- c(1, 2) # Condition 1 has a fold-change of 1 and then for the other condition, specify fold change    

#######################################################################################################################
## precursor parameters
#######################################################################################################################

precursor_detectability <- 0.8

#######################################################################################################################
## Downstream filtering
#######################################################################################################################

minimum_peptide_length <- 7   # minimum peptide length (MaxQuant default)
maximum_peptide_length <- 40  # maximum peptide length (MaxQuant default)

#######################################################################################################################
## Intensity-dependent missingness parameters
#######################################################################################################################

missingness_cdf_mean <- 1
missingness_cdf_sd <- 0

#######################################################################################################################
## Parameter checks
#######################################################################################################################

cat("Checking parameter validity...\n")

if(number_of_cells < 1){
    cat("Number of cells has to be positive!\n")
}

if(number_of_conditions < 1){
    cat("Number of conditions has to be more than 1!\n")
}

if(length(fold_change) != number_of_conditions){
    cat("Fold-change vector has to be the same length as number of conditions!\n")
}

if(number_of_proteins < 1){
    cat("Number of proteins has to be positive!\n")
}

if(number_of_de_proteins > number_of_proteins){
    cat("Numbe of differentially expressed proteins cannot be more than the total number of proteins\n")
}

# if(up_regulated_proteins > number_of_proteins){
#     cat("Numbe of up-regulated proteins cannot be more than the total number of proteins\n")
# }
# 
# if(down_regulated_proteins > number_of_proteins){
#     cat("Numbe of down-regulated proteins cannot be more than the total number of proteins\n")
# }
# 
# if((up_regulated_proteins + down_regulated_proteins) > number_of_proteins){
#     cat("Numbe of differentially expressed proteins cannot be more than the total number of proteins\n")
# }
