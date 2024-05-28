#######################################################################################################################
## set seed
#######################################################################################################################

set.seed(1)

#######################################################################################################################
## source scripts
#######################################################################################################################

source("file_paths.R")          
source("required_packages.R")    ## This file loads all the packages required for running this script
source("helper_functions.R")     ## Load helper functions
source("parameters.R")           ## Load all the simulation parameters
# source("proteome_database.R")  ## Create the proteome database
source("digestion_simulation.R") ## Simulate tryptic digestion

#######################################################################################################################
## operational variables
#######################################################################################################################

# total number of LC-MS/MS runs
total_ms_runs <- number_of_conditions*number_of_technical_replicates*number_of_biological_replicates

#######################################################################################################################
## simulate baseline copy numbers
#######################################################################################################################

# control_copy_numbers <- rnorm(n = , mean = protein_copy_number_mean, sd = protein_copy_number_sd)
control_copy_numbers <- rgamma(n = number_of_proteins, shape = protein_shape, rate = protein_rate) # these are log10
control_copy_numbers <- log2(10^control_copy_numbers) # scale change to log2

#######################################################################################################################
## Create matrix of copy numbers
#######################################################################################################################

## Initialize matrix
copy_number_matrix <- matrix(data = control_copy_numbers,
                             nrow = number_of_proteins,
                             ncol = number_of_conditions)

## Name rows with protein IDs
rownames(copy_number_matrix) <- simulation_protein_IDs

## Create sample names for MS runs
column_names <- paste0("condition_",1:number_of_conditions)

colnames(copy_number_matrix) <- column_names

#######################################################################################################################
## simulate fold-changes in the DE proteins
#######################################################################################################################

for(i in 1:number_of_conditions){
    copy_number_matrix[de_protein_IDs,i] <- copy_number_matrix[de_protein_IDs,i]+log2(fold_change[i])
}

#######################################################################################################################
## simulate extraction efficiencies for all the proteins
#######################################################################################################################

## simulate extraction efficiencies (this is another constant shift)
protein_extraction_efficiencies <- rbeta(number_of_proteins, 10, 2)

protein_extraction_efficiencies <- matrix(log2(protein_extraction_efficiencies), 
                                          nrow = number_of_proteins, ncol = number_of_conditions)

## add extraction errors to proteins
copy_number_matrix_extracted <- copy_number_matrix + protein_extraction_efficiencies

#######################################################################################################################
## simulate biological noise (protein X condition)
#######################################################################################################################

copy_number_matrix_sample <-  matrix(nrow = number_of_proteins,
                                     ncol = number_of_conditions*number_of_biological_replicates)

## populate new matrix
for(i in 1:number_of_conditions){
    copy_number_matrix_sample[,(number_of_biological_replicates*i-(number_of_biological_replicates-1)):
                                  (number_of_biological_replicates*i)] <- copy_number_matrix_extracted[,i]
}

## If there are multiple biological replicates, then create new column names
if(number_of_biological_replicates > 1){
    column_names_sample <- outer(column_names,paste0("_B_",1:number_of_biological_replicates),paste0)
    column_names_sample <- sort(as.character(column_names_sample)) 
} else{
    column_names_sample <- column_names
}

## sample-level matrix
rownames(copy_number_matrix_sample) <- simulation_protein_IDs
colnames(copy_number_matrix_sample) <- column_names_sample

## biological error matrix
biological_error_matrix <- matrix(data = rnorm(n = number_of_proteins*
                                                   number_of_biological_replicates*number_of_conditions,
                                               mean = 0, sd = biological_error), nrow = number_of_proteins)

## add biological error on the log scale
copy_number_matrix_sample_final <- copy_number_matrix_sample + biological_error_matrix

#######################################################################################################################
## scale to sample level to calculate total copy numbers
#######################################################################################################################

copy_number_matrix_sample_final <- copy_number_matrix_sample_final + log2(number_of_cells)

#######################################################################################################################
## Table of tryptic peptide counts in the sample (before LC and ionisation)
#######################################################################################################################

number_of_samples <- dim(copy_number_matrix_sample_final)[2]
peptide_count_matrix <- matrix(nrow = length(simulated_tryptic_digest$Sequence), ncol = number_of_samples)

for(i in 1:number_of_samples){
    for(j in 1:length(simulated_tryptic_digest$Sequence)){
        peptide_count_matrix[j,i] <- 2^(copy_number_matrix_sample_final[simulated_tryptic_digest$Protein[j],i])
    }
    peptide_count_matrix[,i] <- peptide_count_matrix[,i]*simulated_tryptic_digest$Proportion
}

peptide_count_matrix <- log2(peptide_count_matrix)

colnames(peptide_count_matrix) <- colnames(copy_number_matrix_sample_final)
rownames(peptide_count_matrix) <- simulated_tryptic_digest$Sequence

peptide_abundance_data <- bind_cols(simulated_tryptic_digest, peptide_count_matrix)

#######################################################################################################################
## LC efficiency - randomly remove a proportion of peptides as they are undetectable
#######################################################################################################################.

all_peptide_sequences <- unique(peptide_abundance_data$Sequence)
number_of_distinct_peptides <- length(all_peptide_sequences)

detectable_peptide_sequences <- all_peptide_sequences[sample(c(TRUE,FALSE), replace = TRUE,
                                                             number_of_distinct_peptides,
                                                             c(peptide_detectability,1-peptide_detectability))]

peptide_abundance_data <- peptide_abundance_data %>% filter(Sequence %in% detectable_peptide_sequences)

#######################################################################################################################
## Simulate ionisation efficiency
#######################################################################################################################

## simulate ionisation efficiencies of peptides
unique_peptides <- unique(peptide_abundance_data$Sequence)

## Simulate ionisation efficiency
ionisation <- rbeta(length(unique_peptides), 10, 2) # beta distribution for ionisation
# ionisation <- runif(length(unique_peptides), 0, 1) # uniform distribution for ionisation

## Add peptide sequence annotation to ionisation efficiency
names(ionisation) <- unique_peptides

peptide_abundance_data <- peptide_abundance_data %>% mutate(ionisation = ionisation[peptide_abundance_data$Sequence])
peptide_abundance_data <- peptide_abundance_data %>% mutate(condition_1 = log2((2^condition_1)*ionisation)) %>% 
                            mutate(condition_2 = log2((2^condition_2)*ionisation))

#######################################################################################################################
## Length dependent detectability
#######################################################################################################################

## Add a length column
peptide_abundance_data <- peptide_abundance_data %>% mutate(Length = nchar(Sequence))
peptide_abundance_data <- peptide_abundance_data %>% mutate(Normalised_Length = Length/100)

## Simulate detectability with beta PDF

beta_support <- seq(0, 1, length.out = 1000)
# plot(beta_support, dbeta(beta_support,3,16)/range(dbeta(beta_support,3,16))[2], type = "l", xlab = "x",
#      ylab = "f(x)/max(f(x))")

## Add detectability column and filter on it
peptide_abundance_data_detectable <- peptide_abundance_data %>% 
                            mutate(Detectable = dbeta(Normalised_Length,3,16)/range(dbeta(beta_support,3,16))[2])

subtraction_factor <- 1-(range(peptide_abundance_data_detectable$Detectable)[2])

peptide_abundance_data_detectable <- peptide_abundance_data_detectable %>% 
                                        mutate(Detectable = Detectable + subtraction_factor)

## simulate detectability
peptide_detection_vector <- logical(length = length(peptide_abundance_data_detectable$Sequence))

for(i in 1:length(peptide_detection_vector)){
    peptide_detection_vector[i] <- sample(c(TRUE,FALSE),1, replace = TRUE,
                                          prob = c(peptide_abundance_data_detectable$Detectable[i],
                                                   1- peptide_abundance_data_detectable$Detectable[i]))
}

## add annotation for whether or not a peptide is detectable
peptide_abundance_data_detectable <- peptide_abundance_data_detectable %>% 
                                        mutate(Detected = peptide_detection_vector)

peptide_abundance_data_detectable <- peptide_abundance_data_detectable %>% filter(Detected == TRUE)

## remove unnecessary columns
peptide_abundance_data_detectable <- peptide_abundance_data_detectable%>%dplyr::select(-c(Normalised_Length:Detected))

#######################################################################################################################
## Simulate LC runs - add technical noise
#######################################################################################################################

## Extract abundance columns 
peptide_abundance_matrix <- peptide_abundance_data_detectable %>% dplyr::select(-c(Sequence:Missed_Cleavages)) %>% 
                                dplyr::select(-c(ionisation:Length))

peptide_abundance_matrix <- as.matrix(peptide_abundance_matrix)

peptide_abundance_matrix <- matrix(data = c(rep(peptide_abundance_matrix[,1],number_of_technical_replicates),
                                            rep(peptide_abundance_matrix[,2],number_of_technical_replicates)),
                                   nrow = length(peptide_abundance_data_detectable$Sequence),
                                   ncol = 2*number_of_technical_replicates)

rownames(peptide_abundance_matrix) <- peptide_abundance_data_detectable$Sequence

## If there are multiple technical replicates, then create new column names
if(number_of_technical_replicates > 1){
    column_names_technical <- outer(column_names,paste0("_T_",1:number_of_technical_replicates),paste0)
    column_names_technical <- sort(as.character(column_names_technical)) 
} else{
    column_names_technical <- column_names
}

## Add column names to the peptide intensity matrix
colnames(peptide_abundance_matrix) <- column_names_technical

## technical error matrix
technical_error_matrix <- matrix(data = rnorm(n = length(peptide_abundance_data_detectable$Sequence)*
                                                  number_of_conditions*number_of_technical_replicates,
                                               mean = 0, sd = technological_error),
                                 nrow = length(peptide_abundance_data_detectable$Sequence))

## add biological error on the log scale
peptide_abundance_matrix <- peptide_abundance_matrix + technical_error_matrix

#######################################################################################################################
## Add the ion abundance of the same peptide from different proteins
#######################################################################################################################

peptide_abundance_matrix <- tibble(rownames_to_column(as.data.frame(2^peptide_abundance_matrix), "Sequence"))
peptide_abundance_matrix <- peptide_abundance_matrix %>% mutate(Sequence = str_extract(Sequence, "[^.]+")) 
peptide_abundance_matrix <- peptide_abundance_matrix %>% group_by(Sequence) %>% group_by(Sequence) %>% 
                                summarise_all(sum) %>% ungroup()

peptide_abundance_matrix_new <- log2(as.matrix(peptide_abundance_matrix[,-c(1)]))
rownames(peptide_abundance_matrix_new) <- peptide_abundance_matrix$Sequence

#######################################################################################################################
## Generate intensity matrix and simulate intensity-dependent missingness
#######################################################################################################################

## Intialize peptide intensity matrix
Y <- matrix(data = 0, nrow = dim(peptide_abundance_matrix_new)[1], ncol = dim(peptide_abundance_matrix_new)[2])
colnames(Y) <- colnames(peptide_abundance_matrix_new)
rownames(Y) <- rownames(peptide_abundance_matrix_new)

for(i in 1:nrow(Y)){
    for(j in 1:ncol(Y)){
        Y[i,j] <- rlnorm(n = 1, meanlog = log(2^peptide_abundance_matrix_new[i,j]) - ((0.1^2)/2), 0.1)
    }
}

Y <- log2(Y)

missingnes_probability_matrix <- matrix(nrow = dim(Y)[1],
                                        ncol = dim(Y)[2])

for(i in 1:dim(missingnes_probability_matrix)[2]){
    missingnes_probability_matrix[,i] <- pnorm((Y[,i] - mean(Y[,i])) / sd(Y[,i]), 
                                               mean = missingness_cdf_mean,
                                               sd = missingness_cdf_sd)
}

for(i in 1:dim(Y)[1]){
    for(j in 1:dim(Y)[2]){
        if(rbinom(1,1,missingnes_probability_matrix[i,j]) == 0){
            Y[i,j] <- NA
        }
    }
}

#######################################################################################################################
## Construct the mapping matrix
#######################################################################################################################

M <- matrix(data = 0,
            nrow = length(unique(simulated_tryptic_digest$Protein)),
            ncol = length(unique(simulated_tryptic_digest$Sequence)))

rownames(M) <- unique(simulated_tryptic_digest$Protein)
colnames(M) <- unique(simulated_tryptic_digest$Sequence)

for(i in 1:length(simulated_tryptic_digest$Sequence)){
    M[simulated_tryptic_digest$Protein[i],
      simulated_tryptic_digest$Sequence[i]] <- simulated_tryptic_digest$Proportion[i]
}

#######################################################################################################################
## Construct output table
#######################################################################################################################

peptide_table <- peptide_abundance_data_detectable %>% dplyr::select(c(Sequence, Protein, Length, Missed_Cleavages))
peptide_table <- bind_cols(peptide_table, Y)

#######################################################################################################################
## Construct protein table by median summarisation
#######################################################################################################################

## Create median summarised protein table
protein_table <- peptide_table %>% dplyr::select(-c(Length, Missed_Cleavages))
protein_table <- protein_table %>% pivot_longer(-c(Sequence,Protein), names_to = "Experiment", values_to = "Intensity")
protein_table <- protein_table %>% group_by(Protein, Experiment) %>% 
                    summarise(Intensity = median(Intensity, na.rm = T)) %>% 
                        ungroup()

## create annotation map
experiment_names <- column_names_technical
condition <- c(rep(1,3),rep(2,3))
replicate <- c(1:3,1:3)
names(condition) <- experiment_names
names(replicate) <- experiment_names

## add annotations
protein_table <- protein_table %>% mutate(Condition = condition[Experiment]) %>% 
    mutate(Replicate = replicate[Experiment])
    
#######################################################################################################################
## simulation function eventually goes here
#######################################################################################################################

simulate_lfq_data <- function(){
    
}
