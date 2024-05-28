#######################################################################################################################
## filter global tryptic peptide database to have specified min and max peptide lengths 
#######################################################################################################################

tryptic_peptide_database_canonical <- readRDS(tryptic_database_file)

tryptic_peptide_database_canonical <- tryptic_peptide_database_canonical %>% 
    filter(Length <= maximum_peptide_length) %>% 
    filter(Length >= minimum_peptide_length)

#######################################################################################################################
## sample proteins for the simulation 
#######################################################################################################################

## sample indices of proteins to be selected for simlation
simulation_protein_indices <- sample(1:proteome_size, size = number_of_proteins)

## extracted sequences and IDs using sampled indices
simulation_protein_IDs <- global_protein_list$ID[simulation_protein_indices]
simulation_protein_sequences <- global_protein_list$Sequence[simulation_protein_indices]

## Sample selected proteins for proteins that will be differentially expressed in the non-control sample
de_protein_IDs <- sample(simulation_protein_IDs, number_of_de_proteins)

#######################################################################################################################
## simulate tryptic digestion to compute peptide proportions
#######################################################################################################################

## Initialise a tibble
simulated_tryptic_digest <- tibble(Sequence = c(), Protein = c())

for(i in 1:number_of_proteins){
    
    ## extract sequence and UniProtID of the current protein
    current_protein_id <- global_protein_list$ID[simulation_protein_indices][i]
    current_protein_sequence <- global_protein_list$Sequence[simulation_protein_indices][i]
    
    ## Simulate number_of_digestion_simulations iterations of digestion of current protein
    peptides <- unlist(lapply(rep(0,number_of_digestion_simulations),
                              function(x) generate_tryptic_peptides_prob(max_missed_cleavages = x,
                                                                         mc_prob =  missed_cleavage_probability,
                                                                         Protein = current_protein_sequence,
                                                                         remove_proline_sites = FALSE,
                                                                         methionine_peptides = FALSE)))
    ## Store digestion results in a table
    current_digest <- tibble(Sequence = peptides, Protein = current_protein_id)
    
    ## Add to the data from previous protein(s)
    simulated_tryptic_digest <- bind_rows(simulated_tryptic_digest, current_digest)
}

## Only keep the peptides in the database

simulated_tryptic_digest <- simulated_tryptic_digest%>%filter(Sequence %in% tryptic_peptide_database_canonical$Sequence)
simulated_tryptic_digest <- simulated_tryptic_digest %>% group_by(Sequence, Protein) %>% summarise(Count = n()) %>% 
    ungroup()

## prepare data for adding missed cleavages annotation
missed_cleavage_data <- integer(length = length(simulated_tryptic_digest$Sequence))

tryptic_peptide_database_canonical_condensed <- tryptic_peptide_database_canonical %>% 
    filter(Protein %in% simulated_tryptic_digest$Protein) %>% 
    filter(Sequence %in% simulated_tryptic_digest$Sequence) %>% 
    group_by(Sequence, Protein) %>% 
    summarise(MC = mean(MissedCleavages)) %>% ungroup()

## match tables
tryptic_peptide_database_canonical_condensed <-tryptic_peptide_database_canonical_condensed%>%arrange(Sequence,Protein)
simulated_tryptic_digest <- simulated_tryptic_digest %>% arrange(Sequence,Protein)

## Add missed cleavage column
simulated_tryptic_digest <- simulated_tryptic_digest %>% 
    mutate(Missed_Cleavages = tryptic_peptide_database_canonical_condensed$MC)

## Scale counts to per protein molecule 
simulated_tryptic_digest <- simulated_tryptic_digest %>% mutate(Count = Count/number_of_digestion_simulations) %>% 
    dplyr::rename(Proportion = Count)

