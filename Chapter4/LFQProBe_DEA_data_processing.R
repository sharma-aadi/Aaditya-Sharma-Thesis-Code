#######################################################################################################################
## load required packages
#######################################################################################################################

source(file = "required_packages.R")

#######################################################################################################################
## load helper functions
#######################################################################################################################

source(file = "helper_functions.R")

#######################################################################################################################
## load the precursor and protein data tables produced by MaxQuant and DIA-NN for DDA and DIA data, respectively
#######################################################################################################################

## DDA data tables
peptide_table_dda <- read_maxquant_table(file_name = "./Data\ tables/maxquant_peptide_table.txt")

## DIA data tables
peptide_table_dia <- read_maxquant_table(file_name = "./Data\ tables/diann_evidence_table.tsv") 

#######################################################################################################################
## load the database of tryptic peptides (this adds in_silico_human_peptides and in_silico_spike_in_peptides)
#######################################################################################################################

## tryptic peptide databases
load("./Data/tryptic_peptide_databases.RData") 

#######################################################################################################################
## create spike-in protein information summary table: spike_in_protein_statistics
#######################################################################################################################

source("spike_in_protein_statistics.R")

#######################################################################################################################
## everything up until can also be loaded with the following R workspace
#######################################################################################################################

# load("LFQProBe_DEA.RData")

#######################################################################################################################
## spike-in protein letters vector (A-Z and Median and Slack <-> UniProt correspondence for spike-in proteins)
#######################################################################################################################

## vector of protein letter IDs
protein_letters <- spike_in_protein_statistics$ProteinID

names(protein_letters) <- spike_in_protein_statistics$Protein

#######################################################################################################################
## source script to generate latin-square matrix (creates latin_square_matrix object with protein CNs in each mix)
#######################################################################################################################

source("latin_square_matrix.R")

## Ordering of intensity column names
Latin_square_index <- as.integer(sort(as.character(c(1:13))))

#######################################################################################################################
## Standardize DDA peptide table data wrangling (Sequence, Proteins, Leading.razor.protein, Intesities)
#######################################################################################################################

## P13717 is marked as a contaminant by MaxQuant, change this so it doesn't get removed
peptide_table_dda$Potential.contaminant[which(peptide_table_dda$Leading.razor.protein == "P13717")] <- ""

## Remove all potential contaminants and reverse sequences
peptide_table_dda <- peptide_table_dda %>% filter(Potential.contaminant != "+") %>% filter(Reverse != "+")

## P13717 is marked as a contaminant by MaxQuant, change this so it doesn't get removed
peptide_table_dda$Proteins[which(peptide_table_dda$Leading.razor.protein == "P13717")] <- "P13717"

## Keep only the sequence, proteins, gene names, razor protein, and intensity columns for the 39 spike-in runs
peptide_table_dda <- peptide_table_dda[,c(1,35:36,39,211:223,228:240,245:257)]

## Mark if the peptide sequence is unique to the protein in the sample
peptide_table_dda <- peptide_table_dda %>% mutate(Unique = ifelse(str_detect(Proteins, ";"), "No", "Yes"))

## Fetch experiment labels
experiment_labels_dda <- names(peptide_table_dda)[5:43]

## Ordering of intensity column names for MQ output
Latin_square_index <- as.integer(sort(as.character(c(1:13))))

experiment_id_dda <-  c(paste0("Mix",Latin_square_index,"_1"),
                        paste0("Mix", Latin_square_index,"_2"),
                        paste0("Mix", Latin_square_index,"_3"))

## Rename intensity columns with standardized experiment labels
peptide_table_dda <- peptide_table_dda %>% rename_at(vars(experiment_labels_dda), ~ experiment_id_dda)

## Re-arrange intensity columns
peptide_table_dda <- peptide_table_dda[,c(1:4,5,10:17,6:9,18,23:30,19:22,31,36:43,32:35,44)]

#######################################################################################################################
## Standardize DIA - peptide table data wrangling (Sequence, Proteins, Leading.razor.protein, Intesities)
#######################################################################################################################

## Sequence, Leading.razor.protein, genes, and intensity columns for the 39 spike-in runs
peptide_table_dia <- peptide_table_dia[,c(7,1,4,11:23,29:41,47:59)]

peptide_table_dia <- peptide_table_dia %>% dplyr::rename(Sequence = Stripped.Sequence, 
                                                         Leading.razor.protein = Protein.Group)

experiment_labels_dia <- names(peptide_table_dia)[4:42]

experiment_id_dia <-  c(paste0("Mix", 1:13,"_1"),
                        paste0("Mix", 1:13,"_2"),
                        paste0("Mix", 1:13,"_3"))

peptide_table_dia <- peptide_table_dia %>% rename_at(vars(experiment_labels_dia),~ experiment_id_dia)

peptide_table_dia <- peptide_table_dia %>% pivot_longer(-c(Sequence, Leading.razor.protein, Genes),
                                                                          names_to = "Experiment",
                                                                          values_to = "Intensity")

peptide_table_dia <- peptide_table_dia %>% group_by(Sequence, Leading.razor.protein, Genes, Experiment) %>% 
                        summarise(Intensity = sum(Intensity, na.rm = T)) %>% ungroup()

peptide_table_dia <- peptide_table_dia %>% pivot_wider(names_from = Experiment, values_from = Intensity) 

## Create a other proteins column to indicate a peptide is not unique
in_silico_lfqprobe_peptides <- bind_rows(in_silico_human_peptides, in_silico_spike_in_peptides)
in_silico_lfqprobe_peptides <- in_silico_lfqprobe_peptides %>% 
                                filter(Protein %in% peptide_table_dia$Leading.razor.protein)

in_silico_lfqprobe_peptides <- in_silico_lfqprobe_peptides %>% filter(Sequence %in% peptide_table_dia$Sequence)
in_silico_lfqprobe_peptides <- in_silico_lfqprobe_peptides %>% dplyr::select(c(Sequence, Protein))
in_silico_lfqprobe_peptides <- in_silico_lfqprobe_peptides %>% group_by(Sequence) %>% 
                                summarise(Proteins = paste(Protein, collapse = ";")) %>% ungroup()

lfqprobe_proteins <- in_silico_lfqprobe_peptides$Proteins
names(lfqprobe_proteins) <- in_silico_lfqprobe_peptides$Sequence

## Add the Proteins column to the DIA data table
peptide_table_dia <- peptide_table_dia %>% 
                        mutate(Proteins = lfqprobe_proteins[Sequence])

peptide_table_dia <- peptide_table_dia[,c(1,43,2,3,4:42)]

## Rename the genes column
peptide_table_dia <- peptide_table_dia %>% dplyr::rename(Gene.names = Genes)

## Re-arrange intensity columns
peptide_table_dia <- peptide_table_dia[,c(1:4,5,10:17,6:9,18,23:30,19:22,31,36:43,32:35)]

## Mark if the peptide sequence is unique to the protein in the sample
peptide_table_dia <- peptide_table_dia %>% mutate(Unique = ifelse(str_detect(Proteins, ";"), "No", "Yes"))
peptide_table_dia <- peptide_table_dia %>% filter(Leading.razor.protein != "P02769")

#######################################################################################################################
## Intensity column names
#######################################################################################################################

intensity_columns <- names(peptide_table_dia[,c(5:43)])

#######################################################################################################################
## Replace all zeroes with NAs
#######################################################################################################################

peptide_table_dda <- peptide_table_dda %>% mutate(across(intensity_columns, na_if, 0))
peptide_table_dia <- peptide_table_dia %>% mutate(across(intensity_columns, na_if, 0))

#######################################################################################################################
## Replace spike-in protein names with A-Z and Median
#######################################################################################################################

## Rename spike-in proteins with their letter names
peptide_table_dda$Leading.razor.protein[which(peptide_table_dda$Leading.razor.protein %in% protein_uniprot_id)] <- 
protein_letters[peptide_table_dda$Leading.razor.protein[which(peptide_table_dda$Leading.razor.protein%in%protein_uniprot_id)]]

peptide_table_dia$Leading.razor.protein[which(peptide_table_dia$Leading.razor.protein %in% protein_uniprot_id)] <- 
protein_letters[peptide_table_dia$Leading.razor.protein[which(peptide_table_dia$Leading.razor.protein%in%protein_uniprot_id)]]

#######################################################################################################################
## experimental design 
#######################################################################################################################

source("LFQProBe_experimental_design.R")
