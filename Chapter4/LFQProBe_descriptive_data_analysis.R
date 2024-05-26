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
protein_table_dda <- read_maxquant_table(file_name = "./Data\ tables/maxquant_protein_table.txt")

## DIA data tables
precursor_table_dia <- read_maxquant_table(file_name = "./Data\ tables/diann_evidence_table.tsv")
protein_table_dia <- read_maxquant_table(file_name = "./Data\ tables/diann_protein_table.tsv")

#######################################################################################################################
## load the database of tryptic peptides (this adds in_silico_human_peptides and in_silico_spike_in_peptides)
#######################################################################################################################

## tryptic peptide databases
load("./Data/tryptic_peptide_databases.RData") 

#######################################################################################################################
## load the list of proteins in platelet proteome according to the Jingnan Huang et al. 2021 paper 
#######################################################################################################################

## this is the excel file containing the list of platelet proteins from the platelet_proteome paper (Jignan Huang)
platelet_proteome <- read_excel("./Data/platelet_proteome.xlsx")

## extract the uniprot IDs of the platelet proteins
platelet_proteins_uniprot_id <- unique(platelet_proteome$`Human Protein ID UniProtKB`)

#######################################################################################################################
## create a tibble containing only tryptic peptides from platelet proteins
#######################################################################################################################

## filter the human peptide table down to platelet proteins
in_silico_platelet_peptides <- in_silico_human_peptides %>% filter(Protein %in% platelet_proteins_uniprot_id)

#######################################################################################################################
## create a tibble containing tryptic peptides from platelet proteins and spike-in peptides
#######################################################################################################################

## combine platelet and spike-in peptide databases
in_silico_platelet_and_spike_in_peptides <- bind_rows(in_silico_platelet_peptides %>% mutate(Source = "Platelet"), 
                                                      in_silico_spike_in_peptides %>% mutate(Source = "Spike-in"))

## filter out the MyC-tag from the database as it would be shared between some spike-in proteins
in_silico_platelet_and_spike_in_peptides <- in_silico_platelet_and_spike_in_peptides %>% filter(Sequence != "LISEEDL")

#######################################################################################################################
## spike-in protein name and UniProt ID correspondence as defined in Table 3.3
#######################################################################################################################

# This section is to make a map of correspondence between Protein Uniprot ID <-> Protein Name on Google Sheet

## Create vector of Uniprot IDs
protein_uniprot_id <- c("Q99P86", "Q8K426", "Q99P87", "Q9MZR8", "Q07731", "P48540", "P06804",
                        "P16599", "P09056", "P51642", "P20294", "Q65Z15", "P13717", "P00817",
                        "Q9R1E9", "D3ZET1", "Q8R007", "Q8BTJ4", "Q8R4Z1", "O88766", "P11974",
                        "P11980", "P52480", "P02769", "Q9R0B9", "Q5EGZ1", "Q9R1E6", "P00722")

## Create vector of Corresponding protein names (with organism common name) from Google Sheet 
protein_uniprot_name <- c("Mouse RELM beta", "Mouse RELM gamma", "Mouse Resistin", "Rabbit IL4", "Rat GDNF",
                          "Mouse GDNF","Mouse TNF alpha",
                          "Rat TNF alpha", "Mouse LIF", "Mouse CNTF", "Rat CNTF", "Rat Oncostatin M",
                          "S. marcescens NucA", "Yeast PPA1",
                          "Rat CTGF", "Rat Nectin-4", "Mouse Nectin-4", "Mouse ENPP4", "Rat Serpin A12",
                          "Rat MMP8", "Rabbit Pyruvate Kinase",
                          "Rat Pyruvate Kinase", "Mouse Pyruvate Kinase", "Bovine Serum Albumin", "Mouse Plod2",
                          "Rat ACE2", "Mouse ENPP2", "E. Coli Beta galactosidase")

## Create Maps by assigning names to each of the vectors
names(protein_uniprot_id) <- protein_uniprot_name
names(protein_uniprot_name) <- protein_uniprot_id

#######################################################################################################################
## human_protein_statistics table containing UniProtID, Length, Cleavage Site and Weight for the platelet proteins
#######################################################################################################################

## Create a list object from Spike-in proteins' FASTA file
human_proteome <- read_proteome(Proteome_File = "./FASTA/LFQProbeHuman.fasta")

## Create the Statistic table
human_protein_statistics <- protein_statistics(Protein_List = human_proteome)

## Remove irrelevant columns
human_protein_statistics <- human_protein_statistics %>% dplyr::select(-c(Sites_per_AA)) 

## Rename protein length column
human_protein_statistics <- human_protein_statistics %>% dplyr::rename(Length = AA_Length)

## Subset to only platelet proteins
human_protein_statistics <- human_protein_statistics %>% filter(Protein %in% platelet_proteins_uniprot_id)

#######################################################################################################################
## spike_in_protein_statistics table with UniprotID, PairID, Name, Spike-in type, length, cleavage sites, and weight
#######################################################################################################################

## Create a list object from Spike-in proteins' FASTA file
spike_in_proteome <- read_proteome(Proteome_File = "./FASTA/LFQProbeSpike.fasta")

## Create the Statistic table
spike_in_protein_statistics <- protein_statistics(Protein_List = spike_in_proteome)

## Remove irrelevant columns
spike_in_protein_statistics <- spike_in_protein_statistics %>% dplyr::select(-c(Sites_per_AA)) 

## Add protein_uniprot_name

### Initialize
protein_name_list <- character(length = length(spike_in_protein_statistics$Protein))

### Fill 
for(i in 1:length(protein_name_list)){
    protein_name_list[i] <- protein_uniprot_name[spike_in_protein_statistics$Protein[i]]
}

### Append
spike_in_protein_statistics <- spike_in_protein_statistics %>% mutate(ProteinName = protein_name_list)  

## Add Protein_Type column (Latin square, Median, and Slack)

### Initialize
ProteinType <- character(length = length(spike_in_protein_statistics$Protein))

### Fill 
for(i in 1:length(ProteinType)){
    
    if(spike_in_protein_statistics$ProteinName[i] == "Rat MMP8"){
        ProteinType[i] <- "Median"
    } else if(spike_in_protein_statistics$ProteinName[i] == "Bovine Serum Albumin"){
        ProteinType[i] <- "Slack"
    } else {
        ProteinType[i] <- "Latin_Square"
    }
    
}

### Append
spike_in_protein_statistics <- spike_in_protein_statistics %>% mutate(SpikeInProteinType = ProteinType)  

## Re-arrange columns
spike_in_protein_statistics <- spike_in_protein_statistics %>% dplyr::select(c(Protein, ProteinName,
                                                                               SpikeInProteinType, AA_Length,
                                                                               Cleavage_Sites, Molecular_Weight))

## Arrange rows by Type and Weight
spike_in_protein_statistics <- spike_in_protein_statistics %>% arrange(SpikeInProteinType, Molecular_Weight)

## Add Pair_ID column (The above re-arrangement of rows is a pre-requisite for the following simple annotation to work)
### Create list based on above sorting
PairID <- c(1:13,13:1, "Median", "Slack")

## Append
spike_in_protein_statistics <- spike_in_protein_statistics %>% mutate(PairID)

## Re-order according to Pair_ID
spike_in_protein_statistics <- spike_in_protein_statistics %>%  arrange(
    factor(PairID, levels =  c("1","2","3","4","5","6","7","8","9","10","11","12","13","Median", "Slack")),
    desc(PairID), Protein)

## Add Protein_ID annotation
ProteinID <- c(LETTERS, "Median", "Slack")

spike_in_protein_statistics <- spike_in_protein_statistics %>% mutate(ProteinID) %>% 
    dplyr::select(PairID, ProteinID, Protein:Molecular_Weight) %>% 
    dplyr::rename(Length = AA_Length, 
                  CleavageSites = Cleavage_Sites,
                  MolecularWeight = Molecular_Weight)

#######################################################################################################################
## table of counts of theoretical total and theoretical unique tryptic peptides for each spike-in protein (Table 3.3)
#######################################################################################################################

## Retrieve all spike-in peptide sequences 
spike_in_peptides <- unique(in_silico_spike_in_peptides$Sequence)

## Retrieve all proteins (spike-in and platelet) that map to any spike-in peptide sequences
detectable_proteins <- unique((in_silico_platelet_and_spike_in_peptides %>% 
                                   filter(Sequence %in% spike_in_peptides))$Protein)

## Retrieve all peptides for all detectable proteins based on just the spike-in peptides
detectable_peptides_all <- in_silico_platelet_and_spike_in_peptides %>% filter(Protein %in% detectable_proteins)

## Table of how many proteins each detectable peptide can be generated from
detectable_peptide_mapping_counts <- detectable_peptides_all %>% 
    group_by(Sequence) %>% 
    summarise(NumberOfProteins = n()) %>% ungroup()

## Retrieve set of all detectable peptides which are unique to a single protein (spike-in or human)
unique_detectable_peptides <- (detectable_peptide_mapping_counts %>% filter(NumberOfProteins == 1)) #unique peptides

## Retrieve set of all unique spike-in detectable peptides 
unique_spike_in_peptides <- in_silico_spike_in_peptides %>% filter(Sequence %in% unique_detectable_peptides$Sequence)

## Label peptides as unique or not unique in database of all detectable peptides !(this can be vectorized)
is_peptide_unique <- rep("No", length(detectable_peptides_all$Sequence))

for(i in 1:length(is_peptide_unique)){
    if(is.element(detectable_peptides_all$Sequence[i], unique_detectable_peptides$Sequence)){
        is_peptide_unique[i] <- "Yes"
    }
}

## Add uniqueness annotation
detectable_peptides_all <- detectable_peptides_all %>% mutate(Unique = is_peptide_unique)

## Table of total and unique counts for spike-in proteins
spike_in_protein_peptide_counts <- detectable_peptides_all %>% filter(Source == "Spike-in") %>% group_by(Protein) %>% 
    summarise(`Total peptides` = n(), `Unique peptides` = length(which(Unique == "Yes"))) %>%
    mutate(`Non-unique peptides` = `Total peptides` - `Unique peptides`)

#######################################################################################################################
## create a version of spike-in protein statistics with mapping categories
#######################################################################################################################

spike_in_protein_statistics_with_mapping <- spike_in_protein_statistics %>% filter(ProteinID != "Slack") %>% 
    mutate(Mapping = c(1,1,1,1,1,1,1,1,
                       4,2,4,2,4,2,2,1,
                       1,3,2,2,2,2,1,3,
                       1,1,1))

spike_in_protein_statistics_with_mapping <- spike_in_protein_statistics_with_mapping %>% 
    mutate(Category = rep("Category", 27))

spike_in_protein_statistics_with_mapping$Protein <- factor(spike_in_protein_statistics_with_mapping$Protein,
                                                           levels = spike_in_protein_statistics_with_mapping$Protein)


mapping_labels <- c("no other protein", "another spike-in protein", "a platelet protein",
                    "other spike-in proteins & a platelet protein")

names(mapping_labels) <- c(1,2,3,4)

#######################################################################################################################
## generate a comprehensive mapping data structure for Figure 3.7
#######################################################################################################################

## create a table containing all peptide senqueces for the platelet and spike-in proteins
in_silico_lfqprobe_peptides <- dplyr::bind_rows(in_silico_spike_in_peptides %>% mutate(Source = c("Spike-in")),
                                                in_silico_platelet_peptides %>% mutate(Source = c("Platelet")))

protein <- shared_sequence <- sharing_protein <- sharing_source <-  c()

for(i in 1:length(protein_uniprot_id)){
    
    candidate_sequences <- (in_silico_lfqprobe_peptides %>% filter(Protein == protein_uniprot_id[i]))$Sequence
    candidate_table <- in_silico_lfqprobe_peptides %>% filter(Sequence %in% candidate_sequences) %>% 
        filter(Protein != protein_uniprot_id[i])
    
    number_of_shared_sequences <- dim(candidate_table)[1]
    
    if(number_of_shared_sequences != 0){
        protein <- c(protein, rep(protein_uniprot_id[i], number_of_shared_sequences))
        shared_sequence <- c(shared_sequence, candidate_table$Sequence)
        sharing_protein <- c(sharing_protein, candidate_table$Protein)
        sharing_source <- c(sharing_source, candidate_table$Source)
    }
}

sharing_data <- tibble(protein,shared_sequence,sharing_protein,sharing_source)
sharing_data <- sharing_data %>% filter(shared_sequence != "LISEEDL") # remove Myc-tag 

## Create adjacency matrix for peptide sharing

sharing_data_summary <- sharing_data %>% group_by(protein, sharing_protein, sharing_source) %>% summarise(Count = n())
sharing_data_summary <- sharing_data_summary %>% filter(protein != "P02769")

spike_in_proteins <- (spike_in_protein_statistics %>% filter(Protein != "P02769"))$Protein
platelet_proteins <- setdiff(unique(sharing_data_summary$sharing_protein), spike_in_proteins)

all_sharing_proteins <- c(spike_in_proteins, platelet_proteins)

sharing_graph_matrix <- matrix(data = 0, nrow = length(all_sharing_proteins), ncol = length(all_sharing_proteins))
rownames(sharing_graph_matrix) <- colnames(sharing_graph_matrix) <- all_sharing_proteins

for(i in 1:length(sharing_data_summary$protein)){
    sharing_graph_matrix[sharing_data_summary$protein[i],
                         sharing_data_summary$sharing_protein[i]] <- sharing_data_summary$Count[i]
    
    sharing_graph_matrix[sharing_data_summary$sharing_protein[i],
                         sharing_data_summary$protein[i]] <- sharing_data_summary$Count[i]
}

#######################################################################################################################
## source script to generate latin-square matrix (creates latin_square_matrix object with protein CNs in each mix)
#######################################################################################################################

source("latin_square_matrix.R")

#######################################################################################################################
## use latin_square_matrix to create latin_square_long for tileplot (Figure 3.9)
#######################################################################################################################

## Latin-square tibble
latin_square <- (latin_square_matrix)
rownames(latin_square) <- paste0("Sample_",1:13)

latin_square <- data.frame(latin_square) %>% rownames_to_column("Sample") %>% mutate(SampleID = 1:13)

## Long Latin-square tibble
latin_square_long <- latin_square %>% pivot_longer(-c(Sample,SampleID), names_to = "Protein", values_to = "CopyNumber")

latin_square_long$Protein <- factor(latin_square_long$Protein,
                                    levels=rev(unique((spike_in_protein_statistics%>%filter(Protein!="P02769"))$Protein)))

latin_square_long <- latin_square_long %>% mutate(MixID = latin_square_long$SampleID)

#######################################################################################################################
## spike-in protein letters vector (A-Z and Median and Slack <-> UniProt correspondence for spike-in proteins)
#######################################################################################################################

## vector of protein letter IDs
protein_letters <- spike_in_protein_statistics$ProteinID

names(protein_letters) <- spike_in_protein_statistics$Protein

#######################################################################################################################
## unique spike-in peptide sequence (out of platelet+spike-in) <-> protein correspondence vector for spike-in proteins)
#######################################################################################################################

unique_peptide_proteins <- unique_spike_in_peptides$Protein
names(unique_peptide_proteins) <- unique_spike_in_peptides$Sequence

#######################################################################################################################
## create peptide table for data analysis
#######################################################################################################################

## DDA protein tables
peptide_table_dda_analysis <- peptide_table_dda

## DIA protein tables
peptide_table_dia_analysis <- precursor_table_dia # this will be aggregated into a peptide table shortly 

#######################################################################################################################
## create a labelling system for experimens
#######################################################################################################################

experiment_id <-  c(paste0("Platelet_1_Mix_", 1:13),
                    paste0("Platelet_2_Mix_", 1:13),
                    paste0("Platelet_3_Mix_", 1:13))

platelet_id <- c(rep(1,13),
                 rep(2,13),
                 rep(3,13))

mix_id <- c(rep(1:13,3))

names(platelet_id) <- experiment_id
names(mix_id) <- experiment_id

#######################################################################################################################
## DDA data wrangling for the precursor-level data table
#######################################################################################################################

## Remove reverse peptide hits
peptide_table_dda_analysis <- peptide_table_dda_analysis %>% filter(Reverse != "+")

peptide_table_dda_analysis <- peptide_table_dda_analysis[,c(1, 35, 211:223, 228:240, 245:257)] 

## Add letter IDs and rename Proteins column
peptide_table_dda_analysis <- peptide_table_dda_analysis %>% dplyr::rename(Protein = Proteins)

## Fetch experiment labels
experiment_labels_dda <- names(peptide_table_dda_analysis)[3:41]

## Ordering of intensity column names for MQ output
Latin_square_index <- as.integer(sort(as.character(c(1:13))))

experiment_id_dda <-  c(paste0("Platelet_1_Mix_", Latin_square_index),
                        paste0("Platelet_2_Mix_", Latin_square_index),
                        paste0("Platelet_3_Mix_", Latin_square_index))

## Rename intensity columns with standardized experiment labels
peptide_table_dda_analysis <- peptide_table_dda_analysis %>% rename_at(vars(experiment_labels_dda), ~ experiment_id_dda)

## Subset down to only spike-in peptides
peptide_table_dda_analysis <- peptide_table_dda_analysis %>% filter(Sequence %in% in_silico_spike_in_peptides$Sequence)

peptide_table_dda_analysis <- peptide_table_dda_analysis %>% pivot_longer(-c(Sequence, Protein), names_to = "Experiment",
                                              values_to = "Intensity")

peptide_table_dda_analysis <- peptide_table_dda_analysis %>% 
                                mutate(Unique = ifelse(sapply(1:length(peptide_table_dda_analysis$Sequence),
                                    function (x){return(is.element(peptide_table_dda_analysis$Sequence[x],
                                            unique_spike_in_peptides$Sequence))}),
                                                        "yes", "no"))

peptide_table_dda_analysis <- peptide_table_dda_analysis %>% group_by(Sequence, Protein, Experiment, Unique) %>% 
                                summarise(Intensity = sum(Intensity, na.rm = T)) %>% ungroup()


peptide_table_dda_analysis <- peptide_table_dda_analysis %>% 
                                mutate(Platelet = platelet_id[peptide_table_dda_analysis$Experiment],
                                        Mix = mix_id[peptide_table_dda_analysis$Experiment])

for(i in 1:length(peptide_table_dda_analysis$Sequence)){
    if(peptide_table_dda_analysis$Unique[i] == "yes"){
        peptide_table_dda_analysis$Protein[i] <- unique_peptide_proteins[peptide_table_dda_analysis$Sequence[i]]
    }
}

peptide_table_dda_analysis <- peptide_table_dda_analysis %>% filter(Protein != "P02769")

#######################################################################################################################
## Step - 3.1: DIA data wrangling
#######################################################################################################################

peptide_table_dia_analysis <- peptide_table_dia_analysis[,c(7,1,11:23,29:41,47:59)]
peptide_table_dia_analysis <- peptide_table_dia_analysis %>% 
                                dplyr::rename(Sequence = Stripped.Sequence, Proteins = Protein.Group)

experiment_labels_dia <- names(peptide_table_dia_analysis)[3:41]

experiment_id_dia <-  c(paste0("Platelet_1_Mix_", 1:13),
                        paste0("Platelet_2_Mix_", 1:13),
                        paste0("Platelet_3_Mix_", 1:13))

peptide_table_dia_analysis <- peptide_table_dia_analysis %>% rename_at(vars(experiment_labels_dia),~ experiment_id_dia)
peptide_table_dia_analysis <- peptide_table_dia_analysis %>% filter(Sequence %in% in_silico_spike_in_peptides$Sequence)
peptide_table_dia_analysis <- peptide_table_dia_analysis %>% pivot_longer(-c(Sequence, Proteins),
                                                                          names_to = "Experiment",
                                              values_to = "Intensity")

peptide_table_dia_analysis <- peptide_table_dia_analysis %>% 
                                mutate(Unique = ifelse(sapply(1:length(peptide_table_dia_analysis$Sequence),
                                                function (x){return(is.element(peptide_table_dia_analysis$Sequence[x],
                                                                                unique_spike_in_peptides$Sequence))}),
                                                        "yes", "no"))

peptide_table_dia_analysis <- peptide_table_dia_analysis %>% group_by(Sequence, Proteins, Experiment, Unique) %>% 
    summarise(Intensity = sum(Intensity, na.rm = T))

peptide_table_dia_analysis <- peptide_table_dia_analysis %>% ungroup()
peptide_table_dia_analysis <- peptide_table_dia_analysis %>% 
                                mutate(Platelet = platelet_id[peptide_table_dia_analysis$Experiment],
                                        Mix = mix_id[peptide_table_dia_analysis$Experiment])

for(i in 1:length(peptide_table_dia_analysis$Sequence)){
    if(peptide_table_dia_analysis$Unique[i] == "yes"){
        peptide_table_dia_analysis$Proteins[i] <- unique_peptide_proteins[peptide_table_dia_analysis$Sequence[i]]
    }
}

peptide_table_dia_analysis <- peptide_table_dia_analysis %>% filter(Proteins != "P02769")

#######################################################################################################################
## figure 4.4 venn diagram of peptide sequences
#######################################################################################################################

dda_peptide_sequences_quantified <- unique((peptide_table_dda_analysis %>% filter(Intensity != 0))$Sequence)
dia_peptide_sequences_quantified <- unique((peptide_table_dia_analysis %>% filter(Intensity != 0))$Sequence)

sequence_list <- list(DDA = dda_peptide_sequences_quantified, DIA = dia_peptide_sequences_quantified)

sequence_venn_diagram <- ggVennDiagram(sequence_list) + 
    scale_fill_gradient(low="blue",high = "red")+
    theme(legend.position = "none",
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))

## Individual panels
ggsave(plot = sequence_venn_diagram,
       filename = paste0("figure_4_4.pdf"),
       path = paste0("./Figures/"),
       width = 7.99, height = 5.68, units = "in",
       dpi = 320)

#######################################################################################################################
## figure 4.5 umber of samples peptides are quantified in for each dataset and coefficient of variation for spike-in
#######################################################################################################################

## create summary objects for the spike-in peptide quantification in DDA and DIA
dda_peptide_summary <- peptide_table_dda_analysis %>% filter(Intensity != 0) %>%  group_by(Sequence) %>% 
                            summarise(Count = n()) %>% ungroup()

dia_peptide_summary <- peptide_table_dia_analysis %>% filter(Intensity != 0) %>%  group_by(Sequence) %>% 
                            summarise(Count = n()) %>% ungroup()

peptide_summary_all <- bind_rows(dda_peptide_summary %>% mutate(Dataset = "DDA"),
                                 dia_peptide_summary %>% mutate(Dataset = "DIA"))

summary_medians <- peptide_summary_all %>% group_by(Dataset) %>% 
                    summarise(Total = sum(Count), Median = median(Count)) %>% ungroup()

peptide_count_plot <- ggplot(data = peptide_summary_all)+
    geom_bar(fill = "white",color="black", aes(x = Count))+
    facet_wrap(.~Dataset)+
    geom_vline(data = summary_medians, aes(xintercept = Median), color = "red", lty = "dashed")+
    scale_x_continuous(expand = c(0,0), breaks = seq(0,39,13), lim = c(1,40))+
    scale_y_continuous(expand = c(0,0), breaks = seq(0,50,5), lim = c(0,55))+
    theme_linedraw()+
    labs(x = "Number of LC-M/MS runs in which spike-in peptide sequence is quantified",
         y = "Frequency")+
    theme(legend.position = "bottom",
          legend.justification = "right",
          legend.margin = margin(t=-20,b=30),
          plot.margin = margin(5,5,5,5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "white"),
          strip.text = element_text(color = "black"))

peptides_dda_cv <- peptide_table_dda_analysis %>% filter(Intensity != 0)  %>%
                    group_by(Sequence, Mix) %>% summarise(CV = sd(Intensity)/mean(Intensity), Count = n()) %>% 
                        filter(Count == 3) %>% ungroup()

peptides_dia_cv <- peptide_table_dia_analysis %>% filter(Intensity != 0) %>%
                    group_by(Sequence, Mix) %>% summarise(CV = sd(Intensity)/mean(Intensity), Count = n()) %>% 
                        filter(Count == 3) %>% ungroup()

peptides_cv <- bind_rows(peptides_dda_cv %>% mutate(Dataset = "DDA"),
                         peptides_dia_cv %>% mutate(Dataset = "DIA"))

cv_plot <- ggplot(data = peptides_cv, aes(x = Dataset, y = CV))+
    geom_violin()+ geom_boxplot(width=0.1)+ theme_linedraw()+
    labs(y = "Coefficient of variance")+
    theme(legend.position = "bottom",
          legend.justification = "right",
          legend.margin = margin(t=-20,b=30),
          plot.margin = margin(5,5,5,5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "white"),
          strip.text = element_text(color = "black"))

## Individual panels
ggsave(plot = ggarrange(peptide_count_plot, cv_plot, nrow = 1, widths = c(2,1), labels = c("A", "B")),
       filename = paste0("figure_4_5.pdf"),
       path = paste0("./Figures/"),
       width = 11.69, height = 3, units = "in",
       dpi = 320)

#######################################################################################################################
## figure 4.6 regression for re-scaling spike-in peptide intensities for DDA data
#######################################################################################################################

## Keep only peptides which are unique to a spike-in protein
regression_data_dda <- peptide_table_dda_analysis %>% filter(Unique == "yes")
regression_data_dda <- regression_data_dda %>% mutate(Protein = unique_peptide_proteins[Sequence])

#[!Note 436  distinct peptide sequence remain after this this round of filtering] (i.e. there are 436 unique peptides)

## Add letter IDs and rename Proteins column
regression_data_dda <- regression_data_dda %>% mutate(ProteinID = protein_letters[Protein])

## Label missing data with NAs
regression_data_dda <- regression_data_dda %>% mutate(Intensity = na_if(Intensity, 0))

## log2-transform the intensity column
regression_data_dda <- regression_data_dda %>% mutate(Intensity = log2(Intensity))

## Add concentrations
concentration <- numeric(length = length(regression_data_dda$Protein))

for(i in 1:length(concentration)){
    
    concentration[i] <- latin_square_matrix[regression_data_dda$Mix[i], regression_data_dda$Protein[i]]
    
}

regression_data_dda <- regression_data_dda %>% mutate(Concentration = log2(concentration))

## This is the regresstion table for DDA
regression_data_dda <- regression_data_dda %>% filter(!is.na(Intensity))

## Find regression coefficients for each of the experiments
regression_results_dda <- regression_data_dda %>% group_by(Experiment) %>% 
    summarise(Intercept = lm(Intensity~Concentration)$coefficients[1],
              Slope = lm(Intensity~Concentration)$coefficients[2]) %>% ungroup()

## DDA regression results
Intercept_dda <- regression_results_dda$Intercept
Slope_dda <- regression_results_dda$Slope

names(Intercept_dda) <- names(Slope_dda) <- regression_results_dda$Experiment

#######################################################################################################################
## figure 4.6 protein-wise tibble creation
#######################################################################################################################

protein_list<-c("J", "L", "S", "V", "N", "O", "T", "U", "R", "X", "I", "K", "M") #spike-in proteins that share peptides

all_peptide_data_sharing_proteins_dda <- tibble()

for(i in 1:length(protein_list)){
    
    current_protein <- spike_in_protein_statistics$Protein[which(spike_in_protein_statistics$ProteinID == 
                                                                     protein_list[i])]
    
    current_peptides <- peptide_table_dda_analysis %>% filter(str_detect(Protein, pattern = current_protein))
    
    current_peptides <- current_peptides %>% mutate(Intensity = na_if(Intensity, 0))
    
    ## Add concentrations
    concentration <- numeric(length = length(current_peptides$Protein))
    
    for(j in 1:length(concentration)){
        
        concentration[j] <- latin_square_matrix[current_peptides$Mix[j], current_protein]
        
    }
    
    current_peptides <- current_peptides %>% mutate(Concentration = concentration)
    
    current_peptides <- current_peptides %>% mutate(Concentration = log2(Concentration), Intensity = log2(Intensity))
    
    current_peptides <- current_peptides %>% mutate(Concentration = round(Concentration, digits = 2),
                                                    Intensity = round(Intensity, digits = 2))
    
    current_peptides <- current_peptides %>% mutate(Unique = ifelse(Unique == "yes", "Unique", "Shared"))
    
    current_peptides$Unique <- factor(current_peptides$Unique,
                                      levels = c("Unique", "Shared"))
    
    current_peptides <- current_peptides %>% filter(!is.na(Intensity))
    
    current_peptides <- current_peptides %>% 
        mutate(ScaledIntensity = (Intensity - Intercept_dda[Experiment])/Slope_dda[Experiment])
    
    current_peptides <- current_peptides %>% mutate(Protein = current_protein, ProteinID = protein_list[i])
    
    current_peptides <- current_peptides %>% dplyr::select(c(Sequence, Protein, ProteinID, Unique, Experiment,
                                                             Platelet, Mix, 
                                                             Concentration, Intensity, ScaledIntensity))
    
    all_peptide_data_sharing_proteins_dda <- bind_rows(all_peptide_data_sharing_proteins_dda, current_peptides)
}

#######################################################################################################################
## figure 4.6 regression for re-scaling spike-in peptide intensities for DIA data
#######################################################################################################################

## Keep only peptides which are unique to a spike-in protein
regression_data_dia <- peptide_table_dia_analysis %>% filter(Unique == "yes")
regression_data_dia <- regression_data_dia %>% mutate(Proteins = unique_peptide_proteins[Sequence])

#[!Note 390 distinct peptide sequence remain after this this round of filtering] (i.e. there are 390 unique peptides)

## Add letter IDs and rename Proteins column
regression_data_dia <- regression_data_dia %>% mutate(ProteinID = protein_letters[Proteins]) %>% 
                         dplyr::rename(Protein = Proteins)

## Label missing data with NAs
regression_data_dia <- regression_data_dia %>% mutate(Intensity = na_if(Intensity, 0))

## log2-transform the intensity column
regression_data_dia <- regression_data_dia %>% mutate(Intensity = log2(Intensity))

## Add concentrations
concentration <- numeric(length = length(regression_data_dia$Protein))

for(i in 1:length(concentration)){
    
    concentration[i] <- latin_square_matrix[regression_data_dia$Mix[i], regression_data_dia$Protein[i]]
    
}

regression_data_dia <- regression_data_dia %>% mutate(Concentration = log2(concentration))

## This is the regresstion table for DDA
regression_data_dia <- regression_data_dia %>% filter(!is.na(Intensity))

## Find regression coefficients for each of the experiments
regression_results_dia <- regression_data_dia %>% group_by(Experiment) %>% 
    summarise(Intercept = lm(Intensity~Concentration)$coefficients[1],
              Slope = lm(Intensity~Concentration)$coefficients[2]) %>% ungroup()

## DDA regression results
Intercept_dia <- regression_results_dia$Intercept
Slope_dia <- regression_results_dia$Slope

names(Intercept_dia) <- names(Slope_dia) <- regression_results_dia$Experiment

#######################################################################################################################
## figure 4.6 DIA protein-wise tibble creation
#######################################################################################################################

protein_list <- c("J", "L", "S", "V", "N", "O", "T", "U", "R", "X", "I", "K", "M")

all_peptide_data_sharing_proteins_dia <- tibble()

for(i in 1:length(protein_list)){
    
    current_protein <- spike_in_protein_statistics$Protein[which(spike_in_protein_statistics$ProteinID == 
                                                                     protein_list[i])]
    
    current_peptides <- peptide_table_dia_analysis %>% filter(str_detect(Proteins, pattern = current_protein))
    
    current_peptides <- current_peptides%>% dplyr::rename(Protein = Proteins)
    
    current_peptides <- current_peptides %>% mutate(Intensity = na_if(Intensity, 0))
    
    ## Add concentrations
    concentration <- numeric(length = length(current_peptides$Protein))
    
    for(j in 1:length(concentration)){
        
        concentration[j] <- latin_square_matrix[current_peptides$Mix[j], current_protein]
        
    }
    
    current_peptides <- current_peptides %>% mutate(Concentration = concentration)
    
    current_peptides <- current_peptides %>% mutate(Concentration = log2(Concentration), Intensity = log2(Intensity))
    
    current_peptides <- current_peptides %>% mutate(Concentration = round(Concentration, digits = 2),
                                                    Intensity = round(Intensity, digits = 2))
    
    current_peptides <- current_peptides %>% mutate(Unique = ifelse(Unique == "yes", "Unique", "Shared"))
    
    current_peptides$Unique <- factor(current_peptides$Unique,
                                      levels = c("Unique", "Shared"))
    
    current_peptides <- current_peptides %>% filter(!is.na(Intensity))
    
    current_peptides <- current_peptides %>% 
        mutate(ScaledIntensity = (Intensity - Intercept_dia[Experiment])/Slope_dia[Experiment])
    
    current_peptides <- current_peptides %>% mutate(Protein = current_protein, ProteinID = protein_list[i])
    
    current_peptides <- current_peptides %>% dplyr::select(c(Sequence, Protein, ProteinID, Unique, Experiment,
                                                             Platelet, Mix, 
                                                             Concentration, Intensity, ScaledIntensity))
    
    all_peptide_data_sharing_proteins_dia <- bind_rows(all_peptide_data_sharing_proteins_dia, current_peptides)
}

#######################################################################################################################
## figure 4.6 - combine DDA and DIA tables protein wise tables 
#######################################################################################################################

all_peptide_data_sharing_proteins <- bind_rows(all_peptide_data_sharing_proteins_dda %>% mutate(Dataset = "DDA"),
                                               all_peptide_data_sharing_proteins_dia %>% mutate(Dataset = "DIA"))

#######################################################################################################################
## Figure 4.6
#######################################################################################################################

## Create data table for plotting
plot_4_6_data <- all_peptide_data_sharing_proteins %>% filter(!is.na(ScaledIntensity)) %>%  
                    mutate(Concentration_label = "Target concentration")

plot_4_6_data$Concentration_label <- factor(plot_4_6_data$Concentration_label,
                                          labels = c("log[2](expected~protein~concentration)"))

plot_4_6_data <- plot_4_6_data %>% 
            mutate(Unique = ifelse(Unique == "Unique", "Unique peptides", "Shared peptides")) %>% 
                mutate(Concentration = round(Concentration, digits = 1))

## Concentration labels
concentration_labels_1 <- as.character(round(as.numeric((levels(factor(plot_4_6_data$Concentration)))[1:13]),
                                             digits = 1))
concentration_labels_1[c(8,13)] <- c("12.0", "16.0")
names(concentration_labels_1) <- as.character(levels(factor(plot_4_6_data$Concentration)))

## Font sizes
title_size <- 20
subtitle_size <- 18
axis_title_size <- 12 - 2
axis_text_size <- 8.5
strip_text_size <- 5
legend_title_size <- 10
legend_text_size <- 8

## Plot 1
plot_4_6 <- ggplot(data = plot_4_6_data,
                 aes(x = Unique, y = ScaledIntensity))+
    scale_color_manual(values = c("#0072B2", "#D55E00"))+
    geom_hline(aes(yintercept = Concentration, linetype = ""), linewidth = 0.5)+
    geom_boxplot(aes(color = Unique))+
    facet_grid(Dataset~Concentration,
               labeller = labeller(Concentration_label = label_parsed, Concentration = concentration_labels_1))+
    labs(x = "", y = bquote(Scaled~log[2]~peptide~intensity))+
    guides(color = guide_legend(title = "Observations with:", order = 1))+
    guides(linetype = guide_legend(title = expression(log[2]~"spiked-in parent protein concentration")))+
    scale_y_continuous(breaks = seq(from = 0, to = 25, by = 5))+
    theme_linedraw()+
    theme(panel.grid.minor.x=element_blank(),
          panel.grid.minor.y=element_blank(),
          panel.grid.major.x=element_blank(),
          panel.grid.major.y=element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_text(size = axis_text_size+1, margin = margin(r = 5)),
          axis.ticks.x = element_blank(),
          legend.position = "top",
          legend.justification = "right",
          panel.spacing = unit(0, "pt"),
          legend.margin=margin(b=-8, t = 3),
          strip.background = element_rect(color = "black", fill = "white"),
          strip.text = element_text(color = "black"),
          plot.margin =  margin(t = 5, b = -5, r = 5, l = 5, unit = "pt"))

ggsave(plot = ggarrange(plot_4_6, nrow = 1),
       filename = paste0("figure_4_6.pdf"),
       path = paste0("./Figures/"),
       width = 11.69, height = 3, units = "in",
       dpi = 320)

#######################################################################################################################
## figure 4.8: create copies of the protein tables 
#######################################################################################################################

## DDA protein tables
proteins_dda <- protein_table_dda

## DIA protein tables
proteins_dia <- protein_table_dia

#######################################################################################################################
## figure 4.8: create dia protein table with just latin-square samples 
#######################################################################################################################

## Remove unnecessary columns
proteins_dia <- proteins_dia[,-c(2,3,4,5,23,41)]

## remove platelet-only columns
proteins_dia_spike_in <- proteins_dia[,-c(15:18,32:35,49:52)]

#######################################################################################################################
## figure 4.8: create dda protein table with just latin-square samples
#######################################################################################################################

## Remove unnecessary columns
proteins_dda <- proteins_dda[, c(1, 381:431)]

## remove platelet-only columns
proteins_dda_spike_in <- proteins_dda[, c(1, 2:14, 19:31, 36:48)]

#######################################################################################################################
## figure 4.8: keep only spike-in proteins
#######################################################################################################################

## DDA
proteins_dda_intensities <- proteins_dda_spike_in %>% filter(str_detect(Protein.IDs,
                                                                pattern = paste(protein_uniprot_id,collapse = '|')))

## DIA
proteins_dia_intensities <- proteins_dia_spike_in %>% filter(str_detect(Protein.Group,
                                                                pattern = paste(protein_uniprot_id,collapse = '|')))

#######################################################################################################################
## figure 4.8: clean-up protein names column in dda dataset
#######################################################################################################################

protein_id_vector <- c("P02769", "P13717",
                       proteins_dda_intensities$Protein.IDs[3:6],
                       "P06804", "P09056",
                       "P11974", 
                       proteins_dda_intensities$Protein.IDs[10:13], 
                       "P51642", "P52480", "Q07731",
                       proteins_dda_intensities$Protein.IDs[17:18],
                       "Q8BTJ4",
                       proteins_dda_intensities$Protein.IDs[20:28])

proteins_dda_intensities$Protein.IDs <- protein_id_vector

#######################################################################################################################
## figure 4.8: rename protein columns
#######################################################################################################################

## DDA
proteins_dda_intensities <- proteins_dda_intensities %>% dplyr::rename(Protein = Protein.IDs)

## DIA
proteins_dia_intensities <- proteins_dia_intensities %>% dplyr::rename(Protein = Protein.Group)

#######################################################################################################################
## figure 4.8: standardize experiment column names
#######################################################################################################################

experiment_labels_dda <- names(proteins_dda_intensities)[2:40]
experiment_labels_dia <- names(proteins_dia_intensities)[2:40]

## Ordering of intensity column names for MQ output
Latin_square_index <- as.integer(sort(as.character(c(1:13))))

experiment_id_dda <-  c(paste0("Platelet_1_Mix_", Latin_square_index),
                        paste0("Platelet_2_Mix_", Latin_square_index),
                        paste0("Platelet_3_Mix_", Latin_square_index))

experiment_id_dia <-  c(paste0("Platelet_1_Mix_", 1:13),
                        paste0("Platelet_2_Mix_", 1:13),
                        paste0("Platelet_3_Mix_", 1:13))

platelet_id <- c(rep(1,13),
                 rep(2,13),
                 rep(3,13))

mix_id <- c(rep(1:13,3))

names(platelet_id) <- experiment_id_dia
names(mix_id) <- experiment_id_dia

## DDA
proteins_dda_intensities <- proteins_dda_intensities %>% rename_at(vars(experiment_labels_dda),
                                                                   ~ experiment_id_dda)

## DIA
proteins_dia_intensities <- proteins_dia_intensities %>% rename_at(vars(experiment_labels_dia),
                                                                   ~ experiment_id_dia)

#######################################################################################################################
## figure 4.8: long version of protein tables
#######################################################################################################################

## DDA
proteins_dda_intensities_long <- proteins_dda_intensities %>% pivot_longer(-c(Protein), 
                                                                           names_to = "Experiment",
                                                                           values_to = "Intensity")


## DIA
proteins_dia_intensities_long <- proteins_dia_intensities %>% pivot_longer(-c(Protein), 
                                                                           names_to = "Experiment",
                                                                           values_to = "Intensity")

#######################################################################################################################
## figure 4.8: re-level experiment factors
#######################################################################################################################

experiment_levels <- c(paste0("Platelet_1_Mix_", 1:13),
                       paste0("Platelet_2_Mix_", 1:13),
                       paste0("Platelet_3_Mix_", 1:13))

proteins_dda_intensities_long$Experiment <- factor(proteins_dda_intensities_long$Experiment,
                                                   levels = experiment_levels)

proteins_dia_intensities_long$Experiment <- factor(proteins_dia_intensities_long$Experiment,
                                                   levels = experiment_levels)

#######################################################################################################################
## figure 4.8: add platelet, mix, and concentration columns
#######################################################################################################################

# start here
mix_1 <- mix_2 <- integer(length = 1092)

platelet_1 <- platelet_2 <-  integer(length = 1092)

for(i in 1:1092){
    
    mix_1[i] <- mix_id[proteins_dda_intensities_long$Experiment[i]]
    mix_2[i] <- mix_id[proteins_dia_intensities_long$Experiment[i]]
    
    platelet_1[i] <- platelet_id[proteins_dda_intensities_long$Experiment[i]]
    platelet_2[i] <- platelet_id[proteins_dia_intensities_long$Experiment[i]]
    
}

## add the columns
proteins_dda_intensities_long <- proteins_dda_intensities_long %>% mutate(Platelet = platelet_1, Mix = mix_1)
proteins_dia_intensities_long <- proteins_dia_intensities_long %>% mutate(Platelet = platelet_2, Mix = mix_2)

## replace 0 with na
proteins_dda_intensities_long <- proteins_dda_intensities_long %>% mutate(Intensity = na_if(Intensity, 0))
proteins_dia_intensities_long <- proteins_dia_intensities_long %>% mutate(Intensity = na_if(Intensity, 0))

## Remove bovine serum albumin slack protein
proteins_dda_intensities_long <- proteins_dda_intensities_long %>% filter(Protein != "P02769")
proteins_dia_intensities_long <- proteins_dia_intensities_long %>% filter(Protein != "P02769")

## Concentrations
concentration_1 <- concentration_2 <- numeric(length = 1053)

for(i in 1:1053){
    
    concentration_1[i] <- latin_square_matrix[proteins_dda_intensities_long$Mix[i], 
                                              proteins_dda_intensities_long$Protein[i]]
    
    concentration_2[i] <- latin_square_matrix[proteins_dia_intensities_long$Mix[i], 
                                              proteins_dia_intensities_long$Protein[i]]
    
}

## Add concentrations column
proteins_dda_intensities_long <- proteins_dda_intensities_long %>% mutate(Concentration = concentration_1)
proteins_dia_intensities_long <- proteins_dia_intensities_long %>% mutate(Concentration = concentration_2)

## Add dataset tags
proteins_dda_intensities_long <- proteins_dda_intensities_long %>% mutate(Dataset = "DDA")
proteins_dia_intensities_long <- proteins_dia_intensities_long %>% mutate(Dataset = "DIA")

#######################################################################################################################
## figure 4.8: create protein_intensity_data tibble
#######################################################################################################################

## Bind all protein intensity datasets
protein_intensity_data <- rbind(proteins_dda_intensities_long, proteins_dia_intensities_long)

## Re-factor protein column
protein_intensity_data$Protein <- factor(protein_intensity_data$Protein,
                                         levels = spike_in_protein_statistics$Protein[-28])

## log-transform intensity and concentration columns
protein_intensity_data <- protein_intensity_data %>% mutate(Intensity = log2(Intensity)) %>% 
                            mutate(Concentration = log2(Concentration))

#######################################################################################################################
## figure 4.8: summarize protein_intensity_data
#######################################################################################################################

## Remove missing data points
protein_intensity_data_summary <- protein_intensity_data %>% filter(!(is.na(Intensity)))

## Summarize protein_intensity_data
protein_intensity_data_summary <- protein_intensity_data_summary %>% 
    group_by(Experiment, Dataset) %>% 
    summarise(Number_of_proteins = n()) %>% ungroup()

#######################################################################################################################
## figure 4.8: source extra code chunks needed for the figure
#######################################################################################################################

source("extra.R")

#######################################################################################################################
## Normalized protein intensity data 
#######################################################################################################################

## Separate razor and unique datasets
protein_intensity_data_unique <- protein_intensity_data %>% 
                                    mutate(Type = ifelse(Dataset == "DDA", "DDA", "DIA"))

## Collect the median protein intensity measurements
median_data_unique <- protein_intensity_data_unique %>% filter(Protein == "O88766") %>% 
                        group_by(Protein, Type) %>% summarise(MedianIntensity = median(Intensity, na.rm = T)) %>% 
                            ungroup()

## Compute median of median protein intensity measurements
median_gamma_dda_unique <- (median_data_unique %>% filter(Type == "DDA"))$MedianIntensity
median_gamma_dia_unique <- (median_data_unique %>% filter(Type == "DIA"))$MedianIntensity

## Normalized DDA and DIA datasets by median protein intesity
protein_intensity_data_normalized_unique <- bind_rows((protein_intensity_data_unique %>% filter(Type == "DDA") %>% 
                                                           mutate(Intensity = Intensity-median_gamma_dda_unique)),
                                                      (protein_intensity_data_unique %>% filter(Type == "DIA") %>% 
                                                           mutate(Intensity = Intensity-median_gamma_dia_unique)))

## Re-factor levels for non-normalized data
protein_intensity_data_unique$Protein <- factor(protein_intensity_data_unique$Protein,
                                                levels = levels(factor(data_1d$Protein)))

protein_intensity_data_unique$Type <- factor(protein_intensity_data_unique$Type, levels = c("DDA", "DIA"))

## Re-factor levels for normalized data
protein_intensity_data_normalized_unique$Protein <- factor(protein_intensity_data_normalized_unique$Protein,
                                                           levels = levels(factor(data_1d$Protein)))
protein_intensity_data_normalized_unique$Type <- factor(protein_intensity_data_normalized_unique$Type,
                                                        levels = c("DDA", "DIA"))

#######################################################################################################################
## Missing values data
#######################################################################################################################

## Missing values tile plot
missing_values_unique <- protein_intensity_data_unique %>% filter(Protein != "O88766") %>% 
    dplyr::select(Protein,Type, Platelet, Concentration, Intensity)

missing_values_unique <- missing_values_unique %>% mutate(TypeCode = ifelse(Type == "DDA", 1, 2)) 

missing_values_unique$TypeCode[which(is.na(missing_values_unique$Intensity))] <- 0

missing_values_unique <- missing_values_unique %>% group_by(Protein, Concentration, Platelet) %>% 
    summarise(TypeCode = sum(TypeCode)) %>% ungroup()

type_category_labels <- c("Neither",
                          "DDA",
                          "DIA",
                          "Both")

names(type_category_labels) <- c(0,1,2,3)

missing_values_unique$Protein <- factor(missing_values_unique$Protein, 
                                        levels = levels(factor(data_2_top$Protein))[-1])

missing_values_unique <- missing_values_unique %>% 
    mutate(ID = c(sapply(rev(LETTERS), function(x) rep(x,39), simplify = TRUE)))

missing_values_unique$Platelet <- factor(missing_values_unique$Platelet,
                                         levels = rev(levels(factor(missing_values_unique$Platelet))))

missing_values_unique <- missing_values_unique %>% 
    mutate(ConcentrationLevel = ifelse(Concentration < 8.8, "Low",
                                       ifelse(Concentration > 13, "High", "Medium")))

missing_values_unique$ConcentrationLevel <- factor(missing_values_unique$ConcentrationLevel,
                                                   levels = c("Low", "Medium", "High"))

missing_values_unique <- missing_values_unique %>% mutate(Concentration = round(Concentration, digits = 1))

#######################################################################################################################
## Normalized scatterplot
#######################################################################################################################

concentration_level_tibble <- tibble(Concentration = (levels(factor(data_2_top$Concentration)))[1:13],
                                     Level = c(rep("Low", 4),rep("Medium", 5),rep("High", 4)))

scatterplot_normalized_unique <- protein_intensity_data_normalized_unique %>% filter(Protein != "O88766")

scatterplot_normalized_unique$Protein <- factor(scatterplot_normalized_unique$Protein,
                                                levels=rev(levels(factor(scatterplot_normalized_unique$Protein))))

concentration_level_labels_unique <- character(length = length(scatterplot_normalized_unique$Concentration))

for(i in 1:length(concentration_level_labels_unique)){
    
    concentration_level_labels_unique[i] <- concentration_level_tibble$Level[
        which(concentration_level_tibble$Concentration == 
                  as.character(scatterplot_normalized_unique$Concentration[i]))]
}

scatterplot_normalized_unique <- scatterplot_normalized_unique %>% 
    mutate(ConcentrationLevel = concentration_level_labels_unique)

scatterplot_normalized_unique <- scatterplot_normalized_unique %>% 
    mutate(TypeConcentration = paste0(Type," _ ", ConcentrationLevel))

#######################################################################################################################
## Figure - 2.1 (Normalized scatterplot for unique) [newest version - February 1 2023]
#######################################################################################################################

plot_2_scatterplot_color_alpha <- 0.1
plot_2_scatterplot_axis_title <- 9
plot_2_axis_text_size <- 7
plot_legend_title_size <- 8
plot_2_scatterplot_point_stroke <- 0.25
plot_2_scatterplot_expand_y <- 2.5
plot_2_scatterplot_color_alpha <- 0.5

## breaks function
count_1 <- 0
count_2 <- 0
count_3 <- 0

count_1_x <- 0
count_2_x <- 0
count_3_x <- 0

breaks_fun_1 <- function(y) {
    count_1 <<- count_1+1L
    if(count_1 == 1){
        c(-7,5)
    } else {
        c()
    }
}

breaks_fun_2<- function(y) {
    count_2 <<- count_2+1L
    if(count_2 == 1){
        c(-5,4)
    } else {
        c()
    }
}

breaks_fun_3<- function(y) {
    count_3 <<- count_3+1L
    if(count_3 == 1){
        c(-1,8)
    } else {
        c()
    }
}

breaks_fun_1_x <- function(x) {
    count_1_x <<- count_1_x+1L
    if(count_1_x >= 1){
        as.numeric(levels(factor(data_2_top$Concentration)))[c(1,4)]
    } else {
        c()
    }
}

breaks_fun_2_x <- function(x) {
    count_2_x <<- count_2_x+1L
    if(count_2_x >= 1){
        as.numeric(levels(factor(data_2_top$Concentration)))[c(5,9)]
    } else {
        c()
    }
}

breaks_fun_3_x <- function(x) {
    count_3_x <<- count_3_x+1L
    if(count_3_x >= 1){
        as.numeric(levels(factor(data_2_top$Concentration)))[c(10,13)]
    } else {
        c()
    }
}

log_jump <- (as.numeric(levels(factor(data_2_top$Concentration)))[2]-
                 as.numeric(levels(factor(data_2_top$Concentration)))[1])

concentration_labels_2 <- as.character(round(as.numeric((levels(factor(data_2_top$Concentration)))[1:13]), digits = 1))
concentration_labels_2[c(1:5,8,13)] <- c("  6.3", "  7.1", "  7.9", "  8.7", " 9.6", "12.0", "16.0")

plot_2_scatterplot_1 <- ggplot(data = scatterplot_normalized_unique)+
    geom_point(data = scatterplot_normalized_unique %>% filter(Concentration <13) %>% filter(Concentration < 9),
               size = 1.5, stroke = plot_2_scatterplot_point_stroke,
               mapping = aes(x = as.numeric(Concentration),
                             y = Intensity, shape = as_factor(Platelet)),
               inherit.aes = FALSE)+
    scale_shape_manual(values = c(0:2))+
    scale_x_continuous(limits = c(as.numeric(levels(factor(data_2_top$Concentration)))[1] - log_jump/2,
                                  as.numeric(levels(factor(data_2_top$Concentration)))[4] + log_jump/2),
                       breaks = as.numeric(levels(factor(data_2_top$Concentration)))[1:4],
                       labels = concentration_labels_2[1:4])+
    scale_y_continuous(breaks = breaks_fun_1)+
    facet_grid(Protein~Type, labeller = labeller(Protein = protein_labels_3_c))+
    guides(color=guide_legend(override.aes = list(size=3)),
           size=guide_legend(override.aes = list(size=3)))+
    theme_linedraw()+
    labs(x = expression(log[2]~"spiked-in"~protein~concentration),
         y = expression(Standardised~log[2]~protein~intensity),
         color = "Concentration level", shape = "Donor")+
    theme(
        axis.title.x = element_text(size = plot_2_scatterplot_axis_title, color = alpha("white", 0)),
        axis.title.y = element_text(size = plot_2_scatterplot_axis_title, margin = margin(l = 2,r=-10)),
        axis.text.x = element_text(angle = 90, size = plot_2_axis_text_size, 
                                   hjust = 1, vjust = 0.5),
        axis.text.y = element_text(angle = 0, size = plot_2_axis_text_size),
        panel.spacing.x = unit(1, "pt"),
        panel.spacing.y = unit(0, "pt"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_rect(fill = alpha("red",plot_2_scatterplot_color_alpha)),
        strip.text.x = element_text(size = 8, color = "black", margin = margin(1,1,1,1, "pt")),
        strip.text.y = element_text(color = "white", size = 8,
                                    margin = margin(1,1,1,1, "pt"), angle = 0),
        strip.background.x = element_rect(color = "black", fill = "white"),
        strip.background.y = element_blank(),
        legend.position = "bottom",
        legend.justification = "right",
        plot.margin =  margin(t = 0, b = 0, r = -4, l = 0, unit = "pt")
    ) 

plot_2_scatterplot_2 <- ggplot(data = scatterplot_normalized_unique)+
    geom_point(data = scatterplot_normalized_unique %>% filter(Concentration < 13) %>% filter(Concentration > 9),
               size = 1.5, stroke = plot_2_scatterplot_point_stroke,
               mapping = aes(x = as.numeric(Concentration),
                             y = Intensity, shape = as_factor(Platelet)),
               inherit.aes = FALSE)+
    scale_shape_manual(values = c(0:2))+
    scale_x_continuous(limits = c(as.numeric(levels(factor(data_2_top$Concentration)))[5] - log_jump/2,
                                  as.numeric(levels(factor(data_2_top$Concentration)))[9] + log_jump/2),
                       breaks = as.numeric(levels(factor(data_2_top$Concentration)))[5:9],
                       labels = concentration_labels_2[5:9])+
    scale_y_continuous(breaks = breaks_fun_2)+
    facet_grid(Protein~Type, labeller = labeller(Protein = protein_labels_3_c))+
    guides(color=guide_legend(override.aes = list(size=3)),
           size=guide_legend(override.aes = list(size=3)))+
    theme_linedraw()+
    labs(x = expression(log[2]~"spiked-in"~protein~concentration), y = "Measured intensity",
         color = "Concentration level", shape = "Donor")+
    theme(
        axis.title.x = element_text(size = plot_2_scatterplot_axis_title),
        axis.title.y = element_text(size = plot_2_scatterplot_axis_title, color = alpha("white", 0),
                                    margin = margin(l = 2,r=-10)),
        axis.text.x = element_text(angle = 90, size = plot_2_axis_text_size,
                                   hjust = 1, vjust = 0.5),
        axis.text.y = element_text(angle = 0, size = plot_2_axis_text_size),
        panel.spacing.x = unit(1, "pt"),
        panel.spacing.y = unit(0, "pt"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_rect(fill = alpha("blue",plot_2_scatterplot_color_alpha)),
        strip.text.x = element_text(size = 8, color = "black", margin = margin(1,1,1,1, "pt")),
        strip.text.y = element_text(color = "white", size = 8,
                                    margin = margin(1,1,1,1, "pt"), angle = 0),
        strip.background.x = element_rect(color = "black", fill = "white"),
        strip.background.y = element_blank(),
        legend.position = "bottom",
        legend.justification = "right",
        plot.margin =  margin(t = 0, b = 0, r = -2, l = -2, unit = "pt")
    ) 

plot_2_scatterplot_3 <- ggplot(data = scatterplot_normalized_unique)+
    geom_point(data = scatterplot_normalized_unique %>% filter(Concentration > 13),
               size = 1.5, stroke = plot_2_scatterplot_point_stroke,
               mapping = aes(x = as.numeric(Concentration),
                             y = Intensity, shape = as_factor(Platelet)),
               inherit.aes = FALSE)+
    scale_shape_manual(values = c(0:2))+
    scale_x_continuous(limits = c(as.numeric(levels(factor(data_2_top$Concentration)))[10] - log_jump/2,
                                  as.numeric(levels(factor(data_2_top$Concentration)))[13] + log_jump/2),
                       breaks = as.numeric(levels(factor(data_2_top$Concentration)))[10:13],
                       labels = concentration_labels_2[10:13])+
    scale_y_continuous(breaks = breaks_fun_3)+
    facet_grid(Protein~Type, labeller = labeller(Protein = protein_labels_3_c))+
    guides(color=guide_legend(override.aes = list(size=3)),
           size=guide_legend(override.aes = list(size=3)))+
    theme_linedraw()+
    labs(x = expression(log[2]~"spiked-in"~protein~concentration), y = "Measured intensity",
         color = "Concentration level", shape = "Donor")+
    theme(
        axis.title.x = element_text(size = plot_2_scatterplot_axis_title, color = alpha("white", 0)),
        axis.title.y = element_text(size = plot_2_scatterplot_axis_title, color = alpha("white", 0),
                                    margin = margin(l = 2,r=-10)),
        axis.text.x = element_text(angle = 90, size = plot_2_axis_text_size,
                                   hjust = 1, vjust = 0.5),
        axis.text.y = element_text(angle = 0, size = plot_2_axis_text_size),
        panel.spacing.x = unit(1, "pt"),
        panel.spacing.y = unit(0, "pt"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_rect(fill = alpha("green",plot_2_scatterplot_color_alpha)),
        strip.text.x = element_text(size = 8, color = "black", margin = margin(1,1,1,1, "pt")),
        strip.text.y = element_text(color = "black", size = 8,
                                    margin = margin(1,1,1,1, "pt"), angle = 0),
        strip.background.x = element_rect(color = "black", fill = "white"),
        strip.background.y = element_blank(),
        legend.position = "bottom",
        legend.justification = "right",
        plot.margin =  margin(t = 0, b = 0, r = 0, l = -4, unit = "pt"),
        plot.background = element_rect(fill='transparent', color=NA)
    ) 

scatterplot_right_y_axis <- ggplot()+
    labs(y = "Spike-in protein")+
    theme_linedraw()+
    theme(
        panel.grid.minor.x=element_blank(),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.border = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = plot_2_scatterplot_axis_title, margin = margin(t = 0, r = 0, b = 0, l = 0),
                                    angle = 270),
        strip.background = element_rect(fill = "black", colour = NA),
        legend.position = "none",
        plot.margin =  margin(t = 0, b = 0, r = 5, l = 0, unit = "pt")
    ) 

## Assemble bottom part of figure 3
plot_2_scatterplot <- ggarrange(plot_2_scatterplot_1, plot_2_scatterplot_2, plot_2_scatterplot_3,
                                ncol = 3,common.legend = TRUE, legend="top",
                                widths = c(4,5,4), legend.grob = get_legend(plot_2_scatterplot_legend))

scatterplot_unique <- ggarrange(plot_2_scatterplot,scatterplot_right_y_axis, ncol = 2,
                                widths = c(1,0.05))

#######################################################################################################################
## Figure - 1.2 (Missing value plots for unique)
#######################################################################################################################

missing_value_plot_right_y_axis <- ggplot()+
    labs(y = "Spike-in protein")+
    theme_linedraw()+
    theme(
        panel.grid.minor.x=element_blank(),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.border = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = plot_2_scatterplot_axis_title, margin = margin(t = 0, r = 0, b = 0, l = 0),
                                    angle = 270),
        strip.background = element_rect(fill = "black", colour = NA),
        legend.position = "none",
        plot.margin =  margin(t = -29, b = 0, r = 5, l = 0, unit = "pt")
    ) 

## donor labels
donor_id_labels <- paste("Donor",1:3)
names(donor_id_labels) <- c(1,2,3)

## Mix IDs
mix_id_labels <- paste0("Mix ", 1:13)
names(mix_id_labels) <- 1:13

## Unique peptides
count_1_y <- 0

breaks_fun_1_y <- function(y) {
    count_1_y <<- count_1_y+1L
    if(count_1_y == 1){
        c(1:3)
    } else {
        c()
    }
}

concentration_labels <- as.character(round(as.numeric((levels(factor(data_2_top$Concentration)))[1:13]), digits = 1))
concentration_labels[c(8,13)] <- c("12.0", "16.0")

missing_values_unique <- missing_values_unique %>% mutate(Platelet = as.integer(Platelet))

missing_values_plot_legend <- ggplot()+
    geom_point(data = missing_values_unique %>% filter(Platelet == 1),
               size = 1.5, stroke = 0.5,
               mapping = aes(x = as_factor(Protein),
                             y = Concentration, color = as_factor(TypeCode)),
               inherit.aes = FALSE)+
    guides(color=guide_legend(override.aes = list(size=3, shape = 22, color = "black", stroke = 0.25,
                                                  fill = c("white",brewer.pal(12, "Set3")[c(7,3,4)])),
                              title="Quantified in:",
                              nrow = 1, byrow = TRUE , title.position="left"))+
    theme_linedraw()+
    scale_color_manual(labels = c("Neither", "DDA", "DIA", "Both"),
                       values = c("white",brewer.pal(12, "Set3")[c(7,3,4)]))+
    labs(x = "Concentration level", y = "Measured intensity",
         color = "Concentration level group:", shape = "Donor:")+
    theme(
        panel.grid.minor.x = element_blank(),
        strip.text = element_text(size = 8, color = "white", margin = margin(1,0,1,0, "pt")),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, size = 6),
        axis.text.y = element_text(angle = 0, size = 6),
        axis.title.y = element_blank(),
        panel.spacing.x = unit(0, "pt"),
        panel.spacing.y = unit(0, "pt"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_rect(fill = alpha("green",0.5),
                                        colour = "black", 
                                        size = 0.2, linetype = "solid"),
        legend.position = "top",
        legend.justification = "center",
        legend.box="vertical",
        legend.spacing.x = unit(0, 'pt'),
        legend.spacing.y = unit(-2, 'pt'),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.margin=margin(t = 15, b = -2),
        plot.margin =  margin(t = 10, b = 0, r = 2, l = 5, unit = "pt")
    ) 

missing_values_plot_unique <- ggplot(data = missing_values_unique,
                                     aes(x = as_factor(Concentration), y = as_factor(Platelet),
                                         fill = as_factor(TypeCode)))+
    geom_tile(color = "black")+
    facet_grid(ID~ConcentrationLevel, labeller = labeller(Platelet = donor_id_labels), scales = "free")+
    scale_x_discrete(breaks = as.character(round(as.numeric((levels(factor(data_2_top$Concentration)))[1:13]),
                                                 digits = 1)),
                     labels = concentration_labels,
                     expand = c(0,0))+
    scale_fill_manual(labels = type_category_labels, values = c("white",brewer.pal(12, "Set3")[c(7,3,4)]))+
    labs(x = expression(log[2]~"spiked-in"~protein~concentration), y = "Donor")+
    guides(fill=guide_legend(title="Quantified in:",
                             nrow = 1, byrow = TRUE , title.position="left"))+
    theme_linedraw()+
    theme(
        axis.title.x = element_text(size = plot_2_scatterplot_axis_title),
        axis.title.y = element_text(size = plot_2_scatterplot_axis_title, margin = margin(l = 2)),
        axis.text.x = element_text(angle = 90, size = plot_2_axis_text_size,
                                   hjust = 1, vjust = 0.5),
        axis.text.y = element_text(angle = 0, size = plot_2_axis_text_size),
        legend.title = element_text(size = 8, angle = 0),
        legend.text = element_text(size = 8),
        legend.key.size = unit(5, "pt"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.justification = "right",
        strip.text.x = element_text(size = 8, color = "black", margin = margin(1,1,1,1, "pt")),
        strip.background.x = element_rect(color = "black", fill = "white"),
        strip.text.y = element_text(color = "black", size = 8, angle = 0,
                                    margin = margin(1,1,1,1, "pt")),
        strip.background.y = element_blank(),
        panel.spacing = unit(0, "pt"),
        legend.margin = margin(t = 16, b = -2),
        plot.margin =  margin(t = 2, b = 0, r = 2, l = 5, unit = "pt")
    ) 

missing_values_plot_unique <- ggarrange(missing_values_plot_unique,missing_value_plot_right_y_axis, ncol = 2,
                                        legend.grob = get_legend(missing_values_plot_unique),
                                        widths = c(1,0.05))

#######################################################################################################################
## Figure - 2 with normalized scatterplot and missing values plot
#######################################################################################################################

### Combined plot
ggsave(plot = ggarrange(scatterplot_unique, missing_values_plot_unique,
                        nrow = 1, labels = c("A","B")),
       filename = paste0("figure_4_8.pdf"),
       path = paste0("./Figures/"),
       width = 7, height = 7, units = "in",
       dpi = 320)

