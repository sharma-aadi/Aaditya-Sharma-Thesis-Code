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
