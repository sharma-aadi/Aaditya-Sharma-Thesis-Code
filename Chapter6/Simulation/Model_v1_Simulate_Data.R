#######################################################################################################################
## Part - 0.1: Store Current Directory
#######################################################################################################################

Current_Directory <- getwd()

#######################################################################################################################
## Part - 0.2: Define proteome file
#######################################################################################################################

Proteome <- paste0(Proteome_Directory, Simulation_Proteome)

#######################################################################################################################
## Part - 1: Select Proteins
#######################################################################################################################

## Select protein sequences

print(paste0("Reading in Proteome."))

Input_Proteome <- read.fasta(file = Proteome, seqtype = "AA", as.string = TRUE, 
                             set.attributes = TRUE, legacy.mode = TRUE, seqonly = FALSE, strip.desc = FALSE,
                             whole.header = FALSE,
                             bfa = FALSE, sizeof.longlong = .Machine$sizeof.longlong,
                             endian = .Platform$endian, apply.mask = TRUE)

Total_Proteins <- length(Input_Proteome)

## Set Seed and sample
if(Simulation_Seed != 0){
  set.seed(Simulation_Seed)
}

## Check whether sampling needs to be done or not
if(Number_of_Proteins == Total_Proteins){
  Selected_Proteins <- 1:Number_of_Proteins
} else {
  if(Number_of_Proteins > Total_Proteins){
    print(paste0("Proteome file does not contain sufficient proteins for simulation"))
    
  }
  Selected_Proteins <- sample(x = 1:Total_Proteins, size = Number_of_Proteins)
}

print(sort(Selected_Proteins))

print(paste0(Number_of_Proteins, " proteins selected."))

#######################################################################################################################
## Part - 2: Generate tryptic peptides from selected proteins
#######################################################################################################################

## Compile needed functions
source(In_Silico_Digestion)

## Generte Table of Tryptic Peptides
Peptide_Table <- Generate_Peptide_Database(Proteome = Input_Proteome, Maximum_Missed_Cleavages = Max_Missed_Cleavages,
                                           Minimum_Length = 1, 
                                           Protein_Indices = Selected_Proteins)
## Subset peptide length
Peptide_Table<- Peptide_Table %>% filter(Length >= min(Peptide_Length_Range) & Length <= max(Peptide_Length_Range)) %>%
                  arrange(Protein)

#######################################################################################################################
## Part - 2: Generate matrix M from Peptide_Table
#######################################################################################################################

## Generate Matrix
M <- matrix(data = 0, nrow = length(unique(Peptide_Table$Peptide)), ncol = length(unique(Peptide_Table$Protein)))

rownames(M) <- unique(Peptide_Table$Peptide)
colnames(M) <- unique(Peptide_Table$Protein)

## Populate Matrix
for(i in 1:length(Peptide_Table$Peptide)){
  M[Peptide_Table$Peptide[i], Peptide_Table$Protein[i]] <- Peptide_Table$Multiplicity[i]
}

#######################################################################################################################
## Part - 3: Generate Protein concentrations (eta matrix)
#######################################################################################################################

eta <- matrix(data = rgamma(n = length(unique(Peptide_Table$Protein))*Number_of_Samples,
                            shape = shape_protein,
                            rate = rate_protein), 
              ncol = Number_of_Samples,
              nrow = length(unique(Peptide_Table$Protein)))

colnames(eta) <- paste0("Sample_",1:Number_of_Samples)
rownames(eta) <- colnames(M)

#######################################################################################################################
## Part - 4: Generate Peptide XIC's (Y matrix)
#######################################################################################################################

Y <- matrix(data = 0,
            ncol = Number_of_Samples,
            nrow = nrow(M))

rownames(Y) <- rownames(M)
colnames(Y) <- colnames(eta)

for(i in 1:nrow(Y)){
  for(j in 1:Number_of_Samples){
    Y[i,j] <- rlnorm(n = 1, meanlog = log(dot(M[i,], exp(eta[,j]))) - ((Peptide_SD^2)/2), Peptide_SD)
  }
}

#######################################################################################################################
## Part - 5: Print simulation information
#######################################################################################################################

source(Diagnostics)

#######################################################################################################################
## Part - 6: Reset to Original Directory
#######################################################################################################################

setwd(Current_Directory)
