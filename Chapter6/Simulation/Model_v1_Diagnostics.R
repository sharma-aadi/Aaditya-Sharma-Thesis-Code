#########################################################################################################################################################
## Diagnostics 
#########################################################################################################################################################

cat("\n")

print(paste0("Dimension of matrix M [Peptides X Proteins]: ", dim(M)[1], " by ", dim(M)[2]))
print(paste0("Dimension of matrix eta [Proteins X Samples]: ", dim(eta)[1], " by ", dim(eta)[2]))
print(paste0("Dimension of matrix Y [Peptides X Samples]: ", dim(Y)[1], " by ", dim(Y)[2]))

cat("\n")

print(paste0(dim(M)[1], " peptides generated from in-silico digestion of ", dim(M)[2], " proteins."))

cat("\n")

Peptide_Counts <- Peptide_Table %>% group_by(Protein) %>% summarise(Number_of_Peptides = n())

print(Peptide_Counts)

cat("\n")

print(paste0("Simulation finished. Starting inference..."))

