#######################################################################################################################
## Target_CN_Generator()
#######################################################################################################################

Target_CN_Generator <- function(Minimum_Copy_Number_Per_Cell = 80,
                                Ratio = 1.75, 
                                Number_of_Cells = 1,
                                Number_of_Latin_Square_Protein_Pairs = 13){
    
    # This function is used to generate a Copy number profile for the latin square given Min CN (M) and Ratio (R)
    
    # Parameters
    ## Minimum_Copy_Number_Per_Cell is the first non-zero per cell copy number
    ## Ratio is the ratio between consecutive non-zero concentration levels
    ## Number_of_Cells = 1 generates per cell CN profile and Number_of_Cells = 10^7 generates CN_profile for pilot
    ## Number_of_Latin_Square_Protein_Pairs is number of pairs of Latin-square proteins
    
    ## Generate a vector with #LS - 2 M's
    Target_CN <- c(rep(Minimum_Copy_Number_Per_Cell, 2*(Number_of_Latin_Square_Protein_Pairs)))
    
    ## Multiply throughtout by appropriate powers of R
    Target_CN <- Target_CN*(Ratio^c(c(seq(0, (Number_of_Latin_Square_Protein_Pairs - 1), by = 1),
                                      seq(0, (Number_of_Latin_Square_Protein_Pairs - 1), by = 1))))
    
    ## Sort so that the pairs are together and the profile is monotone increasing 
    Target_CN <- sort(Target_CN)
    
    ## Multiply by number of cells if more than per cell profile is required
    Target_CN <- Target_CN*(Number_of_Cells)
    
    ## Add the copy number for the median protein
    Target_CN <- c(Target_CN, median(Target_CN))
    
    return(Target_CN)
}

#######################################################################################################################
## Generate_Rotations()
#######################################################################################################################

Generate_Rotations <- function(Number_of_Pairs = 13, Sample_Number = 1){
    return(c((0:(Number_of_Pairs*2 - 1) + 2*(Sample_Number - 1))%%(Number_of_Pairs*2) + 1, 
             (Number_of_Pairs*2 + 1)))
}

#######################################################################################################################
## create latin square matrix
#######################################################################################################################

## M = 80, R = 1.75
Target_Copy_Number_80 <- Target_CN_Generator(Minimum_Copy_Number_Per_Cell = 80, Ratio = 1.75)

## Initialize latin_square_matrix

latin_square_matrix <- matrix(data = 0, nrow = 13, ncol = 27)
colnames(latin_square_matrix) <- spike_in_protein_statistics$Protein[1:27]
rownames(latin_square_matrix) <- paste0(1:13)

for (i in 1:13) {
    
    latin_square_matrix[i,] <- Target_Copy_Number_80[Generate_Rotations(Sample_Number = i)]   
    
}
