#######################################################################################################################
## Median samples for each pair
#######################################################################################################################

mix_to_median_pair <- c(7:1,13:8)

names(mix_to_median_pair) <- 1:13

spike_in_design_median_mix <- spike_in_protein_statistics %>% mutate(Median = mix_to_median_pair[PairID])

spike_in_design_median_mix <- spike_in_design_median_mix %>% filter(ProteinID %in% LETTERS) %>% 
    dplyr::select(PairID, ProteinID, Protein, Median)


#######################################################################################################################
## Create summarized experiment objects
#######################################################################################################################

## Create experimental design object
experimental_design_matrix <- data.frame(label = names(peptide_table_dda[,c(5:43)]),
                                         condition = rep(paste0("Mix_",1:13),3),  
                                         replicate = rep(1:3,each = 13))

## Re-factor Conditions
experimental_design_matrix$condition <- factor(experimental_design_matrix$condition,
                                               levels = paste0("Mix_",1:13))

## Re-arrange experimental design object
experimental_design_matrix <- experimental_design_matrix %>% arrange(condition)

#######################################################################################################################
## Sample reference matrices
#######################################################################################################################

lower_samples <- matrix(data = 0, nrow  = 13, ncol = 6)
higher_samples <- matrix(data = 0, nrow  = 13, ncol = 6)
rownames(lower_samples) <- rownames(higher_samples) <- 1:13
colnames(lower_samples) <- colnames(higher_samples) <- 1:6

## Fill the lower matrix
lower_samples[1,] <- 13:8
lower_samples[2,] <- c(1,13:9)
lower_samples[3,] <- c(2:1,13:10)
lower_samples[4,] <- c(3:1,13:11)
lower_samples[5,] <- c(4:1,13:12)
lower_samples[6,] <- c(5:1,13:13)
lower_samples[7,] <- c(6:1)
lower_samples[8,] <- c(7:2)
lower_samples[9,] <- c(8:3)
lower_samples[10,] <- c(9:4)
lower_samples[11,] <- c(10:5)
lower_samples[12,] <- c(11:6)
lower_samples[13,] <- c(12:7)

## Fill the higher matrix
higher_samples[1,] <- c(2:7)
higher_samples[2,] <- c(3:8)
higher_samples[3,] <- c(4:9)
higher_samples[4,] <- c(5:10)
higher_samples[5,] <- c(5:10)+1
higher_samples[6,] <- c(5:10)+2
higher_samples[7,] <- c(5:10)+3
higher_samples[8,] <- c(9:13,1)
higher_samples[9,] <- c(10:13,1:2)
higher_samples[10,] <- c(11:13,1:3)
higher_samples[11,] <- c(12:13,1:4)
higher_samples[12,] <- c(13,1:5)
higher_samples[13,] <- c(1:6)
