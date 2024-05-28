#######################################################################################################################
## Test
#######################################################################################################################

tryptic_peptide_database_isoform %>% filter(Protein == current_protein_id) %>% arrange(Start)

generate_tryptic_peptides_prob(max_missed_cleavages = 0,
                               mc_prob = 1,
                               Protein = current_protein_sequence,
                               remove_proline_sites = FALSE,
                               methionine_peptides = FALSE)

# Takes 6 seconds per 1000 simulations of 10 proteins (linear in number of proteins and number of simulations)
# system.time(test <- as.character(sapply(rep(0,1000*10),
#                                         function(x) generate_tryptic_peptides_prob(max_missed_cleavages = x,
#                                                                                    missed_cleavage_probability = 0,
#                                                                                    Protein = current_protein_sequence,
#                                                                                    remove_proline_sites = FALSE,
#                                                                                    methionine_peptides = FALSE))))

## peptide sharing check
test <- tryptic_peptide_database_canonical %>% group_by(Sequence, Length) %>% summarise(Count = n()) %>% ungroup()
test_2 <- test %>% filter(Count == 1)

length(unique((tryptic_peptide_database_canonical %>% filter(Sequence %in% test_2$Sequence) %>% 
                   filter(Sequence %in% peptide_table$Sequence))$Protein))

# only 75 of the proteins generate peptides that are only unique to themselves

## missed cleavage annotation testing
test <- tryptic_peptide_database_isoform %>% 
    filter(Protein %in% simulated_tryptic_digest$Protein) %>% 
    filter(Sequence %in% simulated_tryptic_digest$Sequence)

test <- test %>% arrange(Sequence, Protein, MissedCleavages)
test_2 <- simulated_tryptic_digest%>% arrange(Sequence, Protein)
test[1105:1115,]
test_2[1105:1115,]

## Missed cleavage proportions
test <- simulated_tryptic_digest
test <- test %>% group_by(Missed_Cleavages) %>% summarise(Count = sum(Count)) %>% ungroup()
sum(test$Count)
test
simulated_tryptic_digest
simulated_tryptic_digest
test <- simulated_tryptic_digest %>% filter(Count > 100)
table(simulated_tryptic_digest$Missed_Cleavages)
table(test$Missed_Cleavages)

## abundance distribution
test <- copy_number_matrix_sample_final[,1]
plot(x=1:1000, y = sort(test))

## Model detection with beta distribution
hist(nchar(peptide_abundance_data$Sequence))

## peptide distributions
test <- peptide_abundance_data_detectable %>% filter(Protein == "F8WD26")
test <- peptide_abundance_data_detectable %>% group_by(Protein) %>% summarise(Count = n()) %>% ungroup()
test %>% arrange(desc(Count))

## PCA
peptide_abundance_matrix
test <- tibble(rownames_to_column(data.frame(t(peptide_abundance_matrix)), "Experiment"))

peptides.pca <- prcomp(test[,-1], center = TRUE, scale. = TRUE)

summary(peptides.pca)
autoplot(peptides.pca, data = test, colour = "Experiment")+
    guides(color=guide_legend(title = "LC-MS/MS run:"))+
    scale_color_manual(values = c('#000000', "#E69F00", "#56B4E9",))
    theme_linedraw()+
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "top"
        )
        

## cpndition-level pca

test <- peptide_abundance_data_detectable %>% dplyr::select(-c(Protein:Missed_Cleavages)) %>% 
    dplyr::select(-c(ionisation:Length))
test <- as.matrix(test[,c(2:3)])
rownames(test) <- peptide_abundance_data_detectable$Sequence

test <- tibble(rownames_to_column(data.frame(t(test)), "Experiment"))

peptides.pca <- prcomp(test[,-1], center = TRUE, scale. = TRUE)
peptides.pca
autoplot(peptides.pca)

## 
test <- tibble(rownames_to_column(data.frame(t(peptide_abundance_matrix)), "Experiment"))

## 
hist(colMeans(test[4:6,2:6093]) - colMeans(test[1:3,2:6093]))
hist(peptide_abundance_matrix[,1])
## zscore transformation
mean(a)

# find standard deviation
sd(a)

# calculate z
z_scores <- (peptide_abundance_matrix - mean(peptide_abundance_matrix)) / sd(peptide_abundance_matrix)
hist(z_scores)
plot(peptide_abundance_matrix,z_scores)

## probability of detection as a function of 
plot(peptide_abundance_matrix,)

##
table(is.na(c(peptide_intensity_matrix)))
hist(peptide_table$Missed_Cleavages)

