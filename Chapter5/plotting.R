#######################################################################################################################
## Plotting ideas
#######################################################################################################################

## distribution of protein in sample by condition
test <- as_tibble(copy_number_matrix_sample_final, rownames = "Protein")
test <- test %>% pivot_longer(-c(Protein), names_to = "Condition", values_to = "Abundance")

ggplot(data = test, aes(x = Abundance))+
    geom_histogram(color = "black", bins = 5, alpha = 0.6, aes(fill = Condition))+
    theme_linedraw()

## distribution of contrasts at this point (noise can change this)
hist(copy_number_matrix_sample_final[,2]-copy_number_matrix_sample_final[,1])

## peptide data
plot(x = 1:length(simulated_tryptic_digest$Sequence), y = sort(peptide_abundance_data$condition_2))
hist(copy_number_matrix_sample_final$condition_2-peptide_abundance_data$condition_1)
hist(peptide_abundance_data$condition_2-peptide_abundance_data$condition_1)

ggplot(data = peptide_abundance_data %>% pivot_longer(c(condition_1:condition_2),
                                                      names_to = "Condition", values_to = "Abundance"),
       aes(x = Abundance))+
    facet_grid(Condition~Missed_Cleavages)+
    geom_histogram()+
    theme_linedraw()

## peptide length plot
ggplot(peptide_abundance_data_detectable, aes(x = Length))+
    geom_bar()+
    theme_linedraw()

## peptide data
plot(x = 1:length(peptide_abundance_data_detectable$Sequence),
     y = sort(log10(2^peptide_abundance_data_detectable$condition_2)))
hist(log10(2^peptide_abundance_data_detectable$condition_1))

## After detectability is modeled
hist(peptide_abundance_data_detectable$condition_2-peptide_abundance_data_detectable$condition_1)

## plot rank-intensity plots per run
test <- tibble(rownames_to_column(data.frame(peptide_intensity_matrix), "Sequence"))

test <- test %>% pivot_longer(-c(Sequence), names_to = "Experiment", values_to = "Intensity")

test <- test %>% arrange(Experiment,Intensity) %>% mutate(Rank = rep(1:6092,6))

ggplot(data = test, aes(x = Rank, y = Intensity, colour = Experiment))+
    geom_point()+xlim(c(0,3200))+
        theme_linedraw()

## How many samples are proteins quantified in
test <- protein_table %>% filter(Protein %in% de_protein_IDs)

ggplot(data = test, aes(x = as_factor(Condition), y = Intensity))+
    
    geom_boxplot()

## peptide table boxplot

test <- peptide_table %>% dplyr::select(-c(Length)) %>% pivot_longer(-c(Sequence, Protein, Missed_Cleavages),
                                                                     names_to = "Experiment", values_to = "Intensity")

test <- test %>% mutate(Condition = condition[Experiment]) %>%  mutate(Replicate = replicate[Experiment])
test <- test %>% filter(Protein %in% de_protein_IDs)

ggplot(test, aes(x = as_factor(Condition), y = Intensity))+
    geom_boxplot()+
        facet_wrap(Protein~.)+
            theme_linedraw()

## contrasts in protein table
test <- protein_table %>% group_by(Protein, Condition) %>% summarise(Intensity = median(Intensity, na.rm = TRUE))
test <- test %>% ungroup()        
test <- test %>% filter((Protein %in% de_protein_IDs))
test <- test %>% pivot_wider(names_from = Condition, values_from = Intensity)
test <- test %>% mutate(Contrasts = `2`-`1`)
hist(test$Contrasts)

## compare abundance to intensity
test <- protein_table %>% group_by(Protein, Condition) %>% summarise(Intensity = median(Intensity, na.rm = TRUE))
test <- test %>% ungroup() %>% pivot_wider(names_from = Condition, values_from = Intensity)       

test <- bind_cols(test, copy_number_matrix[test$Protein,])
plot(test$`2`,test$condition_2)
