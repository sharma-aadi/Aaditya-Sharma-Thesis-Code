test <- peptide_table_dda %>% filter(Unique == "Yes")
test <- test %>% filter(Leading.razor.protein %in% protein_uniprot_id)
test <- test %>% dplyr::select(-c(Proteins, Unique))
test <- test %>% dplyr::rename(Protein =Leading.razor.protein)

test <- test %>% pivot_longer(-c(Sequence, Protein), names_to = "Experiment", values_to = "Intensity")

experiment_labels <-  c(paste0("Platelet_1_Mix_", 1:13),
                        paste0("Platelet_2_Mix_", 1:13),
                        paste0("Platelet_3_Mix_", 1:13))

platelet_id <- c(rep(1,13),
                 rep(2,13),
                 rep(3,13))

mix_id <- c(rep(1:13,3))

names(platelet_id) <- experiment_labels
names(mix_id) <- experiment_labels

test <- test %>% mutate(Mix = mix_id[Experiment], Donor = platelet_id[Experiment])

test <- test %>%  mutate(Intensity = ifelse(Intensity == 0, NA, Intensity))

test <- test %>% mutate(Intensity = log2(Intensity))

# ggplot(data = test, aes(x = as_factor(Mix), y = Intensity))+
#     geom_boxplot()+
#     facet_grid(Donor~.)+
#     theme_linedraw()

test <- test %>% group_by(Protein, Experiment, Mix, Donor) %>%
            summarise(MedianIntensity = median(Intensity, na.rm = T)) %>% ungroup()

test <- test %>%  mutate(MedianIntensity = ifelse(MedianIntensity == 0, NA, MedianIntensity))

## Add concentrations
concentration <- numeric(length = length(test$Protein))

for(j in 1:length(concentration)){
    
    concentration[j] <- latin_square_matrix[test$Mix[j], test$Protein[j]]
    
}

test <- test %>% mutate(Concentration = log2(concentration))

test <- test %>% mutate(Protein = protein_letters[Protein])
test <- test %>% filter(Protein != "Median")

test$Protein <- factor(test$Protein, levels = rev(LETTERS))

table(is.na(test$MedianIntensity))

test <- test %>% group_by(Protein, Mix, Concentration) %>% summarise(Intensity = median(MedianIntensity)) %>% ungroup()

ggplot(data = test, aes(x = Concentration, y = Protein, fill = Intensity, colour = "black"))+
    geom_tile()

## 7. Differential test
diff_pep <- test_diff(pe_dda_IN_RS_se, type = "control", control = "Mix1")
dep_pep <- add_rejections(diff_pep, alpha = 0.01)

# plot_pca(dep_pep, x = 1, y = 2, point_size = 3)

### volcano plot on contrast 

plot_volcano(dep_pep, contrast = "Mix7_vs_Mix1", adjusted = F, chooseTolabel = LETTERS)

test <- tibble(get_df_wide(dep_pep))
test$name[which(test$name %in% protein_uniprot_id)]<-protein_letters[test$name[which(test$name%in%protein_uniprot_id)]]
test <- (test %>% filter((name %in% protein_letters)))
test$name <- factor(test$name,
                    levels = c(LETTERS, "Median"))
test <- test %>% arrange(name)
test <- test %>% dplyr::select(c(Mix1_vs_Mix7_significant:Mix9_vs_Mix7_significant))
table(unlist(test))


pe_dda_IN_RS_TP <- pe_dda_IN_RS_TP %>% pivot_longer(c(TruePositives, FalsePositives),
                                                    names_to = "Statistic", values_to = "Count")

ggplot(data = pe_dda_IN_RS_TP, aes( y = Count))+
    facet_grid(.~Statistic)+
    geom_boxplot()
