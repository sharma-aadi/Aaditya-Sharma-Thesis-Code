plot_2_scatterplot_color_alpha <- 0.1
plot_2_scatterplot_axis_title <- 9
plot_2_axis_text_size <- 7
plot_legend_title_size <- 8
plot_2_scatterplot_point_stroke <- 0.25
plot_2_scatterplot_expand_y <- 2.5
plot_2_scatterplot_color_alpha <- 0.5

#######################################################################################################################
## For missing value plot
#######################################################################################################################

data_1b <- spike_in_protein_peptide_counts %>% filter(Protein != "P02769")

data_1b$Protein <- factor(data_1b$Protein,
                          levels = spike_in_protein_statistics_with_mapping$Protein)

data_1b <- data_1b %>% arrange(Protein)

data_1b <- data_1b %>% mutate(Length = spike_in_protein_statistics_with_mapping$Length) %>% 
    mutate(Category = spike_in_protein_statistics_with_mapping$Mapping)

#######################################################################################################################
## For missing value plot
#######################################################################################################################

## protein labels
protein_labels_4 <- paste0(paste0(paste0(spike_in_protein_statistics$ProteinID)))

names(protein_labels_4) <- spike_in_protein_statistics$Protein

## Add a row of white cells for a blank row after Latin-square proteins
data_1d <- bind_rows(latin_square_long %>% mutate(Color = "black") %>% 
                         mutate(Type = ifelse(Protein == "O88766", "Extra", "LS")),
                     tibble(Sample = rep(paste0("Sample_",1:13), 2),
                            SampleID = rep(1:13, 2),
                            Protein = c(rep("Blank", 13), rep("P02769", 13)),
                            CopyNumber = c(rep(NA, 13),
                                           c(82381.69805,83376.1623,81335.92951,74857.92549,63580.60129,55782.46603,
                                             40838.75578,33855.94166,21759.16917,334.214314,2506.17027,4173.098799,
                                             18960.24308)),
                            MixID = rep(1:13, 2),
                            Color = c(rep("white", 13), rep("black", 13)),
                            Type = c("Extra")))

## Modify protein labels 

protein_labels_1_d <- c(protein_labels_4[1:26], expression(gamma), expression(beta))
names(protein_labels_1_d)[27] <- "O88766"
names(protein_labels_1_d)[28] <- "P02769"

## Re-factor protein levels
data_1d$Protein <- factor(data_1d$Protein,
                          levels = c("P02769",
                                     levels(factor(latin_square_long$Protein))[1],
                                     "Blank",
                                     levels(factor(latin_square_long$Protein))[-1]))
## Factor facets
data_1d$Type <- factor(data_1d$Type,
                       levels = c("LS", "Extra"))

## Flip dataset
data_1d_flipped <- data_1d
data_1d_flipped$SampleID <- factor(data_1d_flipped$SampleID,
                                   levels = rev(levels(factor(data_1d$SampleID))))
data_1d_flipped$Protein <- factor(data_1d_flipped$Protein,
                                  levels = rev(levels(factor(data_1d$Protein))))

#######################################################################################################################
## protein-level missing values
#######################################################################################################################

## Load protein intensity data table
data_2_top <- protein_intensity_data 

median_data <- data_2_top %>% filter(Protein == "O88766") %>% 
    group_by(Protein, Dataset) %>% summarise(MedianIntensity = median(Intensity)) %>% ungroup()

median_gamma_dda <- (median_data %>% filter(Dataset == "DDA"))$MedianIntensity
median_gamma_dia <- (median_data %>% filter(Dataset == "DIA"))$MedianIntensity


data_2_top_normalized <- bind_rows((data_2_top %>% filter(Dataset == "DDA") %>% 
                                        mutate(Intensity = Intensity/median_gamma_dda)),
                                   (data_2_top %>% filter(Dataset == "DIA") %>% 
                                        mutate(Intensity = Intensity/median_gamma_dia)))

## re-factor levels
data_2_top$Protein <- factor(data_2_top$Protein, levels = levels(factor(data_1d$Protein)))
data_2_top$Dataset <- factor(data_2_top$Dataset, levels = c("DDA", "DIA"))

data_2_top_normalized$Protein <- factor(data_2_top_normalized$Protein, levels = levels(factor(data_1d$Protein)))
data_2_top_normalized$Dataset <- factor(data_2_top_normalized$Dataset, levels = c("DDA", "DIA"))

## Missing values tile plot
data_2_missing_values <- data_2_top %>% filter(Protein != "O88766") %>% 
    dplyr::select(Protein,Dataset, Platelet, Concentration, Intensity)

data_2_missing_values <- data_2_missing_values %>% mutate(TypeCode = ifelse(Dataset == "DDA", 1, 2)) 

data_2_missing_values$TypeCode[which(is.na(data_2_missing_values$Intensity))] <- 0

data_2_missing_values <- data_2_missing_values %>% group_by(Protein, Concentration, Platelet) %>% 
    summarise(TypeCode = sum(TypeCode)) %>% ungroup()

type_category_labels <- c("Neither",
                          "DDA",
                          "DIA",
                          "Both")

names(type_category_labels) <- c(0,1,2,3)

data_2_missing_values$Protein <- factor(data_2_missing_values$Protein, 
                                        levels = levels(factor(data_2_top$Protein))[-1])

protein_category_labels <- data_1b %>% filter(Protein != "O88766") %>%  arrange(desc(Protein))

data_2_missing_values <- data_2_missing_values %>% arrange(Protein, Platelet, Concentration) %>%
    mutate(Category = c(sapply(protein_category_labels$Category, function(x) rep(x,39), simplify = TRUE)))

data_2_missing_values <- data_2_missing_values %>% 
    mutate(ID = c(sapply(rev(LETTERS), function(x) rep(x,39), simplify = TRUE)))

mapping_labels_new <- c("No sharing", "Another spike-in", "Platelet", "Both")
names(mapping_labels_new) <- 1:4

data_2_missing_values$Platelet <- factor(data_2_missing_values$Platelet,
                                         levels = rev(levels(factor(data_2_missing_values$Platelet))))

#######################################################################################################################
## Protein intensity scatterplot
#######################################################################################################################

data_2_bottom_normalized <- data_2_top_normalized %>% filter(Protein != "O88766")

data_2_bottom_normalized$Protein <- factor(data_2_bottom_normalized$Protein,
                                           levels=rev(levels(factor(data_2_bottom_normalized$Protein))))

concentration_level_tibble <- tibble(Concentration = (levels(factor(data_2_top$Concentration)))[1:13],
                                     Level = c(rep("Low", 4),rep("Medium", 5),rep("High", 4)))

concentration_level_labels <- character(length = length(data_2_bottom_normalized$Concentration))

for(i in 1:length(concentration_level_labels)){
    concentration_level_labels[i] <- concentration_level_tibble$Level[
        which(concentration_level_tibble$Concentration == 
                  as.character(data_2_bottom_normalized$Concentration[i]))]
}

data_2_bottom_normalized <- data_2_bottom_normalized %>% mutate(ConcentrationLevel = concentration_level_labels)
data_2_bottom_normalized <- data_2_bottom_normalized%>%
    mutate(TypeConcentration = paste0(Dataset," _ ", ConcentrationLevel))

#######################################################################################################################
## For protein intensity scatterplot
#######################################################################################################################

protein_labels_3_c <- protein_labels_4[-c(27:28)]

## Create phantom plot for legend
plot_2_scatterplot_legend <- ggplot()+
    geom_point(data = data_2_bottom_normalized,
               size = 1.5, stroke = 0.5,
               mapping = aes(x = as_factor(Concentration),
                             y = Intensity, shape = as_factor(Platelet), color = ConcentrationLevel),
               inherit.aes = FALSE)+
    scale_shape_manual(values = c(0:2))+
    scale_x_discrete(breaks = (levels(factor(data_2_top$Concentration)))[1:13], labels = 1:13,
                     limits = (levels(factor(data_2_top$Concentration)))[1:13])+
    facet_wrap(.~TypeConcentration, scales = "free", labeller = labeller(Protein = protein_labels_3_c))+
    guides(color=guide_legend(override.aes = list(size=2, shape = 22, color = "black", stroke = 0.5,
                                                  fill = c(alpha("red", plot_2_scatterplot_color_alpha),
                                                           alpha("blue", plot_2_scatterplot_color_alpha),
                                                           alpha("green", plot_2_scatterplot_color_alpha)))),
           shape=guide_legend(override.aes = list(size=2)))+
    theme_linedraw()+
    scale_color_manual(breaks=c("Low","Medium","High"), values = c(alpha("red", plot_2_scatterplot_color_alpha),
                                                                   alpha("blue", plot_2_scatterplot_color_alpha),
                                                                   alpha("green", plot_2_scatterplot_color_alpha)))+
    labs(x = "Concentration level", y = "Measured intensity",
         color = "Concentration level group", shape = "Donor")+
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
        # panel.background = element_rect(fill = alpha("green",0.5),
        #                                 colour = "black", 
        #                                 size = 0.2, linetype = "solid"),
        legend.position = "top",
        legend.justification = "center",
        legend.box="vertical",
        legend.spacing.x = unit(0, 'pt'),
        legend.spacing.y = unit(-5, 'pt'),
        legend.margin = margin(0,0,0,0),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        plot.margin =  margin(t = 5, b = 0, r = 5, l = 18, unit = "pt")
    ) 
