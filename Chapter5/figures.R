#######################################################################################################################
## PDF Figures
#######################################################################################################################

gamma_support <- seq(0, 7, length.out = 10000)
plot(gamma_support, dgamma(gamma_support,protein_shape,protein_rate), type = "l", xlab = "x",
     ylab = "f(x)")

beta_support <- seq(0, 1, length.out = 1000)
plot(beta_support, dbeta(beta_support,10,2), type = "l", xlab = "x",
     ylab = "f(x)")

dgamma(x, shape, rate = 1, scale = 1/rate, log = FALSE)

beta_support <- seq(0, 1, length.out = 1000)
plot(beta_support, dbeta(beta_support,3,16)/range(dbeta(beta_support,3,16))[2], type = "l", xlab = "x",
     ylab = "f(x)/max(f(x))")

cdf_support <- seq(-3, 3, length.out = 10000)
plot(cdf_support, pnorm(cdf_support,0,1), type = "l", xlab = "x", ylab = expression(paste(Phi,"(x)")))

#######################################################################################################################
## PCA
#######################################################################################################################

## PCA
pca_data <- tibble(rownames_to_column(data.frame(t(peptide_abundance_matrix)), "Experiment"))
pca_data <- pca_data %>% mutate(Condition = as.factor(c(1,1,1,2,2,2)))

peptides.pca <- prcomp(pca_data[,c(2:3203)], center = T, scale. = T)

pca_plot <- autoplot(peptides.pca, data = pca_data, colour = "Experiment", size = 3, shape = "Condition")+
    guides(color=guide_legend(title = "LC-MS/MS run:", byrow = T, title.position = "top"))+
    guides(shape=guide_legend(title = "Condition:", title.position = "top"))+
    scale_color_manual(values = c('#000000', "#E69F00", "#56B4E9", "#0072B2", "#D55E00", "#CC79A7"))+
theme_linedraw()+
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top"
    )

ggsave(plot = pca_plot,
       filename = paste0("simulation_pca_plot.pdf"),
       path = paste0("./"),
       width = 6, height = 6, units = "in",
       dpi = 320)

#######################################################################################################################
## Peptide abundance data
#######################################################################################################################

peptide_intensity_boxplot_data <- peptide_abundance_data_detectable %>% filter(Protein %in% de_protein_IDs) %>% 
                                    dplyr::select(c(Protein, condition_1, condition_2)) %>% 
                                        pivot_longer(-c(Protein), names_to = "Condition", values_to = "Intensity")

simulation_box_plot <- ggplot(data = peptide_intensity_boxplot_data, aes(x = as_factor(Condition), y = Intensity))+
    facet_wrap(.~Protein, nrow = 4)+
    geom_boxplot()+
    theme_linedraw() +
    scale_x_discrete(labels = c("Control", "DE"))+
    labs(x = "Condition", y = bquote(log[2]-transformed~intensity))+
    theme(
        plot.title = element_text(size=20, face = "bold"),
        plot.subtitle = element_text(size=18, face = "bold"),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12, angle = 0),
        strip.text = element_text(size = 14, color = "black"),
        strip.background = element_rect(fill = "white", color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(plot = simulation_box_plot,
       filename = paste0("simulation_box_plot.pdf"),
       path = paste0("./"),
       width = 11.69, height = 8, units = "in",
       dpi = 320)
