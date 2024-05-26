#######################################################################################################################
## Prototype plots (DIA-NN data) [February 25, 2023] 
#######################################################################################################################

ggplot(data = de_diann %>% mutate(TrueFC = round(TrueFC, digits = 2)),
       aes(x = EstimatedFC, y = -log10(P), label = Protein)) +
    facet_wrap(.~TrueFC, nrow = 2)+
    geom_vline(data = de_diann %>% mutate(TrueFC = round(TrueFC, digits = 2)),  aes(xintercept = TrueFC))+
    geom_text_repel()+
    geom_point()+
    labs(title = "Volcano plot for DIA fold-changes calculated using DEP")+
    theme_linedraw()

ggplot(data = de_dda %>% mutate(TrueFC = round(TrueFC, digits = 2)),
       aes(x = EstimatedFC, y = -log10(P), label = Protein)) +
    facet_wrap(.~TrueFC, nrow = 2)+
    geom_vline(data = de_diann %>% mutate(TrueFC = round(TrueFC, digits = 2)),  aes(xintercept = TrueFC))+
    geom_text_repel()+
    geom_point()+
    labs(title = "Volcano plot for DDA fold-changes calculated using DEP")+
    theme_linedraw()

de_all <- bind_rows(de_dda %>% mutate(Dataset = "DDA"),
                    de_diann %>% mutate(Dataset = "DIA"))

## Scatterplot of differential expression estimates
plot_1_e_alpha <- 1
data_1e_colors <- c(hue_pal()(2)[2], hue_pal()(2)[1], "#7570B3", "#E6AB02",
                    "white",
                    alpha(hue_pal()(2)[2], plot_1_e_alpha),
                    alpha(hue_pal()(2)[1], plot_1_e_alpha),
                    alpha("#7570B3", plot_1_e_alpha),
                    alpha("#E6AB02", plot_1_e_alpha))

fc_plot_colors <- c(alpha(data_1e_colors[c(3)], 0.75),
                    alpha(data_1e_colors[c(4)], 0.75))

names(fc_plot_colors) <- c("DDA", "DIA")

fc_plot <- ggplot(data = de_all,
       aes(x = TrueFC, y = EstimatedFC, color = Dataset))+
    geom_point(shape = 1, stroke = 1)+
    scale_color_manual(values = fc_plot_colors)+
    labs(title = expression(Predicted~log[2]~fold-change~against~actual~log[2]~fold-change~"(DEP pipeline)"))+
    geom_abline(intercept = 0, slope = 1)+
    # geom_errorbar( width=.2) +
    facet_wrap(.~Protein, ncol = 13)+ theme_linedraw()+
    theme(legend.position = "bottom",
          legend.justification = "right",
          legend.margin = margin(t=-20,b=30),
          plot.margin = margin(5,5,5,5),
          axis.text.x = element_text(angle = 90),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "white"),
          strip.text = element_text(color = "black"))

ggsave(plot = fc_plot,
       filename = paste0("FC_plot.pdf"),
       path = paste0("./Figures/"),
       width = 11.69, height = 6, units = "in",
       dpi = 320)
