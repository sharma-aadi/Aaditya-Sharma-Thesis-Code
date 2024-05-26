#######################################################################################################################
## True positives and False Positives
#######################################################################################################################
detection_plot_data <- differential_expression_detection %>% pivot_longer(c(TruePositives, FalsePositives),
                                                                          names_to = "Statistic", values_to = "Count")
detection_plot_data$Statistic <- factor(detection_plot_data$Statistic,
                                        levels = rev(levels(factor(detection_plot_data$Statistic))))

## Numbers
detection_plot_data_numbers <- detection_plot_data %>% group_by(Dataset, Statistic) %>% 
    summarise(Median = median(Count)) %>% ungroup()
detection_plot_data_numbers <- detection_plot_data_numbers %>% arrange(Statistic, Dataset)

## Labels
data_set_labels <- c("DDA-IN-RS", "DDA-IN-MP",
                     "DDA-IQ-RS", "DDA-IQ-MP",
                     "DDA-IM-RS", "DDA-IM-MP",
                     "DIA-IN-RS", "DIA-IN-MP",
                     "DIA-IQ-RS", "DIA-IQ-MP",
                     "DIA-IM-RS", "DIA-IM-MP")
names(data_set_labels) <-  c("DDA_IN_RS", "DDA_IN_MP",
                             "DDA_IQ_RS", "DDA_IQ_MP",
                             "DDA_IM_RS", "DDA_IM_MP",
                             "DIA_IN_RS", "DIA_IN_MP",
                             "DIA_IQ_RS", "DIA_IQ_MP",
                             "DIA_IM_RS", "DIA_IM_MP")

statistics_labels <- c("True Positives", "False Positives")
names(statistics_labels) <- c("TruePositives", "FalsePositives")

##
plot_detection <- ggplot(detection_plot_data, aes(y = Count))+
                    geom_boxplot()+
                        facet_grid(Statistic~Dataset, scales = "free",
                                   labeller = labeller(Dataset = data_set_labels, 
                                                       Statistic = statistics_labels))+
                            theme_linedraw()+
                                theme(panel.grid.minor.x=element_blank(),
                                      panel.grid.minor.y=element_blank(),
                                      panel.grid.major.x=element_blank(),
                                      panel.grid.major.y=element_blank(),
                                      strip.text = element_text(color = "black"),
                                      strip.background = element_rect(color = "black", fill = "white"),
                                      axis.text.x = element_blank(),
                                      axis.ticks.x = element_blank(),
                                      axis.title = element_blank())

plot_detection <- plot_detection + geom_text(data = detection_plot_data_numbers,
                                             aes(label = as.character(Median), y=-5, x=0))

### Combined plot
ggsave(plot = plot_detection,
       filename = paste0("figure_de_detection.pdf"),
       path = paste0("./Figures/"),
       width = 11.69, height = 4, units = "in",
       dpi = 320)

#######################################################################################################################
## FC estimates
#######################################################################################################################

## Labels
data_set_labels_2 <- c("IN-RS", "IN-MP",
                     "IQ-RS", "IQ-MP",
                     "IM-RS", "IM-MP",
                     "IN-RS", "IN-MP",
                     "IQ-RS", "IQ-MP",
                     "IM-RS", "IM-MP")
names(data_set_labels_2) <-  c("DDA_IN_RS", "DDA_IN_MP",
                             "DDA_IQ_RS", "DDA_IQ_MP",
                             "DDA_IM_RS", "DDA_IM_MP",
                             "DIA_IN_RS", "DIA_IN_MP",
                             "DIA_IQ_RS", "DIA_IQ_MP",
                             "DIA_IM_RS", "DIA_IM_MP")


positive_FC <- unique((differential_expression_estimates %>% filter(TrueFC>0))$TrueFC)
positive_FC_labels <- paste0("log2FC = ",round(positive_FC,2))
names(positive_FC_labels) <- positive_FC

negative_FC <- unique((differential_expression_estimates %>% filter(TrueFC<0))$TrueFC)
negative_FC_labels <- paste0("log2FC = ",round(negative_FC,2))
names(negative_FC_labels) <- negative_FC

de_dda <- ggplot(data = differential_expression_estimates %>% filter(Dataset %in% names(data_set_labels)[1:6]),
       aes(x = TrueFC, y = EstimatedFC, color = Dataset))+
    geom_point(shape = 1, stroke = 0.75)+
    scale_color_manual(labels =  data_set_labels_2[1:6], values = c("#332288", "#117733", "#88CCEE",
                                                                  "#DDCC77", "#CC6677","#882255")) +
    guides(colour = guide_legend(nrow = 1, override.aes = list(size = 2)))+
    # labs(title = expression(Predicted~log[2]~fold~change~against~actual~log[2]~fold~change~"(DEP2 pipeline)"))+
    geom_abline(intercept = 0, slope = 1)+
    labs(x = expression(True~log[2]~fold~change), y = expression(Estimated~log[2]~fold~change), color = "Workflow")+
    scale_x_continuous(breaks = c(-5,0,5))+
    # geom_errorbar( width=.2) +
    facet_wrap(.~Protein, ncol = 13)+ theme_linedraw()+
    theme(legend.position = "top",
          legend.justification = "right",
          plot.margin = margin(5,5,5,5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "white"),
          strip.text = element_text(color = "black"))

de_dia <- ggplot(data = differential_expression_estimates %>% filter(Dataset %in% names(data_set_labels)[7:12]),
                 aes(x = TrueFC, y = EstimatedFC, color = Dataset))+
    geom_point(shape = 1, stroke = 0.75)+
    scale_color_manual(labels =  data_set_labels_2[7:12], values = c("#332288", "#117733", "#88CCEE",
                                                                  "#DDCC77", "#CC6677","#882255")) +
    guides(colour = guide_legend(nrow = 1, override.aes = list(size = 2)))+
    # labs(title = expression(Predicted~log[2]~fold~change~against~actual~log[2]~fold~change~"(DEP2 pipeline)"))+
    geom_abline(intercept = 0, slope = 1)+
    labs(x = expression(True~log[2]~fold~change), y = expression(Estimated~log[2]~fold~change), color = "Workflow")+
    scale_x_continuous(breaks = c(-5,0,5))+
    # geom_errorbar( width=.2) +
    facet_wrap(.~Protein, ncol = 13)+ theme_linedraw()+
    theme(legend.position = "top",
          legend.justification = "right",
          plot.margin = margin(5,5,5,5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "white"),
          strip.text = element_text(color = "black"))

### Combined plot
ggsave(plot = ggarrange(de_dda, de_dia, nrow = 2, labels = c("DDA", "DIA"), common.legend = TRUE, vjust = 0),
       filename = paste0("figure_de_estimates.pdf"),
       path = paste0("./Figures/"),
       width = 11.69, height = 6, units = "in",
       dpi = 320)
