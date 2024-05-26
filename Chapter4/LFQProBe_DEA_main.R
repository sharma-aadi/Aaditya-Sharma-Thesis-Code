#######################################################################################################################
## create all 12 datasets for differential expression analysis
#######################################################################################################################

# source(file = "LFQProBe_descriptive_data_analysis.R")
# source(file = "LFQProBe_experimental_design.R")
# source(file = "LFQProbe_dataset_creation.R")

#######################################################################################################################
## Alternatively load all 12 datasets from saved image
#######################################################################################################################

# load("LFQProBe_DEA_datasets.RData")

#######################################################################################################################
## DE analysis - pe_dda_IN_RS_se
#######################################################################################################################

# Specify current dataset
current_dataset <- pe_dda_IN_RS_se

# Intialize tibble to store all the data
pe_dda_IN_RS_de <- tibble() ## This stores point estimates of the FC
pe_dda_IN_RS_TP <- tibble() ## This is to store the number of true positives and false positives

for(i in 1:13){
    
    current_pair <- i
    
    current_proteins <- spike_in_design_median_mix$ProteinID[which(spike_in_design_median_mix$PairID == current_pair)]
    current_median<-unique(spike_in_design_median_mix$Median[which(spike_in_design_median_mix$PairID == current_pair)])
    
    diff_pep <- test_diff(current_dataset, type = "control", control = paste0("Mix",current_median))
    diff_pep <- add_rejections(diff_pep, alpha = 0.01)
    diff_pep <- get_df_wide(diff_pep)
    
    ## DE detection columns
    de_detection_columns <- c(paste0(paste0("Mix",c(1:13)[-current_median]),"_vs_Mix",current_median))
    
    ## Find true positives and false positives in the current contrasts ()
    true_positives_table <- diff_pep %>% filter(name %in% LETTERS) %>% 
                                dplyr::select(paste0(de_detection_columns,"_significant"))
    false_positives_table <- diff_pep %>% filter(!(name %in% LETTERS)) %>% 
                                dplyr::select(paste0(de_detection_columns,"_significant"))
    
    number_of_true_positives <- colSums(true_positives_table == TRUE, na.rm = TRUE)
    number_of_false_positives <- colSums(false_positives_table == TRUE, na.rm = TRUE)
    
    tibble_of_positives <- tibble(TruePositives = number_of_true_positives, FalsePositives = number_of_false_positives)
    
    ## Which contrasts to find the FC estimate for the current proteins
    lower_contrasts <- c(paste0(paste0("Mix",lower_samples[current_median,]),"_vs_Mix",current_median))
    higher_contrasts <- c(paste0(paste0("Mix",higher_samples[current_median,]),"_vs_Mix",current_median))
    
    ## Keep only the current proteins
    diff_pep <- diff_pep %>% filter(name %in% current_proteins)
    
    ## Point estimate
    lower_tibble_point <- diff_pep %>% dplyr::select(c(name, paste0(lower_contrasts,"_diff"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "FC")
    
    higher_tibble_point <- diff_pep %>% dplyr::select(c(name, paste0(higher_contrasts,"_diff"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "FC")
    
    ## Left CI
    lower_tibble_left <- diff_pep %>% dplyr::select(c(name, paste0(lower_contrasts,"_CI.L"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Left")
    
    higher_tibble_left <- diff_pep %>% dplyr::select(c(name, paste0(higher_contrasts,"_CI.L"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Left")
    
    ## Right CI
    lower_tibble_right <- diff_pep %>% dplyr::select(c(name, paste0(lower_contrasts,"_CI.R"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Right")
    
    higher_tibble_right <- diff_pep %>% dplyr::select(c(name, paste0(higher_contrasts,"_CI.R"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Right")
    
    ## Adjusted p value
    lower_tibble_p <- diff_pep %>% dplyr::select(c(name, paste0(lower_contrasts,"_p.adj"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "adjusted_p_value")
    
    higher_tibble_p <- diff_pep %>% dplyr::select(c(name, paste0(higher_contrasts,"_p.adj"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "adjusted_p_value")
    
    ## Create FC tibble
    fc_tibble_current_lower <- tibble(Protein = lower_tibble_point$name,
                                      TrueFC = rep(log2(1.75^-(1:6)),2),
                                      EstimatedFC_Left = lower_tibble_left$Left,
                                      EstimatedFC = lower_tibble_point$FC,
                                      EstimatedFC_Right = lower_tibble_right$Right,
                                      P = lower_tibble_p$adjusted_p_value)
    
    fc_tibble_current_higher <- tibble(Protein = higher_tibble_point$name,
                                       TrueFC = rep(log2(1.75^(1:6)),2),
                                       EstimatedFC_Left = higher_tibble_left$Left,
                                       EstimatedFC = higher_tibble_point$FC,
                                       EstimatedFC_Right = higher_tibble_right$Right,
                                       P = higher_tibble_p$adjusted_p_value) 
    
    fc_tibble_current <- bind_rows(fc_tibble_current_lower,fc_tibble_current_higher)
    
    ## Store the information about the FC estimates and the detection of True and False Positives
    pe_dda_IN_RS_TP <- bind_rows(pe_dda_IN_RS_TP, tibble_of_positives)
    pe_dda_IN_RS_de <- bind_rows(pe_dda_IN_RS_de, fc_tibble_current)
}

#######################################################################################################################
## DE analysis - pe_dia_IN_RS_se
#######################################################################################################################

# Specify current dataset
current_dataset <- pe_dia_IN_RS_se

# Intialize tibble to store all the data
pe_dia_IN_RS_de <- tibble() ## This stores point estimates of the FC
pe_dia_IN_RS_TP <- tibble() ## This is to store the number of true positives and false positives

for(i in 1:13){
    
    current_pair <- i
    
    current_proteins <- spike_in_design_median_mix$ProteinID[which(spike_in_design_median_mix$PairID == current_pair)]
    current_median<-unique(spike_in_design_median_mix$Median[which(spike_in_design_median_mix$PairID == current_pair)])
    
    diff_pep <- test_diff(current_dataset, type = "control", control = paste0("Mix",current_median))
    diff_pep <- add_rejections(diff_pep, alpha = 0.01)
    diff_pep <- get_df_wide(diff_pep)
    
    ## DE detection columns
    de_detection_columns <- c(paste0(paste0("Mix",c(1:13)[-current_median]),"_vs_Mix",current_median))
    
    ## Find true positives and false positives in the current contrasts ()
    true_positives_table <- diff_pep %>% filter(name %in% LETTERS) %>% 
        dplyr::select(paste0(de_detection_columns,"_significant"))
    false_positives_table <- diff_pep %>% filter(!(name %in% LETTERS)) %>% 
        dplyr::select(paste0(de_detection_columns,"_significant"))
    
    number_of_true_positives <- colSums(true_positives_table == TRUE, na.rm = TRUE)
    number_of_false_positives <- colSums(false_positives_table == TRUE, na.rm = TRUE)
    
    tibble_of_positives <- tibble(TruePositives = number_of_true_positives, FalsePositives = number_of_false_positives)
    
    ## Which contrasts to find the FC estimate for the current proteins
    lower_contrasts <- c(paste0(paste0("Mix",lower_samples[current_median,]),"_vs_Mix",current_median))
    higher_contrasts <- c(paste0(paste0("Mix",higher_samples[current_median,]),"_vs_Mix",current_median))
    
    ## Keep only the current proteins
    diff_pep <- diff_pep %>% filter(name %in% current_proteins)
    
    ## Point estimate
    lower_tibble_point <- diff_pep %>% dplyr::select(c(name, paste0(lower_contrasts,"_diff"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "FC")
    
    higher_tibble_point <- diff_pep %>% dplyr::select(c(name, paste0(higher_contrasts,"_diff"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "FC")
    
    ## Left CI
    lower_tibble_left <- diff_pep %>% dplyr::select(c(name, paste0(lower_contrasts,"_CI.L"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Left")
    
    higher_tibble_left <- diff_pep %>% dplyr::select(c(name, paste0(higher_contrasts,"_CI.L"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Left")
    
    ## Right CI
    lower_tibble_right <- diff_pep %>% dplyr::select(c(name, paste0(lower_contrasts,"_CI.R"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Right")
    
    higher_tibble_right <- diff_pep %>% dplyr::select(c(name, paste0(higher_contrasts,"_CI.R"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Right")
    
    ## Adjusted p value
    lower_tibble_p <- diff_pep %>% dplyr::select(c(name, paste0(lower_contrasts,"_p.adj"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "adjusted_p_value")
    
    higher_tibble_p <- diff_pep %>% dplyr::select(c(name, paste0(higher_contrasts,"_p.adj"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "adjusted_p_value")
    
    ## Create FC tibble
    fc_tibble_current_lower <- tibble(Protein = lower_tibble_point$name,
                                      TrueFC = rep(log2(1.75^-(1:6)),2),
                                      EstimatedFC_Left = lower_tibble_left$Left,
                                      EstimatedFC = lower_tibble_point$FC,
                                      EstimatedFC_Right = lower_tibble_right$Right,
                                      P = lower_tibble_p$adjusted_p_value)
    
    fc_tibble_current_higher <- tibble(Protein = higher_tibble_point$name,
                                       TrueFC = rep(log2(1.75^(1:6)),2),
                                       EstimatedFC_Left = higher_tibble_left$Left,
                                       EstimatedFC = higher_tibble_point$FC,
                                       EstimatedFC_Right = higher_tibble_right$Right,
                                       P = higher_tibble_p$adjusted_p_value) 
    
    fc_tibble_current <- bind_rows(fc_tibble_current_lower,fc_tibble_current_higher)
    
    ## Store the information about the FC estimates and the detection of True and False Positives
    pe_dia_IN_RS_TP <- bind_rows(pe_dia_IN_RS_TP, tibble_of_positives)
    pe_dia_IN_RS_de <- bind_rows(pe_dia_IN_RS_de, fc_tibble_current)
}

#######################################################################################################################
## DE analysis - pe_dda_IQ_RS_se
#######################################################################################################################

# Specify current dataset
current_dataset <- pe_dda_IQ_RS_se

# Intialize tibble to store all the data
pe_dda_IQ_RS_de <- tibble() ## This stores point estimates of the FC
pe_dda_IQ_RS_TP <- tibble() ## This is to store the number of true positives and false positives

for(i in 1:13){
    
    current_pair <- i
    
    current_proteins <- spike_in_design_median_mix$ProteinID[which(spike_in_design_median_mix$PairID == current_pair)]
    current_median<-unique(spike_in_design_median_mix$Median[which(spike_in_design_median_mix$PairID == current_pair)])
    
    diff_pep <- test_diff(current_dataset, type = "control", control = paste0("Mix",current_median))
    diff_pep <- add_rejections(diff_pep, alpha = 0.01)
    diff_pep <- get_df_wide(diff_pep)
    
    ## DE detection columns
    de_detection_columns <- c(paste0(paste0("Mix",c(1:13)[-current_median]),"_vs_Mix",current_median))
    
    ## Find true positives and false positives in the current contrasts ()
    true_positives_table <- diff_pep %>% filter(name %in% LETTERS) %>% 
        dplyr::select(paste0(de_detection_columns,"_significant"))
    false_positives_table <- diff_pep %>% filter(!(name %in% LETTERS)) %>% 
        dplyr::select(paste0(de_detection_columns,"_significant"))
    
    number_of_true_positives <- colSums(true_positives_table == TRUE, na.rm = TRUE)
    number_of_false_positives <- colSums(false_positives_table == TRUE, na.rm = TRUE)
    
    tibble_of_positives <- tibble(TruePositives = number_of_true_positives, FalsePositives = number_of_false_positives)
    
    ## Which contrasts to find the FC estimate for the current proteins
    lower_contrasts <- c(paste0(paste0("Mix",lower_samples[current_median,]),"_vs_Mix",current_median))
    higher_contrasts <- c(paste0(paste0("Mix",higher_samples[current_median,]),"_vs_Mix",current_median))
    
    ## Keep only the current proteins
    diff_pep <- diff_pep %>% filter(name %in% current_proteins)
    
    ## Point estimate
    lower_tibble_point <- diff_pep %>% dplyr::select(c(name, paste0(lower_contrasts,"_diff"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "FC")
    
    higher_tibble_point <- diff_pep %>% dplyr::select(c(name, paste0(higher_contrasts,"_diff"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "FC")
    
    ## Left CI
    lower_tibble_left <- diff_pep %>% dplyr::select(c(name, paste0(lower_contrasts,"_CI.L"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Left")
    
    higher_tibble_left <- diff_pep %>% dplyr::select(c(name, paste0(higher_contrasts,"_CI.L"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Left")
    
    ## Right CI
    lower_tibble_right <- diff_pep %>% dplyr::select(c(name, paste0(lower_contrasts,"_CI.R"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Right")
    
    higher_tibble_right <- diff_pep %>% dplyr::select(c(name, paste0(higher_contrasts,"_CI.R"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Right")
    
    ## Adjusted p value
    lower_tibble_p <- diff_pep %>% dplyr::select(c(name, paste0(lower_contrasts,"_p.adj"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "adjusted_p_value")
    
    higher_tibble_p <- diff_pep %>% dplyr::select(c(name, paste0(higher_contrasts,"_p.adj"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "adjusted_p_value")
    
    ## Create FC tibble
    fc_tibble_current_lower <- tibble(Protein = lower_tibble_point$name,
                                      TrueFC = rep(log2(1.75^-(1:6)),2),
                                      EstimatedFC_Left = lower_tibble_left$Left,
                                      EstimatedFC = lower_tibble_point$FC,
                                      EstimatedFC_Right = lower_tibble_right$Right,
                                      P = lower_tibble_p$adjusted_p_value)
    
    fc_tibble_current_higher <- tibble(Protein = higher_tibble_point$name,
                                       TrueFC = rep(log2(1.75^(1:6)),2),
                                       EstimatedFC_Left = higher_tibble_left$Left,
                                       EstimatedFC = higher_tibble_point$FC,
                                       EstimatedFC_Right = higher_tibble_right$Right,
                                       P = higher_tibble_p$adjusted_p_value) 
    
    fc_tibble_current <- bind_rows(fc_tibble_current_lower,fc_tibble_current_higher)
    
    ## Store the information about the FC estimates and the detection of True and False Positives
    pe_dda_IQ_RS_TP <- bind_rows(pe_dda_IQ_RS_TP, tibble_of_positives)
    pe_dda_IQ_RS_de <- bind_rows(pe_dda_IQ_RS_de, fc_tibble_current)
}

#######################################################################################################################
## DE analysis - pe_dia_IQ_RS_se
#######################################################################################################################

# Specify current dataset
current_dataset <- pe_dia_IQ_RS_se

# Intialize tibble to store all the data
pe_dia_IQ_RS_de <- tibble() ## This stores point estimates of the FC
pe_dia_IQ_RS_TP <- tibble() ## This is to store the number of true positives and false positives

for(i in 1:13){
    
    current_pair <- i
    
    current_proteins <- spike_in_design_median_mix$ProteinID[which(spike_in_design_median_mix$PairID == current_pair)]
    current_median<-unique(spike_in_design_median_mix$Median[which(spike_in_design_median_mix$PairID == current_pair)])
    
    diff_pep <- test_diff(current_dataset, type = "control", control = paste0("Mix",current_median))
    diff_pep <- add_rejections(diff_pep, alpha = 0.01)
    diff_pep <- get_df_wide(diff_pep)
    
    ## DE detection columns
    de_detection_columns <- c(paste0(paste0("Mix",c(1:13)[-current_median]),"_vs_Mix",current_median))
    
    ## Find true positives and false positives in the current contrasts ()
    true_positives_table <- diff_pep %>% filter(name %in% LETTERS) %>% 
        dplyr::select(paste0(de_detection_columns,"_significant"))
    false_positives_table <- diff_pep %>% filter(!(name %in% LETTERS)) %>% 
        dplyr::select(paste0(de_detection_columns,"_significant"))
    
    number_of_true_positives <- colSums(true_positives_table == TRUE, na.rm = TRUE)
    number_of_false_positives <- colSums(false_positives_table == TRUE, na.rm = TRUE)
    
    tibble_of_positives <- tibble(TruePositives = number_of_true_positives, FalsePositives = number_of_false_positives)
    
    ## Which contrasts to find the FC estimate for the current proteins
    lower_contrasts <- c(paste0(paste0("Mix",lower_samples[current_median,]),"_vs_Mix",current_median))
    higher_contrasts <- c(paste0(paste0("Mix",higher_samples[current_median,]),"_vs_Mix",current_median))
    
    ## Keep only the current proteins
    diff_pep <- diff_pep %>% filter(name %in% current_proteins)
    
    ## Point estimate
    lower_tibble_point <- diff_pep %>% dplyr::select(c(name, paste0(lower_contrasts,"_diff"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "FC")
    
    higher_tibble_point <- diff_pep %>% dplyr::select(c(name, paste0(higher_contrasts,"_diff"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "FC")
    
    ## Left CI
    lower_tibble_left <- diff_pep %>% dplyr::select(c(name, paste0(lower_contrasts,"_CI.L"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Left")
    
    higher_tibble_left <- diff_pep %>% dplyr::select(c(name, paste0(higher_contrasts,"_CI.L"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Left")
    
    ## Right CI
    lower_tibble_right <- diff_pep %>% dplyr::select(c(name, paste0(lower_contrasts,"_CI.R"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Right")
    
    higher_tibble_right <- diff_pep %>% dplyr::select(c(name, paste0(higher_contrasts,"_CI.R"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Right")
    
    ## Adjusted p value
    lower_tibble_p <- diff_pep %>% dplyr::select(c(name, paste0(lower_contrasts,"_p.adj"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "adjusted_p_value")
    
    higher_tibble_p <- diff_pep %>% dplyr::select(c(name, paste0(higher_contrasts,"_p.adj"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "adjusted_p_value")
    
    ## Create FC tibble
    fc_tibble_current_lower <- tibble(Protein = lower_tibble_point$name,
                                      TrueFC = rep(log2(1.75^-(1:6)),2),
                                      EstimatedFC_Left = lower_tibble_left$Left,
                                      EstimatedFC = lower_tibble_point$FC,
                                      EstimatedFC_Right = lower_tibble_right$Right,
                                      P = lower_tibble_p$adjusted_p_value)
    
    fc_tibble_current_higher <- tibble(Protein = higher_tibble_point$name,
                                       TrueFC = rep(log2(1.75^(1:6)),2),
                                       EstimatedFC_Left = higher_tibble_left$Left,
                                       EstimatedFC = higher_tibble_point$FC,
                                       EstimatedFC_Right = higher_tibble_right$Right,
                                       P = higher_tibble_p$adjusted_p_value) 
    
    fc_tibble_current <- bind_rows(fc_tibble_current_lower,fc_tibble_current_higher)
    
    ## Store the information about the FC estimates and the detection of True and False Positives
    pe_dia_IQ_RS_TP <- bind_rows(pe_dia_IQ_RS_TP, tibble_of_positives)
    pe_dia_IQ_RS_de <- bind_rows(pe_dia_IQ_RS_de, fc_tibble_current)
}

#######################################################################################################################
## DE analysis - pe_dda_IM_RS_se
#######################################################################################################################

# Specify current dataset
current_dataset <- pe_dda_IM_RS_se

# Intialize tibble to store all the data
pe_dda_IM_RS_de <- tibble() ## This stores point estimates of the FC
pe_dda_IM_RS_TP <- tibble() ## This is to store the number of true positives and false positives

for(i in 1:13){
    
    current_pair <- i
    
    current_proteins <- spike_in_design_median_mix$ProteinID[which(spike_in_design_median_mix$PairID == current_pair)]
    current_median<-unique(spike_in_design_median_mix$Median[which(spike_in_design_median_mix$PairID == current_pair)])
    
    diff_pep <- test_diff(current_dataset, type = "control", control = paste0("Mix",current_median))
    diff_pep <- add_rejections(diff_pep, alpha = 0.01)
    diff_pep <- get_df_wide(diff_pep)
    
    ## DE detection columns
    de_detection_columns <- c(paste0(paste0("Mix",c(1:13)[-current_median]),"_vs_Mix",current_median))
    
    ## Find true positives and false positives in the current contrasts ()
    true_positives_table <- diff_pep %>% filter(name %in% LETTERS) %>% 
        dplyr::select(paste0(de_detection_columns,"_significant"))
    false_positives_table <- diff_pep %>% filter(!(name %in% LETTERS)) %>% 
        dplyr::select(paste0(de_detection_columns,"_significant"))
    
    number_of_true_positives <- colSums(true_positives_table == TRUE, na.rm = TRUE)
    number_of_false_positives <- colSums(false_positives_table == TRUE, na.rm = TRUE)
    
    tibble_of_positives <- tibble(TruePositives = number_of_true_positives, FalsePositives = number_of_false_positives)
    
    ## Which contrasts to find the FC estimate for the current proteins
    lower_contrasts <- c(paste0(paste0("Mix",lower_samples[current_median,]),"_vs_Mix",current_median))
    higher_contrasts <- c(paste0(paste0("Mix",higher_samples[current_median,]),"_vs_Mix",current_median))
    
    ## Keep only the current proteins
    diff_pep <- diff_pep %>% filter(name %in% current_proteins)
    
    ## Point estimate
    lower_tibble_point <- diff_pep %>% dplyr::select(c(name, paste0(lower_contrasts,"_diff"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "FC")
    
    higher_tibble_point <- diff_pep %>% dplyr::select(c(name, paste0(higher_contrasts,"_diff"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "FC")
    
    ## Left CI
    lower_tibble_left <- diff_pep %>% dplyr::select(c(name, paste0(lower_contrasts,"_CI.L"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Left")
    
    higher_tibble_left <- diff_pep %>% dplyr::select(c(name, paste0(higher_contrasts,"_CI.L"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Left")
    
    ## Right CI
    lower_tibble_right <- diff_pep %>% dplyr::select(c(name, paste0(lower_contrasts,"_CI.R"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Right")
    
    higher_tibble_right <- diff_pep %>% dplyr::select(c(name, paste0(higher_contrasts,"_CI.R"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Right")
    
    ## Adjusted p value
    lower_tibble_p <- diff_pep %>% dplyr::select(c(name, paste0(lower_contrasts,"_p.adj"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "adjusted_p_value")
    
    higher_tibble_p <- diff_pep %>% dplyr::select(c(name, paste0(higher_contrasts,"_p.adj"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "adjusted_p_value")
    
    ## Create FC tibble
    fc_tibble_current_lower <- tibble(Protein = lower_tibble_point$name,
                                      TrueFC = rep(log2(1.75^-(1:6)),2),
                                      EstimatedFC_Left = lower_tibble_left$Left,
                                      EstimatedFC = lower_tibble_point$FC,
                                      EstimatedFC_Right = lower_tibble_right$Right,
                                      P = lower_tibble_p$adjusted_p_value)
    
    fc_tibble_current_higher <- tibble(Protein = higher_tibble_point$name,
                                       TrueFC = rep(log2(1.75^(1:6)),2),
                                       EstimatedFC_Left = higher_tibble_left$Left,
                                       EstimatedFC = higher_tibble_point$FC,
                                       EstimatedFC_Right = higher_tibble_right$Right,
                                       P = higher_tibble_p$adjusted_p_value) 
    
    fc_tibble_current <- bind_rows(fc_tibble_current_lower,fc_tibble_current_higher)
    
    ## Store the information about the FC estimates and the detection of True and False Positives
    pe_dda_IM_RS_TP <- bind_rows(pe_dda_IM_RS_TP, tibble_of_positives)
    pe_dda_IM_RS_de <- bind_rows(pe_dda_IM_RS_de, fc_tibble_current)
}

#######################################################################################################################
## DE analysis - pe_dia_IM_RS_se
#######################################################################################################################

# Specify current dataset
current_dataset <- pe_dia_IM_RS_se

# Intialize tibble to store all the data
pe_dia_IM_RS_de <- tibble() ## This stores point estimates of the FC
pe_dia_IM_RS_TP <- tibble() ## This is to store the number of true positives and false positives

for(i in 1:13){
    
    current_pair <- i
    
    current_proteins <- spike_in_design_median_mix$ProteinID[which(spike_in_design_median_mix$PairID == current_pair)]
    current_median<-unique(spike_in_design_median_mix$Median[which(spike_in_design_median_mix$PairID == current_pair)])
    
    diff_pep <- test_diff(current_dataset, type = "control", control = paste0("Mix",current_median))
    diff_pep <- add_rejections(diff_pep, alpha = 0.01)
    diff_pep <- get_df_wide(diff_pep)
    
    ## DE detection columns
    de_detection_columns <- c(paste0(paste0("Mix",c(1:13)[-current_median]),"_vs_Mix",current_median))
    
    ## Find true positives and false positives in the current contrasts ()
    true_positives_table <- diff_pep %>% filter(name %in% LETTERS) %>% 
        dplyr::select(paste0(de_detection_columns,"_significant"))
    false_positives_table <- diff_pep %>% filter(!(name %in% LETTERS)) %>% 
        dplyr::select(paste0(de_detection_columns,"_significant"))
    
    number_of_true_positives <- colSums(true_positives_table == TRUE, na.rm = TRUE)
    number_of_false_positives <- colSums(false_positives_table == TRUE, na.rm = TRUE)
    
    tibble_of_positives <- tibble(TruePositives = number_of_true_positives, FalsePositives = number_of_false_positives)
    
    ## Which contrasts to find the FC estimate for the current proteins
    lower_contrasts <- c(paste0(paste0("Mix",lower_samples[current_median,]),"_vs_Mix",current_median))
    higher_contrasts <- c(paste0(paste0("Mix",higher_samples[current_median,]),"_vs_Mix",current_median))
    
    ## Keep only the current proteins
    diff_pep <- diff_pep %>% filter(name %in% current_proteins)
    
    ## Point estimate
    lower_tibble_point <- diff_pep %>% dplyr::select(c(name, paste0(lower_contrasts,"_diff"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "FC")
    
    higher_tibble_point <- diff_pep %>% dplyr::select(c(name, paste0(higher_contrasts,"_diff"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "FC")
    
    ## Left CI
    lower_tibble_left <- diff_pep %>% dplyr::select(c(name, paste0(lower_contrasts,"_CI.L"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Left")
    
    higher_tibble_left <- diff_pep %>% dplyr::select(c(name, paste0(higher_contrasts,"_CI.L"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Left")
    
    ## Right CI
    lower_tibble_right <- diff_pep %>% dplyr::select(c(name, paste0(lower_contrasts,"_CI.R"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Right")
    
    higher_tibble_right <- diff_pep %>% dplyr::select(c(name, paste0(higher_contrasts,"_CI.R"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Right")
    
    ## Adjusted p value
    lower_tibble_p <- diff_pep %>% dplyr::select(c(name, paste0(lower_contrasts,"_p.adj"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "adjusted_p_value")
    
    higher_tibble_p <- diff_pep %>% dplyr::select(c(name, paste0(higher_contrasts,"_p.adj"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "adjusted_p_value")
    
    ## Create FC tibble
    fc_tibble_current_lower <- tibble(Protein = lower_tibble_point$name,
                                      TrueFC = rep(log2(1.75^-(1:6)),2),
                                      EstimatedFC_Left = lower_tibble_left$Left,
                                      EstimatedFC = lower_tibble_point$FC,
                                      EstimatedFC_Right = lower_tibble_right$Right,
                                      P = lower_tibble_p$adjusted_p_value)
    
    fc_tibble_current_higher <- tibble(Protein = higher_tibble_point$name,
                                       TrueFC = rep(log2(1.75^(1:6)),2),
                                       EstimatedFC_Left = higher_tibble_left$Left,
                                       EstimatedFC = higher_tibble_point$FC,
                                       EstimatedFC_Right = higher_tibble_right$Right,
                                       P = higher_tibble_p$adjusted_p_value) 
    
    fc_tibble_current <- bind_rows(fc_tibble_current_lower,fc_tibble_current_higher)
    
    ## Store the information about the FC estimates and the detection of True and False Positives
    pe_dia_IM_RS_TP <- bind_rows(pe_dia_IM_RS_TP, tibble_of_positives)
    pe_dia_IM_RS_de <- bind_rows(pe_dia_IM_RS_de, fc_tibble_current)
}

#######################################################################################################################
## DE analysis - pe_dda_IN_MP_se
#######################################################################################################################

# Specify current dataset
current_dataset <- pe_dda_IN_MP_se

# Intialize tibble to store all the data
pe_dda_IN_MP_de <- tibble() ## This stores point estimates of the FC
pe_dda_IN_MP_TP <- tibble() ## This is to store the number of true positives and false positives

for(i in 1:13){
    
    current_pair <- i
    
    current_proteins <- spike_in_design_median_mix$ProteinID[which(spike_in_design_median_mix$PairID == current_pair)]
    current_median<-unique(spike_in_design_median_mix$Median[which(spike_in_design_median_mix$PairID == current_pair)])
    
    diff_pep <- test_diff(current_dataset, type = "control", control = paste0("Mix",current_median))
    diff_pep <- add_rejections(diff_pep, alpha = 0.01)
    diff_pep <- get_df_wide(diff_pep)
    
    ## DE detection columns
    de_detection_columns <- c(paste0(paste0("Mix",c(1:13)[-current_median]),"_vs_Mix",current_median))
    
    ## Find true positives and false positives in the current contrasts ()
    true_positives_table <- diff_pep %>% filter(name %in% LETTERS) %>% 
        dplyr::select(paste0(de_detection_columns,"_significant"))
    false_positives_table <- diff_pep %>% filter(!(name %in% LETTERS)) %>% 
        dplyr::select(paste0(de_detection_columns,"_significant"))
    
    number_of_true_positives <- colSums(true_positives_table == TRUE, na.rm = TRUE)
    number_of_false_positives <- colSums(false_positives_table == TRUE, na.rm = TRUE)
    
    tibble_of_positives <- tibble(TruePositives = number_of_true_positives, FalsePositives = number_of_false_positives)
    
    ## Which contrasts to find the FC estimate for the current proteins
    lower_contrasts <- c(paste0(paste0("Mix",lower_samples[current_median,]),"_vs_Mix",current_median))
    higher_contrasts <- c(paste0(paste0("Mix",higher_samples[current_median,]),"_vs_Mix",current_median))
    
    ## Keep only the current proteins
    diff_pep <- diff_pep %>% filter(name %in% current_proteins)
    
    ## Point estimate
    lower_tibble_point <- diff_pep %>% dplyr::select(c(name, paste0(lower_contrasts,"_diff"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "FC")
    
    higher_tibble_point <- diff_pep %>% dplyr::select(c(name, paste0(higher_contrasts,"_diff"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "FC")
    
    ## Left CI
    lower_tibble_left <- diff_pep %>% dplyr::select(c(name, paste0(lower_contrasts,"_CI.L"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Left")
    
    higher_tibble_left <- diff_pep %>% dplyr::select(c(name, paste0(higher_contrasts,"_CI.L"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Left")
    
    ## Right CI
    lower_tibble_right <- diff_pep %>% dplyr::select(c(name, paste0(lower_contrasts,"_CI.R"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Right")
    
    higher_tibble_right <- diff_pep %>% dplyr::select(c(name, paste0(higher_contrasts,"_CI.R"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Right")
    
    ## Adjusted p value
    lower_tibble_p <- diff_pep %>% dplyr::select(c(name, paste0(lower_contrasts,"_p.adj"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "adjusted_p_value")
    
    higher_tibble_p <- diff_pep %>% dplyr::select(c(name, paste0(higher_contrasts,"_p.adj"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "adjusted_p_value")
    
    ## Create FC tibble
    fc_tibble_current_lower <- tibble(Protein = lower_tibble_point$name,
                                      TrueFC = rep(log2(1.75^-(1:6)),2),
                                      EstimatedFC_Left = lower_tibble_left$Left,
                                      EstimatedFC = lower_tibble_point$FC,
                                      EstimatedFC_Right = lower_tibble_right$Right,
                                      P = lower_tibble_p$adjusted_p_value)
    
    fc_tibble_current_higher <- tibble(Protein = higher_tibble_point$name,
                                       TrueFC = rep(log2(1.75^(1:6)),2),
                                       EstimatedFC_Left = higher_tibble_left$Left,
                                       EstimatedFC = higher_tibble_point$FC,
                                       EstimatedFC_Right = higher_tibble_right$Right,
                                       P = higher_tibble_p$adjusted_p_value) 
    
    fc_tibble_current <- bind_rows(fc_tibble_current_lower,fc_tibble_current_higher)
    
    ## Store the information about the FC estimates and the detection of True and False Positives
    pe_dda_IN_MP_TP <- bind_rows(pe_dda_IN_MP_TP, tibble_of_positives)
    pe_dda_IN_MP_de <- bind_rows(pe_dda_IN_MP_de, fc_tibble_current)
}

#######################################################################################################################
## DE analysis - pe_dia_IN_MP_se
#######################################################################################################################

# Specify current dataset
current_dataset <- pe_dia_IN_MP_se

# Intialize tibble to store all the data
pe_dia_IN_MP_de <- tibble() ## This stores point estimates of the FC
pe_dia_IN_MP_TP <- tibble() ## This is to store the number of true positives and false positives

for(i in 1:13){
    
    current_pair <- i
    
    current_proteins <- spike_in_design_median_mix$ProteinID[which(spike_in_design_median_mix$PairID == current_pair)]
    current_median<-unique(spike_in_design_median_mix$Median[which(spike_in_design_median_mix$PairID == current_pair)])
    
    diff_pep <- test_diff(current_dataset, type = "control", control = paste0("Mix",current_median))
    diff_pep <- add_rejections(diff_pep, alpha = 0.01)
    diff_pep <- get_df_wide(diff_pep)
    
    ## DE detection columns
    de_detection_columns <- c(paste0(paste0("Mix",c(1:13)[-current_median]),"_vs_Mix",current_median))
    
    ## Find true positives and false positives in the current contrasts ()
    true_positives_table <- diff_pep %>% filter(name %in% LETTERS) %>% 
        dplyr::select(paste0(de_detection_columns,"_significant"))
    false_positives_table <- diff_pep %>% filter(!(name %in% LETTERS)) %>% 
        dplyr::select(paste0(de_detection_columns,"_significant"))
    
    number_of_true_positives <- colSums(true_positives_table == TRUE, na.rm = TRUE)
    number_of_false_positives <- colSums(false_positives_table == TRUE, na.rm = TRUE)
    
    tibble_of_positives <- tibble(TruePositives = number_of_true_positives, FalsePositives = number_of_false_positives)
    
    ## Which contrasts to find the FC estimate for the current proteins
    lower_contrasts <- c(paste0(paste0("Mix",lower_samples[current_median,]),"_vs_Mix",current_median))
    higher_contrasts <- c(paste0(paste0("Mix",higher_samples[current_median,]),"_vs_Mix",current_median))
    
    ## Keep only the current proteins
    diff_pep <- diff_pep %>% filter(name %in% current_proteins)
    
    ## Point estimate
    lower_tibble_point <- diff_pep %>% dplyr::select(c(name, paste0(lower_contrasts,"_diff"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "FC")
    
    higher_tibble_point <- diff_pep %>% dplyr::select(c(name, paste0(higher_contrasts,"_diff"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "FC")
    
    ## Left CI
    lower_tibble_left <- diff_pep %>% dplyr::select(c(name, paste0(lower_contrasts,"_CI.L"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Left")
    
    higher_tibble_left <- diff_pep %>% dplyr::select(c(name, paste0(higher_contrasts,"_CI.L"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Left")
    
    ## Right CI
    lower_tibble_right <- diff_pep %>% dplyr::select(c(name, paste0(lower_contrasts,"_CI.R"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Right")
    
    higher_tibble_right <- diff_pep %>% dplyr::select(c(name, paste0(higher_contrasts,"_CI.R"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Right")
    
    ## Adjusted p value
    lower_tibble_p <- diff_pep %>% dplyr::select(c(name, paste0(lower_contrasts,"_p.adj"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "adjusted_p_value")
    
    higher_tibble_p <- diff_pep %>% dplyr::select(c(name, paste0(higher_contrasts,"_p.adj"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "adjusted_p_value")
    
    ## Create FC tibble
    fc_tibble_current_lower <- tibble(Protein = lower_tibble_point$name,
                                      TrueFC = rep(log2(1.75^-(1:6)),2),
                                      EstimatedFC_Left = lower_tibble_left$Left,
                                      EstimatedFC = lower_tibble_point$FC,
                                      EstimatedFC_Right = lower_tibble_right$Right,
                                      P = lower_tibble_p$adjusted_p_value)
    
    fc_tibble_current_higher <- tibble(Protein = higher_tibble_point$name,
                                       TrueFC = rep(log2(1.75^(1:6)),2),
                                       EstimatedFC_Left = higher_tibble_left$Left,
                                       EstimatedFC = higher_tibble_point$FC,
                                       EstimatedFC_Right = higher_tibble_right$Right,
                                       P = higher_tibble_p$adjusted_p_value) 
    
    fc_tibble_current <- bind_rows(fc_tibble_current_lower,fc_tibble_current_higher)
    
    ## Store the information about the FC estimates and the detection of True and False Positives
    pe_dia_IN_MP_TP <- bind_rows(pe_dia_IN_MP_TP, tibble_of_positives)
    pe_dia_IN_MP_de <- bind_rows(pe_dia_IN_MP_de, fc_tibble_current)
}

#######################################################################################################################
## DE analysis - pe_dda_IQ_MP_se
#######################################################################################################################

# Specify current dataset
current_dataset <- pe_dda_IQ_MP_se

# Intialize tibble to store all the data
pe_dda_IQ_MP_de <- tibble() ## This stores point estimates of the FC
pe_dda_IQ_MP_TP <- tibble() ## This is to store the number of true positives and false positives

for(i in 1:13){
    
    current_pair <- i
    
    current_proteins <- spike_in_design_median_mix$ProteinID[which(spike_in_design_median_mix$PairID == current_pair)]
    current_median<-unique(spike_in_design_median_mix$Median[which(spike_in_design_median_mix$PairID == current_pair)])
    
    diff_pep <- test_diff(current_dataset, type = "control", control = paste0("Mix",current_median))
    diff_pep <- add_rejections(diff_pep, alpha = 0.01)
    diff_pep <- get_df_wide(diff_pep)
    
    ## DE detection columns
    de_detection_columns <- c(paste0(paste0("Mix",c(1:13)[-current_median]),"_vs_Mix",current_median))
    
    ## Find true positives and false positives in the current contrasts ()
    true_positives_table <- diff_pep %>% filter(name %in% LETTERS) %>% 
        dplyr::select(paste0(de_detection_columns,"_significant"))
    false_positives_table <- diff_pep %>% filter(!(name %in% LETTERS)) %>% 
        dplyr::select(paste0(de_detection_columns,"_significant"))
    
    number_of_true_positives <- colSums(true_positives_table == TRUE, na.rm = TRUE)
    number_of_false_positives <- colSums(false_positives_table == TRUE, na.rm = TRUE)
    
    tibble_of_positives <- tibble(TruePositives = number_of_true_positives, FalsePositives = number_of_false_positives)
    
    ## Which contrasts to find the FC estimate for the current proteins
    lower_contrasts <- c(paste0(paste0("Mix",lower_samples[current_median,]),"_vs_Mix",current_median))
    higher_contrasts <- c(paste0(paste0("Mix",higher_samples[current_median,]),"_vs_Mix",current_median))
    
    ## Keep only the current proteins
    diff_pep <- diff_pep %>% filter(name %in% current_proteins)
    
    ## Point estimate
    lower_tibble_point <- diff_pep %>% dplyr::select(c(name, paste0(lower_contrasts,"_diff"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "FC")
    
    higher_tibble_point <- diff_pep %>% dplyr::select(c(name, paste0(higher_contrasts,"_diff"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "FC")
    
    ## Left CI
    lower_tibble_left <- diff_pep %>% dplyr::select(c(name, paste0(lower_contrasts,"_CI.L"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Left")
    
    higher_tibble_left <- diff_pep %>% dplyr::select(c(name, paste0(higher_contrasts,"_CI.L"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Left")
    
    ## Right CI
    lower_tibble_right <- diff_pep %>% dplyr::select(c(name, paste0(lower_contrasts,"_CI.R"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Right")
    
    higher_tibble_right <- diff_pep %>% dplyr::select(c(name, paste0(higher_contrasts,"_CI.R"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Right")
    
    ## Adjusted p value
    lower_tibble_p <- diff_pep %>% dplyr::select(c(name, paste0(lower_contrasts,"_p.adj"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "adjusted_p_value")
    
    higher_tibble_p <- diff_pep %>% dplyr::select(c(name, paste0(higher_contrasts,"_p.adj"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "adjusted_p_value")
    
    ## Create FC tibble
    fc_tibble_current_lower <- tibble(Protein = lower_tibble_point$name,
                                      TrueFC = rep(log2(1.75^-(1:6)),2),
                                      EstimatedFC_Left = lower_tibble_left$Left,
                                      EstimatedFC = lower_tibble_point$FC,
                                      EstimatedFC_Right = lower_tibble_right$Right,
                                      P = lower_tibble_p$adjusted_p_value)
    
    fc_tibble_current_higher <- tibble(Protein = higher_tibble_point$name,
                                       TrueFC = rep(log2(1.75^(1:6)),2),
                                       EstimatedFC_Left = higher_tibble_left$Left,
                                       EstimatedFC = higher_tibble_point$FC,
                                       EstimatedFC_Right = higher_tibble_right$Right,
                                       P = higher_tibble_p$adjusted_p_value) 
    
    fc_tibble_current <- bind_rows(fc_tibble_current_lower,fc_tibble_current_higher)
    
    ## Store the information about the FC estimates and the detection of True and False Positives
    pe_dda_IQ_MP_TP <- bind_rows(pe_dda_IQ_MP_TP, tibble_of_positives)
    pe_dda_IQ_MP_de <- bind_rows(pe_dda_IQ_MP_de, fc_tibble_current)
}

#######################################################################################################################
## DE analysis - pe_dia_IQ_MP_se
#######################################################################################################################

# Specify current dataset
current_dataset <- pe_dia_IQ_MP_se

# Intialize tibble to store all the data
pe_dia_IQ_MP_de <- tibble() ## This stores point estimates of the FC
pe_dia_IQ_MP_TP <- tibble() ## This is to store the number of true positives and false positives

for(i in 1:13){
    
    current_pair <- i
    
    current_proteins <- spike_in_design_median_mix$ProteinID[which(spike_in_design_median_mix$PairID == current_pair)]
    current_median<-unique(spike_in_design_median_mix$Median[which(spike_in_design_median_mix$PairID == current_pair)])
    
    diff_pep <- test_diff(current_dataset, type = "control", control = paste0("Mix",current_median))
    diff_pep <- add_rejections(diff_pep, alpha = 0.01)
    diff_pep <- get_df_wide(diff_pep)
    
    ## DE detection columns
    de_detection_columns <- c(paste0(paste0("Mix",c(1:13)[-current_median]),"_vs_Mix",current_median))
    
    ## Find true positives and false positives in the current contrasts ()
    true_positives_table <- diff_pep %>% filter(name %in% LETTERS) %>% 
        dplyr::select(paste0(de_detection_columns,"_significant"))
    false_positives_table <- diff_pep %>% filter(!(name %in% LETTERS)) %>% 
        dplyr::select(paste0(de_detection_columns,"_significant"))
    
    number_of_true_positives <- colSums(true_positives_table == TRUE, na.rm = TRUE)
    number_of_false_positives <- colSums(false_positives_table == TRUE, na.rm = TRUE)
    
    tibble_of_positives <- tibble(TruePositives = number_of_true_positives, FalsePositives = number_of_false_positives)
    
    ## Which contrasts to find the FC estimate for the current proteins
    lower_contrasts <- c(paste0(paste0("Mix",lower_samples[current_median,]),"_vs_Mix",current_median))
    higher_contrasts <- c(paste0(paste0("Mix",higher_samples[current_median,]),"_vs_Mix",current_median))
    
    ## Keep only the current proteins
    diff_pep <- diff_pep %>% filter(name %in% current_proteins)
    
    ## Point estimate
    lower_tibble_point <- diff_pep %>% dplyr::select(c(name, paste0(lower_contrasts,"_diff"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "FC")
    
    higher_tibble_point <- diff_pep %>% dplyr::select(c(name, paste0(higher_contrasts,"_diff"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "FC")
    
    ## Left CI
    lower_tibble_left <- diff_pep %>% dplyr::select(c(name, paste0(lower_contrasts,"_CI.L"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Left")
    
    higher_tibble_left <- diff_pep %>% dplyr::select(c(name, paste0(higher_contrasts,"_CI.L"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Left")
    
    ## Right CI
    lower_tibble_right <- diff_pep %>% dplyr::select(c(name, paste0(lower_contrasts,"_CI.R"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Right")
    
    higher_tibble_right <- diff_pep %>% dplyr::select(c(name, paste0(higher_contrasts,"_CI.R"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Right")
    
    ## Adjusted p value
    lower_tibble_p <- diff_pep %>% dplyr::select(c(name, paste0(lower_contrasts,"_p.adj"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "adjusted_p_value")
    
    higher_tibble_p <- diff_pep %>% dplyr::select(c(name, paste0(higher_contrasts,"_p.adj"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "adjusted_p_value")
    
    ## Create FC tibble
    fc_tibble_current_lower <- tibble(Protein = lower_tibble_point$name,
                                      TrueFC = rep(log2(1.75^-(1:6)),2),
                                      EstimatedFC_Left = lower_tibble_left$Left,
                                      EstimatedFC = lower_tibble_point$FC,
                                      EstimatedFC_Right = lower_tibble_right$Right,
                                      P = lower_tibble_p$adjusted_p_value)
    
    fc_tibble_current_higher <- tibble(Protein = higher_tibble_point$name,
                                       TrueFC = rep(log2(1.75^(1:6)),2),
                                       EstimatedFC_Left = higher_tibble_left$Left,
                                       EstimatedFC = higher_tibble_point$FC,
                                       EstimatedFC_Right = higher_tibble_right$Right,
                                       P = higher_tibble_p$adjusted_p_value) 
    
    fc_tibble_current <- bind_rows(fc_tibble_current_lower,fc_tibble_current_higher)
    
    ## Store the information about the FC estimates and the detection of True and False Positives
    pe_dia_IQ_MP_TP <- bind_rows(pe_dia_IQ_MP_TP, tibble_of_positives)
    pe_dia_IQ_MP_de <- bind_rows(pe_dia_IQ_MP_de, fc_tibble_current)
}

#######################################################################################################################
## DE analysis - pe_dda_IM_MP_se
#######################################################################################################################

# Specify current dataset
current_dataset <- pe_dda_IM_MP_se

# Intialize tibble to store all the data
pe_dda_IM_MP_de <- tibble() ## This stores point estimates of the FC
pe_dda_IM_MP_TP <- tibble() ## This is to store the number of true positives and false positives

for(i in 1:13){
    
    current_pair <- i
    
    current_proteins <- spike_in_design_median_mix$ProteinID[which(spike_in_design_median_mix$PairID == current_pair)]
    current_median<-unique(spike_in_design_median_mix$Median[which(spike_in_design_median_mix$PairID == current_pair)])
    
    diff_pep <- test_diff(current_dataset, type = "control", control = paste0("Mix",current_median))
    diff_pep <- add_rejections(diff_pep, alpha = 0.01)
    diff_pep <- get_df_wide(diff_pep)
    
    ## DE detection columns
    de_detection_columns <- c(paste0(paste0("Mix",c(1:13)[-current_median]),"_vs_Mix",current_median))
    
    ## Find true positives and false positives in the current contrasts ()
    true_positives_table <- diff_pep %>% filter(name %in% LETTERS) %>% 
        dplyr::select(paste0(de_detection_columns,"_significant"))
    false_positives_table <- diff_pep %>% filter(!(name %in% LETTERS)) %>% 
        dplyr::select(paste0(de_detection_columns,"_significant"))
    
    number_of_true_positives <- colSums(true_positives_table == TRUE, na.rm = TRUE)
    number_of_false_positives <- colSums(false_positives_table == TRUE, na.rm = TRUE)
    
    tibble_of_positives <- tibble(TruePositives = number_of_true_positives, FalsePositives = number_of_false_positives)
    
    ## Which contrasts to find the FC estimate for the current proteins
    lower_contrasts <- c(paste0(paste0("Mix",lower_samples[current_median,]),"_vs_Mix",current_median))
    higher_contrasts <- c(paste0(paste0("Mix",higher_samples[current_median,]),"_vs_Mix",current_median))
    
    ## Keep only the current proteins
    diff_pep <- diff_pep %>% filter(name %in% current_proteins)
    
    ## Point estimate
    lower_tibble_point <- diff_pep %>% dplyr::select(c(name, paste0(lower_contrasts,"_diff"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "FC")
    
    higher_tibble_point <- diff_pep %>% dplyr::select(c(name, paste0(higher_contrasts,"_diff"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "FC")
    
    ## Left CI
    lower_tibble_left <- diff_pep %>% dplyr::select(c(name, paste0(lower_contrasts,"_CI.L"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Left")
    
    higher_tibble_left <- diff_pep %>% dplyr::select(c(name, paste0(higher_contrasts,"_CI.L"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Left")
    
    ## Right CI
    lower_tibble_right <- diff_pep %>% dplyr::select(c(name, paste0(lower_contrasts,"_CI.R"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Right")
    
    higher_tibble_right <- diff_pep %>% dplyr::select(c(name, paste0(higher_contrasts,"_CI.R"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Right")
    
    ## Adjusted p value
    lower_tibble_p <- diff_pep %>% dplyr::select(c(name, paste0(lower_contrasts,"_p.adj"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "adjusted_p_value")
    
    higher_tibble_p <- diff_pep %>% dplyr::select(c(name, paste0(higher_contrasts,"_p.adj"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "adjusted_p_value")
    
    ## Create FC tibble
    fc_tibble_current_lower <- tibble(Protein = lower_tibble_point$name,
                                      TrueFC = rep(log2(1.75^-(1:6)),2),
                                      EstimatedFC_Left = lower_tibble_left$Left,
                                      EstimatedFC = lower_tibble_point$FC,
                                      EstimatedFC_Right = lower_tibble_right$Right,
                                      P = lower_tibble_p$adjusted_p_value)
    
    fc_tibble_current_higher <- tibble(Protein = higher_tibble_point$name,
                                       TrueFC = rep(log2(1.75^(1:6)),2),
                                       EstimatedFC_Left = higher_tibble_left$Left,
                                       EstimatedFC = higher_tibble_point$FC,
                                       EstimatedFC_Right = higher_tibble_right$Right,
                                       P = higher_tibble_p$adjusted_p_value) 
    
    fc_tibble_current <- bind_rows(fc_tibble_current_lower,fc_tibble_current_higher)
    
    ## Store the information about the FC estimates and the detection of True and False Positives
    pe_dda_IM_MP_TP <- bind_rows(pe_dda_IM_MP_TP, tibble_of_positives)
    pe_dda_IM_MP_de <- bind_rows(pe_dda_IM_MP_de, fc_tibble_current)
}

#######################################################################################################################
## DE analysis - pe_dia_IM_MP_se
#######################################################################################################################

# Specify current dataset
current_dataset <- pe_dia_IM_MP_se

# Intialize tibble to store all the data
pe_dia_IM_MP_de <- tibble() ## This stores point estimates of the FC
pe_dia_IM_MP_TP <- tibble() ## This is to store the number of true positives and false positives

for(i in 1:13){
    
    current_pair <- i
    
    current_proteins <- spike_in_design_median_mix$ProteinID[which(spike_in_design_median_mix$PairID == current_pair)]
    current_median<-unique(spike_in_design_median_mix$Median[which(spike_in_design_median_mix$PairID == current_pair)])
    
    diff_pep <- test_diff(current_dataset, type = "control", control = paste0("Mix",current_median))
    diff_pep <- add_rejections(diff_pep, alpha = 0.01)
    diff_pep <- get_df_wide(diff_pep)
    
    ## DE detection columns
    de_detection_columns <- c(paste0(paste0("Mix",c(1:13)[-current_median]),"_vs_Mix",current_median))
    
    ## Find true positives and false positives in the current contrasts ()
    true_positives_table <- diff_pep %>% filter(name %in% LETTERS) %>% 
        dplyr::select(paste0(de_detection_columns,"_significant"))
    false_positives_table <- diff_pep %>% filter(!(name %in% LETTERS)) %>% 
        dplyr::select(paste0(de_detection_columns,"_significant"))
    
    number_of_true_positives <- colSums(true_positives_table == TRUE, na.rm = TRUE)
    number_of_false_positives <- colSums(false_positives_table == TRUE, na.rm = TRUE)
    
    tibble_of_positives <- tibble(TruePositives = number_of_true_positives, FalsePositives = number_of_false_positives)
    
    ## Which contrasts to find the FC estimate for the current proteins
    lower_contrasts <- c(paste0(paste0("Mix",lower_samples[current_median,]),"_vs_Mix",current_median))
    higher_contrasts <- c(paste0(paste0("Mix",higher_samples[current_median,]),"_vs_Mix",current_median))
    
    ## Keep only the current proteins
    diff_pep <- diff_pep %>% filter(name %in% current_proteins)
    
    ## Point estimate
    lower_tibble_point <- diff_pep %>% dplyr::select(c(name, paste0(lower_contrasts,"_diff"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "FC")
    
    higher_tibble_point <- diff_pep %>% dplyr::select(c(name, paste0(higher_contrasts,"_diff"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "FC")
    
    ## Left CI
    lower_tibble_left <- diff_pep %>% dplyr::select(c(name, paste0(lower_contrasts,"_CI.L"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Left")
    
    higher_tibble_left <- diff_pep %>% dplyr::select(c(name, paste0(higher_contrasts,"_CI.L"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Left")
    
    ## Right CI
    lower_tibble_right <- diff_pep %>% dplyr::select(c(name, paste0(lower_contrasts,"_CI.R"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Right")
    
    higher_tibble_right <- diff_pep %>% dplyr::select(c(name, paste0(higher_contrasts,"_CI.R"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Right")
    
    ## Adjusted p value
    lower_tibble_p <- diff_pep %>% dplyr::select(c(name, paste0(lower_contrasts,"_p.adj"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "adjusted_p_value")
    
    higher_tibble_p <- diff_pep %>% dplyr::select(c(name, paste0(higher_contrasts,"_p.adj"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "adjusted_p_value")
    
    ## Create FC tibble
    fc_tibble_current_lower <- tibble(Protein = lower_tibble_point$name,
                                      TrueFC = rep(log2(1.75^-(1:6)),2),
                                      EstimatedFC_Left = lower_tibble_left$Left,
                                      EstimatedFC = lower_tibble_point$FC,
                                      EstimatedFC_Right = lower_tibble_right$Right,
                                      P = lower_tibble_p$adjusted_p_value)
    
    fc_tibble_current_higher <- tibble(Protein = higher_tibble_point$name,
                                       TrueFC = rep(log2(1.75^(1:6)),2),
                                       EstimatedFC_Left = higher_tibble_left$Left,
                                       EstimatedFC = higher_tibble_point$FC,
                                       EstimatedFC_Right = higher_tibble_right$Right,
                                       P = higher_tibble_p$adjusted_p_value) 
    
    fc_tibble_current <- bind_rows(fc_tibble_current_lower,fc_tibble_current_higher)
    
    ## Store the information about the FC estimates and the detection of True and False Positives
    pe_dia_IM_MP_TP <- bind_rows(pe_dia_IM_MP_TP, tibble_of_positives)
    pe_dia_IM_MP_de <- bind_rows(pe_dia_IM_MP_de, fc_tibble_current)
}

#######################################################################################################################
## Store all the information in one big tibble
#######################################################################################################################

differential_expression_estimates <- bind_rows(pe_dda_IN_RS_de %>% mutate(Dataset = "DDA_IN_RS"),
                                               pe_dia_IN_RS_de %>% mutate(Dataset = "DIA_IN_RS"),
                                               pe_dda_IQ_RS_de %>% mutate(Dataset = "DDA_IQ_RS"),
                                               pe_dia_IQ_RS_de %>% mutate(Dataset = "DIA_IQ_RS"),
                                               pe_dda_IM_RS_de %>% mutate(Dataset = "DDA_IM_RS"),
                                               pe_dia_IM_RS_de %>% mutate(Dataset = "DIA_IM_RS"),
                                               pe_dda_IN_MP_de %>% mutate(Dataset = "DDA_IN_MP"),
                                               pe_dia_IN_MP_de %>% mutate(Dataset = "DIA_IN_MP"),
                                               pe_dda_IQ_MP_de %>% mutate(Dataset = "DDA_IQ_MP"),
                                               pe_dia_IQ_MP_de %>% mutate(Dataset = "DIA_IQ_MP"),
                                               pe_dda_IM_MP_de %>% mutate(Dataset = "DDA_IM_MP"),
                                               pe_dia_IM_MP_de %>% mutate(Dataset = "DIA_IM_MP"))

differential_expression_estimates$Dataset <- factor(differential_expression_estimates$Dataset,
                                                    levels = c("DDA_IN_RS", "DDA_IN_MP",
                                                               "DDA_IQ_RS", "DDA_IQ_MP",
                                                               "DDA_IM_RS", "DDA_IM_MP",
                                                               "DIA_IN_RS", "DIA_IN_MP",
                                                               "DIA_IQ_RS", "DIA_IQ_MP",
                                                               "DIA_IM_RS", "DIA_IM_MP"))

differential_expression_detection <- bind_rows(pe_dda_IN_RS_TP %>% mutate(Dataset = "DDA_IN_RS"),
                                               pe_dia_IN_RS_TP %>% mutate(Dataset = "DIA_IN_RS"),
                                               pe_dda_IQ_RS_TP %>% mutate(Dataset = "DDA_IQ_RS"),
                                               pe_dia_IQ_RS_TP %>% mutate(Dataset = "DIA_IQ_RS"),
                                               pe_dda_IM_RS_TP %>% mutate(Dataset = "DDA_IM_RS"),
                                               pe_dia_IM_RS_TP %>% mutate(Dataset = "DIA_IM_RS"),
                                               pe_dda_IN_MP_TP %>% mutate(Dataset = "DDA_IN_MP"),
                                               pe_dia_IN_MP_TP %>% mutate(Dataset = "DIA_IN_MP"),
                                               pe_dda_IQ_MP_TP %>% mutate(Dataset = "DDA_IQ_MP"),
                                               pe_dia_IQ_MP_TP %>% mutate(Dataset = "DIA_IQ_MP"),
                                               pe_dda_IM_MP_TP %>% mutate(Dataset = "DDA_IM_MP"),
                                               pe_dia_IM_MP_TP %>% mutate(Dataset = "DIA_IM_MP"))


differential_expression_detection$Dataset <- factor(differential_expression_detection$Dataset,
                                                    levels = c("DDA_IN_RS", "DDA_IN_MP",
                                                               "DDA_IQ_RS", "DDA_IQ_MP",
                                                               "DDA_IM_RS", "DDA_IM_MP",
                                                               "DIA_IN_RS", "DIA_IN_MP",
                                                               "DIA_IQ_RS", "DIA_IQ_MP",
                                                               "DIA_IM_RS", "DIA_IM_MP"))
