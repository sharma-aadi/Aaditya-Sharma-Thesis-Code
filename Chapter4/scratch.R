#######################################################################################################################
## Step - 1.1: load the protein tables (proteins_dda_51_mbr_on -> proteins_dda, proteins_dia_51_all_lib -> proteins_dia)
#######################################################################################################################

## DDA protein tables
proteins_dda_table <- read_maxquant_table(file_name = "./Data\ tables/maxquant_protein_table.txt") 

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

## Create unique data objects
data_unique_dda <- make_unique(filter(proteins_dda_table, Reverse != "+"), "Gene.names", "Protein.IDs", delim = ";")

## Rename spike-in proteins with their letter names
data_unique_dda$name[which(data_unique_dda$ID %in% protein_uniprot_id)] <- 
    protein_letters[data_unique_dda$ID[which(data_unique_dda$ID %in% protein_uniprot_id)]]

## Fetch column numbers corresponding to the LFQ intensities
lfq_columns_dda <- grep("LFQ.", colnames(proteins_dda_table))

## Column names for each dataset
new_col_names_dda <- str_replace_all(colnames(data_unique_dda[, lfq_columns_dda]),
                                     pattern = "LFQ.intensity.", replacement = "")

## Create experimental design object
exp_design_dda <- data.frame(label = c(new_col_names_dda[c(1:13,18:30,35:47)],new_col_names_dda[-c(1:13,18:30,35:47)]),
                             condition = c(rep(paste0("Mix_",Latin_square_index),3), rep("Control", 12)),  
                             replicate = c(rep(1:3,each = 13), 1:12))

## Re-factor Conditions
exp_design_dda$condition <- factor(exp_design_dda$condition,
                                   levels = c("Control", paste0("Mix_",1:13)))

## Re-arrange experimental design object
exp_design_dda <- exp_design_dda %>% arrange(condition)

## Create summarized experiment object
data_se_dda <- make_se(data_unique_dda, lfq_columns_dda, exp_design_dda)

#######################################################################################################################
## Filtering, normalization, and imputation on the Summarized Experiment object
#######################################################################################################################

## Filter out proteins with too many missing values
missing_values_allowed <- 0

data_filt_dda <- filter_missval(data_se_dda, thr = missing_values_allowed) ## allowed missing value(s) per condition 

## Variance stabilizing normalization
data_norm_dda <- normalize_vsn(data_filt_dda)

## Imputing missing values
data_imp_dda <- impute(data_norm_dda, fun = "MinProb", q = 0.01)

#######################################################################################################################
## Plots - 1
#######################################################################################################################

## Plot normalizations
plot_normalization(data_filt_dda)+theme(axis.text = element_text(size = 10))

plot_normalization(data_norm_dda)+theme(axis.text = element_text(size = 10))

## Average intensity of proteins with missing values
plot_detect(data_filt_dda)

## Plot imputations
plot_imputation(data_norm_dda, data_imp_dda)

#######################################################################################################################
## pca plot
#######################################################################################################################

data_diff_dda <- test_diff(data_imp_dda, type = "control", control = "control")

dep_dda <- add_rejections(data_diff_dda, alpha = 0.00005, lfc = log2(1.75))

plot_pca(dep_dda, x = 1, y = 2, n = 28, point_size = 3)

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

#######################################################################################################################
## DE analysis - DDA (Median contrasts)
#######################################################################################################################

# Intialize tibble to store all the data
de_dda <- tibble()

for(i in 1:13){
    
    current_pair <- i
    
    current_proteins <- spike_in_design_median_mix$ProteinID[which(spike_in_design_median_mix$PairID == current_pair)]
    current_median<-unique(spike_in_design_median_mix$Median[which(spike_in_design_median_mix$PairID == current_pair)])
    
    data_diff <- get_df_wide(test_diff(data_imp_dda, type = "control", control = paste0("Mix_",current_median)))
    data_diff <- data_diff %>% filter(name %in% current_proteins)
    
    lower_contrasts <- c(paste0(paste0("Mix_",lower_samples[current_median,]),"_vs_Mix_",current_median))
    higher_contrasts <- c(paste0(paste0("Mix_",higher_samples[current_median,]),"_vs_Mix_",current_median))
    
    ## Point estimate
    lower_tibble_point <- data_diff %>% dplyr::select(c(name, paste0(lower_contrasts,"_diff"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "FC")
    
    higher_tibble_point <- data_diff %>% dplyr::select(c(name, paste0(higher_contrasts,"_diff"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "FC")
    
    ## Left CI
    lower_tibble_left <- data_diff %>% dplyr::select(c(name, paste0(lower_contrasts,"_CI.L"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Left")
    
    higher_tibble_left <- data_diff %>% dplyr::select(c(name, paste0(higher_contrasts,"_CI.L"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Left")
    
    ## Right CI
    lower_tibble_right <- data_diff %>% dplyr::select(c(name, paste0(lower_contrasts,"_CI.R"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Right")
    
    higher_tibble_right <- data_diff %>% dplyr::select(c(name, paste0(higher_contrasts,"_CI.R"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "Right")
    
    ## Adjusted p value
    lower_tibble_p <- data_diff %>% dplyr::select(c(name, paste0(lower_contrasts,"_p.adj"))) %>% 
        pivot_longer(-c(name), names_to = "Comparison", values_to = "adjusted_p_value")
    
    higher_tibble_p <- data_diff %>% dplyr::select(c(name, paste0(higher_contrasts,"_p.adj"))) %>% 
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
    
    de_dda <- bind_rows(de_dda, fc_tibble_current)
}

#######################################################################################################################
## DE analysis - DDA plots
#######################################################################################################################

positive_FC <- unique((de_dda %>% filter(TrueFC>0))$TrueFC)
positive_FC_labels <- paste0("log2FC = ",round(positive_FC,2))
names(positive_FC_labels) <- positive_FC

negative_FC <- unique((de_dda %>% filter(TrueFC<0))$TrueFC)
negative_FC_labels <- paste0("log2FC = ",round(negative_FC,2))
names(negative_FC_labels) <- negative_FC

ggplot(data = de_dda %>% filter(TrueFC>0), aes(x = as_factor(Protein), y = EstimatedFC))+
    geom_bar(stat="identity")+
    labs(x = "Protein", title = "DDA dataset")+
    geom_hline(data = tibble(TrueFC = positive_FC, FC = positive_FC), aes(yintercept = FC))+
    geom_errorbar(aes(ymin=EstimatedFC_Left, ymax=EstimatedFC_Right), width=.2,) +
    facet_wrap(.~TrueFC, nrow = 6, labeller = labeller(TrueFC = positive_FC_labels)) + theme_linedraw()

ggplot(data = de_dda %>% filter(TrueFC<0), aes(x = as_factor(Protein), y = EstimatedFC))+
    geom_bar(stat="identity")+
    labs(x = "Protein", title = "DDA dataset")+
    geom_hline(data = tibble(TrueFC = negative_FC, FC = negative_FC), aes(yintercept = FC))+
    geom_errorbar(aes(ymin=EstimatedFC_Left, ymax=EstimatedFC_Right), width=.2,) +
    facet_wrap(.~TrueFC, nrow = 6, labeller = labeller(TrueFC = negative_FC_labels)) + theme_linedraw()

ggplot(data = de_dda, aes(x = TrueFC, y = EstimatedFC))+
    geom_point()+
    labs(title = "DDA dataset")+
    geom_abline(intercept = 0, slope = 1)+
    geom_errorbar(aes(ymin=EstimatedFC_Left, ymax=EstimatedFC_Right), width=.2,) +
    facet_wrap(.~Protein)+ theme_linedraw()
