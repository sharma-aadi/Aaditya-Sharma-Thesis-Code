#######################################################################################################################
## list of required packages
#######################################################################################################################
library(devtools, quietly = TRUE)
library(BiocManager, quietly = TRUE)
library(seqinr, quietly = TRUE)          # For reading in fasta files
library(stringr, quietly = TRUE)         # str_sub and str_detect function implementations
library(tidyverse, quietly = TRUE)       # Tidyverse for handling data frames
library(Peptides, quietly = TRUE)        # Add the Peptide Properties
library(viridis, quietly = TRUE)         # viridis color palette
library(ggrepel, quietly = TRUE)         # Repel labels
library(ggpubr, quietly = TRUE)          # Combine multiple ggplots in a single figure
library(UpSetR, quietly = TRUE)          # Base Upset plotter
library(ComplexUpset, quietly = TRUE)    # Fancy Upset plotter
library(ggupset, quietly = TRUE)         # ggplot compatible upset plots
library(gggap, quietly = TRUE)           # break in ggplot axes
library(patchwork, quietly = TRUE)       # Combine multiple ggplots in a single figure
library(msqrob2, quietly = TRUE)         # MSqRob
library(msqrob2, quietly = TRUE)         # MSqRob2
library(reshape2, quietly = TRUE)        # Reshaping and melting tibbles    
library(readxl, quietly = TRUE)          # Read Excel files
library(GGally, quietly = TRUE)          # GGplot extension for pairwise plots 
library(ggVennDiagram, quietly = TRUE)   # GGplot extension for venn diagrams  
library(RColorBrewer, quietly = TRUE)    # Color palettes
library(ComplexHeatmap, quietly = TRUE)
library(rsq, quietly = TRUE)             # Adjusted r squared
library(pheatmap, quietly = TRUE)
library(lme4, quietly = TRUE)
library(sjPlot, quietly = TRUE)
library(sjmisc, quietly = TRUE)
library(scales, quietly = TRUE)
library(ggExtra, quietly = TRUE)
library(ggnewscale, quietly = TRUE)
library(DEP2, quietly = TRUE)
library(VennDiagram)
