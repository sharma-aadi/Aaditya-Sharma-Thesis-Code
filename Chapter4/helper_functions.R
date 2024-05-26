#######################################################################################################################
## Helper function - 1: Read_Proteome
#######################################################################################################################

read_proteome <- function(Proteome_File){
    
    # Parameters
    ## Proteome_File is the file_path to the FASTA file for the spike-in protein sequences
    
    Proteome <- seqinr::read.fasta(file = Proteome_File, seqtype = "AA", as.string = TRUE, 
                                   set.attributes = TRUE, legacy.mode = TRUE, seqonly = FALSE, 
                                   strip.desc = FALSE, whole.header = FALSE, 
                                   bfa = FALSE, sizeof.longlong = .Machine$sizeof.longlong,
                                   endian = .Platform$endian, apply.mask = TRUE)
    
    return(Proteome)
}

#######################################################################################################################
## Helper function - 2: Parse_Proteins()
#######################################################################################################################

parse_proteins <- function(Proteome, Indices = c()){
    
    # Parameters
    ## Proteome is the output list of lists created from a FASTA file by the Read_Proteome function 
    
    ## Indices is the subset of the proteome that is selected for digestion
    ## If no indices are provided, the whole FASTA file is used
    
    if(length(Indices) == 0){
        cat("\nNo indices supplied! Using all available sequences.")
        Indices <- 1:length(Proteome)
    }
    
    Proteins <- vector("list", 2)
    names(Proteins) <- c("Sequence", "ID")
    Protein_Sequences <- Protein_IDs <- character(length = base::length(Indices))
    
    for (i in 1:(base::length(Indices))) {
        Protein_Sequences[i] <- Proteome[Indices[i]][[1]][1]
        Protein_IDs[i] <- stringr::str_replace_all(stringr::str_extract(base::attr(Proteome[Indices[i]],"name"),
                                                                        "\\|(.+)\\|"), "\\|","")
    }
    
    names(Protein_Sequences) <- Protein_IDs
    
    Proteins$Sequence <- Protein_Sequences
    Proteins$ID <- Protein_IDs
    
    return(Proteins)
}

#######################################################################################################################
## Helper function - 3: Substring Position Generator [N-terminal Met and Proline Exceptions included]
#######################################################################################################################

generate_substring_positions <- function(Protein, remove_proline_sites = FALSE){
    
    # Parameters
    ## Protein is a non-empty string containing the entire amino acid sequence of a single protein
    
    Protein_End <- base::nchar(Protein)
    
    if(Protein_End == 0){
        cat("Empty protein sequence!")
        return(invisible(NULL))
    }
    
    ## Get Arginine, Lysine, Proline, and Methionine Sites
    Arginine_Sites <- base::unlist(base::gregexpr(pattern = "R", Protein))
    Lysine_Sites <- base::unlist(base::gregexpr(pattern = "K", Protein))
    Proline_Sites <- base::unlist(base::gregexpr(pattern = "P", Protein))
    Methionine_Sites <- base::unlist(base::gregexpr(pattern = "M", Protein))
    
    ## Remove the "-1" produced by no detection
    Arginine_Sites <- Arginine_Sites[base::which(Arginine_Sites > 0)]
    Lysine_Sites <- Lysine_Sites[base::which(Lysine_Sites > 0)]
    Proline_Sites <- Proline_Sites[base::which(Proline_Sites > 0)]
    Methionine_Sites <- Methionine_Sites[base::which(Methionine_Sites > 0)]
    
    ## Combine to get all cleavage sites
    cleavage_sites <- c(Arginine_Sites, Lysine_Sites)
    
    if(length(Arginine_Sites) > 0){
        names(cleavage_sites)[1:length(Arginine_Sites)] <- "R"
    }
    
    if(length(Lysine_Sites) > 0){
        names(cleavage_sites)[(length(Arginine_Sites) + 1): (length(Lysine_Sites) + length(Arginine_Sites))] <- "K"
    }
    
    ## Add first methionine residue to list of possible cleavage sites (except this should not get counted in missed)
    if(is.element(1,Methionine_Sites)){
        cleavage_sites <- c(1, cleavage_sites)
        names(cleavage_sites)[1] <- "M"
    }
    
    ## mark cleavage sites that are followed by a Proline
    if(base::length(cleavage_sites) > 0 & base::length(Proline_Sites) > 0){
        for(i in 1:base::length(cleavage_sites)){
            if(base::is.element((cleavage_sites[i] + 1), Proline_Sites)){
                names(cleavage_sites)[i] <- "P"
            }
        }
    }
    
    ## remove cleavage sites that are followed by a Proline
    if(remove_proline_sites == TRUE){
        if(length(cleavage_sites) > 0 & length(Proline_Sites) > 0){
            for(i in 1:length(cleavage_sites)){
                if(is.element((cleavage_sites[i] + 1), Proline_Sites)){
                    cleavage_sites <- cleavage_sites[-i]
                }
            }
        }
    }
    
    ## Sort the cleavage sites in ascending order
    if(base::length(cleavage_sites) > 0){
        cleavage_sites <- base::sort(cleavage_sites)
    }
    
    ## A protein terminating R or K is not a cleavage site
    cleavage_sites <- cleavage_sites[which(cleavage_sites != Protein_End)]
    
    ## vector of locations of starts and ends of substrings ("The Tryptic Peptides")
    ## Cleavages "happen" before the protein starts and at the end of the protein
    cleavage_sites <- c(0,cleavage_sites, Protein_End)
    
    names(cleavage_sites)[c(1,length(cleavage_sites))] <- c("Start", "End") 
    
    return(cleavage_sites)
}

#######################################################################################################################
## Helper function - 4: uniprot_statistics()
#######################################################################################################################

protein_statistics <- function(Protein_List){
    
    # Parameters
    ## Proteome_List is the output of the Read_Proteome() function
    
    # Maintain a count of all (non-unique) tryptic peptides generated to truncate list at the end
    count <- 0
    Number_of_Proteins <- length(Protein_List)
    
    print(paste0(Number_of_Proteins, " Proteins in the Database"))
    
    Protein <-  character(length = Number_of_Proteins)
    Sequence_Length <- Cleavage_Sites <- integer(length = Number_of_Proteins)
    Weight <- numeric(length = Number_of_Proteins)
    
    for(i in 1:Number_of_Proteins){
        
        cleavage_sites <- c()
        ## Store Protein Label
        Protein[i] <- str_replace_all(str_extract(attr(Protein_List[i], "names"), "\\|(.+)\\|"), "\\|","") 
        
        ## Store Length of Protein AA Sequence
        Sequence_Length[i] <- nchar(Protein_List[i][[1]][1])
        Weight[i] <- mw(seq = Protein_List[i][[1]][1])
        
        ## Store number of cleavage sites
        Arginine_Sites <- unlist(gregexpr(pattern = "R", Protein_List[i][[1]][1]))
        Lysine_Sites <- unlist(gregexpr(pattern = "K", Protein_List[i][[1]][1]))
        Proline_Sites <- unlist(gregexpr(pattern = "P", Protein_List[i][[1]][1]))
        
        Arginine_Sites <- Arginine_Sites[which(Arginine_Sites > 0)]
        Lysine_Sites <- Lysine_Sites[which(Lysine_Sites > 0)]
        Proline_Sites <- Proline_Sites[which(Proline_Sites > 0)]
        
        ## Combine to get all cleavage sites
        if(length(Arginine_Sites) > 0 || length(Lysine_Sites) > 0){
            cleavage_sites <- sort(c(Arginine_Sites, Lysine_Sites))
        }
        
        Cleavage_Sites[i] <- length(cleavage_sites)
        
    }
    
    Statistics <- tibble(Protein = Protein, AA_Length = Sequence_Length, 
                         Cleavage_Sites,
                         Sites_per_AA = Cleavage_Sites/Sequence_Length, 
                         Molecular_Weight = Weight)
    
    return(Statistics)
}

#######################################################################################################################
## Helper function - 5: generate_tryptic_peptides()
#######################################################################################################################

generate_tryptic_peptides <- function(Protein,
                                      max_missed_cleavages = 0,
                                      remove_proline_sites = FALSE,
                                      methionine_peptides = TRUE,
                                      sequence_context = FALSE){
    
    ## Initialize vector to store tryptic peptides
    Tryptic_Peptides <- character(length = 100000)
    
    MissedCleavages <- integer(length = 100000)
    
    start_position <- integer(length = 100000)
    end_position <- integer(length = 100000)
    
    Number_of_Tryptic_Peptides <- 0
    
    ## Protein is a string of amino acids that is non-empty.
    Protein_End <- as.integer(nchar(Protein))
    
    
    ## vector of locations of starts and ends of substrings ("The Tryptic Peptides")
    substring_positions <- generate_substring_positions(Protein = Protein, remove_proline_sites = remove_proline_sites)
    
    ## Also create a vector of positions starting with M if an N-term Methionine exists ("Methionine peptides)
    if(methionine_peptides == TRUE){
        if(is.element("M", names(substring_positions))){
            m_substring_positions <- substring_positions[2:(min(length(substring_positions),
                                                                (2 + max_missed_cleavages + 1)))]    
        }
    }
    
    ## remove N-term methionine as a cleavage site
    if(is.element("M", names(substring_positions))){
        substring_positions <- substring_positions[-which(names(substring_positions) == "M")]
    } else {
        methionine_peptides <- FALSE
    }
    
    ## Create tryptic peptides
    for(Missed_Cleavages in 0:max_missed_cleavages){
        
        ## Check whether there are enough cleavage sites to "miss"
        if((length(substring_positions) - 2) >= Missed_Cleavages){
            
            ## non-Methionine peptides
            for(j in 1:(length(substring_positions) - Missed_Cleavages - 1)){
                
                ## Define start and end position for the tryptic peptide in the protein AA sequence
                start_pos <- substring_positions[j] + 1
                end_pos <- substring_positions[j + Missed_Cleavages + 1]
                
                ## Update the number of tryptic peptides for the protein
                Number_of_Tryptic_Peptides <- Number_of_Tryptic_Peptides + 1 
                
                Tryptic_Peptides[Number_of_Tryptic_Peptides] <- str_sub(Protein,
                                                                        start = start_pos, 
                                                                        end = end_pos)
                
                ## Add MC annotation
                MissedCleavages[Number_of_Tryptic_Peptides] <- Missed_Cleavages
                
                ## Add sequence context
                start_position[Number_of_Tryptic_Peptides] <- start_pos
                end_position[Number_of_Tryptic_Peptides] <- end_pos
                
            }
        }
        
        ## Methionine peptides
        if(methionine_peptides == TRUE){
            
            if(length(m_substring_positions) > 0){
                
                ## Check whether there are enough cleavage sites to "miss"
                if((length(m_substring_positions) - 2) >= Missed_Cleavages){
                    
                    ## generate Methionine peptides
                    for(j in 1:1){
                        
                        ## Define start and end position for the tryptic peptide in the protein AA sequence
                        start_pos <- m_substring_positions[j] + 1
                        end_pos <- m_substring_positions[j + Missed_Cleavages + 1]
                        
                        ## Update the number of tryptic peptides for the protein
                        Number_of_Tryptic_Peptides <- Number_of_Tryptic_Peptides + 1 
                        
                        Tryptic_Peptides[Number_of_Tryptic_Peptides] <- str_sub(Protein,
                                                                                start = start_pos, 
                                                                                end = end_pos)
                        
                        ## Add MC annotation
                        MissedCleavages[Number_of_Tryptic_Peptides] <- Missed_Cleavages
                        
                        ## Add sequence context
                        start_position[Number_of_Tryptic_Peptides] <- start_pos
                        end_position[Number_of_Tryptic_Peptides] <- end_pos
                        
                    }
                }  
            }
        }
        
        
    }
    
    ## Truncate Length
    length(Tryptic_Peptides) <- Number_of_Tryptic_Peptides
    length(MissedCleavages) <- Number_of_Tryptic_Peptides
    length(start_position) <- Number_of_Tryptic_Peptides
    length(end_position) <- Number_of_Tryptic_Peptides
    
    if(sequence_context == TRUE){
        Tryptic_Peptides <- tibble(Sequence = Tryptic_Peptides, 
                                   Start = start_position,
                                   End = end_position,
                                   MissedCleavages)
        return(Tryptic_Peptides)
    } else {
        return(Tryptic_Peptides)
    }
    
}

#######################################################################################################################
## Helper function - 6: generate_peptide_database()
#######################################################################################################################

## Peptide Database Generator with MC Annotations
generate_peptide_database <- function(Proteome, Protein_Indices = c(),
                                      Maximum_Missed_Cleavages = 0, minimum_length = 1, maximum_length = 100,
                                      remove_proline_sites = FALSE,
                                      methionine_peptides = TRUE){
    
    print(paste0("Obtaining protein sequences."))
    ## Load Sanquin Proteome for input
    Protein_List <- Proteome
    
    ### Protein_List is the output generated by reading in Uniprot Proteome with read.fasta
    ### Maximum_Missed_Cleavages controls the number of missed cleavage sites 
    ### Minimum length is set to 1 Amino Acid
    
    # Maintain a count of all (non-unique) tryptic peptides generated to truncate list at the end
    count <- 0
    
    Sequence <- character(length = 20000000)  #! still hard-coded? 26-01-2022
    Protein <-  character(length = 20000000)
    Start <- integer(length = 20000000)
    End <- integer(length = 20000000)
    MC <- integer(length = 20000000)
    
    Number_of_Proteins <- length(Protein_List)
    Number_of_Tryptic_Peptides <- 0
    
    print(paste0("Generating tryptic peptides from selected proteins..."))
    
    if(length(Protein_Indices) == 0){
        print("No indices supplied! Using all available sequences.")
        Protein_Indices <- 1:length(Protein_List)
        print(paste0("Found ", length(Protein_List), " proteins sequences."))
    }
    
    for(i in Protein_Indices){
        
        tryptic_peptides <- lapply(0:Maximum_Missed_Cleavages,
                                   function(x) generate_tryptic_peptides(Protein = Protein_List[i][[1]][1],
                                                                         max_missed_cleavages = x,
                                                                         remove_proline_sites = remove_proline_sites,
                                                                         methionine_peptides = methionine_peptides,
                                                                         sequence_context = TRUE))
        
        ## Initialize vector containing information about the number of missed cleavages generating the peptide
        mc_vec <- c(rep(0, length(tryptic_peptides[[1]]$Sequence)))
        
        ## Add Missed cleavage annotations for non-zero missed cleavage peptides (if they are generated)
        if(Maximum_Missed_Cleavages > 0){
            for(j in 1:Maximum_Missed_Cleavages){
                mc_vec <- c(mc_vec,
                            rep(j, length(tryptic_peptides[[j+1]]$Sequence) - length(tryptic_peptides[[j]]$Sequence)))
            }
        }
        
        ## Remove Peptides below the minimum length threshold (1 by default)
        mc_vec <- mc_vec[which(nchar(tryptic_peptides[[(Maximum_Missed_Cleavages+1)]]$Sequence) >= minimum_length)]
        
        ## Add sequence context
        start_positions <- tryptic_peptides[[(Maximum_Missed_Cleavages+1)]]$Start[which(
            nchar(tryptic_peptides[[(Maximum_Missed_Cleavages+1)]]$Sequence) >= minimum_length)]
        
        end_positions <- tryptic_peptides[[(Maximum_Missed_Cleavages+1)]]$End[which(
            nchar(tryptic_peptides[[(Maximum_Missed_Cleavages+1)]]$Sequence) >= minimum_length)]
        
        ## Add peptides
        tryptic_peptides <- tryptic_peptides[[(Maximum_Missed_Cleavages+1)]]$Sequence[which(
            nchar(tryptic_peptides[[(Maximum_Missed_Cleavages+1)]]$Sequence) >= minimum_length)]
        
        ## Compute the number of tryptic peptides after length thresholding
        Number_of_Tryptic_Peptides <- length(tryptic_peptides)
        
        if(Number_of_Tryptic_Peptides > 0){
            ## Add Protein Label from the Protein List
            Protein[(count+1):(count + Number_of_Tryptic_Peptides)] <- str_replace_all(
                str_extract(attr(Protein_List[i], "name"), "\\|(.+)\\|"), "\\|","")
            
            ## Add Sequences and corresponding missed cleavages
            Sequence[(count+1):(count + Number_of_Tryptic_Peptides)] <- tryptic_peptides 
            Start[(count+1):(count + Number_of_Tryptic_Peptides)] <- start_positions 
            End[(count+1):(count + Number_of_Tryptic_Peptides)] <- end_positions 
            MC[(count+1):(count + Number_of_Tryptic_Peptides)] <- mc_vec
        }
        
        ## Update the total number of tryptic peptides that have been generated
        count <- count + Number_of_Tryptic_Peptides
    }
    
    ## Truncate the lengths of all the vectors down to the total count of tryptic peptides
    length(Protein) <- length(Sequence) <- length(Start) <- length(End) <- length(MC) <- count
    
    Tryptic_Peptides <- tibble(Sequence, Protein,Start, End, MissedCleavages = MC)
    Tryptic_Peptides <- Tryptic_Peptides %>% group_by(Sequence, Protein, Start, End, MissedCleavages) %>% 
        summarise(Multiplicity = n())
    Tryptic_Peptides <- Tryptic_Peptides %>% dplyr::ungroup()
    Tryptic_Peptides <- Tryptic_Peptides %>% mutate(Length = nchar(Sequence)) %>% filter(Length <= maximum_length)
    
    return(Tryptic_Peptides)
}

#######################################################################################################################
## read_maxquant_table()
#######################################################################################################################

read_maxquant_table <- function(file_name){
    data_table <-  read.delim(file_name, header = TRUE, sep = "\t")
    data_table <- as_tibble(data_table)
    return(data_table)
}
