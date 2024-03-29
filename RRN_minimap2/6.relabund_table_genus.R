library("readr")
library("tidyverse")
library("doBy")
library("reshape2")
library("purrr")

# Genus-level relative abundance tables created from minimap2 outputs for each database, as per Cuscó et al (2019)
# Performed on ONT and PacBio data

###########################
##### 1.FANGORN GTDB ######

GTDB_tax <- read.csv("./minimap_database_tsv_files/taxRep_GTDB_nr.csv", header = F) #path to txaonomy file 
colnames(GTDB_tax) <- c("op", "tax")

#####

# making and defining a function to process each data frame
process_data <- function(df) {
  # Set column names 
  col_names <- c("Query", "Q_length", "Q_start", "Q_end", "Strand", "op", "T_length", "T_start", 
                 "T_end", "N_res_matches", "Align_block", "MapQ", "NM", "ms", "AS", "nn", "P_S") 
  
  df <- df[,1:17] %>%
    setNames(col_names) %>%
    mutate(
      AS = as.numeric(gsub("AS:i:", "", AS)), # getting the alignment scores from AS:i
      per.match = (N_res_matches / Align_block) * 100 # calculating the per match values
    ) %>%
    select(Query, op, N_res_matches, Align_block, AS, MapQ, per.match) # name op as match instead for mirror
  
  # merge df with the taxonomy file #
  df <- left_join(df, GTDB_tax, by ="op") # change to corresponding tax database name
  # changining column names for per.match to Matching and the adding the column name Tax to the merged tax co #
  colnames(df)<-c("Query","op","N_res_matches","Align_block","AS","MapQ","Matching","Tax")
  
  # Filter rows based on Align_block
  df <- subset(df, Align_block > 2999)
  
  # Summarise the data #
  # finds the maximum value of the AS variable for each unique combination of Query, Tax, Matching, and MapQ #
  # essentially picks the highest or the first when ASmax is also the same for Query that have same Matching, Tax and MapQ values #
  df <- summaryBy(AS ~ Query + Tax + Matching + MapQ, data=df, FUN=max)
  
  # selecting the Query-Matches with the highest MapQ value #
  df <- df %>%
    group_by(Query) %>%
    slice_max(MapQ) %>%
    ungroup()
  
  # selecting the Query-Matches with the highest AS.max value #
  df <- df %>%
    group_by(Query) %>%
    slice_max(AS.max) %>%
    ungroup()
  
  # selecting the Query-Matches with the highest Matching value #
    df <- df %>%
    group_by(Query) %>%
    slice_max(Matching) %>%
    ungroup()
    
  # when MapQ, AS.max, and matching are all te same, select the first row of Query with duplicate values #
  df <- df[!duplicated(df$Query), ]
  
  #### for FANGORN GTDB #####
  # Clean and organise the table 
  df <- df %>%
    # Split the Taxon column into separate taxonomic levels
    separate(Tax, into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep="\\|") %>%
    # Remove the prefix from each taxonomic level
    mutate_at(vars(Kingdom:Genus), ~substr(., 4, nchar(.))) %>%
    # getting counts for all the alignments at genus level #
    group_by(Genus) %>%
    summarise(Counts = n())
  
  # Reshape dataframe to wide format
  df <- df %>%
    pivot_wider(names_from = Genus, values_from = Counts, values_fill = 0)
  
  ## changing counts to relative abundance ##
  df <- df/rowSums(df)*100
  rowSums(df) ## to check if each sample adds up to a 100 ##
  
  # Convert tf from wide back to long format
  df <- df %>%
    pivot_longer(cols = everything(), names_to = "Genus", values_to = "Rel_abundance")
  
  return(df)
}

##########################################
### FANGORN_GTDB_RRN - minimap files ###
##########################################

# importing all paf files that were taxonomically classified using FANGORN GTDB_nrRep - minimap2

# Set your working directory to where your files are
setwd("./minimap2/tax_assign/GTDB_Fangorn/")

# going per mock community to reduce computational load #
# Get list of .paf file names per sample #

######### MCAP moc com #######

GTDB_APC_file_list <- list.files(pattern = "APC_.*_nrRep\\.paf$")

# Read all files
GTDB_APC_data_list <- map(GTDB_APC_file_list, ~read.delim(., header = F))

# Apply the function to each data frame in the list
GTDB_APC_processed_data_list <- lapply(GTDB_APC_data_list, process_data)

# making one large dataframe for all the files and adding file names as a column #
# Get file names without the extension
GTDB_APC_file_names <- tools::file_path_sans_ext(GTDB_APC_file_list)

# Apply the function to each data frame and each file name
GTDB_APC_combined_data <- purrr::map2_df(GTDB_APC_processed_data_list, GTDB_APC_file_names, ~cbind(.x, FileName = .y))

######### MCGD moc com #######

GTDB_gDNA_file_list <- list.files(pattern = "gDNA_.*_nrRep\\.paf$")

# Read all files
GTDB_gDNA_data_list <- map(GTDB_gDNA_file_list, ~read.delim(., header = F))

# Apply the function to each data frame in the list
GTDB_gDNA_processed_data_list <- lapply(GTDB_gDNA_data_list, process_data)

# making one large dataframe for all the files and adding file names as a column #
# Get file names without the extension
GTDB_gDNA_file_names <- tools::file_path_sans_ext(GTDB_gDNA_file_list)

# Apply the function to each data frame and each file name
GTDB_gDNA_combined_data <- purrr::map2_df(GTDB_gDNA_processed_data_list, GTDB_gDNA_file_names, ~cbind(.x, FileName = .y))

######### ATCC moc com #######

GTDB_ATCC_file_list <- list.files(pattern = "ATCC_.*_nrRep\\.paf$")

# Read all files
GTDB_ATCC_data_list <- map(GTDB_ATCC_file_list, ~read.delim(., header = F))

# Apply the function to each data frame in the list
GTDB_ATCC_processed_data_list <- lapply(GTDB_ATCC_data_list, process_data)

# making one large dataframe for all the files and adding file names as a column #
# Get file names without the extension
GTDB_ATCC_file_names <- tools::file_path_sans_ext(GTDB_ATCC_file_list)

# Apply the function to each data frame and each file name
GTDB_ATCC_combined_data <- purrr::map2_df(GTDB_ATCC_processed_data_list, GTDB_ATCC_file_names, ~cbind(.x, FileName = .y))

######### Zymo moc com #######

GTDB_Zymo_file_list <- list.files(pattern = "Zymo_.*_nrRep\\.paf$")

# Read all files
GTDB_Zymo_data_list <- map(GTDB_Zymo_file_list, ~read.delim(., header = F))

# Apply the function to each data frame in the list
GTDB_Zymo_processed_data_list <- lapply(GTDB_Zymo_data_list, process_data)

# making one large dataframe for all the files and adding file names as a column #
# Get file names without the extension
GTDB_Zymo_file_names <- tools::file_path_sans_ext(GTDB_Zymo_file_list)

# Apply the function to each data frame and each file name
GTDB_Zymo_combined_data <- purrr::map2_df(GTDB_Zymo_processed_data_list, GTDB_Zymo_file_names, ~cbind(.x, FileName = .y))

######### Neg moc com #######

GTDB_Neg_file_list <- list.files(pattern = "Neg_.*_nrRep\\.paf$")

# Read all files
GTDB_Neg_data_list <- map(GTDB_Neg_file_list, ~read.delim(., header = F))

# Apply the function to each data frame in the list
GTDB_Neg_processed_data_list <- lapply(GTDB_Neg_data_list, process_data)

# making one large dataframe for all the files and adding file names as a column #
# Get file names without the extension
GTDB_Neg_file_names <- tools::file_path_sans_ext(GTDB_Neg_file_list)

# Apply the function to each data frame and each file name
GTDB_Neg_combined_data <- purrr::map2_df(GTDB_Neg_processed_data_list, GTDB_Neg_file_names, ~cbind(.x, FileName = .y))

# rbinding the five sample together for the FANGORN GTDB files #
GTDB_minimap <- rbind(GTDB_APC_combined_data, GTDB_gDNA_combined_data, GTDB_ATCC_combined_data, GTDB_Zymo_combined_data, GTDB_Neg_combined_data)

# export the combined file #
write.csv(GTDB_minimap, "./R_files/output_files_genus/GTDB_minimap_RRN_ONT_genus.csv", row.names = FALSE)

##################################################################################################################################################

#############################
##### 2. FANGORN RefSeq #####

RefSeq_tax <- read.csv("./minimap_database_tsv_files/taxRep_RefSeq_nr.csv", header = F) #path to txaonomy file 
colnames(RefSeq_tax) <- c("op", "tax")

#####

# making and defining a function to process each data frame
process_data <- function(df) {
  # Set column names 
  col_names <- c("Query", "Q_length", "Q_start", "Q_end", "Strand", "op", "T_length", "T_start", 
                 "T_end", "N_res_matches", "Align_block", "MapQ", "NM", "ms", "AS", "nn", "P_S") 
  
  df <- df[,1:17] %>%
    setNames(col_names) %>%
    mutate(
      AS = as.numeric(gsub("AS:i:", "", AS)), # getting the alignment scores from AS:i
      per.match = (N_res_matches / Align_block) * 100 # calculating the per match values
    ) %>%
    select(Query, op, N_res_matches, Align_block, AS, MapQ, per.match) # name op as match instead for mirror
  
  # merge df with the taxonomy file #
  df <- left_join(df, RefSeq_tax, by ="op") # change to corresponding tax database name
  # changining column names for per.match to Matching and the adding the column name Tax to the merged tax co #
  colnames(df)<-c("Query","op","N_res_matches","Align_block","AS","MapQ","Matching","Tax")
  
  # Filter rows based on Align_block
  df <- subset(df, Align_block > 2999)
  
  # Summarise the data #
  # finds the maximum value of the AS variable for each unique combination of Query, Tax, Matching, and MapQ #
  # essentially picks the highest or the first when ASmax is also the same for Query that have same Matching, Tax and MapQ values #
  df <- summaryBy(AS ~ Query + Tax + Matching + MapQ, data=df, FUN=max)
  
  # selecting the Query-Matches with the highest MapQ value #
  df <- df %>%
    group_by(Query) %>%
    slice_max(MapQ) %>%
    ungroup()
  
  # selecting the Query-Matches with the highest AS.max value #
  df <- df %>%
    group_by(Query) %>%
    slice_max(AS.max) %>%
    ungroup()
  
  # selecting the Query-Matches with the highest Matching value #
    df <- df %>%
    group_by(Query) %>%
    slice_max(Matching) %>%
    ungroup()
    
  # when MapQ, AS.max, and matching are all te same, select the first row of Query with duplicate values #
  df <- df[!duplicated(df$Query), ]
  
  #### for FANGORN RefSeq #####
  # Clean and organise the table 
  df <- df %>%
    # Split the Taxon column into separate taxonomic levels
    separate(Tax, into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep="\\|") %>%
    # Remove the prefix from each taxonomic level
    mutate_at(vars(Kingdom:Genus), ~substr(., 4, nchar(.))) %>%
    # getting counts for all the alignments at genus level #
    group_by(Genus) %>%
    summarise(Counts = n())
  
  # Reshape dataframe to wide format
  df <- df %>%
    pivot_wider(names_from = Genus, values_from = Counts, values_fill = 0)
  
  ## changing counts to relative abundance ##
  df <- df/rowSums(df)*100
  rowSums(df) ## to check if each sample adds up to a 100 ##
  
  # Convert tf from wide back to long format
  df <- df %>%
    pivot_longer(cols = everything(), names_to = "Genus", values_to = "Rel_abundance")
  
  return(df)
}

##########################################
### FANGORN_RefSeq_RRN - minimap files ###
##########################################

# importing all paf files that were taxonomically classified using FANGORN GTDB_nrRep - minimap2
# Set your working directory to where your files are
setwd("./minimap2/tax_assign/RefSeq_Fangorn/")

# going per mocl community to reduce computational load #
# Get list of .paf file names per sample #


######### MCAP moc com #######

RefSeq_APC_file_list <- list.files(pattern = "APC_.*_nrRep\\.paf$")

# Read all files
RefSeq_APC_data_list <- map(RefSeq_APC_file_list, ~read.delim(., header = F))

# Apply the function to each data frame in the list
RefSeq_APC_processed_data_list <- lapply(RefSeq_APC_data_list, process_data)

# making one large dataframe for all the files and adding file names as a column #
# Get file names without the extension
RefSeq_APC_file_names <- tools::file_path_sans_ext(RefSeq_APC_file_list)

# Apply the function to each data frame and each file name
RefSeq_APC_combined_data <- purrr::map2_df(RefSeq_APC_processed_data_list, RefSeq_APC_file_names, ~cbind(.x, FileName = .y))

######### MCGD moc com #######

RefSeq_gDNA_file_list <- list.files(pattern = "gDNA_.*_nrRep\\.paf$")

# Read all files
RefSeq_gDNA_data_list <- map(RefSeq_gDNA_file_list, ~read.delim(., header = F))

# Apply the function to each data frame in the list
RefSeq_gDNA_processed_data_list <- lapply(RefSeq_gDNA_data_list, process_data)

# making one large dataframe for all the files and adding file names as a column #
# Get file names without the extension
RefSeq_gDNA_file_names <- tools::file_path_sans_ext(RefSeq_gDNA_file_list)

# Apply the function to each data frame and each file name
RefSeq_gDNA_combined_data <- purrr::map2_df(RefSeq_gDNA_processed_data_list, RefSeq_gDNA_file_names, ~cbind(.x, FileName = .y))

######### ATCC moc com #######

RefSeq_ATCC_file_list <- list.files(pattern = "ATCC_.*_nrRep\\.paf$")

# Read all files
RefSeq_ATCC_data_list <- map(RefSeq_ATCC_file_list, ~read.delim(., header = F))

# Apply the function to each data frame in the list
RefSeq_ATCC_processed_data_list <- lapply(RefSeq_ATCC_data_list, process_data)

# making one large dataframe for all the files and adding file names as a column #
# Get file names without the extension
RefSeq_ATCC_file_names <- tools::file_path_sans_ext(RefSeq_ATCC_file_list)

# Apply the function to each data frame and each file name
RefSeq_ATCC_combined_data <- purrr::map2_df(RefSeq_ATCC_processed_data_list, RefSeq_ATCC_file_names, ~cbind(.x, FileName = .y))

######### Zymo moc com #######

RefSeq_Zymo_file_list <- list.files(pattern = "Zymo_.*_nrRep\\.paf$")

# Read all files
RefSeq_Zymo_data_list <- map(RefSeq_Zymo_file_list, ~read.delim(., header = F))

# Apply the function to each data frame in the list
RefSeq_Zymo_processed_data_list <- lapply(RefSeq_Zymo_data_list, process_data)

# making one large dataframe for all the files and adding file names as a column #
# Get file names without the extension
RefSeq_Zymo_file_names <- tools::file_path_sans_ext(RefSeq_Zymo_file_list)

# Apply the function to each data frame and each file name
RefSeq_Zymo_combined_data <- purrr::map2_df(RefSeq_Zymo_processed_data_list, RefSeq_Zymo_file_names, ~cbind(.x, FileName = .y))

######### Neg moc com #######

RefSeq_Neg_file_list <- list.files(pattern = "Neg_.*_nrRep\\.paf$")

# Read all files
RefSeq_Neg_data_list <- map(RefSeq_Neg_file_list, ~read.delim(., header = F))

# Apply the function to each data frame in the list
RefSeq_Neg_processed_data_list <- lapply(RefSeq_Neg_data_list, process_data)

# making one large dataframe for all the files and adding file names as a column #
# Get file names without the extension
RefSeq_Neg_file_names <- tools::file_path_sans_ext(RefSeq_Neg_file_list)

# Apply the function to each data frame and each file name
RefSeq_Neg_combined_data <- purrr::map2_df(RefSeq_Neg_processed_data_list, RefSeq_Neg_file_names, ~cbind(.x, FileName = .y))

# rbinding the five sample together for the FANGORN RefSeq files #
RefSeq_minimap <- rbind(RefSeq_APC_combined_data, RefSeq_gDNA_combined_data, RefSeq_ATCC_combined_data, RefSeq_Zymo_combined_data, RefSeq_Neg_combined_data)

# export the combined file #
write.csv(RefSeq_minimap, "./R_files/output_files_genus/RefSeq_minimap_RRN_ONT_genus.csv", row.names = FALSE)

#############################################################################################################################################################

##########################
##### 3.rrn_DBv2 #########

rrnDB_tax <- read.csv("./minimap_database_tsv_files/rrn_DBv2_tax.csv", header = F) #path to txaonomy file 
rrnDB_tax <- rrnDB_tax[,-(2:3)]
colnames(rrnDB_tax) <- c("op", "Class", "Order", "Family", "Genus", "Species")
rrnDB_tax <- rrnDB_tax %>% select("op", "Genus")

#####

# making and defining a function to process each data frame
process_data <- function(df) {
  # Set column names
  col_names <- c("Query", "Q_length", "Q_start", "Q_end", "Strand", "op", "T_length", "T_start", 
                 "T_end", "N_res_matches", "Align_block", "MapQ", "NM", "ms", "AS", "nn", "P_S") 
  
  df <- df[,1:17] %>%
    setNames(col_names) %>%
    mutate(
      AS = as.numeric(gsub("AS:i:", "", AS)), # getting the alignment scores from AS:i
      per.match = (N_res_matches / Align_block) * 100 # calculating the per match values
    ) %>%
    select(Query, op, N_res_matches, Align_block, AS, MapQ, per.match) # name op as match instead for mirror
    
  # merge df with the taxonomy file #
  df <- left_join(df, rrnDB_tax, by ="op") # change to corresponding tax database name
  # changining column names for per.match to Matching and the adding the column name Tax to the merged tax co #
  colnames(df)<-c("Query","op","N_res_matches","Align_block","AS","MapQ","Matching","Tax")
  
  # Filter rows based on Align_block
  df <- subset(df, Align_block > 2999)
  
  # Summarise the data #
  # finds the maximum value of the AS variable for each unique combination of Query, Tax, Matching, and MapQ #
  # essentially picks the highest or the first when ASmax is also the same for Query that have same Matching, Tax and MapQ values #
  df <- summaryBy(AS ~ Query + Tax + Matching + MapQ, data=df, FUN=max)
  
  # selecting the Query-Matches with the highest MapQ value #
  df <- df %>%
    group_by(Query) %>%
    slice_max(MapQ) %>%
    ungroup()
  
  # selecting the Query-Matches with the highest AS.max value #
  df <- df %>%
    group_by(Query) %>%
    slice_max(AS.max) %>%
    ungroup()
  
  # selecting the Query-Matches with the highest Matching value #
    df <- df %>%
    group_by(Query) %>%
    slice_max(Matching) %>%
    ungroup()
  
  # when MapQ, AS.max, and matching are all te same, select the first row of Query with duplicate values #
  df <- df[!duplicated(df$Query), ]
  
  #### for rrn_DBv2 #####
  # Clean and organise the table
  df <- df %>%
  # getting counts for all the alignments at genus level #
  group_by(Tax) %>%
    summarise(Counts = n())
  
  # Reshape dataframe to wide format
  df <- df %>%
    pivot_wider(names_from = Tax, values_from = Counts, values_fill = 0)
  
  ## changing counts to relative abundance ##
  df <- df/rowSums(df)*100
  rowSums(df) ## to check if each sample adds up to a 100 ##
  
  # Convert tf from wide back to long format
  df <- df %>%
    pivot_longer(cols = everything(), names_to = "Genus", values_to = "Rel_abundance")
  
  return(df)
}

##########################################
### rrn_DBv2_RRN - minimap files ###
##########################################

# importing all paf files that were taxonomically classified using FANGORN GTDB_nrRep - minimap2
# Set your working directory to where your files are
setwd("./RRN_mc_ONT/minimap2/tax_assign/rrn_DBv2/")

# going per mock community to reduce computational load #
# Get list of .paf file names per sample #

######### MCAP 24 st moc com #######

rrnDB_APC_file_list <- list.files(pattern = "APC_.*_rrn_DBv2\\.paf$")

# Read all files
rrnDB_APC_data_list <- map(rrnDB_APC_file_list, ~read.delim(., header = F))

# Apply the function to each data frame in the list
rrnDB_APC_processed_data_list <- lapply(rrnDB_APC_data_list, process_data)

# making one large dataframe for all the files and adding file names as a column #
# Get file names without the extension
rrnDB_APC_file_names <- tools::file_path_sans_ext(rrnDB_APC_file_list)

# Apply the function to each data frame and each file name
rrnDB_APC_combined_data <- purrr::map2_df(rrnDB_APC_processed_data_list, rrnDB_APC_file_names, ~cbind(.x, FileName = .y))


######### MCGD moc com #######

rrnDB_gDNA_file_list <- list.files(pattern = "gDNA_.*_rrn_DBv2\\.paf$")

# Read all files
rrnDB_gDNA_data_list <- map(rrnDB_gDNA_file_list, ~read.delim(., header = F))

# Apply the function to each data frame in the list
rrnDB_gDNA_processed_data_list <- lapply(rrnDB_gDNA_data_list, process_data)

# making one large dataframe for all the files and adding file names as a column #
# Get file names without the extension
rrnDB_gDNA_file_names <- tools::file_path_sans_ext(rrnDB_gDNA_file_list)

# Apply the function to each data frame and each file name
rrnDB_gDNA_combined_data <- purrr::map2_df(rrnDB_gDNA_processed_data_list, rrnDB_gDNA_file_names, ~cbind(.x, FileName = .y))

######### ATCC moc com #######

rrnDB_ATCC_file_list <- list.files(pattern = "ATCC_.*_rrn_DBv2\\.paf$")

# Read all files
rrnDB_ATCC_data_list <- map(rrnDB_ATCC_file_list, ~read.delim(., header = F))

# Apply the function to each data frame in the list
rrnDB_ATCC_processed_data_list <- lapply(rrnDB_ATCC_data_list, process_data)

# making one large dataframe for all the files and adding file names as a column #
# Get file names without the extension
rrnDB_ATCC_file_names <- tools::file_path_sans_ext(rrnDB_ATCC_file_list)

# Apply the function to each data frame and each file name
rrnDB_ATCC_combined_data <- purrr::map2_df(rrnDB_ATCC_processed_data_list, rrnDB_ATCC_file_names, ~cbind(.x, FileName = .y))

######### Zymo moc com #######

rrnDB_Zymo_file_list <- list.files(pattern = "Zymo_.*_rrn_DBv2\\.paf$")

# Read all files
rrnDB_Zymo_data_list <- map(rrnDB_Zymo_file_list, ~read.delim(., header = F))

# Apply the function to each data frame in the list
rrnDB_Zymo_processed_data_list <- lapply(rrnDB_Zymo_data_list, process_data)

# making one large dataframe for all the files and adding file names as a column #
# Get file names without the extension
rrnDB_Zymo_file_names <- tools::file_path_sans_ext(rrnDB_Zymo_file_list)

# Apply the function to each data frame and each file name
rrnDB_Zymo_combined_data <- purrr::map2_df(rrnDB_Zymo_processed_data_list, rrnDB_Zymo_file_names, ~cbind(.x, FileName = .y))

######### Neg moc com #######

rrnDB_Neg_file_list <- list.files(pattern = "Neg_.*_rrn_DBv2\\.paf$")

# Read all files
rrnDB_Neg_data_list <- map(rrnDB_Neg_file_list, ~read.delim(., header = F))

# Apply the function to each data frame in the list
rrnDB_Neg_processed_data_list <- lapply(rrnDB_Neg_data_list, process_data)

# making one large dataframe for all the files and adding file names as a column #
# Get file names without the extension
rrnDB_Neg_file_names <- tools::file_path_sans_ext(rrnDB_Neg_file_list)

# Apply the function to each data frame and each file name
rrnDB_Neg_combined_data <- purrr::map2_df(rrnDB_Neg_processed_data_list, rrnDB_Neg_file_names, ~cbind(.x, FileName = .y))


# rbinding the five sample together for the rrnDB files #
rrnDB_minimap <- rbind(rrnDB_APC_combined_data, rrnDB_gDNA_combined_data, rrnDB_ATCC_combined_data, rrnDB_Zymo_combined_data, rrnDB_Neg_combined_data)

# export the combined file #
write.csv(rrnDB_minimap, "./R_files/output_files_genus/rrnDB_minimap_RRN_ONT_genus.csv", row.names = FALSE)

#######################################################################################################################################################

##########################
##### 4.MIrROR ###########

mirror_tax <- read.csv("./minimap_database_tsv_files/MIrROR_DB_r01.csv", header = T) #path to txaonomy file 
mirror_tax <- mirror_tax %>% select(X.Accession, gtdb_taxonomy_g)
colnames(mirror_tax) <- c("op", "Genus")

#####

# making and defining a function to process each data frame
process_data <- function(df) {
  # Set column names
  col_names <- c("Query", "Q_length", "Q_start", "Q_end", "Strand", "match", "T_length", "T_start", 
                 "T_end", "N_res_matches", "Align_block", "MapQ", "NM", "ms", "AS", "nn", "P_S") 
  
  df <- df[,1:17] %>%
    setNames(col_names) %>%
    mutate(
      AS = as.numeric(gsub("AS:i:", "", AS)), # getting the alignment scores from AS:i
      per.match = (N_res_matches / Align_block) * 100 # calculating the per match values
    ) %>%
    select(Query, match, N_res_matches, Align_block, AS, MapQ, per.match) # name op as match instead for mirror
  
  ## only for mirror database (this followed by the rrnDB amd mirror bit below) ##
  ## getting only the first and second undrescore parts of the matched sequence/accession number to merge later ##
  df$op1 <- sapply(strsplit(df$match, "_"), "[", 1)
  df$op2 <- sapply(strsplit(df$match, "_"), "[", 2)
  ## now joining the two columns op_1 and op_2 ##
  df$op <- paste(df$op1, df$op2, sep = "_")
  ## removing the columns match, op1 and op 2 ##
  df <- df %>% select(-match, -op1, -op2)
  ## rearraging the order of the columns since we changed it to bring op to the right of Query ##
  df <- select(df, "Query", "op", "N_res_matches","Align_block","AS","MapQ","per.match")
  #############################################
    
  # merge df with the taxonomy file #
  df <- left_join(df, mirror_tax, by ="op") # change to corresponding tax database name
  # changining column names for per.match to Matching and the adding the column name Tax to the merged tax co #
  colnames(df)<-c("Query","op","N_res_matches","Align_block","AS","MapQ","Matching","Tax")
  
  # Filter rows based on Align_block
  df <- subset(df, Align_block > 2999)
  
  # Summarise the data #
  # finds the maximum value of the AS variable for each unique combination of Query, Tax, Matching, and MapQ #
  # essentially picks the highest or the first when ASmax is also the same for Query that have same Matching, Tax and MapQ values #
  df <- summaryBy(AS ~ Query + Tax + Matching + MapQ, data=df, FUN=max)
  
  # selecting the Query-Matches with the highest MapQ value #
  df <- df %>%
    group_by(Query) %>%
    slice_max(MapQ) %>%
    ungroup()
  
  # selecting the Query-Matches with the highest AS.max value #
  df <- df %>%
    group_by(Query) %>%
    slice_max(AS.max) %>%
    ungroup()
  
  # selecting the Query-Matches with the highest Matching value #
    df <- df %>%
    group_by(Query) %>%
    slice_max(Matching) %>%
    ungroup()
  
  # when MapQ, AS.max, and matching are all te same, select the first row of Query with duplicate values #
  df <- df[!duplicated(df$Query), ]
  
  #### for mirror #####
  # Clean and organise the table
  df <- df %>%
  # getting counts for all the alignments at genus level #
  group_by(Tax) %>%
    summarise(Counts = n())
  
  # Reshape dataframe to wide format
  df <- df %>%
    pivot_wider(names_from = Tax, values_from = Counts, values_fill = 0)
  
  ## changing counts to relative abundance ##
  df <- df/rowSums(df)*100
  rowSums(df) ## to check if each sample adds up to a 100 ##
  
  # Convert tf from wide back to long format
  df <- df %>%
    pivot_longer(cols = everything(), names_to = "Genus", values_to = "Rel_abundance")
  
  return(df)
}

##########################################
### MIrROR_RRN - minimap files ###
##########################################

# importing all paf files that were taxonomically classified using FANGORN GTDB_nrRep - minimap2
# Set your working directory to where your files are
setwd("./RRN_mc_ONT/minimap2/tax_assign/mirror/")

# going per mock community to reduce computational load #
# Get list of .paf file names per sample #

######### MCAP 24 st moc com #######
mirror_APC_file_list <- list.files(pattern = "APC_.*_mirror\\.paf$")

# Read all files
mirror_APC_data_list <- map(mirror_APC_file_list, ~read.delim(., header = F))

# Apply the function to each data frame in the list
mirror_APC_processed_data_list <- lapply(mirror_APC_data_list, process_data)

# making one large dataframe for all the files and adding file names as a column #
# Get file names without the extension
mirror_APC_file_names <- tools::file_path_sans_ext(mirror_APC_file_list)

# Apply the function to each data frame and each file name
mirror_APC_combined_data <- purrr::map2_df(mirror_APC_processed_data_list, mirror_APC_file_names, ~cbind(.x, FileName = .y))

######### MCGD moc com #######
mirror_gDNA_file_list <- list.files(pattern = "gDNA_.*_mirror\\.paf$")

# Read all files
mirror_gDNA_data_list <- map(mirror_gDNA_file_list, ~read.delim(., header = F))

# Apply the function to each data frame in the list
mirror_gDNA_processed_data_list <- lapply(mirror_gDNA_data_list, process_data)

# making one large dataframe for all the files and adding file names as a column #
# Get file names without the extension
mirror_gDNA_file_names <- tools::file_path_sans_ext(mirror_gDNA_file_list)

# Apply the function to each data frame and each file name
mirror_gDNA_combined_data <- purrr::map2_df(mirror_gDNA_processed_data_list, mirror_gDNA_file_names, ~cbind(.x, FileName = .y))

######### ATCC moc com #######
mirror_ATCC_file_list <- list.files(pattern = "ATCC_.*_mirror\\.paf$")

# Read all files
mirror_ATCC_data_list <- map(mirror_ATCC_file_list, ~read.delim(., header = F))

# Apply the function to each data frame in the list
mirror_ATCC_processed_data_list <- lapply(mirror_ATCC_data_list, process_data)

# making one large dataframe for all the files and adding file names as a column #
# Get file names without the extension
mirror_ATCC_file_names <- tools::file_path_sans_ext(mirror_ATCC_file_list)

# Apply the function to each data frame and each file name
mirror_ATCC_combined_data <- purrr::map2_df(mirror_ATCC_processed_data_list, mirror_ATCC_file_names, ~cbind(.x, FileName = .y))

######### Zymo moc com #######
mirror_Zymo_file_list <- list.files(pattern = "Zymo_.*_mirror\\.paf$")

# Read all files
mirror_Zymo_data_list <- map(mirror_Zymo_file_list, ~read.delim(., header = F))

# Apply the function to each data frame in the list
mirror_Zymo_processed_data_list <- lapply(mirror_Zymo_data_list, process_data)

# making one large dataframe for all the files and adding file names as a column #
# Get file names without the extension
mirror_Zymo_file_names <- tools::file_path_sans_ext(mirror_Zymo_file_list)

# Apply the function to each data frame and each file name
mirror_Zymo_combined_data <- purrr::map2_df(mirror_Zymo_processed_data_list, mirror_Zymo_file_names, ~cbind(.x, FileName = .y))

######### Neg moc com #######
mirror_Neg_file_list <- list.files(pattern = "Neg_.*_mirror\\.paf$")

# Read all files
mirror_Neg_data_list <- map(mirror_Neg_file_list, ~read.delim(., header = F))

# Apply the function to each data frame in the list
mirror_Neg_processed_data_list <- lapply(mirror_Neg_data_list, process_data)

# making one large dataframe for all the files and adding file names as a column #
# Get file names without the extension
mirror_Neg_file_names <- tools::file_path_sans_ext(mirror_Neg_file_list)

# Apply the function to each data frame and each file name
mirror_Neg_combined_data <- purrr::map2_df(mirror_Neg_processed_data_list, mirror_Neg_file_names, ~cbind(.x, FileName = .y))

# rbinding the five sample together for the mirror files #
mirror_minimap <- rbind(mirror_APC_combined_data, mirror_gDNA_combined_data, mirror_ATCC_combined_data, mirror_Zymo_combined_data, mirror_Neg_combined_data)

# export the combined file #
write.csv(mirror_minimap, "./R_files/output_files_genus/mirror_minimap_RRN_ONT_genus.csv", row.names = FALSE)
