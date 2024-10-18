library("readr")
library("tidyverse")
library("doBy")
library("reshape2")
library("purrr")

# Species-level relative abundance tables created from minimap2 outputs for each database, as per Cuscó et al (2019)
# Performed on ONT and PacBio data

############################
##### 1. FANGORN GTDB ######

GTDB_tax <- read.csv("./minimap_database_tsv_files/taxRep_GTDB_nr.csv", header = F) #path to taxonomy files
colnames(GTDB_tax) <- c("op", "tax")

#####

# making and defining a function to process each data frame
process_data <- function(df) {
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
  
  # getting all the unique hits out for a query sequence
  df_uni <- df %>% filter(MapQ > 0)
  df_uni_ids <- df_uni$Query
  
  # filtering out all the query sequneces with unique hits from df that alos have non unique hits for them
  df_uni_mapG0 <- df %>% filter(!(Query %in% df_uni_ids & MapQ == 0))
  # in the non unique hits first filter by AS the ones with the highest alignment score is selected
  # first getting out all the query that have MapQ = 0
  df_map0 <- df_uni_mapG0 %>% filter(MapQ == 0)
  
  ## now looking at the the queries with multiple hits so MapQ = 0 with the AS score are different i.e. there is one max score comapred to the others,     #selecting based on that max 
  df_map0_dAS <- df_map0 %>% 
  group_by(Query) %>%
  slice_max(AS, n = 1) %>%
  ungroup() # This will leave the rows that have identical AS for a query
  # get only the the one hist when mapQ = 0 and AS is different for one query
  df_map0_dAS2 <- df_map0_dAS %>%
  group_by(Query, AS) %>%
  filter(!(n() > 1)) %>%
  ungroup()
  # rbind the unique hit where mapq > 0 (uni) and for the hits where mapQ = 0 but AS was diff with one higher than the other (dAS2)
  df_uni_dAS2 <- rbind(df_uni, df_map0_dAS2)
  # cleaning up the tax column 
  df_uni_dAS2 <- df_uni_dAS2 %>%
    # Split the Taxon column into separate taxonomic levels
    separate(Tax, into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep="\\|") %>%
    # Remove the prefix from each taxonomic level
    mutate_at(vars(Kingdom:Species), ~substr(., 4, nchar(.))) 
    # removing the _letters in species names
    df_uni_dAS2$Species <- gsub("_[A-Z]", "", df_uni_dAS2$Species) 
    # adding an lca column 
    df_uni_dAS2$lca <- df_uni_dAS2$Species 
  
  ## pull out the reads that have equal AS scores, i.e. multiple rows for the same AS
  df_map0_sAS <- df_map0_dAS %>%
  group_by(Query, AS) %>%
  filter(n() > 1) %>%
  ungroup()
  # cleaning up the tax column 
  df_map0_sAS <- df_map0_sAS %>%
    # Split the Taxon column into separate taxonomic levels
    separate(Tax, into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep="\\|") %>%
    # Remove the prefix from each taxonomic level
    mutate_at(vars(Kingdom:Species), ~substr(., 4, nchar(.))) 
  # if the species is the same for the mutiple hits keep only the first entry
    df_map0_sAS_sp <- df_map0_sAS %>%
    group_by(Query, Species) %>%
    slice_head(n = 1) %>%
    ungroup()  # so queries with (i) mapQ = 0 and same AS that all have hits to the same species only the first hit is taken, and (ii) all hits for mapQ = 0 and same AS that have different hits at species level as also kept
    ## identifying (i) and (ii) from df_min0_sAS_sp and then adding lca 
    # (i) pulling out the hits that have mapQ = 0 and same AS that all have hits to the same species only the first hit has been taken
    df_map0_sAS_1sp <- df_map0_sAS_sp %>%
    group_by(Query) %>%
    filter(!(n() > 1)) %>%
    ungroup()
    # adding an lca column to it 
    df_map0_sAS_1sp$lca <- df_map0_sAS_1sp$Species
    # (ii) pulling out hits for mapQ = 0 and same AS that have different hits at species level as also kept
    df_map0_sAS_msp <- df_map0_sAS_sp %>%
    group_by(Query) %>%
    filter(n() > 1) %>%
    ungroup()
    # removing the _letters in species names
    df_map0_sAS_msp$Species <- gsub("_[A-Z]", "", df_map0_sAS_msp$Species) 
    # now lets add the lowest common acestor for them 
    df_map0_sAS_msp <- df_map0_sAS_msp %>%
    group_by(Query) %>%
    mutate(
    lca = case_when(
      n_distinct(Genus) == 1 ~ Genus,
      n_distinct(Family) == 1 ~ Family,
      n_distinct(Order) == 1 ~ Order,
      n_distinct(Class) == 1 ~ Class,
      n_distinct(Phylum) == 1 ~ Phylum,
      n_distinct(Kingdom) == 1 ~ Kingdom,
      TRUE ~ NA_character_  # If no common taxonomic level is found
    )
  ) %>%
  ungroup() %>%
  # keep only one hit now that lca has been added for muyltiple hits for one query
  group_by(Query) %>%
  slice_head(n = 1) %>%
  ungroup()
    # rbinding the map0_sAS dups nad no dups 
    df_m0_sAS <- rbind(df_map0_sAS_1sp, df_map0_sAS_msp)
    
    ## rbinding it all now 
    df2 <- rbind(df_uni_dAS2, df_m0_sAS)
    
# getting counts for all the alignments at species level #
    df2 <- df2 %>% group_by(lca) %>%
          summarise(Counts = n())
  
  # Reshape dataframe to wide format
  df2 <- df2 %>%
    pivot_wider(names_from = lca, values_from = Counts, values_fill = 0)
  
  ## changing counts to relative abundance ##
  df2 <- df2/rowSums(df2)*100
  rowSums(df2) ## to check if each sample adds up to a 100 ##
  
  # Convert tf from wide back to long format
  df2 <- df2 %>%
    pivot_longer(cols = everything(), names_to = "lca_sp", values_to = "Rel_abundance")
    
return(df2)
}

##########################################
### FANGORN_GTDB_RRN - minimap files ###
##########################################

# importing all paf files that were taxonomically classified using FANGORN GTDB_nrRep - minimap2

# Set your working directory to where your files are
setwd("./minimap2/tax_assign/GTDB_Fangorn/")

# going per mock community to reduce computational load #
# Get list of .paf file names per sample #

######### MCAP st moc com #######

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
write.csv(GTDB_minimap, "./R_files/output_files/GTDB_minimap_RRN_ONT.csv", row.names = FALSE)

###################################################################################################################################################

############################
##### 2.FANGORN RefSeq #####

RefSeq_tax <- read.csv("./minimap_database_tsv_files/taxRep_RefSeq_nr.csv", header = F) #path to taxonomy files
colnames(RefSeq_tax) <- c("op", "tax")

#####

# making and defining a function to process each data frame
process_data <- function(df) {
  # Set column names ## NOTE: name op as match instead for mirror ##
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
  
# getting all the unique hits out for a query sequence
  df_uni <- df %>% filter(MapQ > 0)
  df_uni_ids <- df_uni$Query
  
  # filtering out all the query sequneces with unique hits from df that alos have non unique hits for them
  df_uni_mapG0 <- df %>% filter(!(Query %in% df_uni_ids & MapQ == 0))
  # in the non unique hits first filter by AS the ones with the highest alignment score is selected
  # first getting out all the query that have MapQ = 0
  df_map0 <- df_uni_mapG0 %>% filter(MapQ == 0)
  
  ## now looking at the the queries with multiple hits so MapQ = 0 with the AS score are different i.e. there is one max score comapred to the others,     #selecting based on that max 
  df_map0_dAS <- df_map0 %>% 
  group_by(Query) %>%
  slice_max(AS, n = 1) %>%
  ungroup() # This will leave the rows that have identical AS for a query
  # get only the the one hist when mapQ = 0 and AS is different for one query
  df_map0_dAS2 <- df_map0_dAS %>%
  group_by(Query, AS) %>%
  filter(!(n() > 1)) %>%
  ungroup()
  # rbind the unique hit where mapq > 0 (uni) and for the hits where mapQ = 0 but AS was diff with one higher than the other (dAS2)
  df_uni_dAS2 <- rbind(df_uni, df_map0_dAS2)
  # cleaning up the tax column 
  df_uni_dAS2 <- df_uni_dAS2 %>%
    # Split the Taxon column into separate taxonomic levels
    separate(Tax, into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep="\\|") %>%
    # Remove the prefix from each taxonomic level
    mutate_at(vars(Kingdom:Species), ~substr(., 4, nchar(.))) 
    # removing the _letters in species names
    df_uni_dAS2$Species <- gsub("_[A-Z]", "", df_uni_dAS2$Species) 
    # adding an lca column 
    df_uni_dAS2$lca <- df_uni_dAS2$Species 
  
  ## pull out the reads that have equal AS scores, i.e. multiple rows for the same AS
  df_map0_sAS <- df_map0_dAS %>%
  group_by(Query, AS) %>%
  filter(n() > 1) %>%
  ungroup()
  # cleaning up the tax column 
  df_map0_sAS <- df_map0_sAS %>%
    # Split the Taxon column into separate taxonomic levels
    separate(Tax, into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep="\\|") %>%
    # Remove the prefix from each taxonomic level
    mutate_at(vars(Kingdom:Species), ~substr(., 4, nchar(.))) 
  # if the species is the same for the mutiple hits keep only the first entry
    df_map0_sAS_sp <- df_map0_sAS %>%
    group_by(Query, Species) %>%
    slice_head(n = 1) %>%
    ungroup()  # so queries with (i) mapQ = 0 and same AS that all have hits to the same species only the first hit is taken, and (ii) all hits for mapQ = 0 and same AS that have different hits at species level as also kept
    ## identifying (i) and (ii) from df_min0_sAS_sp and then adding lca 
    # (i) pulling out the hits that have mapQ = 0 and same AS that all have hits to the same species only the first hit has been taken
    df_map0_sAS_1sp <- df_map0_sAS_sp %>%
    group_by(Query) %>%
    filter(!(n() > 1)) %>%
    ungroup()
    # adding an lca column to it 
    df_map0_sAS_1sp$lca <- df_map0_sAS_1sp$Species
    # (ii) pulling out hits for mapQ = 0 and same AS that have different hits at species level as also kept
    df_map0_sAS_msp <- df_map0_sAS_sp %>%
    group_by(Query) %>%
    filter(n() > 1) %>%
    ungroup()
    # removing the _letters in species names
    df_map0_sAS_msp$Species <- gsub("_[A-Z]", "", df_map0_sAS_msp$Species) 
    # now lets add the lowest common acestor for them 
    df_map0_sAS_msp <- df_map0_sAS_msp %>%
    group_by(Query) %>%
    mutate(
    lca = case_when(
      n_distinct(Genus) == 1 ~ Genus,
      n_distinct(Family) == 1 ~ Family,
      n_distinct(Order) == 1 ~ Order,
      n_distinct(Class) == 1 ~ Class,
      n_distinct(Phylum) == 1 ~ Phylum,
      n_distinct(Kingdom) == 1 ~ Kingdom,
      TRUE ~ NA_character_  # If no common taxonomic level is found
    )
  ) %>%
  ungroup() %>%
  # keep only one hit now that lca has been added for muyltiple hits for one query
  group_by(Query) %>%
  slice_head(n = 1) %>%
  ungroup()
    # rbinding the map0_sAS dups nad no dups 
    df_m0_sAS <- rbind(df_map0_sAS_1sp, df_map0_sAS_msp)
    
    ## rbinding it all now 
    df2 <- rbind(df_uni_dAS2, df_m0_sAS)
    
# getting counts for all the alignments at species level #
    df2 <- df2 %>% group_by(lca) %>%
          summarise(Counts = n())
  
  # Reshape dataframe to wide format
  df2 <- df2 %>%
    pivot_wider(names_from = lca, values_from = Counts, values_fill = 0)
  
  ## changing counts to relative abundance ##
  df2 <- df2/rowSums(df2)*100
  rowSums(df2) ## to check if each sample adds up to a 100 ##
  
  # Convert tf from wide back to long format
  df2 <- df2 %>%
    pivot_longer(cols = everything(), names_to = "lca_sp", values_to = "Rel_abundance")
    
return(df2)
}
##########################################
### FANGORN_RefSeq_RRN - minimap files ###
##########################################

# importing all paf files that were taxonomically classified using FANGORN GTDB_nrRep - minimap2
# Set your working directory to where your files are
setwd("./minimap2/tax_assign/RefSeq_Fangorn/")

# going per mock community to reduce computational load #
# Get list of .paf file names per sample #

######### MCAP 24 st moc com #######

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
write.csv(RefSeq_minimap, "./R_files/output_files/RefSeq_minimap_RRN_ONT.csv", row.names = FALSE)  

################################################################################################################################################################

##########################
##### 3.rrn_DBv2 #########

rrnDB_tax <- read.csv("/data/Food/analysis/R1150_biotransformation/minimap_database_tsv_files/rrn_DBv2_tax.csv", header = F) # path to the taxonomy file
rrnDB_tax <- rrnDB_tax[,-(2:3)]
colnames(rrnDB_tax) <- c("op", "Class", "Order", "Family", "Genus", "Species")
rrnDB_tax <- rrnDB_tax %>% select("op", "Species")

#####

# making and defining a function to process each data frame
process_data <- function(df) {
  # Set column names ## NOTE: name op as match instead for mirror ##
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
  
  # getting all the unique hits out for a query sequence
  df_uni <- df %>% filter(MapQ > 0)
  df_uni_ids <- df_uni$Query
  
  # filtering out all the query sequneces with unique hits from df that alos have non unique hits for them
  df_uni_mapG0 <- df %>% filter(!(Query %in% df_uni_ids & MapQ == 0))
  # in the non unique hits first filter by AS the ones with the highest alignment score is selected
  # first getting out all the query that have MapQ = 0
  df_map0 <- df_uni_mapG0 %>% filter(MapQ == 0)
  
  ## now looking at the the queries with multiple hits so MapQ = 0 with the AS score are different i.e. there is one max score comapred to the others,     #selecting based on that max 
  df_map0_dAS <- df_map0 %>% 
  group_by(Query) %>%
  slice_max(AS, n = 1) %>%
  ungroup() # This will leave the rows that have identical AS for a query
  # get only the the one hist when mapQ = 0 and AS is different for one query
  df_map0_dAS2 <- df_map0_dAS %>%
  group_by(Query, AS) %>%
  filter(!(n() > 1)) %>%
  ungroup()
  # rbind the unique hit where mapq > 0 (uni) and for the hits where mapQ = 0 but AS was diff with one higher than the other (dAS2)
  df_uni_dAS2 <- rbind(df_uni, df_map0_dAS2)
  # adding an lca column 
  df_uni_dAS2$lca <- df_uni_dAS2$Tax 
  
  ## pull out the reads that have equal AS scores, i.e. multiple rows for the same AS
  df_map0_sAS <- df_map0_dAS %>%
  group_by(Query, AS) %>%
  filter(n() > 1) %>%
  ungroup()
  # if the species is the same for the mutiple hits keep only the first entry
    df_map0_sAS_sp <- df_map0_sAS %>%
    group_by(Query, Tax) %>%
    slice_head(n = 1) %>%
    ungroup()  # so queries with (i) mapQ = 0 and same AS that all have hits to the same species only the first hit is taken, and (ii) all hits for mapQ = 0 and same AS that have different hits at species level as also kept
    ## identifying (i) and (ii) from df_min0_sAS_sp and then adding lca 
    # (i) pulling out the hits that have mapQ = 0 and same AS that all have hits to the same species only the first hit has been taken
    df_map0_sAS_1sp <- df_map0_sAS_sp %>%
    group_by(Query) %>%
    filter(!(n() > 1)) %>%
    ungroup()
    # formating tax a bit
    df_map0_sAS_1sp <- df_map0_sAS_1sp %>%
    separate(Tax, into = c("Genus", "Species"), sep = "_", remove = F) %>%  # Separate the 'Tax' column into 'Genus' and 'Species'
    mutate(Tax = gsub("_", " ", Tax)) %>%  # Remove underscores from the 'Species' column if there are any left
    select(-Species)
    # adding an lca column to it 
    df_map0_sAS_1sp$lca <- df_map0_sAS_1sp$Tax
    # (ii) pulling out hits for mapQ = 0 and same AS that have different hits at species level as also kept
    df_map0_sAS_msp <- df_map0_sAS_sp %>%
    group_by(Query) %>%
    filter(n() > 1) %>%
    ungroup() 
    # now lets add the lowest common acestor for them 
    # first separting tax column to genus and species 
    df_map0_sAS_msp <- df_map0_sAS_msp %>%
    separate(Tax, into = c("Genus", "Species"), sep = "_", remove = F) %>%  # Separate the 'Tax' column into 'Genus' and 'Species'
    mutate(Tax = gsub("_", " ", Tax)) %>%  # Remove underscores from the 'Species' column if there are any left
    select(-Species)
    # low adding at LCA
    df_map0_sAS_msp <- df_map0_sAS_msp %>%
    group_by(Query) %>%
    mutate(
    lca = case_when(
      n_distinct(Genus) == 1 ~ Genus,
      TRUE ~ NA_character_  # If no common taxonomic level is found
    )
  ) %>%
  ungroup() %>%
  # keep only one hit now that lca has been added for muyltiple hits for one query
  group_by(Query) %>%
  slice_head(n = 1) %>%
  ungroup()
    # rbinding the map0_sAS dups and no dups 
    df_m0_sAS <- rbind(df_map0_sAS_1sp, df_map0_sAS_msp)
  
  ## rbinding it all now 
  # before that
  df_uni_dAS2 <- df_uni_dAS2 %>%
    separate(Tax, into = c("Genus", "Species"), sep = "_", remove = F) %>%  # Separate the 'Tax' column into 'Genus' and 'Species'
    mutate(Tax = gsub("_", " ", Tax)) %>%  # Remove underscores from the 'Species' column if there are any left
    mutate(lca = gsub("_", " ", lca)) %>%
    select(-Species)
    # no binding it all
    df2 <- rbind(df_uni_dAS2, df_m0_sAS)
  
  #### for rrn_DBv2 and mirror #####
  # Clean and organise the table
  df2 <- df2 %>%
  # getting counts for all the alignments at species level #
  group_by(lca) %>%
    summarise(Counts = n())
  
  # Reshape dataframe to wide format
  df2 <- df2 %>%
    pivot_wider(names_from = lca, values_from = Counts, values_fill = 0)
  
  ## changing counts to relative abundance ##
  df2 <- df2/rowSums(df2)*100
  rowSums(df2) ## to check if each sample adds up to a 100 ##
  
  # Convert tf from wide back to long format
  df2 <- df2 %>%
    pivot_longer(cols = everything(), names_to = "lca_sp", values_to = "Rel_abundance")
  
  return(df2)
}

##########################################
### rrn_DBv2_RRN - minimap files ###
##########################################

# importing all paf files that were taxonomically classified using FANGORN GTDB_nrRep - minimap2
# Set your working directory to where your files are
setwd("./minimap2/tax_assign/rrn_DBv2/")

# going per mock community to reduce computational load #
# Get list of .paf file names per sample #

#########  MCAP moc com #######

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
write.csv(rrnDB_minimap, "./R_files/output_files/rrnDB_minimap_RRN_ONT.csv", row.names = FALSE)

#######################################################################################################################################################

##########################
##### 4.MIrROR ###########

mirror_tax <- read.csv("/data/Food/analysis/R1150_biotransformation/minimap_database_tsv_files/MIrROR_DB_r01.csv", header = T)
mirror_tax <- mirror_tax %>% select(X.Accession, gtdb_taxonomy_s)
colnames(mirror_tax) <- c("op", "Species")

#####

# making and defining a function to process each data frame
process_data <- function(df) {
  # Set column names ## NOTE: name op as match instead for mirror ##
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
  
  # getting all the unique hits out for a query sequence
  df_uni <- df %>% filter(MapQ > 0)
  df_uni_ids <- df_uni$Query
  
  # filtering out all the query sequneces with unique hits from df that alos have non unique hits for them
  df_uni_mapG0 <- df %>% filter(!(Query %in% df_uni_ids & MapQ == 0))
  # in the non unique hits first filter by AS the ones with the highest alignment score is selected
  # first getting out all the query that have MapQ = 0
  df_map0 <- df_uni_mapG0 %>% filter(MapQ == 0)
  
  ## now looking at the the queries with multiple hits so MapQ = 0 with the AS score are different i.e. there is one max score comapred to the others,     #selecting based on that max 
  df_map0_dAS <- df_map0 %>% 
  group_by(Query) %>%
  slice_max(AS, n = 1) %>%
  ungroup() # This will leave the rows that have identical AS for a query
  # get only the the one hist when mapQ = 0 and AS is different for one query
  df_map0_dAS2 <- df_map0_dAS %>%
  group_by(Query, AS) %>%
  filter(!(n() > 1)) %>%
  ungroup()
  # rbind the unique hit where mapq > 0 (uni) and for the hits where mapQ = 0 but AS was diff with one higher than the other (dAS2)
  df_uni_dAS2 <- rbind(df_uni, df_map0_dAS2)
  # adding an lca column 
  df_uni_dAS2$lca <- df_uni_dAS2$Tax 
  
  ## pull out the reads that have equal AS scores, i.e. multiple rows for the same AS
  df_map0_sAS <- df_map0_dAS %>%
  group_by(Query, AS) %>%
  filter(n() > 1) %>%
  ungroup()
  # if the species is the same for the mutiple hits keep only the first entry
    df_map0_sAS_sp <- df_map0_sAS %>%
    group_by(Query, Tax) %>%
    slice_head(n = 1) %>%
    ungroup()  # so queries with (i) mapQ = 0 and same AS that all have hits to the same species only the first hit is taken, and (ii) all hits for mapQ = 0 and same AS that have different hits at species level as also kept
    ## identifying (i) and (ii) from df_min0_sAS_sp and then adding lca 
    # (i) pulling out the hits that have mapQ = 0 and same AS that all have hits to the same species only the first hit has been taken
    df_map0_sAS_1sp <- df_map0_sAS_sp %>%
    group_by(Query) %>%
    filter(!(n() > 1)) %>%
    ungroup()
    # formating tax a bit
    df_map0_sAS_1sp <- df_map0_sAS_1sp %>%
    separate(Tax, into = c("Genus", "Species"), sep = "_", remove = F) %>%  # Separate the 'Tax' column into 'Genus' and 'Species'
    mutate(Tax = gsub("_", " ", Tax)) %>%  # Remove underscores from the 'Species' column if there are any left
    select(-Species)
    # adding an lca column to it 
    df_map0_sAS_1sp$lca <- df_map0_sAS_1sp$Tax
    # (ii) pulling out hits for mapQ = 0 and same AS that have different hits at species level as also kept
    df_map0_sAS_msp <- df_map0_sAS_sp %>%
    group_by(Query) %>%
    filter(n() > 1) %>%
    ungroup() 
    # now lets add the lowest common acestor for them 
    # first separting tax column to genus and species 
    df_map0_sAS_msp <- df_map0_sAS_msp %>%
    separate(Tax, into = c("Genus", "Species"), sep = "_", remove = F) %>%  # Separate the 'Tax' column into 'Genus' and 'Species'
    mutate(Tax = gsub("_", " ", Tax)) %>%  # Remove underscores from the 'Species' column if there are any left
    select(-Species)
    # low adding at LCA
    df_map0_sAS_msp <- df_map0_sAS_msp %>%
    group_by(Query) %>%
    mutate(
    lca = case_when(
      n_distinct(Genus) == 1 ~ Genus,
      TRUE ~ NA_character_  # If no common taxonomic level is found
    )
  ) %>%
  ungroup() %>%
  # keep only one hit now that lca has been added for muyltiple hits for one query
  group_by(Query) %>%
  slice_head(n = 1) %>%
  ungroup()
    # rbinding the map0_sAS dups and no dups 
    df_m0_sAS <- rbind(df_map0_sAS_1sp, df_map0_sAS_msp)
  
  ## rbinding it all now 
  # before that
  df_uni_dAS2 <- df_uni_dAS2 %>%
    separate(Tax, into = c("Genus", "Species"), sep = "_", remove = F) %>%  # Separate the 'Tax' column into 'Genus' and 'Species'
    mutate(Tax = gsub("_", " ", Tax)) %>%  # Remove underscores from the 'Species' column if there are any left
    mutate(lca = gsub("_", " ", lca)) %>%
    select(-Species)
    # no binding it all
    df2 <- rbind(df_uni_dAS2, df_m0_sAS)
    
  #### for rrn_DBv2 and mirror #####
  # Clean and organise the table
  df2 <- df2 %>%
  # getting counts for all the alignments at species level #
  group_by(lca) %>%
    summarise(Counts = n())
  
  # Reshape dataframe to wide format
  df2 <- df2 %>%
    pivot_wider(names_from = lca, values_from = Counts, values_fill = 0)
  
  ## changing counts to relative abundance ##
  df2 <- df2/rowSums(df2)*100
  rowSums(df2) ## to check if each sample adds up to a 100 ##
  
  # Convert tf from wide back to long format
  df2 <- df2 %>%
    pivot_longer(cols = everything(), names_to = "lca_sp", values_to = "Rel_abundance")
  
  return(df2)
}

##########################################
### MIrROR_RRN - minimap files ###
##########################################

# importing all paf files that were taxonomically classified using FANGORN GTDB_nrRep - minimap2
# Set your working directory to where your files are
setwd("./minimap2/tax_assign/mirror/")

# going per mock community to reduce computational load #
# Get list of .paf file names per sample #

######### MCAP moc com #######
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
write.csv(mirror_minimap, "./RRN_mc_ONT/R_files/output_files/mirror_minimap_RRN_ONT.csv", row.names = FALSE)
