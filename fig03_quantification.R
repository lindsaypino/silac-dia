
# Setup -------------------------------------------------------------------

options(stringsAsFactors = FALSE)
library(tidyverse)
library(rio)
library(ggthemes)
library(DBI)
library(RSQLite)
library(patchwork)


# Formatting DDA data -----------------------------------------------------

process_dda_pd <- function(data_dda, input_files){
  
  # Input file setup
  # The input file contains the raw file names to match to the file_id
  input_files <- input_files %>% 
    rename_all(tolower) %>% 
    rename_all(~str_replace_all(., "\\W", "_")) %>% # replaces all non-word characters to an underscore
    rename(file_id = "study_file_id") %>% 
    mutate(file_id = tolower(file_id)) %>% 
    select(file_id, file_name)
  
  # Removing that one row which has a blank space
  input_files[input_files == ""] <- NA
  input_files <- input_files %>% 
    filter(!is.na(file_id))
  
  # Parsing the raw file names
  input_files$file_name <- unlist(lapply(input_files$file_name, function(x){
    unlist(strsplit(x, split = "\\", fixed = TRUE))[6]
  }))
  
  input_files$file_name <- unlist(lapply(input_files$file_name, function(x){
    unlist(strsplit(x, split = ".", fixed = TRUE))[1]
  }))
  
  # Formatting the raw file name
  # This is to match the file names in the DIA dataset
  input_files$file_name <- gsub("DDA_", "", input_files$file_name)
  
  # formatting names and selecting columns for use
  data_dda <- data_dda %>% 
    rename_all(tolower) %>% 
    rename_all(~str_replace_all(., "\\W", "_")) %>% # replaces all non-word characters to an underscore
    rename_all(~str_replace_all(., "_{2,}", "_")) %>% # replaces multiple underscores to single underscores
    select(-contains("found_in_"), 
           -contains("ratio"),
           -contains("abundances_"))
  
  # Cleaning up the modifications 
  # Using regex to clean SILAC label
  data_dda$modifications <- gsub("; ", ";", data_dda$modifications)
  data_dda$modifications <- gsub("K\\d{,2};", "", data_dda$modifications)
  data_dda$modifications <- gsub("R\\d{,2};", "", data_dda$modifications)
  data_dda$modifications <- gsub("K\\d{,2}", "", data_dda$modifications)
  data_dda$modifications <- gsub("R\\d{,2}", "", data_dda$modifications)
  data_dda$modifications <- gsub("\\[]", "", data_dda$modifications)
  data_dda$modifications <- gsub("[1-9]xLabel:13C\\(6)15N\\(4)", "", data_dda$modifications)
  data_dda$modifications <- gsub("[1-9]xLabel:13C\\(6)15N\\(2)", "", data_dda$modifications)
  data_dda$modifications <- gsub("; $", "", data_dda$modifications)
  data_dda$modifications <- gsub(" $", "", data_dda$modifications)
  
  # Cleaning up Carbamidomethylation
  data_dda$modifications <- gsub("[1-9]xCarbamidomethyl ", "", data_dda$modifications)
  
  # Cleaning up Oxidation
  data_dda$modifications <- gsub("[1-9]xOxidation ", "", data_dda$modifications)
  
  # Removing miscellaneous characters
  data_dda$modifications <- gsub("];\\[", ";", data_dda$modifications)
  data_dda$modifications <- gsub("^;", "", data_dda$modifications)
  data_dda$modifications <- gsub(";$", "", data_dda$modifications)
  data_dda$modifications <- gsub("^\\[", "", data_dda$modifications)
  data_dda$modifications <- gsub("]$", "", data_dda$modifications)
  
  # Counting the number of cysteines and carbamidomethylation
  # Identifying N-term acetyl mods -- used to remove
  data_dda <- data_dda %>% 
    mutate(c_pep_count = str_count(sequence, "C"),
           c_mod_count = str_count(modifications, "C"),
           m_ox_count = str_detect(modifications, "M"),
           nterm_ac = str_detect(modifications, "1xAcetyl"))
  
  # Data wrangling to create long form of data
  tidy_dda <- data_dda %>% 
    filter(contaminant == FALSE,
           nterm_ac == FALSE, # Removing N-term acetylation
           m_ox_count == FALSE, # Removing methionine oxidation
           c_pep_count == c_mod_count) %>% 
    select(sequence, modifications, master_protein_accessions,  
           contains("abundance")) %>% 
    gather(temp, abundance, contains("abundance")) %>% 
    filter(abundance > 0) %>% 
    separate(temp, c("temp2", "file_id", "isotope", "temp3", "temp4", "temp5")) %>% 
    right_join(input_files, .) %>% 
    mutate(dataset = "DDA") %>% 
    select(-contains("temp"),
           -file_id) %>% 
    rename(protein = "master_protein_accessions",
           peptide = "sequence",
           peptide_modified_sequence = "modifications") %>% 
    mutate(peptide_modified_sequence = ifelse(peptide_modified_sequence == "", 
                                              peptide, peptide_modified_sequence))

  return(tidy_dda)
}

# Peptide group data 
# Output from Proteome Discoverer
#data_dda <- rio::import("data/curve/data_DDA_PD23/20191104_QEHFX_lkp_pSILAC-DDA_PeptideGroups.txt", setclass = "tibble")

# Input file table
# Output table from Proteome Discoverer
# This table contains the raw file names to match to the file id
#input_files <- rio::import("data/curve/data_DDA_PD23/20191104_QEHFX_lkp_pSILAC-DDA_InputFiles.txt", setclass = "tibble")

process_dda_mm <- function(quant_file){
  
  # PSM data 
  # Output from Metamorpheus (.psmtsv file)
  data_dda <- read.table(file=quant_file, sep="\t", header=TRUE)
  
  # formatting names 
  data_dda <- data_dda %>% 
    rename_all(tolower) %>% 
    rename_all(~str_replace_all(., "\\W", "_")) %>% # replaces all non-word characters to an underscore
    rename_all(~str_replace_all(., "_{2,}", "_")) # replaces multiple underscores to single underscores
  
  # filter for FDR? is this how you do it?
  #data_dda <- data_dda[data_dda$qvalue < 0.01,]
  
  # Data wrangling to create long form of data
  tidy_dda <- data_dda %>% 
    select('peak_intensity', 'protein_group', 'file_name', 'base_sequence', 'full_sequence') %>% 
    rename(peptideseq = "base_sequence")
  
  # Clean up and fix raw file names
  tidy_dda$file_name <- gsub(".raw", "", tidy_dda$file_name)
  tidy_dda$file_name <- gsub("_20200316162818", "", tidy_dda$file_name)  # fix the file that was rerun
  tidy_dda$file_name <- gsub("dia_007", "dda_007", tidy_dda$file_name)  # fix the file that was mislabeled
  
  # label peptides as light or heavy
  tidy_dda$isotope <- gsub(".*(\\+8.014).*|.*(\\+10.008).*", "heavy", tidy_dda$peptideseq)
  
  # remove the modification from the peptideseq variable
  tidy_dda$peptideseq <- gsub("\\(\\+8.014\\)|\\(\\+10.008\\)", "", tidy_dda$peptideseq)
  
  # label peptides as light or heavy
  tidy_dda$isotope[which(tidy_dda$isotope != "heavy")] = "light"
  
  # Adding peptidemodseq column to match DIA convention
  tidy_dda$peptidemodseq <- paste(tidy_dda$peptideseq, tidy_dda$isotope, sep = "_")
  
  # Parsing the raw file name
  # Creating new variables to indicate sample conditions
  tidy_dda <- tidy_dda %>%
    separate(file_name, c("temp1", "temp2", "temp3", "temp4", "ratio_id", "heavy_ug", "dataset", "temp5"),
             sep = "_", remove = FALSE, convert = TRUE) %>%
    select(-contains("temp")) %>%
    mutate(heavy_ug = str_replace_all(heavy_ug, "ug", ""),
           heavy_ug = as.numeric(str_replace_all(heavy_ug, "-", ".")))
  
  tidy_dda <- tidy_dda %>%
    rename(abundance = "peak_intensity", peptide = "peptideseq", 
           peptide_modified_sequence = "full_sequence") %>%
    select(file_name, peptide, peptide_modified_sequence, 
           isotope, abundance, heavy_ug, ratio_id) %>%
    distinct(file_name, peptide, peptide_modified_sequence, 
             isotope, abundance, .keep_all = TRUE)
  
  tidy_dda$dataset <- "DDA"
  tidy_dda$abundance <- tidy_dda$abundance + 1
  
  return(tidy_dda)
  
}


# Rolfs: When you search as "turnover", the program automatically splits the intensity of 
# peptides with missed cleavages that have both a heavy and a light label. If you want to see 
# if any such peptides were identified/quantified, you need to look int the 
# "AllQuantifiedPeaks.psmtsv" output file instead of the "AllQuantifiedPeptides.psmtsv" file.
dda_df <- process_dda_mm("data/curve/2020-10-23-00-11-30/Task1-turnover curve/AllQuantifiedPeaks.tsv")


# Formatting DIA data -----------------------------------------------------

process_dia <- function(data_dia){
  
  # formatting names 
  data_dia <- data_dia %>% 
    rename_all(tolower) %>% 
    rename_all(~str_replace_all(., "\\W", "_"))
  
  # Counting the number of cysteines and carbamidomethylation
  data_dia <- data_dia %>% 
    mutate(c_pep_count = str_count(peptide_modified_sequence, "C"),
           c_mod_count = str_count(peptide_modified_sequence, "57"))
  
  data_dia$peptide <- gsub("[+57]", "", data_dia$peptide_modified_sequence)
  
  # Data wrangling to create a long form of the dia data
  tidy_dia <- data_dia %>% 
    filter(c_pep_count == c_mod_count,
           total_area_fragment > 0) %>% 
    select(replicate, peptide, peptide_modified_sequence, 
           isotope_label_type, total_area_fragment) %>% 
    distinct(replicate, peptide, peptide_modified_sequence, 
             isotope_label_type, total_area_fragment) %>% 
    rename(isotope = "isotope_label_type",
           abundance = "total_area_fragment",
           file_name = "replicate")
  
  tidy_dia$dataset <- "MS2"
  
  tidy_dia <- tidy_dia %>%
    separate(file_name, c("temp1", "temp2", "temp3", "temp4", "ratio_id", "heavy_ug", "temp5", "temp6"),
             sep = "_", remove = FALSE, convert = TRUE) %>%
    select(-contains("temp")) %>%
    mutate(heavy_ug = str_replace_all(heavy_ug, "ug", ""),
           heavy_ug = as.numeric(str_replace_all(heavy_ug, "-", ".")))
  
  tidy_dia$dataset <- "DIA"
  
  # Data wrangling again for the MS1 data
  tidy_ms1 <- data_dia %>% 
    filter(c_pep_count == c_mod_count,
           total_area_fragment > 0) %>% 
    select(replicate, peptide, peptide_modified_sequence, 
           isotope_label_type, area, isotope_dist_rank) %>% 
    filter(isotope_dist_rank == 1) %>% 
    distinct(replicate, peptide, peptide_modified_sequence, 
             isotope_label_type, area) %>% 
    rename(isotope = "isotope_label_type",
           abundance = "area",
           file_name = "replicate")
  
  tidy_ms1 <- tidy_ms1 %>%
    separate(file_name, c("temp1", "temp2", "temp3", "temp4", "ratio_id", "heavy_ug", "temp5", "temp6"),
             sep = "_", remove = FALSE, convert = TRUE) %>%
    select(-contains("temp")) %>%
    mutate(heavy_ug = str_replace_all(heavy_ug, "ug", ""),
           heavy_ug = as.numeric(str_replace_all(heavy_ug, "-", ".")))
  
  tidy_ms1$dataset <- "MS1"
  
  # 
  data <- rbind(tidy_dia, tidy_ms1)
  
  # There are duplicate entries in the dataset. This is likely due to protein grouping
  # Removing the duplicate entry
  data <- data %>% 
    distinct(file_name, peptide, peptide_modified_sequence, 
             isotope, abundance, .keep_all = TRUE)
  
  return(data)
}

# DIA Peptide quant data
# Output from Skyline post Encyclopedia
dia_df <- process_dia(import("data/curve/QUANTIFICATION - DIA MS1 and MS2.csv"))


# Combining datasets ------------------------------------------------------


# Merging DIA and DDA datasets
data <- bind_rows(dia_df, dda_df)

# removing datasets
rm(dia_df, dda_df);gc()
  
##
## Q: WHY dda_df HAS SO MANY NA QUANTS? LEADS TO NO RATIO...
## A: MANY ZERO QUANTS WHICH GET LOG TRANSFORMED TO INF
##

# Calculating fractional abundance 
# Heavy / (Heavy + Light)
data_fraction <- data %>%
  group_by(file_name, ratio_id, heavy_ug, peptide, isotope, dataset) %>% 
  summarise(abundance = mean(abundance)) %>% 
  spread(isotope, abundance) %>%
  mutate(hl_ratio_log10 = log10(heavy / light),
         fraction = heavy / (heavy + light)) %>% 
  filter(!is.na(fraction))


# High Conf peptides ------------------------------------------------------

# Selecting for peptides with high heavy label incorporation
# The HeLa cells were grown in culture for 6 cell passages. 
# 6 passages is the "recommended" number of passages for full isotopic labeling 
# The assumption is that the heavy labeled sample should be ~100% heavy labeled.
# Peptides not meeting this criteria are filtered out
# I used the median value of heavy incorporation for filtering criteria
high_conf_pep <- data_fraction %>% 
  filter(heavy_ug == 1) %>% 
  select(-heavy, -light, -hl_ratio_log10) %>% 
  group_by(peptide, dataset) %>% 
  summarize(fraction_mean = mean(fraction, na.rm = TRUE),
            n = n()) %>% 
  ungroup()

temp <- high_conf_pep %>% 
  filter(n == 3, 
         dataset == "DIA",
         fraction_mean > .90) %>%   # is this necessary?
  select(peptide)


high_conf_pep <- high_conf_pep %>% 
  # filter(n == 3) %>%
  select(-n) %>% 
  spread(dataset, fraction_mean) %>% 
  filter(!is.na(DIA) & !is.na(DDA))

# filtered data
# found in three technical replicates
# heavy incorporation in 100% heavy sample > 95.7% (median value)
data_high_conf <- inner_join(temp, data_fraction)

rm(temp);gc()


# Wide data ---------------------------------------------------------------


data_wide <- data_fraction %>% 
  filter(heavy_ug > 0) %>% 
  group_by(peptide, heavy_ug, dataset) %>% 
  summarize(fraction_mean = mean(fraction, na.rm = TRUE),
            fraction_n = n()) %>% 
  ungroup() %>% 
  spread(dataset, fraction_mean)


# TIC-based normalization factor ------------------------------------------

dia_elib <- "data/curve/20200828_QEHFX_lkp_pSILAC-DIA_curve_2020-10-13_09-33-04/20200828_QEHFX_lkp_pSILAC-DIA_curve_QUANT.elib"
con <- dbConnect(drv = RSQLite::SQLite(), 
                 dbname = dia_elib)
elib_df <- dbGetQuery(conn = con, statement = paste("SELECT * FROM metadata"))
dbDisconnect(con)

elib_df <- elib_df[grep("TIC_", elib_df$Key), ]
heavy_tic <- mean(elib_df[grep("A_1-0ug", elib_df$Key), ]$Value)
light_tic <- mean(elib_df[grep("N_0-0ug", elib_df$Key), ]$Value)
adj_ticratio <- heavy_tic/light_tic

rm(elib_df)

  
# MA plot -----------------------------------------------------------------

# selecting the heavy log2 abundance in the 100% heavy sample of the DDA and DIA
# Calculating the mean heavy abundance per peptide
# Calculating the sd of heavy abundance per peptide
# Calculating the CV of heavy abundance per peptide
heavy_log2 <- data_fraction %>% 
  filter(heavy_ug == 1,
         !is.na(heavy)) %>% 
  rename(heavy_100 = "heavy") %>% 
  select(peptide, dataset, heavy_100) %>% 
  group_by(peptide, dataset) %>% 
  summarize(heavy_mean = mean(heavy_100, na.rm = TRUE),
            heavy_sd = sd(heavy_100, na.rm = TRUE),
            heavy_n = n()) %>% 
  ungroup() %>% 
  mutate(heavy_cv = heavy_sd / heavy_mean * 100)

# Joining the mean heavy abundance with the complete dataset
data_ma <- data_fraction %>% 
  filter(!is.na(hl_ratio_log10)) %>% 
  inner_join(heavy_log2, .)

# Subsetting data for the DDA plots
ma_dda <- data_ma %>% 
  filter(dataset == "DDA", 
         heavy_ug == c(0.7, 0.5, 0.3 , 0.100, 0.010, 0.001)) #, 
         #heavy_cv < 20,
         #heavy_mean > 0.90,
         #heavy_n == 3) 

# Subsetting data for the DIA plots
ma_dia <- data_ma %>% 
  filter(dataset == "DIA", 
         heavy_ug == c(0.7, 0.5, 0.3 , 0.100, 0.010, 0.001)) #, 
         #heavy_cv < 20,
         #heavy_mean > 0.90,
         #heavy_n == 3) 


# Plots ---------------------------------------------------------------

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")

# DDA MA plot
dda1 <- ggplot(ma_dda) +
  geom_point(aes(x = log2(heavy_mean), y = hl_ratio_log10, color = factor(heavy_ug)), 
             size = 2, alpha = 0.1) +
  geom_smooth(aes(x = log2(heavy_mean), y = hl_ratio_log10, 
                  group = factor(heavy_ug)), se = FALSE, color = "darkgray", 
              method = "loess", linetype = "dashed") + 
  geom_hline(aes(yintercept = log10((70/30)*adj_ticratio), color = factor("0.7"))) +
  geom_hline(aes(yintercept = log10((50/50)*adj_ticratio), color = factor("0.5"))) +
  geom_hline(aes(yintercept = log10((30/70)*adj_ticratio), color = factor("0.3"))) +
  geom_hline(aes(yintercept = log10((10/90)*adj_ticratio), color = factor("0.1"))) +
  geom_hline(aes(yintercept = log10((1/99)*adj_ticratio), color = factor("0.01"))) +
  geom_hline(aes(yintercept = log10((.1/99.9)*adj_ticratio), color = factor("0.001"))) +
  geom_text(data=data.frame(x=35,y=log10((70/30)*adj_ticratio)), aes(x, y, color = factor("0.7")), 
            label="70%", vjust=-1) +
  geom_text(data=data.frame(x=35,y=log10((50/50)*adj_ticratio)), aes(x, y, color = factor("0.5")), 
            label="50%", vjust=-1) +
  geom_text(data=data.frame(x=35,y=log10((30/70)*adj_ticratio)), aes(x, y, color = factor("0.3")), 
            label="30%", vjust=-1) +
  geom_text(data=data.frame(x=35,y=log10((10/90)*adj_ticratio)), aes(x, y, color = factor("0.1")), 
            label="10%", vjust=-1) +
  geom_text(data=data.frame(x=35,y=log10((1/99)*adj_ticratio)), aes(x, y, color = factor("0.01")), 
            label="1%", vjust=-1) +
  geom_text(data=data.frame(x=35,y=log10((.1/99.9)*adj_ticratio)), aes(x, y, color = factor("0.001")), 
            label="0.1%", vjust=-1) +
  scale_y_continuous(limits = c(-5,1)) +
  scale_x_continuous(limits = c(10, 35)) +
  theme_classic(base_size = 12) +
  guides(color = FALSE) +
  labs(y = expression(Log[10]~(Heavy/Light)),
       x = expression(Log[2]~DDA~Heavy~abundance),
       title = "DDA") +
  scale_color_manual(values=cbPalette)

# DDA Boxplot
dda2 <- ggplot(ma_dda) +
  geom_boxplot(aes(x = factor(heavy_ug), y = hl_ratio_log10, color = factor(heavy_ug)), alpha = 0.1) +
  geom_hline(aes(yintercept = log10((70/30)*adj_ticratio), color = factor("0.7"))) +
  geom_hline(aes(yintercept = log10((50/50)*adj_ticratio), color = factor("0.5"))) +
  geom_hline(aes(yintercept = log10((30/70)*adj_ticratio), color = factor("0.3"))) +
  geom_hline(aes(yintercept = log10((10/90)*adj_ticratio), color = factor("0.1"))) +
  geom_hline(aes(yintercept = log10((1/99)*adj_ticratio), color = factor("0.01"))) +
  geom_hline(aes(yintercept = log10((.1/99.9)*adj_ticratio), color = factor("0.001"))) +
  scale_color_manual(values=cbPalette) +
  scale_y_continuous(limits = c(-5,1)) +
  theme_void() +
  guides(color = FALSE) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank())

# DIA MA plot
dia1 <- ggplot(ma_dia) +
  geom_point(aes(x = log2(heavy_mean), y = hl_ratio_log10, color = factor(heavy_ug)), 
             size = 2, alpha = 0.1) +
  geom_smooth(aes(x = log2(heavy_mean), y = hl_ratio_log10, 
                  group = factor(heavy_ug)), se = FALSE, color = "darkgray", 
              method = "loess", linetype = "dashed") + 
  geom_hline(aes(yintercept = log10((70/30)*adj_ticratio), color = factor("0.7"))) +
  geom_hline(aes(yintercept = log10((50/50)*adj_ticratio), color = factor("0.5"))) +
  geom_hline(aes(yintercept = log10((30/70)*adj_ticratio), color = factor("0.3"))) +
  geom_hline(aes(yintercept = log10((10/90)*adj_ticratio), color = factor("0.1"))) +
  geom_hline(aes(yintercept = log10((1/99)*adj_ticratio), color = factor("0.01"))) +
  geom_hline(aes(yintercept = log10((.1/99.9)*adj_ticratio), color = factor("0.001"))) +
  geom_text(data=data.frame(x=35,y=log10((70/30)*adj_ticratio)), aes(x, y, color = factor("0.7")), 
            label="70%", vjust=-1) +
  geom_text(data=data.frame(x=35,y=log10((50/50)*adj_ticratio)), aes(x, y, color = factor("0.5")), 
            label="50%", vjust=-1) +
  geom_text(data=data.frame(x=35,y=log10((30/70)*adj_ticratio)), aes(x, y, color = factor("0.3")), 
            label="30%", vjust=-1) +
  geom_text(data=data.frame(x=35,y=log10((10/90)*adj_ticratio)), aes(x, y, color = factor("0.1")), 
            label="10%", vjust=-1) +
  geom_text(data=data.frame(x=35,y=log10((1/99)*adj_ticratio)), aes(x, y, color = factor("0.01")), 
            label="1%", vjust=-1) +
  geom_text(data=data.frame(x=35,y=log10((.1/99.9)*adj_ticratio)), aes(x, y, color = factor("0.001")), 
            label="0.1%", vjust=-1) +
  scale_y_continuous(limits = c(-5,1)) +
  scale_x_continuous(limits = c(10, 35)) +
  theme_classic(base_size = 12) +
  labs(y = expression(Log[10]~(Heavy/Light)),
       x = expression(Log[2]~DIA~Heavy~abundance),
       title = "DIA") +
  guides(color = FALSE) +
  scale_color_manual(values=cbPalette)

# DIA boxplot
dia2 <- ggplot(ma_dia) +
  geom_boxplot(aes(x = factor(heavy_ug), y = hl_ratio_log10, color = factor(heavy_ug)), alpha = 0.1) +
  geom_hline(aes(yintercept = log10((70/30)*adj_ticratio), color = factor("0.7"))) +
  geom_hline(aes(yintercept = log10((50/50)*adj_ticratio), color = factor("0.5"))) +
  geom_hline(aes(yintercept = log10((30/70)*adj_ticratio), color = factor("0.3"))) +
  geom_hline(aes(yintercept = log10((10/90)*adj_ticratio), color = factor("0.1"))) +
  geom_hline(aes(yintercept = log10((1/99)*adj_ticratio), color = factor("0.01"))) +
  geom_hline(aes(yintercept = log10((.1/99.9)*adj_ticratio), color = factor("0.001"))) +
  scale_color_manual(values=cbPalette) +
  scale_y_continuous(limits = c(-5,1)) +
  theme_void() +
  guides(color = FALSE) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank())

# Plotting the figure
dda1 + dda2 + dia1 + dia2 +
  plot_layout(widths = c(1,0.5),
              heights = c(4,4)) +
  plot_annotation(tag_levels = "A")

ggsave(filename = "figures/fig03_lfqbench.svg", width = 8, height = 10)


## summary stats 
summary(subset(ma_dia, heavy_ug=="0.7"))

## calculate stats for n
length(unique(ma_dda[ma_dda$heavy_ug == '1',]$peptide))
length(unique(ma_dda[ma_dda$heavy_ug == '0.1',]$peptide))
length(unique(ma_dda[ma_dda$heavy_ug == '0.01',]$peptide))
length(unique(ma_dda[ma_dda$heavy_ug == '0.001',]$peptide))

length(unique(ma_dia[ma_dia$heavy_ug == '1',]$peptide))
length(unique(ma_dia[ma_dia$heavy_ug == '0.1',]$peptide))
length(unique(ma_dia[ma_dia$heavy_ug == '0.01',]$peptide))
length(unique(ma_dia[ma_dia$heavy_ug == '0.001',]$peptide))
