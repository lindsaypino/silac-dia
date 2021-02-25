
# Setup -------------------------------------------------------------------

options(stringsAsFactors = FALSE)
library(tidyverse)
library(rio)
library(ggthemes)
library(DBI)
library(RSQLite)
library(patchwork)
library(svglite)

# Formatting DDA data -----------------------------------------------------

process_dda <- function(psm_file){
  
  # PSM data 
  # Output from Proteome Discoverer
  data_dda <- rio::import(psm_file, setclass = "tibble")
  
  # formatting names 
  data_dda <- data_dda %>% 
    rename_all(tolower) %>% 
    rename_all(~str_replace_all(., "\\W", "_")) %>% # replaces all non-word characters to an underscore
    rename_all(~str_replace_all(., "_{2,}", "_")) # replaces multiple underscores to single underscores
  
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
  
  
  # Data wrangling to create long form of data
  tidy_dda <- data_dda %>% 
    select(annotated_sequence, modifications, spectrum_file) %>% 
    rename(file_name = "spectrum_file",
           peptideseq = "annotated_sequence") %>% 
    mutate(peptideseq = toupper(peptideseq)) ### modified AA denoted in lowercase affecting how duplicate peptides would be removed.
  
  # Clean up and fix raw file names
  tidy_dda$file_name <- gsub(".raw", "", tidy_dda$file_name)
  tidy_dda$file_name <- gsub("_20200316162818", "", tidy_dda$file_name)  # fix the file that was rerun
  tidy_dda$file_name <- gsub("dia_007", "dda_007", tidy_dda$file_name)  # fix the file that was mislabeled
  
  # label the detected peptides as light or heavy
  tidy_dda$isotope <- gsub(".*[(Label:13C(6)15N(2))].*|.*[(Label:13C(6)15N(4))].*", "heavy", tidy_dda$modifications)
  tidy_dda$isotope <- sub("^$", "light", tidy_dda$isotope)
  
  # Adding peptidemodseq column to match DIA convention
  tidy_dda$peptidemodseq <- paste(tidy_dda$peptideseq, tidy_dda$isotope, sep = "_")
  
  
  tidy_dda <- tidy_dda %>%
    select(file_name, peptideseq, isotope) %>%
    distinct(file_name, peptideseq, isotope, .keep_all = TRUE)

  return(tidy_dda)
}

process_dda_metamorpheus <- function(psm_file){
  
  # PSM data 
  # Output from Metamorpheus (.psmtsv file)
  data_dda <- read.table(file=psm_file, sep="\t", header=TRUE)
  
  # formatting names 
  data_dda <- data_dda %>% 
    rename_all(tolower) %>% 
    rename_all(~str_replace_all(., "\\W", "_")) %>% # replaces all non-word characters to an underscore
    rename_all(~str_replace_all(., "_{2,}", "_")) # replaces multiple underscores to single underscores
  
  # Cleaning up the modifications 
  # Using gsub to clean SILAC label
  #data_dda$base_sequence <- gsub("\\(", "[", data_dda$base_sequence)
  #data_dda$base_sequence <- gsub(")", "]", data_dda$base_sequence)
  
  # filter for FDR? is this how you do it?
  #data_dda <- data_dda[data_dda$qvalue < 0.01,]
  
  # Data wrangling to create long form of data
  tidy_dda <- data_dda %>% 
    select(base_sequence, file_name) %>% 
    rename(peptideseq = "base_sequence")
  
  # Clean up and fix raw file names
  tidy_dda$file_name <- gsub(".raw", "", tidy_dda$file_name)
  tidy_dda$file_name <- gsub("_20200316162818", "", tidy_dda$file_name)  # fix the file that was rerun
  tidy_dda$file_name <- gsub("dia_007", "dda_007", tidy_dda$file_name)  # fix the file that was mislabeled
  
  
  tidy_dda$isotope <- gsub(".*(\\+8.014).*|.*(\\+10.008).*", "heavy", tidy_dda$peptideseq)
  
  
  # label peptides as light or heavy
  tidy_dda$isotope <- gsub(".*(\\+8.014).*|.*(\\+10.008).*", "heavy", tidy_dda$peptideseq)
  tidy_dda$isotope[which(tidy_dda$isotope != "heavy")] = "light"
  
  # Adding peptidemodseq column to match DIA convention
  tidy_dda$peptidemodseq <- paste(tidy_dda$peptideseq, tidy_dda$isotope, sep = "_")
  
  tidy_dda <- tidy_dda %>%
    select(file_name, peptideseq, isotope) %>%
    distinct(file_name, peptideseq, isotope, .keep_all = TRUE)
  
  return(tidy_dda)
}



# Formatting DIA data -----------------------------------------------------

process_dia <- function(elib_files){
  # DIA Peptide quant data
  # Output from EncyclopeDIA (Quant Report)
  
  # function to retrieve detection table from elib database
  data_dia <- data.frame()
  for (infile in elib_files) {
    con <- dbConnect(drv = RSQLite::SQLite(), 
                     dbname=infile)
    
    temp_df <- dbGetQuery(conn=con, statement=paste("SELECT * FROM peptidescores"))
    
    dbDisconnect(con)
    
    data_dia <- rbind(data_dia, temp_df)
  }
  
  # formatting names 
  data_dia <- data_dia %>% 
    rename_all(tolower) %>% 
    rename_all(~str_replace_all(., "\\W", "_")) %>% 
    rename(file_name = "sourcefile")
  
  # Cleaning the raw file name
  data_dia$file_name <- gsub(".mzML", "", data_dia$file_name)
  
  # label peptides as light or heavy
  data_dia <- data_dia %>% 
    mutate(peptidemodseq = str_remove_all(peptidemodseq, "\\[\\+57.0214635]"), # Removing carbamidomethylation
           heavy = str_detect(peptidemodseq, "\\[\\+")) %>% # using [+ to denote heavy peptides
    mutate(peptidemodseq = ifelse(heavy == TRUE, paste0(peptideseq, "_heavy"), paste0(peptideseq, "_light")))
  
  # # label detected peptides as light or heavy
  # data_dia$peptidemodseq <- paste(data_dia$peptidemodseq, "light", sep="_")
  # data_dia$peptidemodseq <- gsub(c("\\[\\+10.008269\\]_light|\\[\\+8.014199\\]_light$"), "_heavy", data_dia$peptidemodseq)
  
  # Data wrangling to create a long form of the dia data
  tidy_dia <- data_dia %>% 
    filter(isdecoy == 0) %>% 
    select(file_name, peptidemodseq) %>%
    distinct(file_name, peptidemodseq, .keep_all = TRUE) %>%
    separate(peptidemodseq, c("peptideseq", "isotope"), "_")
  
  return(tidy_dia)
}



# Single-shot detections ------------------------------------------------------

dia_inputlist <- dir("../../data/dda_vs_dia/", pattern = "*.mzML.elib", 
                     full.names = TRUE, ignore.case = TRUE)

# Read in the single-shot DDA and DIA data
data_detect <- bind_rows(process_dda_metamorpheus("../../data/dda_vs_dia/results_metamorpheus/Task1-silac_multiplex/AllPeptides.psmtsv"), 
                  process_dia(dia_inputlist))

data_detect <- data_detect %>%
  select(file_name, peptideseq, isotope) %>%
  distinct(file_name, peptideseq, isotope, .keep_all = TRUE) %>%
  separate(file_name, c("date", "inst", "operator", "exp", 
                        "sample", "acquisition", "run"), "_") %>% 
  mutate(isotope = factor(isotope, levels = c("light", "heavy")),
         sample = factor(sample, levels = c("light", "mix", "heavy")))

data_detect$acquisition <- toupper(data_detect$acquisition)


## COMPARE PROTEOME DISCOVERER TO METAMORPHEUS
#dda_pd <- process_dda("data/dda_vs_dia/20200315_QEHFX_lkp_silacdia_detections_dda_PSMs.txt")
#dda_mm <- process_dda_metamorpheus("data/dda_vs_dia/results_metamorpheus/Task1-silac_multiplex/AllPeptides.psmtsv")
#dda_pd$acquisition <- "DDA (PD)"
#dda_mm$acquisition <- "DDA (MM)"

#data_detect <- bind_rows(dda_pd, dda_mm)

#data_detect <- data_detect %>%
#  select(file_name, peptideseq, isotope, acquisition) %>%
#  distinct(file_name, peptideseq, isotope, .keep_all = TRUE) %>%
#  separate(file_name, c("date", "inst", "operator", "exp", "sample", "acq", "run"), "_") %>% 
#  mutate(isotope = factor(isotope, levels = c("light", "heavy")),
#         sample = factor(sample, levels = c("light", "mix", "heavy")))

data_detect <- data_detect[complete.cases(data_detect), ]  # should all be complete, but just in case

# suggestion to color columns by rep1, rep2, rep3 rather than acquisition file
# requires new column with rep label
run <- unique(data_detect$run)
reps <- c("rep3", "rep1", "rep2", "rep1", "rep1", "rep2", "rep2", 
          "rep3", "rep3", "rep1", "rep2", "rep3", "rep1", "rep2", 
          "rep3", "rep1", "rep2", "rep3")
runrepmap <- as.data.frame(cbind(run, reps))
data_detect <- merge(data_detect, runrepmap, by="run")


# Fractionation detections ------------------------------------------------------

gpf_inputlist <- dir("../../data/fractionation/", pattern = "*.elib", 
                     full.names = TRUE, ignore.case = TRUE)

# Read in the fractionated DDA and DIA data
data_frx <- bind_rows(process_dda("../../data/fractionation/20200315_QEHFX_lkp_silacdia_hprp__PSMs.txt"), 
                         process_dia(gpf_inputlist))

data_frx <- data_frx %>%
  select(file_name, peptideseq, isotope) %>%
  distinct(file_name, peptideseq, isotope, .keep_all = TRUE) %>%
  separate(file_name, c("date", "inst", "operator", "exp", 
                        "sample", "acquisition", "run"), "_") %>% 
  mutate(isotope = factor(isotope, levels = c("light", "heavy")),
         sample = factor(sample, levels = c("hprp", "gpf")))

data_frx$sample <- gsub("gpf", "GPF-DIA", data_frx$sample)
data_frx$sample <- gsub("hprp", "HPRP-DDA", data_frx$sample)

data_frx <- within(data_frx, 
                   sample = factor(sample, 
                                   levels = names(sort(table(sample), decreasing = TRUE))))



# Windowing detections ----------------------------------------------------

windowing_inputlist <- dir("../../data/windowing/", pattern = "*.elib", 
                     full.names = TRUE, ignore.case = TRUE)

# Read in the window scheme comparison DIA data
data_windows <- process_dia(windowing_inputlist)

data_windows <- data_windows %>%
  select(file_name, peptideseq, isotope) %>%
  distinct(file_name, peptideseq, isotope, .keep_all = TRUE) %>%
  separate(file_name, c("date", "inst", "operator", "exp", 
                        "sample", "acquisition", "method", "run"), "_") %>% 
  mutate(isotope = factor(isotope, levels = c("light", "heavy")),
         method = factor(method, levels = c("25x24mzol", "50x12mz",
                                            "75x8mz", "75x8mzol"))) %>%
  filter(method %in% c("75x8mz", "75x8mzol"))  # remove incorrect method acquisitions

data_windows <- within(data_windows, 
                       method = factor(method, 
                                        levels = names(sort(table(method), decreasing = TRUE))))

# suggestion to color columns by rep1, rep2, rep3 rather than acquisition file
# requires new column with rep label
run <- unique(data_windows$run)
reps <- c("rep1", "rep2", "rep3", "rep1", "rep2", "rep3", "rep1", "rep2", "rep3", "rep1", "rep2", "rep3")
runrepmap <- as.data.frame(cbind(run, reps))
data_windows <- merge(data_windows, runrepmap, by="run")


# Plots -------------------------------------------------------------------

# color blind palette
cbPalette <- c("#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")

# bar graph counting the number of peptides measured in single-shot samples
singleshot_detect <- ggplot(data_detect) +
  facet_grid(rows = vars(sample), cols = vars(acquisition)) +
  geom_bar(aes(x = isotope, fill = reps), 
           position = "dodge") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none") +
  scale_fill_manual(values=cbPalette) +
  labs(x = "Isotope type",
       y = "Detections at 1% FDR",
       fill = NULL)


# bar graph counting the number of peptides measured in fractionated samples
frx_detect <- ggplot(data_frx) +
  geom_bar(aes(x = sample, fill = isotope), 
           position = "dodge") +
  theme_bw(base_size = 14) +
  theme(legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6) ) +
  labs(x = "",
    y = "Detections at 1% FDR",
    fill = NULL)


# cumulative number of peptides measured in each fraction
frx_cumsum <- ggplot(data_frx %>% 
         group_by(sample, isotope, run) %>%
         mutate(count_pep = length(!duplicated(peptideseq))) %>%
         select(sample, isotope, run, count_pep) %>%
         distinct(.keep_all = TRUE) %>%
         group_by(sample, isotope) %>%
         mutate(frx_ranked = order(order(count_pep, decreasing = TRUE))) %>%
         group_by(sample, isotope) %>%
         arrange(frx_ranked) %>%
         mutate(cumsumpep = cumsum(count_pep)),
       aes(x = frx_ranked, y = cumsumpep, color = factor(isotope))) + 
  facet_grid(cols = vars(sample)) +
  geom_point() +
  geom_line() + 
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_blank(),
        legend.position = "none") +
  labs(x = "Fraction",
    y = "Cumulative detections \nat 1% FDR across fractions",
    fill = NULL)


# Patchwork to piece each element into one figure
singleshot_detect + (frx_detect / frx_cumsum) +
  plot_annotation(tag_levels = "A")


# bar graph counting the number of peptides measured in window scheme samples
window_detect <- ggplot(data_windows) +
  facet_grid(cols = vars(method)) +
  geom_bar(aes(x = isotope, fill = reps), 
           position = "dodge") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none") +
  scale_fill_manual(values=cbPalette) +
  labs(x = "Isotope type",
       y = "Detections at 1% FDR",
       fill = NULL)

# Patchwork to piece each element into one figure
singleshot_detect / window_detect +
  plot_layout(heights = c(3, 1)) +
  plot_annotation(tag_levels = "A")

ggsave(filename = "../../figures/fig02_detections.svg", width = 7, height =14)
