# Document Setup ----------------------------------------------------------
set.seed(1234)

options(stringsAsFactors = FALSE)
library(tidyverse)
library(dplyr)
library(rio)
library(lubridate)
library(broom)
library(modelr)
library(ggthemes)
library(pander)
library(ggridges)
library(patchwork)
library(ggrepel)
library(reshape2)
library(UniprotR)


# Functions for data import and processing  -------------------------------

process_dda_pd <- function(input_files, data){
  
  # Formatting Input files
  # Formatting names and subsetting columns
  input_files <- input_files %>% 
    rename_all(tolower) %>% 
    rename_all(~str_replace_all(., "\\W", "_")) %>% 
    select(file_id, file_name, creation_date) %>% 
    mutate(file_id = tolower(file_id),
           creation_date = mdy_hms(creation_date),
           file_name = str_replace_all(file_name, "pm", "pM"))
  
  
  # Selecting raw file name from file path
  input_files$file_name <- unlist(lapply(input_files$file_name, function(x){
    unlist(strsplit(x, split = "\\", fixed = TRUE))[8]
  }))
  
  # Cleaning raw file name
  input_files$file_name <- unlist(lapply(input_files$file_name, function(x){
    unlist(strsplit(x, split = "\\."))[1]
  }))
  
  input_files <- input_files %>% 
    separate(file_name, c("temp1", "temp2", "btz_conc", "time", "sample_id"), sep = "_", convert = TRUE, remove = FALSE) %>% 
    select(-contains("temp"))
  
  
  # Formatting data
  # Formatting names and subsetting columns
  data <- data %>% 
    rename_all(tolower) %>% # converts to lowercase
    rename_all(~str_replace_all(., "\\W", "_")) %>% # substitutes spaces to underscore
    select(-contains("abundance_ratio"),
           -contains("abundances"),
           -contains("found_in_")) # removes columns containing these names
  
  
  # Cleaning protein/peptide modifications ----------------------------------
  
  
  # Cleaning up the modifications 
  # Using regex to clean SILAC label
  data$modifications <- gsub("; ", ";", data$modifications)
  data$modifications <- gsub("K\\d{,2};", "", data$modifications)
  data$modifications <- gsub("R\\d{,2};", "", data$modifications)
  data$modifications <- gsub("K\\d{,2}", "", data$modifications)
  data$modifications <- gsub("R\\d{,2}", "", data$modifications)
  data$modifications <- gsub("\\[]", "", data$modifications)
  data$modifications <- gsub("[1-9]xLabel:13C\\(6)15N\\(4)", "", data$modifications)
  data$modifications <- gsub("[1-9]xLabel:13C\\(6)15N\\(2)", "", data$modifications)
  data$modifications <- gsub("; $", "", data$modifications)
  data$modifications <- gsub(" $", "", data$modifications)
  
  # Cleaning up Carbamidomethylation
  data$modifications <- gsub("[1-9]xCarbamidomethyl ", "", data$modifications)
  
  # Cleaning up Oxidation
  data$modifications <- gsub("[1-9]xOxidation ", "", data$modifications)
  
  # N-term acetyl
  data$modifications <- gsub("[1-9]xAcetyl \\[N-Term]", "N_ac", data$modifications)
  
  # N-term Gln to pyro-Glu
  data$modifications <- gsub("[1-9]xGln->pyro-Glu \\[N-Term]", "N_QpE", data$modifications)
  
  # N-term pyro-Carbamidomethylation
  data$modifications <- gsub("[1-9]xPyro-carbamidomethyl \\[N-Term]", "N_pCarb", data$modifications)
  
  # N-term Glu to pyro-Glu
  data$modifications <- gsub("[1-9]xGlu->pyro-Glu \\[N-Term]", "N_EpE", data$modifications)
  
  
  # Removing miscellaneous characters
  data$modifications <- gsub("];\\[", ";", data$modifications)
  data$modifications <- gsub("^;", "", data$modifications)
  data$modifications <- gsub(";$", "", data$modifications)
  data$modifications <- gsub("^\\[", "", data$modifications)
  data$modifications <- gsub("]$", "", data$modifications)
  data$modifications <- gsub("\\[", "", data$modifications)
  data$modifications <- gsub("]", "", data$modifications)
  
  
  # Counting the number of cysteines and carbamidomethylation
  # Identifying N-term acetyl mods -- used to remove
  data <- data %>% 
    mutate(c_pep_count = str_count(sequence, "C"),
           c_mod_count = str_count(modifications, "C"),
           m_ox = str_detect(modifications, "M"),
           nterm_mod = str_detect(modifications, "N_"))
  
  # verifying that all cysteines are carbamidomethylated
  all(data$c_pep_count == data$c_mod_count)
  
  
  # Data Wrangling ----------------------------------------------------------
  
  
  # Creating a long form of the data
  data_tidy <- data %>% 
    filter(contaminant == FALSE, # filtering out contaminants
           m_ox == FALSE, # filtering out peptides with methionine oxidation
           nterm_mod == FALSE # filtering out peptides with n-terminal modifications
    ) %>%  
    select(sequence, master_protein_accessions, 
           contains("abundance")) %>% 
    gather(temp, abundance, contains("abundance")) %>%
    separate(temp, c("temp2", "file_id", "isotope", "temp3")) %>% 
    right_join(input_files, .) %>% 
    select(-contains("temp")) %>% 
    filter(!is.na(abundance))
  
  
  # Generating model input data ---------------------------------------------
  
  
  # Calculating fraction Heavy / Total
  data_ria <- data_tidy %>% 
    filter(!is.na(abundance)) %>%
    filter(btz_conc == "DMSO" | btz_conc == "1000pM") %>%
    select(-file_id, -file_name, -creation_date) %>% 
    spread(isotope, abundance) %>% 
    mutate(total = ifelse(is.na(light) & !is.na(heavy), heavy, (heavy + light)),  # switched for light degrad
           fraction = ifelse(is.na(light), NA, (light / total)),  # switched for light degrad
           has_fraction = ifelse(is.na(fraction), FALSE, TRUE)) %>%
    
    #data_ria$fraction[data_ria$time == 0] <- 0
    #data_ria$fraction[data_ria$time == "Inf"] <- 1
    #data_ria$time[data_ria$time == "Inf"] <- 500
    
    return(data_ria)
  
}

process_dda_mm <- function(quant_file, annotations){
  
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
  
  # label peptides as light or heavy
  tidy_dda$isotope <- gsub(".*(\\+8.014).*|.*(\\+10.008).*", "heavy", tidy_dda$peptideseq)
  
  # remove the modification from the peptideseq variable
  tidy_dda$peptideseq <- gsub("\\(\\+8.014\\)|\\(\\+10.008\\)", "", tidy_dda$peptideseq)
  
  # remove the "missed cleavage" annotation
  tidy_dda$peptideseq <- gsub("[[:punct:]]", "", tidy_dda$peptideseq)
  
  
  # label peptides as light or heavy
  tidy_dda$isotope[which(tidy_dda$isotope != "heavy")] = "light"
  
  # Adding peptidemodseq column to match DIA convention
  tidy_dda$peptidemodseq <- paste(tidy_dda$peptideseq, tidy_dda$isotope, sep = "_")
  
  tidy_dda <- tidy_dda %>%
    rename(abundance = "peak_intensity", sequence = "peptideseq", 
           #peptide_modified_sequence = "full_sequence",
           `File Name` = "file_name", master_protein_accessions = "protein_group") %>%
    mutate(abundance = if_else(is.na(abundance), 0, abundance)) %>% 
    select(`File Name`, sequence, peptidemodseq, 
           isotope, abundance, master_protein_accessions) %>%
    distinct(`File Name`, sequence, peptidemodseq, 
             isotope, abundance, master_protein_accessions, .keep_all = TRUE)
  
  #tidy_dda$abundance <- tidy_dda$abundance + 1
  
  # Calculating fraction Heavy / Total
  data_ria <- tidy_dda %>% 
    #filter(!is.na(abundance)) %>%
    select(`File Name`, `sequence`, `isotope`, `abundance`, `master_protein_accessions`) %>%
    distinct(`File Name`, `sequence`, `isotope`, `abundance`, `master_protein_accessions`, .keep_all = TRUE) %>%
    group_by(`File Name`, `sequence`, `isotope`, `master_protein_accessions`) %>%  # TODO: why are there duplicate peptides? maybe modified+unmodified forms?
    summarise(abundance = sum(abundance)) %>%  # collapse duplicate rows
    spread(`isotope`, `abundance`) %>%  # reshape isotope label quants into same row
    mutate(total = ifelse(is.na(light) & !is.na(heavy), heavy, (heavy + light)),  # switched for light degrad
           fraction = ifelse(is.na(light), NA, (light / total)),  # switched for light degrad
           has_fraction = ifelse(is.na(fraction), FALSE, TRUE)) %>%
    inner_join(.,annotations)  # map file names to metadata annotations
  
  
  return(data_ria)
  
  #return(tidy_dda)
  
}

process_dia <- function(skyline_df, annotations, peptoprot_df){
  
  # build peptide to protein group lookup table using Encyclopedia results
  peptoprot_df <- peptoprot_df %>%
    select(Peptide, Protein) %>%
    mutate(Peptide = str_remove_all(Peptide, "\\[\\+57.021464]")) %>%
    distinct(Peptide, Protein, .keep_all = TRUE) %>% 
    rename(master_protein_accessions = Protein,
           `Peptide Modified Sequence` = Peptide)
  
  peptoprot_df$master_protein_accessions <- gsub('sp\\|', '', peptoprot_df$master_protein_accessions)  # remove swissprot identifier from proteins
  peptoprot_df$master_protein_accessions <- gsub('\\|.*_HUMAN;', ';', peptoprot_df$master_protein_accessions)  # remove gene names for all 'internal' proteins in group
  peptoprot_df$master_protein_accessions <- gsub('\\|.*_HUMAN', '', peptoprot_df$master_protein_accessions)  # remove gene name from only/last protein in group
  peptoprot_df$master_protein_accessions <- gsub(';', '; ', peptoprot_df$master_protein_accessions)  # format to match conventions

  temp <- skyline_df %>%
    rename(`File Name` = Replicate) %>%
    select(`File Name`, `Peptide Modified Sequence`, `Isotope Label Type`, `Total Area Fragment`) %>%
    mutate(`Total Area Fragment` = if_else(is.na(`Total Area Fragment`), 0, `Total Area Fragment`)) %>% 
    mutate(`Peptide Modified Sequence` = str_remove_all(`Peptide Modified Sequence`, "\\[\\+57]")) %>%  # Removing carbamidomethylation
    distinct(`File Name`, `Peptide Modified Sequence`, `Isotope Label Type`, `Total Area Fragment`, .keep_all = TRUE) %>%
    spread(`Isotope Label Type`, `Total Area Fragment`)  %>% # reshape isotope label quants into same row
    mutate(total = ifelse(is.na(light) & !is.na(heavy), heavy, (heavy + light)),  # switched for light degrad
           fraction = ifelse(is.na(light), NA, (light / total)),  # switched for light degrad
           has_fraction = ifelse(is.na(fraction), FALSE, TRUE)) %>%
    inner_join(.,annotations) %>%  # map file names to "time" and "btz_conc"
    inner_join(.,peptoprot_df) %>%  # map peptides to protein groups
    rename(sequence = `Peptide Modified Sequence`)
  
  return(temp)
}


# Data Import -------------------------------------------------------------

# metadata file mapping File Name (w/o extension) to required annotations
metadata_df <- import("../../data/curve/20200828_QEHFX_lkp_pSILAC-DIA_curve_annotations.csv", 
                      setclass = "tibble")

# data from PD 2.2
# MSPepSearch, Sequest and MSAmanda
#data_file <- import("data/bortezomib/20190115_HFF_DDA_Bortezomib_Turnover-(1)_PeptideGroups.txt", 
#                    setclass = "tibble")
#file_names <- import("data/bortezomib/20190115_HFF_DDA_Bortezomib_Turnover-(1)_InputFiles.txt", 
#                     setclass = "tibble")

# Rolfs: When you search as "turnover", the program automatically splits the intensity of 
# peptides with missed cleavages that have both a heavy and a light label. If you want to see 
# if any such peptides were identified/quantified, you need to look int the 
# "AllQuantifiedPeaks.psmtsv" output file instead of the "AllQuantifiedPeptides.psmtsv" file.
dda_input <- process_dda_mm("../../data/curve/2020-10-23-00-11-30/Task1-turnover curve/AllQuantifiedPeaks.tsv", metadata_df)

# data from Skyline-daily (20.1.xx) quantifications
# Encyclopedia (0.9.0) matrix for peptide-to-protein group mapping
skyline_quant <- import("../../data/curve/QUANTIFICATION - DIA MS1 and MS2.csv", 
                        setclass = "tibble")
encyc_matrix <- import("../../data/curve/20200828_QEHFX_lkp_pSILAC-DIA_curve_QUANT.elib.peptides.txt", 
                       setclass = "tibble")

dia_input <- process_dia(skyline_quant, metadata_df, encyc_matrix)
rm(skyline_quant, metadata_df, encyc_matrix);gc()



## PLOT: Check the fractional incorporation stacked density plot
plot_incorporation <- function(data_input){
  ggplot(data_input) +
    geom_density_ridges(aes(y = factor(time_hrs), x = (1 - fraction), fill = factor(time_hrs)), scale = 4) +
    theme_pander() +
    scale_fill_brewer(palette = "Spectral") +
    guides(fill = FALSE) +
    labs(title = "Fraction heavy incorporation over time",
         y = "Time (hours)",
         x = "Fraction Heavy")
}

plot5_dda <- plot_incorporation(dda_input) + labs(subtitle = "DDA")
plot5_dia <- plot_incorporation(dia_input) + labs(subtitle = "DIA")

plot5_dda / plot5_dia +
  plot_annotation(tag_levels = "A")
ggsave(filename = "../../figures/ratios_model_plots/exploratory_frxincorporation.png")

rm(plot5_dda, plot5_dia);gc()


#############################################
## FILTER FOR ABUNDANCE ASSUMPTIONS
## UNION OF LIST FOR DDA AND DIA TO BE FAIR
#dda_keep <- subset(dda_input, (time_hrs == 500) & (fraction > 0.45) & (fraction < 0.55))
#dia_keep <- subset(dia_input, (time_hrs == 500) & (fraction > 0.45) & (fraction < 0.55))
#dda_keep <- subset(dda_input, (time_hrs == 1000) & (fraction > 0.9))
#dia_keep <- subset(dia_input, (time_hrs == 1000) & (fraction > 0.9))
#newlist <- c(dda_keep, dia_keep)
#pep_keep <- unique(newlist[['sequence']])
#dda_input <- subset(dda_input, sequence %in% unique(pep_keep))
#dia_input <- subset(dia_input, sequence %in% unique(pep_keep))

#############################################
## SPLIT EACH DDA AND DIA INTO REPLICATE 1 AND REPLICATE 2
# randomly pick 1 file from each concentration point

rep1_dda <- dda_input %>% 
  group_by(time_hrs) %>% 
  sample_n(1) %>% 
  select(`File Name`) %>% 
  distinct(`File Name`)
rep1_dda <- rep1_dda$`File Name`
rep2_dda <- dda_input %>% 
  filter(!`File Name` %in% rep1_dda) %>%
  group_by(time_hrs) %>% 
  sample_n(1) %>% 
  select(`File Name`) %>% 
  distinct(`File Name`)
rep2_dda <- rep2_dda$`File Name`

rep1_dia <- dia_input %>% 
  group_by(time_hrs) %>% 
  sample_n(1) %>% 
  select(`File Name`) %>% 
  distinct(`File Name`)
rep1_dia <- rep1_dia$`File Name`
rep2_dia <- dia_input %>% 
  filter(!`File Name` %in% rep1_dia) %>%
  group_by(time_hrs) %>% 
  sample_n(1) %>% 
  select(`File Name`) %>% 
  distinct(`File Name`)
rep2_dia <- rep2_dia$`File Name`


# Fitting the model -------------------------------------------------------

# Function to count the number of time points for each protein 
time_count <- function(x){
  length(unique(x$time_hrs))
}

# Function to count the number of peptides for each protein
peptide_count <- function(x){
  length(unique(x$peptide))
}

# Function used to model turnover
non_linear_least_squares <- function(data){
  tryCatch(
    nls(
      #fraction ~ (1 - exp(-k_deg * time)),  # frx labeled (heavy)
      fraction ~ (exp((-1*k_deg)*time_hrs)),  # frx unlabeled
      data = data,
      control = nls.control(maxiter = 1000, minFactor = 1/1024, warnOnly = TRUE),
      start = list(
        k_deg = abs(getInitial(fraction ~ SSasymp(time_hrs, Asym, lrc, R0 = 0), data = data)[[2]])
      )
    ), 
    error = function(e) NA, 
    warning = function(w) NA
  )
}

# Function that creates a dataframe with model predictions
best_fit_line <- function(x){
  data.frame(time_hrs = seq(0, 1000, 2)) %>% 
    add_predictions(x)
}

# Function to fit peptide models to get peptide degradation time
fit_nonlinear <- function(data_input){
  
  # nest and filter the dataframe
  data_model <- data_input %>% 
    filter(!is.na(fraction)) %>% 
    filter(fraction > 0) %>%  # new! filter finite fractions
    group_by(sequence, master_protein_accessions) %>% 
    nest() %>% 
    mutate(n = map_dbl(data, nrow),
           time_n = map_dbl(data, time_count)) #%>%
    #filter(time_n > 8)  # ADDED TIME_N FILTER <- REMOVE FILTER, DROPS TOO MANY PROTEINS
  
  # run the nonlinear least squares model across the dataframe
  data_model <- data_model %>%
    mutate(nls_model = map(data, non_linear_least_squares)) %>% 
    filter(nls_model != "NA")
  
  return(data_model)
}

model_dda_rep1 <- fit_nonlinear((dda_input %>% filter(`File Name` %in% rep1_dda)))
model_dda_rep2 <- fit_nonlinear((dda_input %>% filter(`File Name` %in% rep2_dda)))
model_dia_rep1 <- fit_nonlinear((dia_input %>% filter(`File Name` %in% rep1_dia)))
model_dia_rep2 <- fit_nonlinear((dia_input %>% filter(`File Name` %in% rep2_dia)))

# Function to unnest the model predictions
extract_pred <- function(data_model){
  data_best_fit_line <- data_model %>% 
    mutate(best_fit = map(nls_model, best_fit_line)) %>% 
    unnest(best_fit) %>% 
    select(-data, -nls_model)
  
  return(data_best_fit_line)
}

pred_dda_rep1 <- extract_pred(model_dda_rep1)
pred_dda_rep2 <- extract_pred(model_dda_rep2)
pred_dia_rep1 <- extract_pred(model_dia_rep1)
pred_dia_rep2 <- extract_pred(model_dia_rep2)


# color blind palette
cbPalette <- c("#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")

## PLOT: Line plot of the nls model predictions
plot_predictions <- function(data_best_fit_line){
  ggplot(data_best_fit_line, aes(x = time_hrs, y = pred)) +
    geom_line(aes(group = sequence), alpha = 0.1) +
    scale_colour_manual(values=cbPalette) +
    theme_pander() +
    theme(text = element_text(size=26),
          axis.text.x = element_text(size=20),
          axis.text.y = element_text(size=20),
          axis.title.x = element_text(size=20),
          axis.title.y = element_text(size=20),
          axis.text=element_text(size=20),
          axis.title=element_text(size=20)) +
    guides(color = FALSE) +
    labs(#title = "HFF turnover rates",
      y = "Fraction light/total",
      x = "Time (hours)")
}

#plot4_dda <- plot_predictions(pred_dda) + 
#  labs(title = "DDA") + 
#  theme(plot.title = element_text(hjust = 0.5, size=36))
#plot4_dia <- plot_predictions(pred_dia) + 
#  labs(title = "DIA") +
#  theme(plot.title = element_text(hjust = 0.5, size=36))

#plot4_dda / plot4_dia +
#  plot_annotation(tag_levels = "A")
#ggsave(filename = "../../figures/ratios_model_plots/model_turnover_rate_profile.svg", width = 10, height = 14)
#rm(pred_dda, pred_dia);gc()


# Function to unnest the model results
extract_tidydata <- function(data_model){
  model_tidy <- data_model %>% 
    mutate(tidy_model = map(nls_model, tidy)) %>% 
    unnest(tidy_model) %>% 
    select(-data, -nls_model) %>% 
    # estimate is the value of k_deg, this converts k_deg to hours
    mutate(half_life = log(2)/estimate) %>%
    mutate(half_life_stderr = log(2)/std.error)
  
  return(model_tidy)
}

tidy_dda_rep1 <- extract_tidydata(model_dda_rep1)
tidy_dda_rep2 <- extract_tidydata(model_dda_rep2)
tidy_dia_rep1 <- extract_tidydata(model_dia_rep1)
tidy_dia_rep2 <- extract_tidydata(model_dia_rep2)
rm(model_dda, model_dia);gc()

# label the half life columns to merge the replicates into a single df
names(tidy_dda_rep1)[names(tidy_dda_rep1) == 'half_life'] <- 'half_life_rep1'
names(tidy_dda_rep2)[names(tidy_dda_rep2) == 'half_life'] <- 'half_life_rep2'
names(tidy_dia_rep1)[names(tidy_dia_rep1) == 'half_life'] <- 'half_life_rep1'
names(tidy_dia_rep2)[names(tidy_dia_rep2) == 'half_life'] <- 'half_life_rep2'
tidy_dda <- merge(tidy_dda_rep1,tidy_dda_rep2, by='sequence')
tidy_dia <- merge(tidy_dia_rep1,tidy_dia_rep2, by='sequence')

# calculate a min and max stderr for error bars in scatterplot below
# INCORRECT! stderr from lmfit isn't what I think it is
tidy_dda$min_halflife_x <- tidy_dda$half_life_rep1 - tidy_dda$half_life_stderr.x
tidy_dda$min_halflife_y <- tidy_dda$half_life_rep2 - tidy_dda$half_life_stderr.y
tidy_dia$min_halflife_x <- tidy_dia$half_life_rep1 - tidy_dia$half_life_stderr.x
tidy_dia$min_halflife_y <- tidy_dia$half_life_rep2 - tidy_dia$half_life_stderr.y
tidy_dda$max_halflife_x <- tidy_dda$half_life_rep1 + tidy_dda$half_life_stderr.x
tidy_dda$max_halflife_y <- tidy_dda$half_life_rep2 + tidy_dda$half_life_stderr.y
tidy_dia$max_halflife_x <- tidy_dia$half_life_rep1 + tidy_dia$half_life_stderr.x
tidy_dia$max_halflife_y <- tidy_dia$half_life_rep2 + tidy_dia$half_life_stderr.y

#####################################################
## PLOT: Half life correlations between replicates
# to add x and y std err bars:
# https://stackoverflow.com/questions/9231702/ggplot2-adding-two-errorbars-to-each-point-in-scatterplot
plot_halfliferep_corr <- function(tidy_data){
  ggplot(data=tidy_data, aes(x = log2(half_life_rep1/24), y=log2(half_life_rep2/24))) +
    geom_point() +
    geom_smooth(method="lm",aes(x = log2(half_life_rep1/24), y=log2(half_life_rep2/24)),se=F) +
    scale_color_brewer(palette = "Dark2") +
    theme_light(base_size = 14) +
    geom_abline(intercept=0, slope=1, color="red") +  
    labs(x = "Rep 1, log2 Half Life (days)",
         y = "Rep 2, log2 Half Life (days)") #+
    #geom_errorbar(aes(ymin = min_halflife_y, ymax = max_halflife_y)) + 
    #geom_errorbarh(aes(xmin = min_halflife_x, xmax = max_halflife_x))
}

lm_eqn <- function(x,y){
  # browser()
  m <- lm(y ~ x)
  a <- coef(m)[1]
  a <- ifelse(sign(a) >= 0, 
              paste0(" + ", format(a, digits = 4)), 
              paste0(" - ", format(-a, digits = 4))  )
  eq1 <- substitute( paste( italic(y) == b, italic(x), a ), 
                     list(a = a, 
                          b = format(coef(m)[2], digits = 4)))
  eq2 <- substitute( paste( italic(R)^2 == r2 ), 
                     list(r2 = format(summary(m)$r.squared, digits = 3)))
  c( as.character(as.expression(eq1)), as.character(as.expression(eq2)))
}

# calculate R^2 for the correlation
labels_dda <- lm_eqn(log2(tidy_dda$half_life_rep2/24), log2(tidy_dda$half_life_rep1/24))
labels_dia <- lm_eqn(log2(tidy_dia$half_life_rep2/24), log2(tidy_dia$half_life_rep1/24))

# make each correlation plot, plus a few specific parameters to make it prettier
rep_corr_dda <- plot_halfliferep_corr(tidy_dda) + labs(subtitle = "DDA") +
  xlim(min=-15,max=15) + 
  ylim(min=-15,max=15) +
  geom_text(x = 12, y = -9, label = labels_dda[1], parse = TRUE,  check_overlap = TRUE ) +
  geom_text(x = 12, y = -11, label = labels_dda[2], parse = TRUE, check_overlap = TRUE )
rep_corr_dia <- plot_halfliferep_corr(tidy_dia) + labs(subtitle = "DIA") +
  xlim(min=-15,max=15) +
  ylim(min=-15,max=15) +
  geom_text(x = 12, y = -9, label = labels_dia[1], parse = TRUE,  check_overlap = TRUE ) +
  geom_text(x = 12, y = -11, label = labels_dia[2], parse = TRUE, check_overlap = TRUE )

rep_corr_dda / rep_corr_dia +
  plot_annotation(tag_levels = "A", title="Correlation between rep1 and rep2 half lives")
ggsave(filename = "../../figures/ratios_model_plots/rep_correlation_half_life.png", width = 7, height =14)

#####################################################


## PLOT: Half life distribution
plot_halflife_dist <- function(model_tidy){
  ggplot(model_tidy %>% filter(p.value < 0.05)) +
    geom_density(aes(x = log2(half_life/24))) +
    scale_color_brewer(palette = "Dark2") +
    theme_light(base_size = 14) +
    geom_vline(xintercept=log2(500/24)) +  # "known" half life is the 50/50 sample
    labs(title = "Half life distribution", 
         x = "log2 Half Life (days)")
}

plot2_dda <- plot_halflife_dist(tidy_dda) + labs(subtitle = "DDA")
plot2_dia <- plot_halflife_dist(tidy_dia) + labs(subtitle = "DIA")

plot2_dda / plot2_dia +
  plot_annotation(tag_levels = "A")
ggsave(filename = "../../figures/ratios_model_plots/density_half_life_days.png")
rm(plot2_dda, plot2_dia);gc()

## PLOT: Turnover rate distribution
plot_turnover <- function(model_tidy){
  ggplot(model_tidy %>% filter(p.value < 0.05)) +
    geom_density(aes(x = estimate)) +
    theme_light(base_size = 14) +
    scale_color_brewer(palette = "Dark2") +
    scale_x_log10() +
    labs(title = "Turnover Rate distribution", 
         x = "k (1/hour)")
}

plot3_dda <- plot_turnover(tidy_dda) + labs(subtitle = "DDA")
plot3_dia <- plot_turnover(tidy_dia) + labs(subtitle = "DIA")

plot3_dda / plot3_dia +
  plot_annotation(tag_levels = "A")
ggsave(filename = "../../figures/ratios_model_plots/density_turnover_rate_k.png")
rm(plot3_dda, plot3_dia);gc()


