# Document Setup ----------------------------------------------------------


options(stringsAsFactors = FALSE)
library(tidyverse)
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
    select(`File Name`, sequence, peptidemodseq, 
           isotope, abundance, master_protein_accessions) %>%
    distinct(`File Name`, sequence, peptidemodseq, 
             isotope, abundance, master_protein_accessions, .keep_all = TRUE)

  tidy_dda$abundance <- tidy_dda$abundance + 1
  
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
  
  return(tidy_dda)
  
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
metadata_df <- import("data/bortezomib/20200828_QEHFX_lkp_pSILAC-DIA_btz_annotations.csv", 
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
dda_input <- process_dda_mm("data/bortezomib/dda_mm/Task1-btz silac turnover/AllQuantifiedPeaks.tsv", metadata_df)

# data from Skyline-daily (20.1.xx) quantifications
# Encyclopedia (0.9.0) matrix for peptide-to-protein group mapping
skyline_quant <- import("data/bortezomib/dia/20200828_QEHFX_lkp_pSILAC-DIA_btz_SKYLINEQUANTEXPORT.csv", 
                        setclass = "tibble")
encyc_matrix <- import("data/bortezomib/dia/20200828_QEHFX_lkp_pSILAC-DIA_btz_QUANT.elib.peptides.txt", 
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
ggsave(filename = "figures/btz_model_plots/exploratory_frxincorporation.png")

rm(plot5_dda, plot5_dia);gc()

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
      fraction ~ (exp(-(k_deg * time_hrs))),  # frx unlabeled
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
  data.frame(time_hrs = seq(0, 500, 2)) %>% 
    add_predictions(x)
}

# Function to fit peptide models to get peptide degradation time
fit_nonlinear <- function(data_input){
  
  # nest and filter the dataframe
  data_model <- data_input %>% 
    filter(!is.na(fraction)) %>% 
    group_by(sequence, btz_conc, master_protein_accessions) %>% 
    nest() %>% 
    mutate(n = map_dbl(data, nrow),
           time_n = map_dbl(data, time_count)) %>%
    filter(time_n > 8)  # ADDED TIME_N FILTER
  
  # run the nonlinear least squares model across the dataframe
  data_model <- data_model %>%
    mutate(nls_model = map(data, non_linear_least_squares)) %>% 
    filter(nls_model != "NA")
  
  return(data_model)
}

model_dda <- fit_nonlinear(dda_input)
model_dia <- fit_nonlinear(dia_input)

# Function to unnest the model predictions
extract_pred <- function(data_model){
  data_best_fit_line <- data_model %>% 
    mutate(best_fit = map(nls_model, best_fit_line)) %>% 
    unnest(best_fit) %>% 
    select(-data, -nls_model) %>% 
    ungroup(btz_conc) %>%
    mutate(btz_conc = str_replace_all(btz_conc, "1000pM", "bortezomib"),
           btz_conc = factor(btz_conc, levels = c("DMSO", "bortezomib")))
  
  return(data_best_fit_line)
}

pred_dda <- extract_pred(model_dda)
pred_dia <- extract_pred(model_dia)

# color blind palette
cbPalette <- c("#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")

## PLOT: Line plot of the nls model predictions
plot_predictions <- function(data_best_fit_line){
  ggplot(data_best_fit_line, aes(x = time_hrs, y = pred)) +
    geom_line(aes(group = sequence, color = btz_conc), alpha = 0.1) +
    scale_colour_manual(values=cbPalette) +
    theme_pander() +
    theme(text = element_text(size=26),
          axis.text.x = element_text(size=20),
          axis.text.y = element_text(size=20),
          axis.title.x = element_text(size=20),
          axis.title.y = element_text(size=20),
          axis.text=element_text(size=20),
          axis.title=element_text(size=20)) +
    facet_wrap(~btz_conc, nrow=2, ncol=1) +
    guides(color = FALSE) +
    labs(#title = "HFF turnover rates",
      y = "Fraction light/total",
      x = "Time (hours)")
}

plot4_dda <- plot_predictions(pred_dda) + 
  labs(title = "DDA") + 
  theme(plot.title = element_text(hjust = 0.5, size=36))
plot4_dia <- plot_predictions(pred_dia) + 
  labs(title = "DIA") +
  theme(plot.title = element_text(hjust = 0.5, size=36))

plot4_dda / plot4_dia +
  plot_annotation(tag_levels = "A")
ggsave(filename = "figures/btz_model_plots/btz_model_turnover_rate_profile.svg", width = 10, height = 14)
#rm(pred_dda, pred_dia);gc()


# Function to unnest the model results
extract_tidydata <- function(data_model){
  model_tidy <- data_model %>% 
    mutate(tidy_model = map(nls_model, tidy)) %>% 
    unnest(tidy_model) %>% 
    select(-data, -nls_model) %>% 
    # estimate is the value of k_deg, this converts k_deg to hours
    mutate(half_life = log(2)/estimate) %>%
    ungroup(btz_conc) %>%
    mutate(btz_conc = str_replace_all(btz_conc, "pM", " pM"))
  
  return(model_tidy)
}

tidy_dda <- extract_tidydata(model_dda)
tidy_dia <- extract_tidydata(model_dia)
rm(model_dda, model_dia);gc()

## PLOT: Density plot of delta half life
plot_deltahalf_density <- function(data_delta){
  ggplot(data_delta) +
    geom_density(aes(x = delta_half)) +
    scale_color_brewer(palette = "Set1") +
    theme_bw(base_size = 14) +
    expand_limits(x = c(-106, 106)) +
    labs(x = expression(Delta~half~life~("BTZ - DMSO")),
         color = "Bortezomib\nconcentration",
         title = "Global half life changes",
         subtitle = "Peptide-centric")
}

# calculate the change (delta) in half life
calculate_halflife <- function(model_tidy){
  data_delta <- model_tidy %>% 
    select(btz_conc, sequence, master_protein_accessions, half_life) %>% 
    spread(btz_conc, half_life) %>% 
    gather(btz_conc, half_life, contains("pM")) %>% 
    mutate(delta_half = half_life - DMSO)
  
  return(data_delta)
  
}

plot_dda <- plot_deltahalf_density(calculate_halflife(tidy_dda)) + labs(subtitle = "DDA")
plot_dia <- plot_deltahalf_density(calculate_halflife(tidy_dia)) + labs(subtitle = "DIA")

plot_dda / plot_dia +
  plot_annotation(tag_levels = "A")
ggsave(filename = "figures/btz_model_plots/delta_half_life_peptide_centric.png")
rm(plot_dda, plot_dia);gc()

## PLOT: Half life distribution
plot_halflife_dist <- function(model_tidy){
  ggplot(model_tidy %>% filter(p.value < 0.05)) +
    geom_density(aes(x = half_life/24)) +
    scale_color_brewer(palette = "Dark2") +
    theme_light(base_size = 14) +
    labs(title = "Half life distribution", 
         x = "Half Life (days)",
         color = "Bortezomib\nconcentration")
}

plot2_dda <- plot_halflife_dist(tidy_dda) + labs(subtitle = "DDA")
plot2_dia <- plot_halflife_dist(tidy_dia) + labs(subtitle = "DIA")

plot2_dda / plot2_dia +
  plot_annotation(tag_levels = "A")
ggsave(filename = "figures/btz_model_plots/btz_density_half_life_days.png")
rm(plot2_dda, plot2_dia);gc()

## PLOT: Turnover rate distribution
plot_turnover <- function(model_tidy){
  ggplot(model_tidy %>% filter(p.value < 0.05)) +
    geom_density(aes(x = estimate)) +
    theme_light(base_size = 14) +
    scale_color_brewer(palette = "Dark2") +
    scale_x_log10() +
    labs(title = "Turnover Rate distribution", 
         x = "k (1/hour)",
         color = "Bortezomib\nconcentration")
}

plot3_dda <- plot_turnover(tidy_dda) + labs(subtitle = "DDA")
plot3_dia <- plot_turnover(tidy_dia) + labs(subtitle = "DIA")

plot3_dda / plot3_dia +
  plot_annotation(tag_levels = "A")
ggsave(filename = "figures/btz_model_plots/btz_density_turnover_rate_k.png")
rm(plot3_dda, plot3_dia);gc()

# Statistics --------------------------------------------------------------

## use the tidy_d*a data to find proteins with differential halflives 
## use the distribution of *peptide* halflives for each protein
## and compare between the two conditions using linear model

diff_halflife <- function(tidy_data, uniprot_map){
  
  # first, melt the dataframe to long form 
  tidy_long <- melt(tidy_data, c("master_protein_accessions", "sequence", "btz_conc", "half_life"))
  
  tidy_long <- tidy_long %>% 
    select("master_protein_accessions", "sequence", "btz_conc", "half_life") %>%
    distinct(.) %>%
    mutate(log2_halflife = ifelse(is.na(half_life), log2(.Machine$double.eps), log2(half_life))) %>%  # TODO what does it mean if half_life is NA? means that k_deg was NA?
    mutate(btz_conc = factor(btz_conc, levels = c("bortezomib", "DMSO")))
  
  # define the linear model performed on categorical data (DMSO/btz conditions), e.g. an ANOVA
  tidy_anova <- function(data){
    
    # needs exception check bc some proteins don't have measurements for both conditions
    tryCatch(
      lm(log2_halflife ~ btz_conc,  # linear model itself
         data = data),
      error = function(e) NA,
      warning = function(w) NA
    )
  }
  
  # This step performs the model
  quant_nest <- tidy_long %>% 
    group_by(master_protein_accessions) %>% 
    nest() %>% 
    mutate(model_data = map(data, tidy_anova)) %>%
    filter(!is.na(model_data))
  
  # This chunk cleans up the output of the model we just performed
  model_data <- quant_nest %>% 
    mutate(model_tidy = map(model_data, tidy)) %>% 
    select(-model_data) %>% 
    unnest(model_tidy) %>% 
    unnest(data) %>%
    filter(term != "(Intercept)") %>%
    # now need to group by protein, calculate average delta half life
    group_by(master_protein_accessions, btz_conc,estimate, std.error, statistic, p.value) %>% 
    summarise(log2_prot_halflife = mean(log2_halflife)) %>% 
    spread(btz_conc,log2_prot_halflife)  %>% 
    mutate(log2fc_halflife = bortezomib - DMSO) %>%
    ungroup() %>%
    mutate(adj.p.value = p.adjust(p.value, method = "BH")) %>%   
    ungroup()

  # map master_protein_accessions from uniprot to gene names
  #model_data$gene_name <- ConvertID(model_data$master_protein_accessions, 
  #                                  ID_from = "ACC+ID", 
  #                                  ID_to = "GENENAME", 
  #                                  directorypath = NULL)[,2]
  
  model_data <- model_data %>%
    inner_join(.,uniprot_map, by="master_protein_accessions") %>%
    rename(gene_name = "Entry name")
  
  return(model_data)
}

##
## VOLCANO PLOTS
##

uniprot_annotations <- read_tsv(file = "data/uniprot-proteome_UP000005640.tab")
uniprot_annotations$`Entry name` <- gsub("_HUMAN", "", uniprot_annotations$`Entry name`)
uniprot_annotations <- uniprot_annotations %>%
  rename(master_protein_accessions = Entry)

diff_dda <- diff_halflife(tidy_dda, uniprot_annotations)
diff_dia <- diff_halflife(tidy_dia, uniprot_annotations)

extract_model_factors <- function(data_in, shared_proteins){
  data_out <- data_in %>%
    #filter(master_protein_accessions %in% shared_proteins)  %>% 
    select(master_protein_accessions, btz_conc, estimate, half_life, time_n) %>%
    group_by(master_protein_accessions, btz_conc) %>%
    spread(btz_conc, half_life) %>%
    mutate(peptide_n = n(), mean_estimate = mean(estimate),
           min_time_n = min(time_n), max_time_n = max(time_n)) %>%
    #ungroup() %>%
    distinct(master_protein_accessions, min_time_n, max_time_n, peptide_n) %>%
    distinct(master_protein_accessions, min_time_n, max_time_n, peptide_n)
  
  return(data_out)
}

modeled_proteins <- intersect(diff_dda$master_protein_accessions, 
                              diff_dia$master_protein_accessions)

modelfactors_dda <- extract_model_factors(tidy_dda, modeled_proteins)
modelfactors_dia <- extract_model_factors(tidy_dia, modeled_proteins)




diff_dia <- diff_dia %>%
  left_join(modelfactors_dia, by="master_protein_accessions")

diff_dda <- diff_dda %>%
  left_join(modelfactors_dda, by="master_protein_accessions")


write.csv(diff_dda, file = "results/btz_results_dda.csv", row.names = FALSE)
write.csv(diff_dia, file = "results/btz_results_dia.csv", row.names = FALSE)



plot_pvalues <- function(diff_data){
  
  # Make a p value histogram of the model results
  pvalue_plot <- ggplot(diff_data) +  # using the "model_data" data set...
    geom_histogram(aes(x = p.value), bins = 30) +  # ... make a histogram of the pvalues
    geom_vline(xintercept = 0.05, linetype = "dashed", color = "red", size = 1)  # add a red line at p=0.05
  
  return(pvalue_plot)
  
}

plot_volcano_fc <- function(diff_data){


  volcano_fc <- ggplot(diff_data, aes(x = (log2fc_halflife), y = -log10(adj.p.value))) +  
    geom_point() +
    theme_bw(base_size = 14) +
    geom_text_repel(aes(label = ifelse(adj.p.value < 0.01, gene_name, ""))) +
    labs(x = expression(Log[2]~(hours[btz]/hours[DMSO])),
         y = expression(-Log[10]~p.value~(fdr~corrected))) +
    xlim(-17,17)
    
  
  return(volcano_fc)
}

plot_volcano_delta <- function(diff_data){
  
  volcano_delta <- ggplot(diff_data, aes(x = (2^(bortezomib) - 2^(DMSO)), y = -log10(adj.p.value))) +  
    geom_point() +
    theme_bw()  +
    geom_text_repel(aes(label = ifelse(adj.p.value < 0.01, gene_name, ""), 
                        size=18),
                    box.padding = 1,
                    show.legend = FALSE) + 
    theme(text = element_text(size=20),
          axis.text.x = element_text(size=18),
          axis.text.y = element_text(size=18)) +
    labs(x = "change in half life (btz - DMSO)",
         y = expression(-Log[10]~p.value~(fdr~corrected))) +
    xlim(-750,750)
  
}


volcano_dda <- plot_volcano_delta(diff_dda)
volcano_dia <- plot_volcano_delta(diff_dia)

volcano2_dda <- plot_volcano_fc(diff_dda)
volcano2_dia <- plot_volcano_fc(diff_dia)

# Patchwork to piece each element into one figure
volcano_dda / volcano_dia +
  plot_annotation(tag_levels = "A")

ggsave(filename = "figures/btz_model_plots/volcano.svg", width = 7, height = 10)

# Patchwork to piece each element into one figure
(plot4_dda + plot4_dia) / (volcano2_dda + volcano2_dia) +
  plot_layout() + 
  plot_annotation(tag_levels = "A")

ggsave(filename = "figures/fig04_bortezomib.svg", width = 14, height = 14)



## halflife histograms per Birgit's suggestion
# Make a histogram of the half lives
plot_halflife <- function(diff_data){
  
  diff_data <- subset(diff_data, adj.p.value < 0.05)
  
  halflife_hist <- ggplot(diff_data) +  
    geom_histogram(aes(x = DMSO),
                   fill=cbPalette[1], color=cbPalette[1], alpha=0.5, bins=100) +  
    geom_histogram(aes(x = bortezomib),
                   fill=cbPalette[2], color=cbPalette[2], alpha=0.5, bins=100) +
    theme_bw(base_size = 14) +
    xlim(2, 16)+
    labs(#title = "Turnover Rate distribution", 
      x = "Halflife (hours)",
      labels = c('DMSO', 'Bortezomib'))

  return(halflife_hist)
  
}

halflives_dda <- plot_halflife(diff_dda) + labs(title = "DDA")
halflives_dia <- plot_halflife(diff_dia) + labs(title = "DIA")

halflives_dda / halflives_dia +
  plot_layout() + 
  plot_annotation(tag_levels = "A")
ggsave(filename = "figures/btz_model_plots/halflives_hist.png", width = 10, height = 7)


# plot_predictions pred_dia, pred_dda
sig_models_dda <- plot_predictions(subset(pred_dda, 
                                          master_protein_accessions %in% subset(diff_dda, adj.p.value < 0.5)$master_protein_accessions))
sig_models_dia <- plot_predictions(subset(pred_dia, 
                        master_protein_accessions %in% subset(diff_dia, adj.p.value < 0.5)$master_protein_accessions))

sig_models_dda + sig_models_dia

# Comparison of DDA vs DIA----------------------------------------------------------

#
# protein level comparisons
#
plot_halflifecorr <- function(halflife_df){
  ggplot(halflife_df, aes(x = DMSO_dia, y = DMSO_dda)) +
    geom_point(alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, color = 'red') +
    #xlim(-50, 50) + ylim(-50, 50) +
    theme_pander(base_size = 14) +
    scale_color_brewer(palette = "Dark2") +
    theme_pander(base_size = 14) +
    #geom_text_repel(aes(label = ifelse(adj.p.value_dia < 0.01, gene_name_dia, ""), color="blue")) +
    #geom_text_repel(aes(label = ifelse(adj.p.value_dda < 0.01, gene_name_dda, ""), color="green")) +
    #geom_point(aes(color = ifelse(adj.p.value_dda < 0.01, "green", "black"), alpha = 0.5)) +
    #geom_point(aes(color = ifelse(adj.p.value_dia < 0.01, "blue", "black"), alpha = 0.5)) +
    #facet_wrap(~btz_conc) +
    guides(color = FALSE) +
    labs(#title = "Delta (bortezomib/DMSO) Half-life",
         y = "halflife, DMSO (SILAC-DDA)",
         x = "halflife, DMSO (SILAC-DIA)")
}

plot_pvalcorr <- function(halflife_df){
  ggplot(halflife_df, aes(x = adj.p.value_dia, y = adj.p.value_dda)) +
    geom_point(alpha = 0.5) +
    #geom_abline(slope = 1, intercept = 0, color = 'red') +
    geom_hline(yintercept = 0.05, color = 'red') +
    geom_vline(xintercept = 0.05, color = 'red') +
    #xlim(-50, 50) + ylim(-50, 50) +
    theme_pander(base_size = 14) +
    scale_color_brewer(palette = "Dark2") +
    geom_text_repel(aes(label = ifelse(adj.p.value_dia < 0.01, gene_name_dia, ""), color="blue")) +
    geom_text_repel(aes(label = ifelse(adj.p.value_dda < 0.01, gene_name_dda, ""), color="green")) +
    #geom_rect(aes(xmin=0, xmax=0.05, ymin=0, ymax=0.05, color="black", alpha=0.5)) +
    #geom_point(aes(color = ifelse(adj.p.value_dda < 0.01, "green", "black"), alpha = 0.5)) +
    #geom_point(aes(color = ifelse(adj.p.value_dia < 0.01, "blue", "black"), alpha = 0.5)) +
    #facet_wrap(~btz_conc) +
    guides(color = FALSE) +
    labs(#title = "Delta (bortezomib/DMSO) Half-life",
      y = "adjusted p value (SILAC-DDA)",
      x = "adjusted p value (SILAC-DIA)")
}


diff_all <- merge(diff_dda, diff_dia, 
                  by="master_protein_accessions", 
                  suffixes=c("_dda", "_dia"))

plot_halflifecorr(diff_all)
ggsave(filename = "figures/btz_model_plots/halflife_correlation.png", width = 10, height = 7)
plot_pvalcorr(diff_all)
ggsave(filename = "figures/btz_model_plots/pvalue_correlation.png", width = 10, height = 7)


#
# peptide level comparisons
#

plot_corr_peptide <- function(tidy_df){
  ggplot(tidy_df, aes(x = half_life_dia, y = half_life_dda)) +
    geom_point(alpha = 0.25) +
    #scale_color_gradient("Spectral", name = "light abundance at time 0 (decile)") +
    geom_abline(slope = 1, intercept = 0, color = 'red') +
    geom_smooth(method = "lm", se = FALSE, color = 'gray') +
    xlim(0, 150) + ylim(0, 150) +
    theme_pander(base_size = 14) +
    #scale_color_brewer(palette = "Dark2") +
    #geom_text_repel(aes(label = ifelse(adj.p.value_dia < 0.01, gene_name_dia, ""), color="blue")) +
    #geom_text_repel(aes(label = ifelse(adj.p.value_dda < 0.01, gene_name_dda, ""), color="green")) +
    #geom_point(aes(color = ifelse(adj.p.value_dda < 0.01, "green", "black"), alpha = 0.5)) +
    #geom_point(aes(color = ifelse(adj.p.value_dia < 0.01, "blue", "black"), alpha = 0.5)) +
    #facet_wrap(~btz_conc) +
    #guides(color = FALSE) +
    labs(#title = "Delta (bortezomib/DMSO) Half-life",
      y = "peptide halflife, DMSO (SILAC-DDA)",
      x = "peptide halflife, DMSO (SILAC-DIA)")
}

plot_abund_halflife <- function(tidy_df){
  temp1 <- ggplot(tidy_df) +
    geom_boxplot(aes(x = factor(bin_dda), y = log2(half_life_dda), 
                     fill="green", alpha=0.5)) +
    scale_fill_identity(name = 'the fill', guide = 'legend',labels = c('DIA', 'DDA')) +
    labs(#title = "Delta (bortezomib/DMSO) Half-life",
      x = "abundance of light peptide at time zero (deciles)",
      y = "log2 peptide halflife, DMSO")
  temp2 <- ggplot(tidy_df) +
    geom_boxplot(aes(x = factor(bin_dia), y = log2(half_life_dia), 
                     fill="blue", alpha=0.5)) +
    scale_fill_identity(name = 'the fill', guide = 'legend',labels = c('DIA', 'DDA')) +
    labs(#title = "Delta (bortezomib/DMSO) Half-life",
      x = "abundance of light peptide at time zero (deciles)",
      y = "log2 peptide halflife, DMSO")
  temp <- temp1 / temp2
  return(temp)
}

plot_ntime_halflife <- function(tidy_df){
  temp1 <- ggplot(tidy_df) +
    geom_boxplot(aes(x = factor(n_dda), y = log2(half_life_dda), 
                     fill="green", alpha=0.5)) +
    geom_jitter(aes(x = factor(n_dia), y = log2(half_life_dia)), 
                    position=position_jitter(width=.1, height=0), alpha=0.1) +
    scale_fill_identity(name = 'the fill', guide = 'legend',labels = c('DDA')) +
    labs(#title = "Delta (bortezomib/DMSO) Half-life",
      x = "timepoints available for halflife modeling",
      y = "log2 peptide halflife, DMSO")
  temp2 <- ggplot(tidy_df) +
    geom_boxplot(aes(x = factor(n_dia), y = log2(half_life_dia), 
                     fill="blue", alpha=0.5)) +
    scale_fill_identity(name = 'the fill', guide = 'legend',labels = c('DIA')) +
    geom_jitter(aes(x = factor(n_dia), y = log2(half_life_dia)), 
                    position=position_jitter(width=.1, height=0), alpha=0.1) +
    labs(#title = "Delta (bortezomib/DMSO) Half-life",
      x = "timepoints available for halflife modeling",
      y = "log2 peptide halflife, DMSO")
  temp <- temp1 / temp2
  return(temp)
}

tidy_all <- merge(tidy_dda, tidy_dia, 
                  by="sequence", 
                  suffixes=c("_dda", "_dia"))
tidy_all <- merge(tidy_all,
                  subset(dda_input, dda_input$`File Name` == "20200828_QEHFX_lkp_pSILAC-DIA_btz_01_DDA_155", select=c("sequence", "light")), 
                  by = "sequence")
tidy_all <- rename(tidy_all, c("light_dda"="light"))
tidy_all <- merge(tidy_all,
                  subset(dia_input, dia_input$`File Name` == "20200828_QEHFX_lkp_pSILAC-DIA_btz_01_DIA_156", select=c("sequence", "light")), 
                  by = "sequence")
tidy_all <- rename(tidy_all, c("light_dia"="light"))
tidy_all <- tidy_all %>%
  mutate(bin_dda = floor(rank(light_dda) * 10 / (length(light_dda) + 1)),
         bin_dia = floor(rank(light_dia) * 10 / (length(light_dia) + 1)))

plot_corr_peptide(tidy_all)
ggsave(filename = "figures/btz_model_plots/halflife_correlation_peptide.png", width = 10, height = 7)
plot_abund_halflife(tidy_all)
ggsave(filename = "figures/btz_model_plots/halflife_abundance_correlation.png", width = 10, height = 7)
plot_ntime_halflife(tidy_all)
ggsave(filename = "figures/btz_model_plots/halflife_timepoint_correlation.png", width = 10, height = 7)


density_halflife_dda <- ggplot(tidy_all) + 
  geom_density(aes(x=half_life_dda, fill=factor(btz_conc_dda)), alpha=0.5) +
  theme_pander(base_size = 14) +
  labs(title = "DDA",
    x = "Protein half life (hours)") +
  xlim(-50, 200)
density_halflife_dia <- ggplot(tidy_all) + 
  geom_density(aes(x=half_life_dia, fill=factor(btz_conc_dia)), alpha=0.5) +
  labs(title = "DIA",
       x = "Protein half life (hours)") +
  theme_pander(base_size = 14) +
  xlim(-50, 200)

density_halflife_dda_zoom <- ggplot(tidy_all) + 
  geom_density(aes(x=half_life_dda, fill=factor(btz_conc_dda)), alpha=0.5) +
  theme_pander(base_size = 14) +
  labs(title = "DDA",
       x = "Protein half life (hours)") +
  xlim(100, 200)
density_halflife_dia_zoom <- ggplot(tidy_all) + 
  geom_density(aes(x=half_life_dia, fill=factor(btz_conc_dia)), alpha=0.5) +
  labs(title = "DIA",
       x = "Protein half life (hours)") +
  theme_pander(base_size = 14) +
  xlim(100, 200)
  
temp1 <- density_halflife_dda / density_halflife_dia  +
  plot_layout() + 
  plot_annotation(tag_levels = "A")
temp2 <- density_halflife_dda_zoom / density_halflife_dia_zoom +
  plot_layout() + 
  plot_annotation(tag_levels = "A")
temp1 + temp2
ggsave(filename = "figures/btz_model_plots/density_halflife.png", width = 10, height = 14)


## stat checks
length(unique(dda_input$master_protein_accessions))
length(unique(diff_dda$master_protein_accessions))
length(unique(dia_input$master_protein_accessions))
length(unique(diff_dia$master_protein_accessions))

(sum(is.na(dda_input$fraction))/nrow(dda_input))
(sum(is.na(dia_input$fraction))/nrow(dia_input))


