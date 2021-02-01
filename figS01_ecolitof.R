
# Setup -------------------------------------------------------------------

options(stringsAsFactors = FALSE)
library(tidyverse)
library(rio)
library(ggthemes)
library(DBI)
library(RSQLite)
library(patchwork)

# Formatting DIA data -----------------------------------------------------

# DIA Peptide quant data
# Output from Skyline post Encyclopedia
dia_df <- import("../../data/ecoli_tof/schilling_swath_silac-ecoli_6600.csv")

# formatting names 
data_dia <- dia_df %>% 
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

tidy_dia <- tidy_dia %>%
  separate(file_name, c("temp1", "temp2", "temp3", "temp4", "ratio_l", "ratio_h", 
                        "temp7", "temp8", "temp9", "temp10", "temp11", "temp12", "temp13"),
           sep = "_", remove = FALSE, convert = TRUE) %>%
  select(-contains("temp")) %>%
  mutate(ratio_id = ratio_l/ratio_h)

tidy_dia$dataset <- "MS2"

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
  separate(file_name, c("temp1", "temp2", "temp3", "temp4", "ratio_l", "ratio_h", 
                        "temp7", "temp8", "temp9", "temp10", "temp11", "temp12", "temp13"),
           sep = "_", remove = FALSE, convert = TRUE) %>%
  select(-contains("temp")) %>%
  mutate(ratio_id = ratio_l/ratio_h)

tidy_ms1$dataset <- "MS1"



# Combining datasets ------------------------------------------------------

# merge the MS2 and MS1 datasets together
data <- rbind(tidy_dia, tidy_ms1)

# There are duplicate entries in the dataset. This is likely due to protein grouping
# Removing the duplicate entry
data <- data %>% 
  select(file_name, peptide, peptide_modified_sequence, ratio_id,
         isotope, abundance,dataset) %>%
  distinct(file_name, peptide, peptide_modified_sequence, ratio_id,
           isotope, abundance, .keep_all = TRUE)

# Calculating fractional abundance 
# Light / (Heavy + Light)
data_fraction <- data %>%
  group_by(file_name, ratio_id, peptide, isotope, dataset) %>% 
  summarise(abundance = mean(abundance)) %>% 
  spread(isotope, abundance) %>%
  mutate(hl_ratio_log10 = log10(light / heavy),
         fraction = light / (heavy + light)) %>% 
  filter(!is.na(fraction))


# High Conf peptides ------------------------------------------------------

# Selecting for peptides with high heavy label incorporation
# The HeLa cells were grown in culture for 6 cell passages. 
# 6 passages is the "recommended" number of passages for full isotopic labeling 
# The assumption is that the heavy labeled sample should be ~100% heavy labeled.
# Peptides not meeting this criteria are filtered out
# I used the median value of heavy incorporation for filtering criteria
high_conf_pep <- data_fraction %>% 
  filter(ratio_id == 20) %>% 
  select(-heavy, -light, -hl_ratio_log10) %>% 
  group_by(peptide, dataset) %>% 
  summarize(fraction_mean = mean(fraction, na.rm = TRUE),
            n = n()) %>% 
  ungroup()

temp <- high_conf_pep %>% 
  filter(n == 3, 
         dataset == "MS2",
         fraction_mean > .90) %>%   # is this necessary?
  select(peptide)


high_conf_pep <- high_conf_pep %>% 
  # filter(n == 3) %>%
  select(-n) %>% 
  spread(dataset, fraction_mean) %>% 
  filter(!is.na(MS2))

# filtered data
# found in three technical replicates
# heavy incorporation in 100% heavy sample > 95.7% (median value)
data_high_conf <- inner_join(temp, data_fraction)

rm(temp);gc()


# Wide data ---------------------------------------------------------------


data_wide <- data_fraction %>% 
  group_by(peptide, ratio_id, dataset) %>% 
  summarize(fraction_mean = mean(fraction, na.rm = TRUE),
            fraction_n = n()) %>% 
  ungroup() %>% 
  spread(dataset, fraction_mean)


# MA plot -----------------------------------------------------------------

# selecting the heavy log2 abundance in the 100% heavy sample of the DDA and DIA
# Calculating the mean heavy abundance per peptide
# Calculating the sd of heavy abundance per peptide
# Calculating the CV of heavy abundance per peptide
light_log2 <- data_fraction %>% 
  filter(ratio_id == 20,
         !is.na(light)) %>% 
  rename(light_100 = "light") %>% 
  select(peptide, dataset, light_100) %>% 
  group_by(peptide, dataset) %>% 
  summarize(light_mean = mean(light_100, na.rm = TRUE),
            light_sd = sd(light_100, na.rm = TRUE),
            light_n = n()) %>% 
  ungroup() %>% 
  mutate(light_cv = light_sd / light_mean * 100)

# Joining the mean heavy abundance with the complete dataset
data_ma <- data_fraction %>% 
  filter(!is.na(hl_ratio_log10)) %>% 
  inner_join(light_log2, .)

data_ma$ratio_id <- factor(round(data_ma$ratio_id, digits = 3))

# Subsetting data for the DIA plots
ma_ms1 <- data_ma %>% 
  filter(dataset == "MS1") 
         #heavy_ug == c(1.000, 0.100, 0.010, 0.001)) #, 
         #heavy_cv < 20,
         #heavy_mean > 0.90,
         #heavy_n == 3) 

ma_ms2 <- data_ma %>% 
  filter(dataset == "MS2") 
#heavy_ug == c(1.000, 0.100, 0.010, 0.001)) #, 
#heavy_cv < 20,
#heavy_mean > 0.90,
#heavy_n == 3) 





# Plots ---------------------------------------------------------------

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "darkseagreen", "darksalmon")

# DIA MS1 plot
ms1_1 <- ggplot(ma_ms1) +
  geom_point(aes(x = log2(light_mean), y = hl_ratio_log10, color = factor(ratio_id)), 
             size = 2, alpha = 0.1) +
  geom_smooth(aes(x = log2(light_mean), y = hl_ratio_log10, 
                  group = factor(ratio_id)), se = FALSE, color = "darkgray", 
              method = "loess", linetype = "dashed") + 
  geom_hline(aes(yintercept = log10(0.05), color = factor("0.05"))) +
  geom_hline(aes(yintercept = log10(0.1), color = factor("0.1"))) +
  geom_hline(aes(yintercept = log10(0.2), color = factor("0.2"))) +
  geom_hline(aes(yintercept = log10(0.333333333333333), color = factor("0.333"))) +
  geom_hline(aes(yintercept = log10(1), color = factor("1"))) +
  geom_hline(aes(yintercept = log10(3), color = factor("3"))) +
  geom_hline(aes(yintercept = log10(5), color = factor("5"))) +
  geom_hline(aes(yintercept = log10(10), color = factor("10"))) +
  geom_hline(aes(yintercept = log10(20), color = factor("20"))) +
  geom_text(data=data.frame(x=20,y=log10(0.05)), aes(x, y, color = factor("0.05")), 
            label="1:20", vjust=-1) +
  geom_text(data=data.frame(x=20,y=log10(0.1)), aes(x, y, color = factor("0.1")), 
            label="1:10", vjust=-1) +
  geom_text(data=data.frame(x=20,y=log10(0.2)), aes(x, y, color = factor("0.2")), 
            label="1:5", vjust=-1) +
  geom_text(data=data.frame(x=20,y=log10(0.333333333333333)), aes(x, y, color = factor("0.333")), 
            label="1:3", vjust=-1) +
  geom_text(data=data.frame(x=20,y=log10(1)), aes(x, y, color = factor("1")), 
            label="1:1", vjust=-1) +
  geom_text(data=data.frame(x=20,y=log10(3)), aes(x, y, color = factor("3")), 
            label="3:1", vjust=-1) +
  geom_text(data=data.frame(x=20,y=log10(5)), aes(x, y, color = factor("5")), 
            label="5:1", vjust=-1) +
  geom_text(data=data.frame(x=20,y=log10(10)), aes(x, y, color = factor("10")), 
            label="10:1", vjust=-1) +
  geom_text(data=data.frame(x=20,y=log10(20)), aes(x, y, color = factor("20")), 
            label="20:1", vjust=-1) +
  scale_y_continuous(limits = c(-2.5,2.5)) +
  scale_x_continuous(limits = c(1, 20)) +
  theme_classic(base_size = 12) +
  labs(y = expression(Log[10]~(Light/Heavy)),
       x = expression(Log[2]~DIA~Light~abundance)) +
  guides(color = FALSE) +
  scale_color_manual(values=cbPalette)

# DIA boxplot
ms1_2 <- ggplot(ma_ms1) +
  geom_boxplot(aes(x = factor(ratio_id), y = hl_ratio_log10, color = factor(ratio_id)), alpha = 0.1) +
  geom_hline(aes(yintercept = log10(0.05), color = factor("0.05"))) +
  geom_hline(aes(yintercept = log10(0.1), color = factor("0.1"))) +
  geom_hline(aes(yintercept = log10(0.2), color = factor("0.2"))) +
  geom_hline(aes(yintercept = log10(0.333333333333333), color = factor("0.333"))) +
  geom_hline(aes(yintercept = log10(1), color = factor("1"))) +
  geom_hline(aes(yintercept = log10(3), color = factor("3"))) +
  geom_hline(aes(yintercept = log10(5), color = factor("5"))) +
  geom_hline(aes(yintercept = log10(10), color = factor("10"))) +
  geom_hline(aes(yintercept = log10(20), color = factor("20"))) +
  scale_color_manual(values=cbPalette) +
  scale_y_continuous(limits = c(-2.5,2.5)) +
  theme_void() +
  guides(color = FALSE) +
  labs(color = "Ratio (Light/Heavy)") + 
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank())

# DIA MA plot
ms2_1 <- ggplot(ma_ms2) +
  geom_point(aes(x = log2(light_mean), y = hl_ratio_log10, color = factor(ratio_id)), 
             size = 2, alpha = 0.1) +
  geom_smooth(aes(x = log2(light_mean), y = hl_ratio_log10, 
                  group = factor(ratio_id)), se = FALSE, color = "darkgray", 
              method = "loess", linetype = "dashed") + 
  geom_hline(aes(yintercept = log10(0.05), color = factor("0.05"))) +
  geom_hline(aes(yintercept = log10(0.1), color = factor("0.1"))) +
  geom_hline(aes(yintercept = log10(0.2), color = factor("0.2"))) +
  geom_hline(aes(yintercept = log10(0.333333333333333), color = factor("0.333"))) +
  geom_hline(aes(yintercept = log10(1), color = factor("1"))) +
  geom_hline(aes(yintercept = log10(3), color = factor("3"))) +
  geom_hline(aes(yintercept = log10(5), color = factor("5"))) +
  geom_hline(aes(yintercept = log10(10), color = factor("10"))) +
  geom_hline(aes(yintercept = log10(20), color = factor("20"))) +
  geom_text(data=data.frame(x=20,y=log10(0.05)), aes(x, y, color = factor("0.05")), 
            label="1:20", vjust=-1) +
  geom_text(data=data.frame(x=20,y=log10(0.1)), aes(x, y, color = factor("0.1")), 
            label="1:10", vjust=-1) +
  geom_text(data=data.frame(x=20,y=log10(0.2)), aes(x, y, color = factor("0.2")), 
            label="1:5", vjust=-1) +
  geom_text(data=data.frame(x=20,y=log10(0.333333333333333)), aes(x, y, color = factor("0.333")), 
            label="1:3", vjust=-1) +
  geom_text(data=data.frame(x=20,y=log10(1)), aes(x, y, color = factor("1")), 
            label="1:1", vjust=-1) +
  geom_text(data=data.frame(x=20,y=log10(3)), aes(x, y, color = factor("3")), 
            label="3:1", vjust=-1) +
  geom_text(data=data.frame(x=20,y=log10(5)), aes(x, y, color = factor("5")), 
            label="5:1", vjust=-1) +
  geom_text(data=data.frame(x=20,y=log10(10)), aes(x, y, color = factor("10")), 
            label="10:1", vjust=-1) +
  geom_text(data=data.frame(x=20,y=log10(20)), aes(x, y, color = factor("20")), 
            label="20:1", vjust=-1) +
  scale_y_continuous(limits = c(-2.5,2.5)) +
  scale_x_continuous(limits = c(1, 20)) +
  theme_classic(base_size = 12) +
  labs(y = expression(Log[10]~(Light/Heavy)),
       x = expression(Log[2]~DIA~Light~abundance)) +
  guides(color = FALSE) +
  scale_color_manual(values=cbPalette)

# DIA boxplot
ms2_2 <- ggplot(ma_ms2) +
  geom_boxplot(aes(x = factor(ratio_id), y = hl_ratio_log10, color = factor(ratio_id)), alpha = 0.1) +
  geom_hline(aes(yintercept = log10(0.05), color = factor("0.05"))) +
  geom_hline(aes(yintercept = log10(0.1), color = factor("0.1"))) +
  geom_hline(aes(yintercept = log10(0.2), color = factor("0.2"))) +
  geom_hline(aes(yintercept = log10(0.333333333333333), color = factor("0.333"))) +
  geom_hline(aes(yintercept = log10(1), color = factor("1"))) +
  geom_hline(aes(yintercept = log10(3), color = factor("3"))) +
  geom_hline(aes(yintercept = log10(5), color = factor("5"))) +
  geom_hline(aes(yintercept = log10(10), color = factor("10"))) +
  geom_hline(aes(yintercept = log10(20), color = factor("20"))) +
  scale_color_manual(values=cbPalette) +
  scale_y_continuous(limits = c(-2.5,2.5)) +
  theme_void() +
  guides(color = FALSE) +
  labs(color = "Ratio (Light/Heavy)") + 
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank())

# Plotting the figure
ms1_1 + ms1_2 + ms2_1 + ms2_2 +
  plot_layout(widths = c(1,0.5),
              heights = c(4,4)) +
  plot_annotation(tag_levels = "A")


ggsave(filename = "../../figures/figS01_ecolitof.svg", width = 8, height = 10)


## calculate stats for n
#length(unique(ma_dia[ma_dia$heavy_ug == '1',]$peptide))
#length(unique(ma_dia[ma_dia$heavy_ug == '0.1',]$peptide))
#length(unique(ma_dia[ma_dia$heavy_ug == '0.01',]$peptide))
#length(unique(ma_dia[ma_dia$heavy_ug == '0.001',]$peptide))


##
## REVIEWER SUGGESTIONS
##

# make a plot comparing the MS1 versus MS2 coefficient of variation (CV)
cv_df <- data_fraction %>%
  group_by(dataset, ratio_id, peptide) %>%
  mutate(fraction_cv = (sd(fraction)/mean(fraction))) %>%
  select(dataset, ratio_id, peptide, fraction_cv) %>%
  distinct(dataset, ratio_id, peptide, fraction_cv)

ggplot(cv_df, aes(x=fraction_cv, fill=dataset, color=dataset)) +
  geom_histogram(alpha = 0.5, bins = 100) +
  facet_grid(rows=vars(ratio_id)) +
  scale_x_continuous(limits = c(-0.1, 1)) +
  theme_classic(base_size = 12)

ggsave(filename = "../../figures/figS0x_ecolitof_cvs.svg", width = 8, height = 10)
  