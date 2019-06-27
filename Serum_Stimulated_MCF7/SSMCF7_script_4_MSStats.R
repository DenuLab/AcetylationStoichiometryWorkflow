
# Loading Libraries -------------------------------------------------------
# source("http://bioconductor.org/biocLite.R")
# biocLite("MSstats")
# library(MSstats)

options(stringsAsFactors = FALSE)
library(tidyverse)
library(rio)
library(stringr)
library(ggthemes)
library(MSstats)
#source("../Functions/Function_3_MSstats_dataProcess_3_13_2.R")
library(corrplot) # For correlation matrix plot
library(viridisLite)

# Data import -------------------------------------------------------------


# Spectronaut output
ms_data <- read.csv("./Serum_Stimulated_MCF7/Raw/20181031_SSMCF7_stoich_Spectronaut10_reprocess_correctFASTA_Report.csv")


# Annotation file ---------------------------------------------------------

# Creating dataframe
annotation_file <- data.frame(Run = unique(ms_data$R.FileName),
                              Condition = NA,
                              BioReplicate = NA,
                              Fraction = NA)

annotation_file$Condition <- unlist(lapply(annotation_file$Run, function(x){
  unlist(strsplit(x, split = "_"))[6]
}))

annotation_file$BioReplicate <- as.numeric(unlist(lapply(annotation_file$Run, function(x){
  unlist(strsplit(x, split = "_"))[7]
})))

annotation_file$Fraction <- as.numeric(unlist(lapply(annotation_file$Run, function(x){
  unlist(strsplit(x, split = "_"))[8]
})))


# Preprocessing data ------------------------------------------------------


# Adding required column
ms_data$F.ExcludedFromQuantification <- ifelse(ms_data$F.PossibleInterference == TRUE, "True", "False")

# Filtering contaminants
data <- ms_data %>% 
  mutate(remove = str_detect(.$PG.ProteinGroups, "CON_")) %>% 
  filter(remove == FALSE) %>% 
  select(-remove)


# MSstats converter -------------------------------------------------------


input_data <- SpectronauttoMSstatsFormat(input = data,
                                         annotation = annotation_file,
                                         intensity = 'PeakArea',
                                         filter_with_Qvalue = TRUE,
                                         qvalue_cutoff = 0.01,
                                         useUniquePeptide = TRUE,
                                         fewMeasurements = "remove",
                                         removeProtein_with1Feature = FALSE,
                                         summaryforMultipleRows = max)


# Adding fraction column for the dataProcess function to run properly
input_data$Run <- as.character(input_data$Run)
input_data$Fraction <- as.numeric(unlist(lapply(input_data$Run, function(x){
  unlist(strsplit(x, split = "_"))[8]
})))
input_data$Run <- as.factor(input_data$Run)

# Clearing unnecessary data
rm(ms_data, data, annotation_file)
gc()

# Data Quantification -----------------------------------------------------

quant_data <- dataProcess(input_data,
                          logTrans = 2,
                          normalization = "equalizeMedians",
                          nameStandards = NULL,
                          address = "",
                          fillIncompleteRows = TRUE,
                          featureSubset = "all",
                          remove_noninformative_feature_outlier = FALSE,
                          n_top_feature = 3,
                          summaryMethod = "TMP",
                          equalFeatureVar = TRUE,
                          censoredInt = "NA",
                          cutoffCensored = "minFeature",
                          MBimpute = TRUE,
                          remove50missing = FALSE,
                          maxQuantileforCensored = 0.999,
                          clusters = NULL)


# Extracting the quant data as dataframe
processed_data <- quant_data$ProcessedData
rl_data <- quant_data$RunlevelData


# Plots -------------------------------------------------------------------

dataProcessPlots(quant_data, type = "QCPlot", which.Protein = "allonly",
                 address = "figures/")
dataProcessPlots(quant_data, type = "ProfilePlot", address = "figures/")
dataProcessPlots(quant_data, type = "ConditionPlot", address = "figures/")


# Post processing of data -------------------------------------------------


# Formatting names
names(processed_data) <- tolower(names(processed_data))
names(rl_data) <- tolower(names(rl_data))

## Processed data
# time point 
processed_data$time <- unlist(lapply(processed_data$originalrun, function(x){
  as.numeric(unlist(strsplit(x, split = "_"))[1])
}))

# Biorep
processed_data$biorep <- unlist(lapply(processed_data$originalrun, function(x){
  as.numeric(unlist(strsplit(x, split = "_"))[2])
}))

## Run level data
# time point
rl_data$time <- unlist(lapply(rl_data$originalrun, function(x){
  as.numeric(unlist(strsplit(x, split = "_"))[1])
}))

# Biorep
rl_data$biorep <- unlist(lapply(rl_data$originalrun, function(x){
  as.numeric(unlist(strsplit(x, split = "_"))[2])
}))




# Writing tables ----------------------------------------------------------


# write.csv(processed_data, file = "./Serum_Stimulated_MCF7/Cleaned/msstats/20181101_processed_data.csv", row.names = FALSE)
# write.csv(rl_data, file = "./Serum_Stimulated_MCF7/Cleaned/msstats/20181101_run_level_data.csv", row.names = FALSE)


# Summary  ----------------------------------------------------------------


# Summary table
rl_summary <- rl_data %>% 
  select(protein, logintensities, time, biorep) %>% 
  group_by(protein, time) %>% 
  summarize(mean_abundance = mean(logintensities),
            sd_abundance = sd(logintensities),
            count = n(),
            cv = (sd_abundance/mean_abundance)*100)



# Plots -------------------------------------------------------------------


# Boxplot - Processed data
ggplot(processed_data %>% filter(censored == FALSE)) +
  geom_boxplot(aes(x = factor(biorep), y = abundance, fill = factor(time))) +
  scale_fill_brewer(palette = "Set1") +
  theme_pander(base_size = 14) +
  labs(x = "Replicate",
       y = expression(Log[2]~Feature~abundance),
       fill = "time point")


# Boxplot - run level data
ggplot(rl_data) +
  geom_boxplot(aes(x = factor(biorep), y = logintensities, fill = factor(time))) +
  scale_fill_brewer(palette = "Set1") +
  theme_pander(base_size = 14) +
  labs(x = "Replicate",
       y = expression(Log[2]~Protein~abundance),
       fill = "time point")
#ggsave(filename = "figures/Boxplot_prot_abundance.png")

# CV density plot
ggplot(rl_summary) +
  geom_density(aes(x = cv, color = factor(time))) +
  scale_color_brewer(palette = "Set1") +
  theme_bw(base_size = 14) +
  labs(y = "Density",
       x = "Coefficient of Variation",
       color = "time point")
#ggsave(filename = "figures/DensityPlot_CoefVaration.png")


# Fuzzy c-means data ------------------------------------------------------
### Input for cluster analysis

rl_wide <- rl_summary %>% 
  filter(cv <= 10) %>% 
  select(protein, time, mean_abundance) %>% 
  spread(time, mean_abundance)


#write.csv(rl_wide, file = "./Serum_Stimulated_MCF7/Fuzzy_cm_data/20181101_MSstats_protein_quant_fcm_input_data.csv", row.names = FALSE)
