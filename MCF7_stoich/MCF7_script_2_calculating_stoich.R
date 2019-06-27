
# Loading Libraries -------------------------------------------------------

options(stringsAsFactors = FALSE)
library(tidyverse)
library(stringr)
library(seqinr) # used to read FASTA files
library(ggthemes)
source("./Functions/function_1_Isotopic_distribution_BRAIN.R")
library(BRAIN) # Required for the isotopic distribution function

# Reading in data ---------------------------------------------------------


# Spectronaut output
data <- read.csv("./MCF7_stoich/Raw/20181031_MCF7_stoich_spectronaut10_correctFASTA_Report.csv")

# Spectral library with isotope labels
spec_lib_isotope <- read.csv("./Spectral_Libraries/20181031_tidy_mcf7_spec_lib_spectronaut10_correctFASTA_isotope.csv")

# Human Fasta file
hs_fasta <- read.fasta(file = "./Fasta_files/human_swissprot_canonical_iRT_20171212.fasta", seqtype = "AA")


# Functions ---------------------------------------------------------------


## Function to extract the first id
extract_first_item <- function(x, a = ";", b = 1){
  unlist(lapply(x, function(y){
    unlist(strsplit(y, split = a))[b]
  }))
}


# Formatting data ---------------------------------------------------------


# Annotating the contaminating proteins
data$contaminant <- FALSE
data$contaminant[grep("CON", data$PG.ProteinGroups)] <- TRUE

# changing these columns to logicals
data$F.PossibleInterference[which(data$F.PossibleInterference == "False")] <- FALSE
data$F.PossibleInterference[which(data$F.PossibleInterference == "True")] <- TRUE

# Sample ID
data$sample_id <- substring(data$R.FileName, 17, 22)

# Fraction
data$fraction <- substring(data$R.FileName, 24, 24)

# K count
data$k_count <- str_count(data$EG.StrippedSequence, "K")

# Acetyl count
data$acetyl_count <- str_count(data$EG.ModifiedSequence, "K\\[")

# Fragment Ion Sequence
data$fragment_ion_sequence <- 
  substring(data$EG.StrippedSequence, 
            ifelse(data$F.FrgType == "b", 1, (nchar(data$EG.StrippedSequence) - (data$F.FrgNum - 1))), 
            ifelse(data$F.FrgType == "y", nchar(data$EG.StrippedSequence), data$F.FrgNum))

# Number of lysines on the fragment ion sequence
data$frag_k_count <- str_count(data$fragment_ion_sequence, "K")

# adding an methionine oxidation column
data$oxidation_count <- str_count(data$EG.ModifiedSequence, "M\\[")


# Adding the Isotope Label ------------------------------------------------


spec_lib_isotope2 <- spec_lib_isotope %>% unite(col = "FG.Id",
                                                LabeledPeptide, FG.Charge, 
                                                sep = ".", remove = FALSE)
spec_lib_merge <- spec_lib_isotope2 %>% 
  select(FG.Charge, EG.StrippedSequence, FG.Id, k_count, 
         EG.ModifiedPeptide, fragment_ion_sequence, 
         F.FrgLossType, F.FrgZ, isotope_label)

# Removing duplicate rows
spec_lib_merge <- spec_lib_merge %>% 
  distinct()

# Merging the data and merge dataframes
data <- merge(data, spec_lib_merge)
data <- data %>% 
  arrange(EG.StrippedSequence, fragment_ion_sequence, sample_id)

rm(spec_lib_merge, spec_lib_isotope2)
gc()


# Selecting acetyl fragment ions ------------------------------------------


# selecting for 1K fragment ions
data_1k <- data %>% 
  filter(acetyl_count > 0, 
         frag_k_count == 1,
         EG.DatapointsPerPeak >= 4,
         contaminant == FALSE, 
         oxidation_count == 0,
         F.PredictedRelativeIntensity < 95,
         F.PeakArea > 10)

# Subsetting the columns
data_1k <- data_1k %>% 
  select(sample_id, fraction, PG.Genes, PG.ProteinDescriptions, PG.ProteinGroups, 
         EG.DatapointsPerPeak, PEP.NrOfMissedCleavages, FG.ShapeQualityScore, 
         EG.StrippedSequence, EG.ModifiedSequence, 
         FG.PrecMz, k_count, fragment_ion_sequence, 
         F.FrgIon, F.FrgType, F.FrgNum, F.FrgMz, 
         F.FrgLossType, F.FrgZ, isotope_label, EG.Qvalue, F.PeakArea) %>%
  arrange(sample_id, EG.StrippedSequence, -FG.ShapeQualityScore)

# Selecting the top id for protein groups - placed this function here to save time
data_1k[3:5] <- apply(data_1k[3:5], 2, extract_first_item)

# Removing duplicate rows
data_1k <- data_1k %>% 
  distinct(sample_id, fraction, PG.Genes, PG.ProteinDescriptions, PG.ProteinGroups, 
           PEP.NrOfMissedCleavages, k_count, EG.StrippedSequence, EG.ModifiedSequence,  
           fragment_ion_sequence, F.FrgIon, F.FrgType, F.FrgNum, 
           F.FrgLossType, F.FrgIon, F.FrgZ, isotope_label, .keep_all = TRUE)


# Annotating K-site -------------------------------------------------------


# Generates a dataframe with the column header names
hs_fasta_df <- data.frame("Fasta.header" = names(hs_fasta))

# Adds the Uniprot accession number
hs_fasta_df$PG.ProteinGroups <- unlist(lapply(names(hs_fasta), function(x){
  unlist(strsplit(x, split = "|", fixed = TRUE))[2]
}))

# Adds the sequence of the protein as a new column
hs_fasta_df$ProteinSequence <- unlist(lapply(hs_fasta, function(x){
  paste(x, collapse = "")
}))

# Merging the fasta file with the stoich file
data_1k <- merge(hs_fasta_df, data_1k)


# Indexing peptide start position on protein
peptide_index <- str_locate(data_1k$ProteinSequence, 
                            data_1k$EG.StrippedSequence)[,1]

# Indexing fragment ion sequence position on peptide
fragment_index <- str_locate(data_1k$EG.StrippedSequence, 
                             data_1k$fragment_ion_sequence)[,1]

# Indexing K site on fragment ion sequence
k_index <- str_locate(data_1k$fragment_ion_sequence, "K")[,1]

# Annotating K site on the protein
data_1k$k_site <- paste("K", 
                        peptide_index + fragment_index + k_index - 2, 
                        sep = "")

# Removing files
rm(hs_fasta_df, hs_fasta, fragment_index, k_index, peptide_index);gc()

# Reorganizing columns
data_1k <- data_1k %>% 
  select(sample_id, PG.ProteinGroups, PG.ProteinDescriptions, PG.Genes, 
         EG.StrippedSequence, EG.ModifiedSequence, k_count, k_site, 
         everything()) %>%
  select(-Fasta.header, -ProteinSequence) %>%
  arrange(EG.StrippedSequence, sample_id)

# Spreading isotope labels ------------------------------------------------


frg_1k <- data_1k %>% 
  select(-F.FrgMz, -FG.PrecMz, -EG.DatapointsPerPeak, 
         -FG.ShapeQualityScore) %>% 
  gather(name, value, 18:19) %>% 
  unite(temp, name, isotope_label, sep = "_") %>% 
  spread(temp, value)


# Correcting for Isotopic Distribution ------------------------------------


# fragment ion sequences used for correction
frg_seq <- frg_1k %>% 
  select(PG.ProteinGroups, k_site, fragment_ion_sequence) %>% 
  distinct()

# Quantifying the relative abundance of the isotopic envelope
frg_seq$M0 <- NA
frg_seq$M0 <- unlist(lapply(frg_seq$fragment_ion_sequence, function(x){
  unlist(Isotopic_dist_BRAIN(x, nrPeaks = 50))[1]
}))

frg_seq$M3 <- NA
frg_seq$M3 <- unlist(lapply(frg_seq$fragment_ion_sequence, function(x){
  unlist(Isotopic_dist_BRAIN(x, nrPeaks = 50))[4]
}))

# Merging correction factors with data
frg_1k <- merge(frg_1k, frg_seq)
rm(frg_seq);gc()

# Corrected heavy peak area
frg_1k$H_corr <- frg_1k$F.PeakArea_H - 
  (frg_1k$F.PeakArea_L * (1 / (frg_1k$M0 / frg_1k$M3)))


# Calculating Stoichiometry -----------------------------------------------


# Fragment level stoichiometry
frg_stoich <- frg_1k %>% 
  filter(EG.Qvalue_H <= 0.01 | EG.Qvalue_L <= 0.01) %>% 
  mutate(stoich_obs = F.PeakArea_L / (F.PeakArea_L + F.PeakArea_H),
         stoich_corr = F.PeakArea_L / (F.PeakArea_L + H_corr)) %>% 
  filter(!is.na(stoich_obs)) %>% 
  select(-M0, -M3)


# Fragment level stoichiometry --------------------------------------------


# Summing fragment ion species across fractions
frg_sum <- frg_stoich %>%
  group_by(sample_id, PG.ProteinGroups, PG.ProteinDescriptions, 
           PG.Genes, EG.StrippedSequence, EG.ModifiedSequence, k_count, 
           k_site, fragment_ion_sequence, F.FrgIon, F.FrgType, 
           F.FrgNum, F.FrgZ, F.FrgLossType) %>%
  summarize(L_sum = sum(F.PeakArea_L, na.rm = TRUE),
            H_sum = sum(F.PeakArea_H, na.rm = TRUE),
            H_corr_sum = sum(H_corr, na.rm = TRUE),
            EG.Qvalue_L_rms = sqrt(mean((EG.Qvalue_L)^2, na.rm = TRUE)),
            EG.Qvalue_H_rms = sqrt(mean((EG.Qvalue_H)^2, na.rm = TRUE))) %>%
  ungroup() %>%
  mutate(stoich_obs = L_sum / (L_sum + H_sum),
         stoich_corr = L_sum / (L_sum + H_corr_sum)) %>%
  arrange(PG.ProteinGroups, sample_id)


# Peptide level stoichiometry ---------------------------------------------


# Summing the fragment ion areas across all the fractions
pep_sum <- frg_stoich %>% 
  group_by(sample_id, PG.ProteinGroups, PG.ProteinDescriptions, PG.Genes, 
           EG.StrippedSequence, EG.ModifiedSequence, k_count, k_site) %>% 
  summarize(L_sum = sum(F.PeakArea_L, na.rm = TRUE),
            H_sum = sum(F.PeakArea_H, na.rm = TRUE),
            H_corr_sum = sum(H_corr, na.rm = TRUE),
            EG.Qvalue_L_rms = sqrt(mean((EG.Qvalue_L)^2, na.rm = TRUE)),
            EG.Qvalue_H_rms = sqrt(mean((EG.Qvalue_H)^2, na.rm = TRUE))) %>% 
  mutate(stoich_obs = L_sum / (L_sum + H_sum),
         stoich_corr = L_sum / (L_sum + H_corr_sum)) %>% 
  ungroup() %>% 
  arrange(PG.ProteinGroups, k_site, sample_id)


# Writing tables ----------------------------------------------------------


# write.table(frg_sum, file = "./MCF7_stoich/Cleaned/20181031_tidy_frg_stoich_output_correctFASTA.csv", sep = ",", row.names = FALSE)
# write.table(pep_sum, file = "./MCF7_stoich/Cleaned/20181031_tidy_peptide_stoich_output_correctFASTA.csv", sep = ",", row.names = FALSE)


# Quality check of the data -----------------------------------------------


## How does the raw data look?

# Peak Area
ggplot(data = data) +
  scale_color_brewer(palette = "Dark2") +
  theme_bw(base_size = 15) +
  scale_fill_brewer(palette = "Dark2") +
  geom_density(aes(x = log10(F.PeakArea), color = sample_id)) +
  facet_grid(~isotope_label)

# Q-Value
ggplot(data = data) +
  scale_color_brewer(palette = "Dark2") +
  theme_bw(base_size = 15) +
  scale_fill_brewer(palette = "Dark2") +
  geom_density(aes(x = EG.Qvalue, color = sample_id)) +
  xlim(c(0,0.1)) +
  facet_grid(~isotope_label)

# Data Points per Peaks
ggplot(data = data) +
  scale_color_brewer(palette = "Dark2") +
  theme_bw(base_size = 15) +
  scale_fill_brewer(palette = "Dark2") +
  geom_histogram(aes(x = EG.DatapointsPerPeak, fill = sample_id), position = "dodge") +
  xlim(c(0,10)) +
  facet_grid(~isotope_label)

# Shape Quality Score
ggplot(data = data) +
  scale_color_brewer(palette = "Dark2") +
  theme_bw(base_size = 15) +
  scale_fill_brewer(palette = "Dark2") +
  geom_density(aes(x = FG.ShapeQualityScore, color = sample_id)) +
  facet_grid(~isotope_label)

# Interference Score
ggplot(data = data) +
  scale_color_brewer(palette = "Dark2") +
  theme_bw(base_size = 15) +
  scale_fill_brewer(palette = "Dark2") +
  geom_density(aes(x = F.InterferenceScore, color = sample_id)) +
  xlim(c(0,1)) +
  facet_grid(~isotope_label)

# Predicted Relative Intensity
ggplot(data = data) +
  scale_color_brewer(palette = "Dark2") +
  theme_bw(base_size = 15) +
  scale_fill_brewer(palette = "Dark2") +
  geom_density(aes(x = F.PredictedRelativeIntensity, color = sample_id)) +
  facet_grid(~isotope_label)

# Missed Cleavage
ggplot(data = data) +
  scale_color_brewer(palette = "Dark2") +
  theme_bw(base_size = 15) +
  scale_fill_brewer(palette = "Dark2") +
  geom_histogram(aes(x = PEP.NrOfMissedCleavages, fill = sample_id), stat = "count", position = "dodge")


# Quality check of fragment ions ------------------------------------------


# How does the raw data look?

# Peak Area
ggplot(data = data_1k) +
  scale_color_brewer(palette = "Dark2") +
  theme_bw(base_size = 15) +
  scale_fill_brewer(palette = "Dark2") +
  geom_density(aes(x = log10(F.PeakArea), color = sample_id)) +
  facet_grid(~isotope_label)

# Q-Value
ggplot(data = data_1k) +
  scale_color_brewer(palette = "Dark2") +
  theme_bw(base_size = 15) +
  scale_fill_brewer(palette = "Dark2") +
  geom_density(aes(x = EG.Qvalue, color = sample_id)) + xlim(c(0,0.1)) +
  facet_wrap(~isotope_label)

# Data Points per Peaks
ggplot(data = data_1k) +
  scale_color_brewer(palette = "Dark2") +
  theme_bw(base_size = 15) +
  scale_fill_brewer(palette = "Dark2") +
  geom_histogram(aes(x = EG.DatapointsPerPeak, fill = sample_id), position = "dodge") +
  xlim(c(0,10)) +
  facet_grid(~isotope_label)

# Shape Quality Score
ggplot(data = data_1k) +
  scale_color_brewer(palette = "Dark2") +
  theme_bw(base_size = 15) +
  scale_fill_brewer(palette = "Dark2") +
  geom_density(aes(x = FG.ShapeQualityScore, color = sample_id)) +
  facet_grid(~isotope_label)

# Missed Cleavage
ggplot(data = data_1k) +
  scale_color_brewer(palette = "Dark2") +
  theme_bw(base_size = 15) +
  scale_fill_brewer(palette = "Dark2") +
  geom_histogram(aes(x = PEP.NrOfMissedCleavages, fill = sample_id), stat = "count", position = "dodge")


# Stoichiometry figures ---------------------------------------------------


# Rank order plot of fragment stoichiometry
ggplot(subset(frg_stoich, !is.na(frg_stoich$stoich_corr))) +
  geom_point(aes(x = rank(stoich_corr), y = stoich_corr)) +
  facet_wrap(~sample_id)

# Density plot of fragment stoichiometry
ggplot(subset(frg_sum, stoich_corr < 1 & !is.na(frg_stoich$stoich_corr))) +
  geom_density(aes(x = stoich_corr, color = sample_id))

