
# Loading Libraries -------------------------------------------------------


options(stringsAsFactors = FALSE)
library(tidyverse)
library(stringr)
library(ggthemes)
library(BRAIN)
library(seqinr)
library(GGally)
source("./Functions/function_1_Isotopic_distribution_BRAIN.R")
library(rio)


# Data import -------------------------------------------------------------


data <- read.csv("./Serum_Stimulated_HCT116/Raw/20190219_SSHCT116_AcK_Pulsar_Report.csv")

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

# Time point 
data$time_point1 <- as.character(unlist(lapply(data$R.FileName, function(x){
  unlist(strsplit(x, split = "_"))[6]
})))
data$time_point <- NA
data$time_point[which(data$time_point1 == "0")] <- 0
data$time_point[which(data$time_point1 == "1")] <- 1
data$time_point[which(data$time_point1 == "15m")] <- 0.25
data$time_point[which(data$time_point1 == "2")] <- 2
data$time_point[which(data$time_point1 == "4")] <- 4
data$time_point1 <- NULL

# Sample ID
data$sample_id <- as.numeric(unlist(lapply(data$R.FileName, function(x){
  unlist(strsplit(x, split = "_"))[7]
})))

# K count
data$k_count <- str_count(data$PEP.StrippedSequence, "K")

# Acetyl count
data$acetyl_count <- str_count(data$EG.ModifiedPeptide, "K\\[")

# Fragment Ion Sequence
data$fragment_ion_sequence <- 
  substring(data$PEP.StrippedSequence, 
            ifelse(data$F.FrgType == "b", 1, (nchar(data$PEP.StrippedSequence) - (data$F.FrgNum - 1))), 
            ifelse(data$F.FrgType == "y", nchar(data$PEP.StrippedSequence), data$F.FrgNum))

# Number of lysines on the fragment ion sequence
data$frag_k_count <- str_count(data$fragment_ion_sequence, "K")

# adding an methionine oxidation column
data$oxidation_count <- str_count(data$EG.ModifiedPeptide, "M\\[")



# Adding the Isotope Label ------------------------------------------------


# Need Modified peptide column, not IntModifiedPeptide
spec_lib_isotope <- spec_lib_isotope %>% 
  rename("EG.ModifiedPeptide" = "IntModifiedPeptide",
         "ModifiedPeptide" = "EG.ModifiedPeptide")

spec_lib_isotope2 <- spec_lib_isotope %>% unite(col = "FG.Id",
                                                LabeledPeptide, FG.Charge, 
                                                sep = ".", remove = FALSE)

spec_lib_merge <- spec_lib_isotope2 %>% 
  select(FG.Charge, EG.StrippedSequence, FG.Id, k_count, 
         EG.ModifiedPeptide, fragment_ion_sequence, 
         F.FrgLossType, F.FrgZ, isotope_label) %>% 
  rename("F.FrgZ" = "F.Charge",
         "EG.StrippedSequence" = "PEP.StrippedSequence")

# Removing duplicate rows
spec_lib_merge <- spec_lib_merge %>% 
  distinct()

# Merging the data and merge dataframes
data <- merge(data, spec_lib_merge)
data <- data %>% 
  arrange(PEP.StrippedSequence, fragment_ion_sequence, sample_id)

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
#data_1k2 <- data_1k

# Subsetting the columns
data_1k <- data_1k %>% 
  select(sample_id, time_point, PG.Genes, PG.ProteinDescriptions, PG.ProteinGroups, 
         EG.DatapointsPerPeak, PEP.NrOfMissedCleavages, FG.ShapeQualityScore, 
         PEP.StrippedSequence, EG.ModifiedPeptide, 
         FG.PrecMz, k_count, fragment_ion_sequence, 
         F.FrgIon, F.FrgType, F.FrgNum, F.FrgMz, 
         F.FrgLossType, F.Charge, isotope_label, EG.Qvalue, F.PeakArea) %>%
  arrange(sample_id, PEP.StrippedSequence, -FG.ShapeQualityScore)


# Selecting the top id for protein groups - placed this function here to save time
data_1k[3:5] <- apply(data_1k[3:5], 2, extract_first_item)

# Removing duplicate rows
data_1k <- data_1k %>% 
  distinct(sample_id, time_point, PG.Genes, PG.ProteinDescriptions, PG.ProteinGroups, 
           PEP.NrOfMissedCleavages, k_count, PEP.StrippedSequence, EG.ModifiedPeptide,  
           fragment_ion_sequence, F.FrgIon, F.FrgType, F.FrgNum, 
           F.FrgLossType, F.FrgIon, F.Charge, isotope_label, .keep_all = TRUE)


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
                            data_1k$PEP.StrippedSequence)[,1]

# Indexing fragment ion sequence position on peptide
fragment_index <- str_locate(data_1k$PEP.StrippedSequence, 
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
  select(sample_id, time_point, PG.ProteinGroups, PG.ProteinDescriptions, PG.Genes, 
         PEP.StrippedSequence, EG.ModifiedPeptide, k_count, k_site, 
         everything()) %>%
  select(-Fasta.header, -ProteinSequence) %>%
  arrange(PEP.StrippedSequence, sample_id)


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


# Peptide level stoichiometry ---------------------------------------------


# Summing the fragment ion areas across all the fractions
pep_sum <- frg_stoich %>% 
  group_by(sample_id, time_point, PG.ProteinGroups, PG.ProteinDescriptions, PG.Genes, 
           PEP.StrippedSequence, EG.ModifiedPeptide, k_count, k_site) %>% 
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

# write.table(frg_stoich, file = "./Serum_Stimulated_HCT116/Cleaned/20190220_SSHCT116_tc_tidy_frag_stoich_output.csv", sep = ",", row.names = FALSE)
# write.table(pep_sum, file = "./Serum_Stimulated_HCT116/Cleaned/20190220_SSHCT116_tc_tidy_peptide_stoich_output.csv", sep = ",", row.names = FALSE)


