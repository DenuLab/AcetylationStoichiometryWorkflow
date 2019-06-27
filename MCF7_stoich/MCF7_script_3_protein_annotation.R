
# Loading Libraries -------------------------------------------------------

options(stringsAsFactors = FALSE)
library(tidyverse)
library(seqinr)
library(GGally)
library(ggpubr)
library(PMCMR) # for stats


# Importing data ----------------------------------------------------------

# Fragment stoichiometry
frg_stoich <- read.csv("./MCF7_stoich/Cleaned/20181031_tidy_frg_stoich_output_correctFASTA.csv")

# Peptide stoichiometry
pep_stoich <- read.csv("./MCF7_stoich/Cleaned/20181031_tidy_peptide_stoich_output_correctFASTA.csv")

# Human Fasta file
hs_fasta <- read.fasta(file = "./Fasta_files/human_swissprot_canonical_iRT_20171212.fasta", seqtype = "AA")

# Mitocarta fasta file
hs_mitocarta <- read.fasta(file = "./Fasta_files/Human_MitoCarta2.fasta", seqtype = "AA")

# Uniprot subcellular location
uniprot_id <- read.csv("./Fasta_files/uniprot_subcellular_location_library_only.csv")
names(uniprot_id)[1] <- "PG.ProteinGroups"
uniprot_id$id <- NULL


# Formatting Fasta files --------------------------------------------------


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

# Removes column
hs_fasta_df$Fasta.header <- NULL

# Generates a dataframe with the column header names
hs_mitocarta_df <- data.frame("Fasta.header" = names(hs_mitocarta))

# Adds the sequence of the protein as a new column
hs_mitocarta_df$ProteinSequence <- unlist(lapply(hs_mitocarta, function(x){
  paste(x, collapse = "")
}))

# Removes column
hs_mitocarta_df$Fasta.header <- NULL

# merging based on protein sequence and adding mitocarta label
hs_mitocarta_df$mitocarta <- TRUE
hs_fasta_df <- merge(hs_fasta_df, hs_mitocarta_df, all.x = TRUE)
hs_fasta_df$mitocarta[which(is.na(hs_fasta_df$mitocarta))] <- FALSE

# Removing column
hs_fasta_df$ProteinSequence <- NULL

# Removing duplicates
hs_fasta_df <- distinct(hs_fasta_df)


# Annotating subcellular location -----------------------------------------


# Adds the mitocarta label to the df
frg_stoich <- merge(hs_fasta_df, frg_stoich, all.y = TRUE)

# Adds the Uniprot subcellular localization to the df
frg_stoich <- merge(uniprot_id, frg_stoich, all.y = TRUE)

# Removing files
rm(hs_fasta, hs_fasta_df, hs_mitocarta, hs_mitocarta_df, uniprot_id);gc()

# Subcellular Localization
frg_stoich$subcellular_location1[which(frg_stoich$subcellular_location1 == "")] <- NA
frg_stoich$subcellular_location2[which(frg_stoich$subcellular_location2 == "")] <- NA


# where I left off --------------------------------------------------------


frg_stoich$localization <- NA
frg_stoich$localization[grep(TRUE, frg_stoich$mitocarta)] <- "Cytoplasm" #"Mitochondria"
frg_stoich$localization[grep("Mito", frg_stoich$subcellular_location1)] <- "Cytoplasm" #"Mitochondria"
frg_stoich$localization[grep("Lyso", frg_stoich$subcellular_location1)] <- "Cytoplasm"
frg_stoich$localization[grep("mitochondrial", frg_stoich$PG.ProteinDescriptions)] <- "Cytoplasm" #"Mitochondria"

frg_stoich$localization[grep("Nucl", frg_stoich$subcellular_location1)] <- "Nucleus"
frg_stoich$localization[grep("Endoplas", frg_stoich$subcellular_location1)] <- "Nucleus"
frg_stoich$localization[grep("Golgi", frg_stoich$subcellular_location1)] <- "Nucleus"
frg_stoich$localization[grep("Cell projec", frg_stoich$subcellular_location1)] <- "Nucleus"
frg_stoich$localization[grep("cell membrane", frg_stoich$subcellular_location1)] <- "Nucleus"
frg_stoich$localization[grep("Cell membrane", frg_stoich$subcellular_location1)] <- "Nucleus"
frg_stoich$localization[grep("Secreted", frg_stoich$subcellular_location1)] <- "Nucleus"
frg_stoich$localization[grep("junction", frg_stoich$subcellular_location1)] <- "Nucleus"
frg_stoich$localization[grep("Membrane", frg_stoich$subcellular_location1)] <- "Nucleus"
frg_stoich$localization[grep("Cell surface", frg_stoich$subcellular_location1)] <- "Nucleus"

frg_stoich$localization[grep("Cytopla", frg_stoich$subcellular_location1)] <- "Cytoplasm"
frg_stoich$localization[grep("Early endo", frg_stoich$subcellular_location1)] <- "Cytoplasm"
frg_stoich$localization[grep("Endosome", frg_stoich$subcellular_location1)] <- "Cytoplasm"
frg_stoich$localization[grep("Prevacuolar", frg_stoich$subcellular_location1)] <- "Cytoplasm"
frg_stoich$localization[grep("ribosomal protein", frg_stoich$protein_name)] <- "Cytoplasm"
frg_stoich$localization[grep("Eukaryotic initiation factor", frg_stoich$protein_name)] <- "Cytoplasm"
frg_stoich$localization[grep("Eukaryotic translation", frg_stoich$protein_name)] <- "Cytoplasm"
frg_stoich$localization[grep("Glutathione synth", frg_stoich$protein_name)] <- "Cytoplasm"
frg_stoich$localization[grep("Phosphoglycerate mutase", frg_stoich$protein_name)] <- "Cytoplasm"
frg_stoich$localization[grep("Glucose-6", frg_stoich$protein_name)] <- "Cytoplasm"
frg_stoich$localization[grep("Alcohol dehydrogenase", frg_stoich$protein_name)] <- "Cytoplasm"

frg_stoich$localization[grep("Histone", frg_stoich$protein_name)] <- "Nucleus"
frg_stoich$localization[grep("chromatin", frg_stoich$protein_name)] <- "Nucleus"
frg_stoich$localization[grep("Chromatin", frg_stoich$protein_name)] <- "Nucleus"
frg_stoich$localization[grep("chromosom", frg_stoich$protein_name)] <- "Nucleus"

frg_stoich$localization[grep("Keratin", frg_stoich$PG.ProteinDescriptions)] <- NA

######################################################################
# Localization 2
frg_stoich$localization2 <- NA

frg_stoich$localization2[grep(TRUE, frg_stoich$mitocarta)] <- "Mitochondria" 
frg_stoich$localization2[grep("Mito", frg_stoich$subcellular_location1)] <- "Mitochondria" 
frg_stoich$localization2[grep(", mitochondrial", frg_stoich$PG.ProteinDescriptions)] <- "Mitochondria"

frg_stoich$localization2[grep("Nucl", frg_stoich$subcellular_location1)] <- "Nucleus"
frg_stoich$localization2[grep("Cell projec", frg_stoich$subcellular_location1)] <- "Nucleus"
frg_stoich$localization2[grep("Secreted", frg_stoich$subcellular_location1)] <- "Nucleus"
frg_stoich$localization2[grep("junction", frg_stoich$subcellular_location1)] <- "Nucleus"
frg_stoich$localization2[grep("Cell surface", frg_stoich$subcellular_location1)] <- "Nucleus"

frg_stoich$localization2[grep("Endoplas", frg_stoich$subcellular_location1)] <- "ER"
frg_stoich$localization2[grep("Golgi", frg_stoich$subcellular_location1)] <- "Golgi"
frg_stoich$localization2[grep("cell membrane", frg_stoich$subcellular_location1)] <- "Membrane"
frg_stoich$localization2[grep("Cell membrane", frg_stoich$subcellular_location1)] <- "Membrane"
frg_stoich$localization2[grep("Membrane", frg_stoich$subcellular_location1)] <- "Membrane"

frg_stoich$localization2[grep("Cytopla", frg_stoich$subcellular_location1)] <- "Cytoplasm"
frg_stoich$localization2[grep("Early endo", frg_stoich$subcellular_location1)] <- "Cytoplasm"
frg_stoich$localization2[grep("Endosome", frg_stoich$subcellular_location1)] <- "Cytoplasm"
frg_stoich$localization2[grep("Prevacuolar", frg_stoich$subcellular_location1)] <- "Cytoplasm"
frg_stoich$localization2[grep("ribosomal protein", frg_stoich$protein_name)] <- "Cytoplasm"
frg_stoich$localization2[grep("Eukaryotic initiation factor", frg_stoich$protein_name)] <- "Cytoplasm"
frg_stoich$localization2[grep("Eukaryotic translation", frg_stoich$protein_name)] <- "Cytoplasm"
frg_stoich$localization2[grep("Glutathione synth", frg_stoich$protein_name)] <- "Cytoplasm"
frg_stoich$localization2[grep("Phosphoglycerate mutase", frg_stoich$protein_name)] <- "Cytoplasm"
frg_stoich$localization2[grep("Glucose-6", frg_stoich$protein_name)] <- "Cytoplasm"
frg_stoich$localization2[grep("Alcohol dehydrogenase", frg_stoich$protein_name)] <- "Cytoplasm"
frg_stoich$localization2[grep("Lyso", frg_stoich$subcellular_location1)] <- "Lysosome"

frg_stoich$localization2[grep("Histone H", frg_stoich$protein_name)] <- "Histone"
frg_stoich$localization2[grep("chromatin", frg_stoich$protein_name)] <- "Chromatin"
frg_stoich$localization2[grep("Chromatin", frg_stoich$protein_name)] <- "Chromatin"
frg_stoich$localization2[grep("chromosom", frg_stoich$protein_name)] <- "Chromatin"
frg_stoich$localization2[grep(", mitochondrial", frg_stoich$PG.ProteinDescriptions)] <- "Mitochondria"

table(frg_stoich$localization2)

frg_stoich$localization2 <- factor(frg_stoich$localization2, levels = 
                                     c("Cytoplasm", "Nucleus", "Chromatin", "Histone", "Lysosome", "Membrane", 
                                       "Golgi", "ER", "Mitochondria"))


# Fragment ion stoichiometry ----------------------------------------------


# Fragment stoichiometry data
frg_wide <- frg_stoich %>%
  select(sample_id, localization, localization2,
         PG.ProteinGroups, EG.StrippedSequence, EG.ModifiedSequence, 
         fragment_ion_sequence,
         k_site, F.FrgZ, F.FrgLossType, stoich_corr) %>%
  subset(stoich_corr <= 1) %>% 
  spread(sample_id, stoich_corr) %>% 
  mutate(quant = "fragment")

frg_wide2 <- frg_wide %>% 
  group_by(localization, localization2,
           PG.ProteinGroups, EG.StrippedSequence, EG.ModifiedSequence, k_site) %>% 
  summarise(mean_MCF71 = mean(MCF7_1, na.rm = TRUE),
            mean_MCF72 = mean(MCF7_2, na.rm = TRUE), 
            mean_MCF73 = mean(MCF7_3, na.rm = TRUE))

# Mean
frg_wide$stoich_mean <- rowMeans(frg_wide[10:12], na.rm = TRUE)

# Standard deviation
frg_wide$stoich_sd <- apply(frg_wide[10:12], 1, function(x){
  sd(x, na.rm = TRUE)
})

# n obs
frg_wide$stoich_n <- apply(frg_wide[10:12], 1, function(x){
  length(which(!is.na(x)))
})

# Standard error
frg_wide$stoich_se <- frg_wide$stoich_sd / sqrt(frg_wide$stoich_n)


# fragment and peptide scatterplot ----------------------------------------

pep_wide <- pep_stoich %>% 
  select(sample_id, PG.ProteinGroups, EG.StrippedSequence, 
         EG.ModifiedSequence, k_site, stoich_corr) %>%
  subset(stoich_corr <= 1) %>% 
  spread(sample_id, stoich_corr) %>% 
  mutate(quant = "peptide")

# Mean
pep_wide$stoich_mean <- rowMeans(pep_wide[5:7], na.rm = TRUE)

# Standard deviation
pep_wide$stoich_sd <- apply(pep_wide[5:7], 1, function(x){
  sd(x, na.rm = TRUE)
})

# n obs
pep_wide$stoich_n <- apply(pep_wide[5:7], 1, function(x){
  length(which(!is.na(x)))
})

# Standard error
pep_wide$stoich_se <- pep_wide$stoich_sd / sqrt(pep_wide$stoich_n)

frg_wide_temp <- frg_wide %>% 
  select(PG.ProteinGroups, EG.StrippedSequence, EG.ModifiedSequence, k_site,
         MCF7_1:stoich_se)

# combining fragment and peptide stoichiometry
frg_pep_wide <- rbind(frg_wide_temp, pep_wide)
rm(frg_wide_temp);gc()

frg_pep_wide <- frg_pep_wide %>% 
  filter(stoich_n == 3)


# Plots -------------------------------------------------------------------

frg_wide$singleaa_loc <- NA
frg_wide$singleaa_loc[which(frg_wide$localization == "Cytoplasm")] <- "Cytoplasm"
frg_wide$singleaa_loc[which(frg_wide$localization == "Nucleus")] <- "Nucleus"
frg_wide$singleaa_loc[which(frg_wide$localization2 == "Mitochondria")] <- "Mitochondria"

# Unique fragment ions
ggplot(frg_wide %>% filter(!is.na(localization))) +
  geom_bar(aes(x = localization, fill = localization)) +
  scale_fill_brewer(palette = "Set1") +
  ggthemes::theme_pander(base_size = 18) +
  guides(fill = FALSE) +
  labs(x = NULL)
#ggsave(filename = "figures/20181107_fragment_barplot.pdf", height = 5, width = 7)

# Scatterplot matrix fragment and peptide level stoichiometry
### Figure 4A
ggscatmat(subset(frg_pep_wide), 
          columns = 5:7, alpha = 0.1, color = "quant", corMethod = "spearman") +
  scale_color_brewer(palette = "Set1", direction = 1) +
  theme_light(base_size = 14) +
  labs(y = "Stoichiometry",
       x = "Stoichiometry")
#ggsave(filename = "figures/20190625_fragment_peptide_scatmat_spearman.pdf", height = 4, width = 5.6, useDingbats = FALSE)

# Density plot
ggplot(frg_wide %>% filter(stoich_sd < 0.075, stoich_n == 3, !is.na(localization))) + 
  geom_density(aes(x = stoich_mean, color = localization)) +
  theme_light(base_size = 18) +
  scale_color_brewer(palette = "Set1")

# Dotplot
### Figure 4B
ggplot(frg_wide %>% filter(stoich_sd < 0.075, !is.na(localization), stoich_n == 3),
       aes(x = localization, y = stoich_mean)) +
  geom_dotplot(aes(fill = localization2, color = localization2),
               stackdir = "center", 
               binaxis = "y", 
               method = "dotdensity", 
               binwidth = 0.01) +
  geom_boxplot(aes(), alpha = 0) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  annotate("segment", x = 1, xend = 2, y = 1.01, yend = 1.01) +
  annotate("text", x = 1.5, y = 1.05, label = "p.value = 0.00027", size = 4) +
  labs(fill = "Subcellular\nlocation",
       color = "Subcellular\nlocation",
       x = NULL,
       y = "Corrected Stoichiometry") +
  ggthemes::theme_pander(base_size = 16)
#ggsave(filename = "figures/20181107_Cyto_nuclear_dotplot.pdf", height = 5, width = 10, useDingbats = FALSE)

# sd less than 0.053 (median of standard deviation)
posthoc.kruskal.nemenyi.test(stoich_mean ~ as.factor(localization),
                             data = subset(frg_wide, stoich_sd < 0.075 & 
                                             !is.na(localization) &
                                             stoich_n == 3),
                             dist = "Tukey")
#pvalue = 0.00027


fivenum(frg_wide$stoich_sd[frg_wide$stoich_n == 3])[4]
# 0.075
