
# Loading packages --------------------------------------------------------


options(stringsAsFactors = FALSE)
library(tidyverse)
library(dplyr)
library(seqinr)
library(stringr)
source("./Functions/function_2_QSSA.R")
library(viridisLite)
library(gplots)



# Data Import -------------------------------------------------------------


# load data sets
pep_stoich <- read.csv("./MCF7_stoich/Cleaned/20181031_tidy_peptide_stoich_output_correctFASTA.csv")
kegg <- load_kegg("./Fasta_files/Human_KEGG_October_01_2017_UniProt.gmt")
sequences <- load_fasta("./Fasta_files/human_swissprot_canonical_iRT_20171212.fasta")
background <- read.csv("./Spectral_Libraries/20181029_MCF7_Library_Remake_correctFASTA_inflated.csv")["UniProtIds"]


# Cleaning peptide data ---------------------------------------------------


# peptide summary
pep_wide <- pep_stoich %>% 
  select(sample_id, PG.ProteinGroups, EG.StrippedSequence, 
         EG.ModifiedSequence, k_site, stoich_corr) %>% 
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

stoich <- pep_wide %>% 
  filter(stoich_n == 3)


# Formatting data ---------------------------------------------------------


# Cleaning the background
background <- unique(background)
background$UniProtIds <- unlist(lapply(background$UniProtIds, function(x){
  unlist(strsplit(x, split = ";"))[1]
}))

# Adding stoichiometry bins
stoich$bin <- cut(rank(stoich$stoich_mean), breaks = 4, labels = FALSE)

# trim kegg map to detected uniprots
kegg <- trim_kegg(kegg, unique(stoich$PG.ProteinGroups))

# Truncates fasta to just those found in our data
sequences <- sequences[names(sequences) %in% background$UniProtIds]


# Formatting QSSA ---------------------------------------------------------


# Subsets the binned stoichiometry data
stoich_1 <- stoich %>% filter(bin == 1)
stoich_2 <- stoich %>% filter(bin == 2)
stoich_3 <- stoich %>% filter(bin == 3)
stoich_4 <- stoich %>% filter(bin == 4)


# 1. Calculate the coverage ratio for the pathway
#    i.e. What percent of the lysine sites in the pathway
#    have we detected as acetylated?
pathway_summary <- data.frame("kegg_pathway" = names(kegg))
pathway_summary$n_K <- lysine_sites_per_pathway(kegg, sequences)
pathway_summary$n_acK_1 <- detected_sites_per_pathway(kegg, stoich_1)
pathway_summary$n_acK_2 <- detected_sites_per_pathway(kegg, stoich_2)
pathway_summary$n_acK_3 <- detected_sites_per_pathway(kegg, stoich_3)
pathway_summary$n_acK_4 <- detected_sites_per_pathway(kegg, stoich_4)
pathway_summary <- pathway_summary %>%
  mutate(coverage_ratio_1 = n_acK_1 / n_K,
         coverage_ratio_2 = n_acK_2 / n_K,
         coverage_ratio_3 = n_acK_3 / n_K,
         coverage_ratio_4 = n_acK_4 / n_K)

pathway_summary$zcoverage_1 <- calc_z(pathway_summary$coverage_ratio_1)
pathway_summary$zcoverage_2 <- calc_z(pathway_summary$coverage_ratio_2)
pathway_summary$zcoverage_3 <- calc_z(pathway_summary$coverage_ratio_3)
pathway_summary$zcoverage_4 <- calc_z(pathway_summary$coverage_ratio_4)


# 2. Calculate the sum of the stoichiometries 
#    across all conditions and peptides in each pathway.
pathway_summary$sum_stoich_1 <- sum_per_pathway(kegg, stoich_1)
pathway_summary$sum_stoich_2 <- sum_per_pathway(kegg, stoich_2)
pathway_summary$sum_stoich_3 <- sum_per_pathway(kegg, stoich_3)
pathway_summary$sum_stoich_4 <- sum_per_pathway(kegg, stoich_4)

pathway_summary$zsum_stoich_1 <- calc_z(pathway_summary$sum_stoich_1)
pathway_summary$zsum_stoich_2 <- calc_z(pathway_summary$sum_stoich_2)
pathway_summary$zsum_stoich_3 <- calc_z(pathway_summary$sum_stoich_3)
pathway_summary$zsum_stoich_4 <- calc_z(pathway_summary$sum_stoich_4)


# 3. Take the sum of the two z-scores
pathway_summary$qssa_1 <- pathway_summary$zcoverage_1 + pathway_summary$zsum_stoich_1
pathway_summary$qssa_2 <- pathway_summary$zcoverage_2 + pathway_summary$zsum_stoich_2
pathway_summary$qssa_3 <- pathway_summary$zcoverage_3 + pathway_summary$zsum_stoich_3
pathway_summary$qssa_4 <- pathway_summary$zcoverage_4 + pathway_summary$zsum_stoich_4

rm(stoich_1, stoich_2, stoich_3, stoich_4);gc()

# Cleaning QSSA -----------------------------------------------------------


# pathway_summary <- pathway_summary %>% arrange(-qssa)
pathway_summary$kegg <- sapply(pathway_summary$kegg_pathway, function(x){
  unlist(strsplit(x, "%", fixed = T))[1]
})

# Rounding values
pathway_summary[7:26] <- mapply(pathway_summary[7:26], FUN = round, digits = 3)

# selecting qssa data
qssa_results <- pathway_summary %>% 
  select(kegg, qssa_1:qssa_4)

# Arranging by row mean
qssa_results$mean <- rowMeans(qssa_results[2:5])

# Long format of data for ggplot
qssa_results <- qssa_results %>% 
  arrange(-mean) %>% 
  mutate(rank = 1:nrow(.)) %>% 
  select(-mean) %>% 
  gather(condition, qssa, qssa_1:qssa_4)

# ## Matrix for heatmap
# qssa_matrix <- pathway_summary[28:33]
# qssa_matrix$mean <- rowMeans(qssa_matrix[1:5])
# qssa_matrix <- qssa_matrix %>% arrange(-mean)
# qssa_mat <- as.matrix(qssa_matrix[1:5])
# row.names(qssa_mat) <- qssa_matrix$kegg


# Plots -------------------------------------------------------------------


# basic heatmap
#pdf("figures/test.pdf", width = 8, height = 6)
ggplot(subset(qssa_results, rank <= 20)) +
  geom_tile(aes(x = condition, y = reorder(kegg, -rank), fill = qssa)) +
  scale_fill_distiller(palette = "RdBu", direction = -1) +
  theme_light(base_size = 12) +
  theme(axis.ticks = element_blank(),
        axis.title.y = element_blank()) +
  labs(title = "Quantitative Site Set functional score Analysis",
       x = "Stoichiometry bins")
#dev.off()

#pdf("figures/MCF7_QSSA_20190129_blues.pdf", width = 8, height = 6)
ggplot(subset(qssa_results, rank <= 20)) +
  geom_tile(aes(x = condition, y = reorder(kegg, -rank), fill = qssa)) +
  scale_fill_distiller(palette = "Blues", direction = -1, limits = c(-2.5, 0)) +
  # scale_fill_gradient2(high = "white", low = "blue",
  #                     midpoint = 0, na.value = "grey50", guide = "colourbar",
  #                     limits = c(-3, 0)) +
  theme_light(base_size = 12) +
  theme(axis.ticks = element_blank(),
        axis.title.y = element_blank()) +
  labs(title = "Quantitative Site Set functional score Analysis",
       x = "Stoichiometry bins")
#dev.off()

# pdf("figures/MCF7_QSSA_20190129_reds.pdf", width = 8, height = 6)
ggplot(subset(qssa_results, rank <= 20)) +
  geom_tile(aes(x = condition, y = reorder(kegg, -rank), fill = qssa)) +
  scale_fill_gradient2(high = "red3", low = "blue4", mid = "white",
                       midpoint = 0, na.value = "grey50", guide = "colourbar",
                       limits = c(0, NA)) +
  theme_light(base_size = 12) +
  theme(axis.ticks = element_blank(),
        axis.title.y = element_blank()) +
  labs(title = "Quantitative Site Set functional score Analysis",
       x = "Stoichiometry bins")
# dev.off()

fivenum( stoich$stoich_mean[stoich$bin == 1] )

# # Top n hits using Heatmap.2 
# heatmap.2(qssa_mat[1:20,],
#           # Dendogram control
#           Rowv = TRUE,
#           Colv = FALSE,
#           distfun = dist,
#           hclustfun = hclust,
#           dendrogram = "row",
#           # data scaling
#           scale = "none",
#           # Colors
#           col = "viridis",
#           # Level trace
#           trace = "none",
#           # Row/Column labeling
#           cexRow = 0.75,
#           cexCol = 1,
#           margins = c(5,20)
#           )

