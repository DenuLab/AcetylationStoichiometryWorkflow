
# Loading packages --------------------------------------------------------

options(stringsAsFactors = FALSE)
library(tidyverse)
library(dplyr)
library(seqinr)
library(stringr)
source("./Functions/function_2_QSSA.R")
library(viridisLite)
library(gplots)
library(data.table)
library(stringi)
library(RColorBrewer)


# Data Import -------------------------------------------------------------


# load data sets
cluster <- read.csv("./Serum_Stimulated_MCF7/Fuzzy_cm_data/20181102_Protein_stoich_FCM_clusters_reprocess_correctFASTA.csv")
cluster <- cluster %>% select(PG.ProteinGroups, PG.Genes, PG.ProteinDescriptions, everything())
extract_first_item <- function(x, a = ";", b = 1){
  unlist(lapply(x, function(y){
    unlist(strsplit(y, split = a))[b]
  }))
}
cluster[1:3] <- apply(cluster[1:3], 2, extract_first_item)

kegg <- load_kegg("./Fasta_files/Human_KEGG_October_01_2017_UniProt.gmt")
sequences <- load_fasta("./Fasta_files/human_swissprot_canonical_iRT_20171212.fasta")
background <- read.csv("./Spectral_Libraries/20181029_MCF7_Library_Remake_correctFASTA_inflated.csv")["UniProtIds"]

### Focusing on the stoichiometry clusters
cluster_stoich <- cluster %>% 
  select(PG.ProteinGroups, PG.Genes, PG.ProteinDescriptions, time,
         EG.StrippedSequence, EG.ModifiedSequence, k_site, stoich_mean,
         stoich_sd, stoich_n, stoich_cluster, stoich_isClusterMember,
         stoich_maxMembership) %>% 
  subset(!is.na(stoich_mean)) %>% 
  mutate(quant = "peptide")


ggplot(subset(cluster, stoich_cluster >= 1), aes(time, stoich_zscore, group = PG.ProteinGroups)) +
  geom_line(aes(group = k_site, color = stoich_maxMembership)) +
  facet_wrap(~ stoich_cluster)


# Formatting data ---------------------------------------------------------


# Cleaning the background
background <- unique(background)
background$UniProtIds <- unlist(lapply(background$UniProtIds, function(x){
  unlist(strsplit(x, split = ";"))[1]
}))

# trim kegg map to detected uniprots
kegg <- trim_kegg(kegg, unique(cluster_stoich$PG.ProteinGroups))

# Truncates fasta to just those found in our data
sequences <- sequences[names(sequences) %in% background$UniProtIds]

kegg_df <- melt(kegg)
kegg_df <- kegg_df %>%
  group_by(L1) %>%
  summarise(PG.ProteinGroups = toString(value)) %>%
  ungroup()

#### What I want is a list of the proteins that are in each kegg cluster to see
#### why certain clusters may have been enriched. Therefore, I am going to make
#### separate datafiles of each cluster to make the same kegg files...
cluster_1 <- cluster_stoich %>% subset(stoich_cluster == 1)
cluster_2 <- cluster_stoich %>% subset(stoich_cluster == 2)
cluster_3 <- cluster_stoich %>% subset(stoich_cluster == 3)
cluster_4 <- cluster_stoich %>% subset(stoich_cluster == 4)

kegg_1 <- trim_kegg(kegg, unique(cluster_1$PG.ProteinGroups))
kegg_2 <- trim_kegg(kegg, unique(cluster_2$PG.ProteinGroups))
kegg_3 <- trim_kegg(kegg, unique(cluster_3$PG.ProteinGroups))
kegg_4 <- trim_kegg(kegg, unique(cluster_4$PG.ProteinGroups))

kegg_1_df <- melt(kegg_1)
kegg_1_df <- kegg_1_df %>%
  group_by(L1) %>%
  summarise(PG.ProteinGroups = toString(value),
            n_protein = n()) %>%
  ungroup()

kegg_2_df <- melt(kegg_2)
kegg_2_df <- kegg_2_df %>%
  group_by(L1) %>%
  summarise(PG.ProteinGroups = toString(value),
            n_protein = n()) %>%
  ungroup()
kegg_3_df <- melt(kegg_3)
kegg_3_df <- kegg_3_df %>%
  group_by(L1) %>%
  summarise(PG.ProteinGroups = toString(value),
            n_protein = n()) %>%
  ungroup()
kegg_4_df <- melt(kegg_4)
kegg_4_df <- kegg_4_df %>%
  group_by(L1) %>%
  summarise(PG.ProteinGroups = toString(value),
            n_protein = n()) %>%
  ungroup()

rm(kegg_1, kegg_2, kegg_3, kegg_4, cluster_1, cluster_2, cluster_3, cluster_4);gc()


# Formatting QSSA ---------------------------------------------------------

# Subsets the binned stoichiometry data
cluster_stoich_1 <- cluster_stoich %>% filter(stoich_cluster == 1)
cluster_stoich_2 <- cluster_stoich %>% filter(stoich_cluster == 2)
cluster_stoich_3 <- cluster_stoich %>% filter(stoich_cluster == 3)
cluster_stoich_4 <- cluster_stoich %>% filter(stoich_cluster == 4)


# 1. Calculate the coverage ratio for the pathway
#    i.e. What percent of the lysine sites in the pathway
#    have we detected as acetylated?
pathway_summary <- data.frame("kegg_pathway" = names(kegg))
pathway_summary$n_K <- lysine_sites_per_pathway(kegg, sequences)
pathway_summary$n_acK_1 <- detected_sites_per_pathway(kegg, cluster_stoich_1)
pathway_summary$n_acK_2 <- detected_sites_per_pathway(kegg, cluster_stoich_2)
pathway_summary$n_acK_3 <- detected_sites_per_pathway(kegg, cluster_stoich_3)
pathway_summary$n_acK_4 <- detected_sites_per_pathway(kegg, cluster_stoich_4)
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
pathway_summary$sum_stoich_1 <- sum_per_pathway(kegg, cluster_stoich_1)
pathway_summary$sum_stoich_2 <- sum_per_pathway(kegg, cluster_stoich_2)
pathway_summary$sum_stoich_3 <- sum_per_pathway(kegg, cluster_stoich_3)
pathway_summary$sum_stoich_4 <- sum_per_pathway(kegg, cluster_stoich_4)

# pathway_summary$mean_stoich_1 <- mean_per_pathway(kegg, cluster_stoich_1)
# pathway_summary$mean_stoich_2 <- mean_per_pathway(kegg, cluster_stoich_2)
# pathway_summary$mean_stoich_3 <- mean_per_pathway(kegg, cluster_stoich_3)
# pathway_summary$mean_stoich_4 <- mean_per_pathway(kegg, cluster_stoich_4)
# 
# pathway_summary$mean_stoich_1[is.nan(pathway_summary$mean_stoich_1)] <- 0
# pathway_summary$mean_stoich_2[is.nan(pathway_summary$mean_stoich_2)] <- 0
# pathway_summary$mean_stoich_3[is.nan(pathway_summary$mean_stoich_3)] <- 0
# pathway_summary$mean_stoich_4[is.nan(pathway_summary$mean_stoich_4)] <- 0

pathway_summary$zsum_stoich_1 <- calc_z(pathway_summary$sum_stoich_1)
pathway_summary$zsum_stoich_2 <- calc_z(pathway_summary$sum_stoich_2)
pathway_summary$zsum_stoich_3 <- calc_z(pathway_summary$sum_stoich_3)
pathway_summary$zsum_stoich_4 <- calc_z(pathway_summary$sum_stoich_4)

# pathway_summary$zmean_stoich_1 <- calc_z(pathway_summary$mean_stoich_1)
# pathway_summary$zmean_stoich_2 <- calc_z(pathway_summary$mean_stoich_2)
# pathway_summary$zmean_stoich_3 <- calc_z(pathway_summary$mean_stoich_3)
# pathway_summary$zmean_stoich_4 <- calc_z(pathway_summary$mean_stoich_4)

# 3. Take the sum of the two z-scores
pathway_summary$qssa_1 <- pathway_summary$zcoverage_1 + pathway_summary$zsum_stoich_1
pathway_summary$qssa_2 <- pathway_summary$zcoverage_2 + pathway_summary$zsum_stoich_2
pathway_summary$qssa_3 <- pathway_summary$zcoverage_3 + pathway_summary$zsum_stoich_3
pathway_summary$qssa_4 <- pathway_summary$zcoverage_4 + pathway_summary$zsum_stoich_4

# pathway_summary$qssa_1 <- pathway_summary$zcoverage_1 + pathway_summary$zmean_stoich_1
# pathway_summary$qssa_2 <- pathway_summary$zcoverage_2 + pathway_summary$zmean_stoich_2
# pathway_summary$qssa_3 <- pathway_summary$zcoverage_3 + pathway_summary$zmean_stoich_3
# pathway_summary$qssa_4 <- pathway_summary$zcoverage_4 + pathway_summary$zmean_stoich_4

rm(cluster_stoich_1, cluster_stoich_2, cluster_stoich_3, cluster_stoich_4);gc()

# Cleaning QSSA -----------------------------------------------------------

# pathway_summary <- pathway_summary %>% arrange(-qssa)
pathway_summary$kegg <- sapply(pathway_summary$kegg_pathway, function(x){
  unlist(strsplit(x, "%", fixed = T))[1]
})

# Rounding values
pathway_summary[7:26] <- mapply(pathway_summary[7:26], FUN = round, digits = 3)

# selecting qssa data
qssa_results <- pathway_summary %>% 
  subset(n_acK_1 >= 2 | n_acK_2 >= 2 | n_acK_3 >= 2 | n_acK_4 >= 2) %>% 
  select(kegg, qssa_1:qssa_4) # qssa_1:qssa_5

# Arranging by row mean
qssa_results$mean <- rowMeans(qssa_results[2:5]) # qssa_results[2:6]


#############
test <- qssa_results %>% 
  arrange(-mean) %>% 
  mutate(rank_mean = 1:nrow(.)) %>% 
  arrange(-qssa_1) %>% 
  mutate(rank_1 = 1:nrow(.)) %>% 
  arrange(-qssa_2) %>% 
  mutate(rank_2 = 1:nrow(.)) %>% 
  arrange(-qssa_3) %>% 
  mutate(rank_3 = 1:nrow(.)) %>% 
  arrange(-qssa_4) %>% 
  mutate(rank_4 = 1:nrow(.))

test <- test %>% subset(rank_1 <= 10 | rank_2 <= 10 | rank_3 <= 10 | rank_4 <= 10)
test <- test %>% 
  arrange(rank_mean) %>% 
  select(-mean, rank_1, rank_2, rank_3, rank_4) %>% 
  gather(condition, qssa, qssa_1:qssa_4) # qssa_1:qssa_5

ggplot(subset(test)) +
  geom_tile(aes(x = condition, y = reorder(kegg, -rank_mean), fill = qssa)) +
  #scale_fill_distiller(palette = "RdBu", direction = -1) +
  scale_fill_gradient2(high = "red3", mid = "white",
                       low = "blue4", midpoint = 0, space = "Lab",
                       na.value = "grey50", guide = "colourbar", aesthetics = "fill") +
  theme_light(base_size = 12) +
  theme(axis.ticks = element_blank(),
        axis.title.y = element_blank()) +
  labs(title = "Quantitative Site Set functional score Analysis",
       x = "Stoichiometry bins")

############
# Long format of data for ggplot
qssa_results <- qssa_results %>% 
  arrange(-mean) %>% 
  mutate(rank = 1:nrow(.)) %>% 
  select(-mean) %>% 
  gather(condition, qssa, qssa_1:qssa_4) # qssa_1:qssa_5


# Plots -------------------------------------------------------------------

ggplot(subset(qssa_results, rank <= 20)) +
  geom_tile(aes(x = condition, y = reorder(kegg, -rank), fill = qssa)) +
  scale_fill_distiller(palette = "RdBu", direction = -1) +
  theme_light(base_size = 12) +
  theme(axis.ticks = element_blank(),
        axis.title.y = element_blank()) +
  labs(title = "Quantitative Site Set functional score Analysis",
       x = "Stoichiometry bins")

# basic heatmap blue less than 0
#pdf("stoich_cluster_sig_blue3_QSSA.pdf", width = 8, height = 6)
ggplot(subset(qssa_results, rank <= 20)) +
  geom_tile(aes(x = condition, y = reorder(kegg, -rank), fill = qssa)) +
  scale_fill_distiller(palette = "Blues", direction = -1, limits = c(-2.5, 0)) +
  theme_light(base_size = 12) +
  theme(axis.ticks = element_blank(),
        axis.title.y = element_blank()) +
  labs(title = "Quantitative Site Set functional score Analysis",
       x = "Stoichiometry bins")
#dev.off()

# basic heatmap red greater than 0
#pdf("stoich_cluster_sig_red_QSSA.pdf", width = 8, height = 6)
ggplot(subset(qssa_results, rank <= 20)) +
  geom_tile(aes(x = condition, y = reorder(kegg, -rank), fill = qssa)) +
  scale_fill_gradient2(high = "red3", low = "white",
                       midpoint = 0, na.value = "grey50", guide = "colourbar",
                       limits = c(0, NA)) +
  theme_light(base_size = 12) +
  theme(axis.ticks = element_blank(),
        axis.title.y = element_blank()) +
  labs(title = "Quantitative Site Set functional score Analysis",
       x = "Stoichiometry bins")
#dev.off()

