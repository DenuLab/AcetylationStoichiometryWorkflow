
# Libraries ---------------------------------------------------------------

options(stringsAsFactors = FALSE)
library(tidyverse)
library(rio)


# Importing data ----------------------------------------------------------

# Peptide data
pep1 <- import("./Serum_Stimulated_MCF7/Cleaned/20181101_SS_tc_tidy_peptide_stoich_output_REPROCESS_1sthalf.csv")
pep2 <- import("./Serum_Stimulated_MCF7/Cleaned/20181101_SS_tc_tidy_peptide_stoich_output_REPROCESS_2ndhalf.csv")
pep <- rbind(pep1, pep2)
rm(pep1, pep2);gc()

pep_fcm <- import("./Serum_Stimulated_MCF7/Fuzzy_cm_data/20181101_FCMResults_pep_reprocess_correctFASTA.csv")

#Protein data
prot <- import("./Serum_Stimulated_MCF7/Raw/20181031_SSMCF7_stoich_Spectronaut10_reprocess_correctFASTA_Report.csv")
prot_fcm <- import("./Serum_Stimulated_MCF7/Fuzzy_cm_data/20181101_FCMResults_prot_reprocess_correctFASTA.csv")
prot_msstats <- import("./Serum_Stimulated_MCF7/Cleaned/msstats/20181101_run_level_data.csv")


# Formatting peptide data -------------------------------------------------


# labeling the replicates
pep$rep <- NA
pep$rep[which(pep$sample_id == 1)] <- 1
pep$rep[which(pep$sample_id == 2)] <- 2
pep$rep[which(pep$sample_id == 3)] <- 3
pep$rep[which(pep$sample_id == 4)] <- 1
pep$rep[which(pep$sample_id == 5)] <- 2
pep$rep[which(pep$sample_id == 6)] <- 3
pep$rep[which(pep$sample_id == 7)] <- 1
pep$rep[which(pep$sample_id == 8)] <- 2
pep$rep[which(pep$sample_id == 9)] <- 3
pep$rep[which(pep$sample_id == 10)] <- 1
pep$rep[which(pep$sample_id == 11)] <- 2
pep$rep[which(pep$sample_id == 12)] <- 3

# Rearranging columns and deleting sample_id
pep <- pep %>% 
  select(time_point, rep, everything()) %>% 
  select(-sample_id) %>% 
  arrange(PG.ProteinGroups, k_site, time_point, rep)

# Cleaning modified sequence
pep$EG.ModifiedSequence <- pep$EG.ModifiedSequence %>% 
  str_replace_all("(_)$", "") %>% 
  str_replace_all("^_", "")

# Summarizing stoich data
pep_sum <- pep %>% 
  group_by(time_point, PG.ProteinGroups, PG.ProteinDescriptions, PG.Genes, 
           EG.StrippedSequence, EG.ModifiedSequence, k_site) %>% 
  summarise(stoich_mean = mean(stoich_corr, na.rm = TRUE),
            stoich_sd = sd(stoich_corr, na.rm = TRUE),
            stoich_n = n()) %>% 
  ungroup() %>% 
  rename("time" = "time_point") %>% 
  filter(stoich_n == 3) %>%
  arrange(PG.ProteinGroups, k_site, time)

# Wide format of data
pep_wide <- pep_sum %>% 
  select(-stoich_sd, -stoich_n) %>% 
  spread(time, stoich_mean)

# Counting time point observations
pep_wide$tp_count <- apply(pep_wide[7:10], 1, function(x){
  length(which(!is.na(x)))
})

# Filtering for sites with all four observations
pep_wide <- pep_wide %>% 
  filter(tp_count == 4) %>% 
  select(-tp_count, -PG.ProteinDescriptions) %>% 
  unite(Identifier, PG.ProteinGroups, PG.Genes, 
        EG.StrippedSequence, EG.ModifiedSequence, k_site, sep = "_")


# Formatting protein data -------------------------------------------------


# Removing duplicates
prot <- prot %>% 
  select(R.FileName:FG.ShapeQualityScore) %>% 
  distinct()

# Contaminating proteins
prot$contaminant <- FALSE
prot$contaminant[grep("CON", prot$PG.ProteinGroups)] <- TRUE

# Time point 
prot$time_point <- as.numeric(unlist(lapply(prot$R.FileName, function(x){
  unlist(strsplit(x, split = "_"))[6]
})))

# Sample_id 
prot$sample_id <- as.numeric(unlist(lapply(prot$R.FileName, function(x){
  unlist(strsplit(x, split = "_"))[7]
})))

# Fraction
prot$fraction <- as.numeric(unlist(lapply(prot$R.FileName, function(x){
  unlist(strsplit(x, split = "_"))[8]
})))

# Removing instances where string starts with ;
prot$PG.Genes <- gsub(pattern = "^;", replacement = "", prot$PG.Genes)
prot$PG.ProteinDescriptions <- gsub(pattern = "^;", replacement = "", prot$PG.ProteinDescriptions)
prot$PG.ProteinGroups <- gsub(pattern = "^;", replacement = "", prot$PG.ProteinGroups)

# labeling the biological replicate
prot$rep <- NA
prot$rep[which(prot$sample_id == 1)] <- 1
prot$rep[which(prot$sample_id == 2)] <- 2
prot$rep[which(prot$sample_id == 3)] <- 3
prot$rep[which(prot$sample_id == 4)] <- 1
prot$rep[which(prot$sample_id == 5)] <- 2
prot$rep[which(prot$sample_id == 6)] <- 3
prot$rep[which(prot$sample_id == 7)] <- 1
prot$rep[which(prot$sample_id == 8)] <- 2
prot$rep[which(prot$sample_id == 9)] <- 3
prot$rep[which(prot$sample_id == 10)] <- 1
prot$rep[which(prot$sample_id == 11)] <- 2
prot$rep[which(prot$sample_id == 12)] <- 3

# Ordering the protframe
prot <- prot %>% 
  arrange(PG.ProteinGroups, time_point, rep)

merge_data <- prot %>% 
  filter(contaminant == FALSE) %>% 
  select(PG.ProteinGroups, PG.Genes, PG.ProteinDescriptions) %>% 
  distinct() %>% 
  rename("protein" = "PG.ProteinGroups")

prot_msstats <- merge(merge_data, prot_msstats, all.y = TRUE)

# Removing data
rm(merge_data);gc()

# Summary table
prot_summary <- prot_msstats %>% 
  select(protein, PG.Genes, PG.ProteinDescriptions, logintensities, time, biorep) %>% 
  group_by(protein, PG.Genes, PG.ProteinDescriptions, time) %>% 
  summarize(mean_abundance = mean(logintensities),
            sd_abundance = sd(logintensities),
            count = n(),
            cv = (sd_abundance/mean_abundance)*100) %>% 
  ungroup() %>% 
  rename("PG.ProteinGroups" = "protein")


# Formatting FCM peptide data ---------------------------------------------


# Formatting names
names(pep_fcm)[1] <- "Identifier"
names(pep_fcm)[3] <- "0"
names(pep_fcm)[4] <- "1"
names(pep_fcm)[5] <- "2"
names(pep_fcm)[6] <- "4"

# Selecting only needed columns
pep_fcm <- pep_fcm %>% 
  select(Identifier:maxMembership) %>% 
  separate(Identifier, c("PG.ProteinGroups", "PG.Genes", "EG.StrippedSequence", 
                         "EG.ModifiedSequence", "k_site"), sep = "_")

# Adding the non-cluster group
pep_fcm$cluster[pep_fcm$isClusterMember == FALSE] <- 0

# Long format of data
pep_fcm_long <- pep_fcm %>% 
  gather(time, value, 7:10) %>% 
  arrange(PG.ProteinGroups) %>% 
  rename("stoich_cluster" = cluster,
         "stoich_isClusterMember" = isClusterMember,
         "stoich_maxMembership" = maxMembership,
         "stoich_zscore" = value)

# Adding cluster info with peptide summary data
pep_sum <- merge(pep_sum, pep_fcm_long, all.x = TRUE)
rm(pep_fcm_long);gc()


# Formatting FCM protein data ---------------------------------------------


# Formatting names
names(prot_fcm)[1] <- "Identifier"
names(prot_fcm)[3] <- "0"
names(prot_fcm)[4] <- "1"
names(prot_fcm)[5] <- "2"
names(prot_fcm)[6] <- "4"

# Selecting only needed columns
prot_fcm <- prot_fcm %>% 
  select(Identifier:maxMembership) %>%
  rename("PG.ProteinGroups" = "Identifier")


# Adding the non-cluster group
prot_fcm$cluster[prot_fcm$isClusterMember == FALSE] <- 0

# Long format of data
prot_fcm_long <- prot_fcm %>% 
  gather(time, value, 3:6) %>% 
  arrange(PG.ProteinGroups) %>% 
  rename("prot_cluster" = "cluster",
         "prot_isClusterMember" = "isClusterMember",
         "prot_maxMembership" = "maxMembership",
         "prot_zscore" = "value")

# Combining protein with cluster data
prot_summary <- merge(prot_summary, prot_fcm_long, all = TRUE)
rm(prot_fcm_long);gc()


# Merging Protein & Peptide data ------------------------------------------


protein_peptide_fcm <- merge(prot_summary, pep_sum, all = TRUE)

# Plots -------------------------------------------------------------------


# data for bar graph
pep_fcm2 <- pep_fcm %>% filter(isClusterMember == TRUE)
fcm_pep_bar <- as.data.frame(table(pep_fcm2$cluster))

# Bar graph of cluster members
ggplot(fcm_pep_bar %>% filter(Var1 != 0), 
       aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", fill = "grey", color = "black") +
  ggthemes::theme_pander(base_size = 14) +
  labs(x = "Cluster", y = "Count") +
  geom_text(aes(label = Freq), vjust = -0.5) +
  expand_limits(y = c(0,NA))
#ggsave(filename = "figures/20181105_Peptide_cluster_bar_graph.pdf", height = 6, width = 8)

# Peptide clustering
### Figure 5D
ggplot(pep_sum %>% filter(stoich_cluster != 0 & stoich_maxMembership >= 0.5), 
       aes(x = time, y = stoich_zscore)) +
  geom_line(aes(group = reorder(EG.StrippedSequence, stoich_maxMembership), color = stoich_maxMembership), 
            size = 0.75) +
  scale_color_distiller(palette = "YlGnBu", direction = -1, limits = c(0.5, 1)) +
  theme_light(base_size = 14) +
  facet_wrap(~stoich_cluster, nrow = 2) +
  labs(color = "Cluster\nMembership\n",
       y = "z-score",
       x = "Time (hours)",
       title = "Stoichiometry dynamics")
#ggsave(filename = "figures/20181105_Peptide_cluster_profiles.pdf", height = 6, width = 8)


# data for bar graph
prot_fcm2 <- prot_fcm %>% subset(isClusterMember == TRUE)
fcm_prot_bar <- as.data.frame(table(prot_fcm2$cluster))

# Bar graph of cluster members
ggplot(fcm_prot_bar %>% filter(Var1 != 0), aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", fill = "grey", color = "black") +
  ggthemes::theme_pander(base_size = 14) +
  labs(x = "Cluster", y = "Count") +
  geom_text(aes(label = Freq), vjust = -0.5) +
  expand_limits(y = c(0, 330))
#ggsave(filename = "figures/20181105_Protein_cluster_bar_graph.pdf", height = 6, width = 8)

# Protein clustering
### Figure 5E
prot_summary$time <- as.numeric(prot_summary$time)
ggplot(prot_summary %>% subset(prot_maxMembership >= 0.5), 
       aes(x = time, y = prot_zscore)) +
  geom_line(aes(group = reorder(PG.ProteinGroups, prot_maxMembership), color = prot_maxMembership), 
            size = 0.75) +
  scale_color_distiller(palette = "YlGnBu", direction = -1, limits = c(0.5,1)) +
  theme_light(base_size = 14) +
  facet_wrap(~prot_cluster, nrow = 2) +
  labs(color = "Cluster\nMembership\n",
       y = "z-score",
       x = "Time (hours)", 
       title = "Protein dynamics")
#ggsave(filename = "figures/20181105_Protein_cluster_profiles.pdf", height = 6, width = 8)


# Combined protein & stoich bar graph
prot_pep_fcm2 <- protein_peptide_fcm %>% 
  select(PG.ProteinGroups, PG.Genes, k_site, EG.ModifiedSequence, EG.StrippedSequence,
         stoich_cluster, stoich_isClusterMember, prot_cluster, prot_isClusterMember) %>%
  distinct()

### Figure 5G
ggplot(prot_pep_fcm2 %>% filter(stoich_isClusterMember == TRUE, prot_isClusterMember == TRUE)) +
  geom_bar(aes(x = factor(prot_cluster), fill = factor(stoich_cluster)), 
           stat = "count") +
  scale_fill_brewer(palette = "Set1") +
  ggthemes::theme_pander(base_size = 14) +
  labs(title = "Cluster comparisons", 
       y = "Count", 
       x = "Protein Cluster", 
       fill = "Stoichiometry\nCluster")
#ggsave(filename = "figures/20181105_protein_stoichiometry_cluster_bargraph2.pdf", height = 6, width = 8)

# Writing table for FCM clustering ----------------------------------------

# For peptide Fuzzy c-means clustering
#write.table(pep_wide, file = "./Serum_Stimulated_MCF7/Fuzzy_cm_data/20181101_SS_tc_peptide_fuzzy_clustering_data_REPROCESS_correctFASTA.csv", sep = ",", row.names = FALSE)

# Combined protein and stoichiometry cluster data
#write.table(protein_peptide_fcm, file = "./Serum_Stimulated_MCF7/Fuzzy_cm_data/20181102_Protein_stoich_FCM_clusters_reprocess_correctFASTA.csv", sep = ",", row.names = FALSE)

