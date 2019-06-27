
# Libraries ---------------------------------------------------------------

options(stringsAsFactors = FALSE)
library(tidyverse)
library(rio)
library(stringr)



# Importing data ----------------------------------------------------------

# Peptide data
pep <- import("./Serum_Stimulated_HCT116/Cleaned/20190220_SSHCT116_tc_tidy_peptide_stoich_output.csv")

pep_fcm <- import("./Serum_Stimulated_HCT116/Fuzzy_clust/20190310_FCMResults_SerumStim_HCT116.csv")


# Formatting peptide data -------------------------------------------------


# labeling the replicates
pep$rep <- pep$sample_id
pep$sample_id <- NA
pep$sample_id[which(pep$sample_id == 1 & pep$time_point == 0)] <- 1
pep$sample_id[which(pep$sample_id == 2 & pep$time_point == 0)] <- 2
pep$sample_id[which(pep$sample_id == 3 & pep$time_point == 0)] <- 3
pep$sample_id[which(pep$sample_id == 1 & pep$time_point == 0.25)] <- 4
pep$sample_id[which(pep$sample_id == 2 & pep$time_point == 0.25)] <- 5
pep$sample_id[which(pep$sample_id == 3 & pep$time_point == 0.25)] <- 6
pep$sample_id[which(pep$sample_id == 1 & pep$time_point == 1)] <- 7
pep$sample_id[which(pep$sample_id == 2 & pep$time_point == 1)] <- 8
pep$sample_id[which(pep$sample_id == 3 & pep$time_point == 1)] <- 9
pep$sample_id[which(pep$sample_id == 1 & pep$time_point == 2)] <- 10
pep$sample_id[which(pep$sample_id == 2 & pep$time_point == 2)] <- 11
pep$sample_id[which(pep$sample_id == 3 & pep$time_point == 2)] <- 12
pep$sample_id[which(pep$sample_id == 1 & pep$time_point == 4)] <- 13
pep$sample_id[which(pep$sample_id == 2 & pep$time_point == 4)] <- 14
pep$sample_id[which(pep$sample_id == 3 & pep$time_point == 4)] <- 15

# Rearranging columns and deleting sample_id
pep <- pep %>% 
  select(time_point, rep, everything()) %>% 
  select(-sample_id) %>% 
  arrange(PG.ProteinGroups, k_site, time_point, rep)

# Cleaning modified sequence
pep$EG.ModifiedPeptide <- pep$EG.ModifiedPeptide %>% 
  str_replace_all("(_)$", "") %>% 
  str_replace_all("^_", "")

# Summarizing stoich data
pep_sum <- pep %>% 
  group_by(time_point, PG.ProteinGroups, PG.ProteinDescriptions, PG.Genes, 
           PEP.StrippedSequence, EG.ModifiedPeptide, k_site) %>% 
  summarise(stoich_mean = mean(stoich_corr, na.rm = TRUE),
            stoich_sd = sd(stoich_corr, na.rm = TRUE),
            stoich_n = n()) %>% 
  ungroup() %>% 
  rename("time" = "time_point") %>% 
  arrange(PG.ProteinGroups, k_site, time)

# Wide format of data
pep_wide <- pep_sum %>% 
  select(-stoich_sd, -stoich_n) %>% 
  spread(time, stoich_mean)

# Counting time point observations
pep_wide$tp_count <- apply(pep_wide[7:11], 1, function(x){
  length(which(!is.na(x)))
})

# Filtering for sites with all four observations
pep_wide <- pep_wide %>% 
  filter(tp_count == 5) %>% 
  select(-tp_count, -PG.ProteinDescriptions) %>% 
  unite(Identifier, PG.ProteinGroups, PG.Genes, 
        PEP.StrippedSequence, EG.ModifiedPeptide, k_site, sep = "_")


# Formatting FCM peptide data ---------------------------------------------


# Formatting names
names(pep_fcm)[1] <- "Identifier"
names(pep_fcm)[3] <- "0"
names(pep_fcm)[4] <- "0.25"
names(pep_fcm)[5] <- "1"
names(pep_fcm)[6] <- "2"
names(pep_fcm)[7] <- "4"

# Selecting only needed columns
pep_fcm <- pep_fcm %>% 
  select(Identifier:maxMembership) %>% 
  separate(Identifier, c("PG.ProteinGroups", "PG.Genes", "PEP.StrippedSequence", 
                         "EG.ModifiedPeptide", "k_site"), sep = "_")

# Adding the non-cluster group
pep_fcm$cluster[pep_fcm$isClusterMember == FALSE] <- 0

# Long format of data
pep_fcm_long <- pep_fcm %>% 
  gather(time, value, 7:11) %>% 
  arrange(PG.ProteinGroups) %>% 
  rename("stoich_cluster" = "cluster",
         "stoich_isClusterMember" = "isClusterMember",
         "stoich_maxMembership" = "maxMembership",
         "stoich_zscore" = "value")

# Adding cluster info with peptide summary data
pep_sum <- merge(pep_sum, pep_fcm_long, all.x = TRUE)
rm(pep_fcm_long);gc()


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
ggplot(pep_sum %>% filter(stoich_cluster != 0 & stoich_maxMembership >= 0.5), 
       aes(x = time, y = stoich_zscore)) +
  geom_line(aes(group = reorder(PEP.StrippedSequence, stoich_maxMembership), color = stoich_maxMembership), 
            size = 0.6) +
  scale_color_distiller(palette = "YlGnBu", direction = -1, limits = c(0.5, 1)) +
  theme_light(base_size = 14) +
  facet_wrap(~stoich_cluster, nrow = 2) +
  labs(color = "Cluster\nMembership\n",
       y = "z-score",
       x = "Time (hours)",
       title = "Stoichiometry dynamics")
#ggsave(filename = "figures/20190310_Peptide_cluster_profiles_HCT116.pdf", height = 6, width = 8)

### Writing Tables
# write.table(pep_sum, "20190310_SSHCT116_Clusters.csv", sep = ",", row.names = FALSE)



# Writing table for FCM clustering ----------------------------------------

# For peptide Fuzzy c-means clustering
#write.table(pep_wide, file = "./Serum_Stimulated_HCT116/Fuzzy_clust/20190310_SSHCT116_fuzzycluster_input_larger.csv", sep = ",", row.names = FALSE)


# Peptide Cluster Information
#write.table(pep_sum, file = "./Serum_Stimulated_HCT116/Fuzzy_clust/20190305_SSHCT116_FuzzyClust_n3.csv", sep = ",", row.names = FALSE)

