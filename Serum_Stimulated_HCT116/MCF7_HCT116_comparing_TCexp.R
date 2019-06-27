options(stringsAsFactors = FALSE)
library(tidyverse)
library(rio)
library(stringr)
library(ggthemes)
library(BRAIN)
library(seqinr)
library(GGally)
library(ggExtra)

# Importing data ----------------------------------------------------------

# Peptide data
pep1 <- import("./Serum_Stimulated_MCF7/Cleaned/20181101_SS_tc_tidy_peptide_stoich_output_REPROCESS_1sthalf.csv")
pep2 <- import("./Serum_Stimulated_MCF7/Cleaned/20181101_SS_tc_tidy_peptide_stoich_output_REPROCESS_2ndhalf.csv")
pep_MCF7 <- rbind(pep1, pep2)
rm(pep1, pep2);gc()

pep_HCT <- read.csv("./Serum_Stimulated_HCT116/Cleaned/20190220_SSHCT116_tc_tidy_peptide_stoich_output.csv")


# labeling the replicates
pep_MCF7$rep <- NA
pep_MCF7$rep[which(pep_MCF7$sample_id == 1)] <- 1
pep_MCF7$rep[which(pep_MCF7$sample_id == 2)] <- 2
pep_MCF7$rep[which(pep_MCF7$sample_id == 3)] <- 3
pep_MCF7$rep[which(pep_MCF7$sample_id == 4)] <- 1
pep_MCF7$rep[which(pep_MCF7$sample_id == 5)] <- 2
pep_MCF7$rep[which(pep_MCF7$sample_id == 6)] <- 3
pep_MCF7$rep[which(pep_MCF7$sample_id == 7)] <- 1
pep_MCF7$rep[which(pep_MCF7$sample_id == 8)] <- 2
pep_MCF7$rep[which(pep_MCF7$sample_id == 9)] <- 3
pep_MCF7$rep[which(pep_MCF7$sample_id == 10)] <- 1
pep_MCF7$rep[which(pep_MCF7$sample_id == 11)] <- 2
pep_MCF7$rep[which(pep_MCF7$sample_id == 12)] <- 3

pep_HCT$rep <- pep_HCT$sample_id
pep_HCT$sample_id <- NA
pep_HCT$sample_id[which(pep_HCT$rep == 1 & pep_HCT$time_point == 0)] <- 1
pep_HCT$sample_id[which(pep_HCT$rep == 2 & pep_HCT$time_point == 0)] <- 2
pep_HCT$sample_id[which(pep_HCT$rep == 3 & pep_HCT$time_point == 0)] <- 3
pep_HCT$sample_id[which(pep_HCT$rep == 1 & pep_HCT$time_point == 0.25)] <- 4
pep_HCT$sample_id[which(pep_HCT$rep == 2 & pep_HCT$time_point == 0.25)] <- 5
pep_HCT$sample_id[which(pep_HCT$rep == 3 & pep_HCT$time_point == 0.25)] <- 6
pep_HCT$sample_id[which(pep_HCT$rep == 1 & pep_HCT$time_point == 1)] <- 7
pep_HCT$sample_id[which(pep_HCT$rep == 2 & pep_HCT$time_point == 1)] <- 8
pep_HCT$sample_id[which(pep_HCT$rep == 3 & pep_HCT$time_point == 1)] <- 9
pep_HCT$sample_id[which(pep_HCT$rep == 1 & pep_HCT$time_point == 2)] <- 10
pep_HCT$sample_id[which(pep_HCT$rep == 2 & pep_HCT$time_point == 2)] <- 11
pep_HCT$sample_id[which(pep_HCT$rep == 3 & pep_HCT$time_point == 2)] <- 12
pep_HCT$sample_id[which(pep_HCT$rep == 1 & pep_HCT$time_point == 4)] <- 13
pep_HCT$sample_id[which(pep_HCT$rep == 2 & pep_HCT$time_point == 4)] <- 14
pep_HCT$sample_id[which(pep_HCT$rep == 3 & pep_HCT$time_point == 4)] <- 15

MCF7_sum <- pep_MCF7 %>%
  filter(stoich_corr <= 1) %>%
  group_by(time_point, PG.ProteinGroups, PG.ProteinDescriptions, PG.Genes,
           EG.StrippedSequence, EG.ModifiedSequence, k_site) %>%
  summarise(stoich_mean = mean(stoich_corr, na.rm = TRUE),
            stoich_sd = sd(stoich_corr, na.rm = TRUE),
            stoich_n = n()) %>%
  ungroup() %>%
  rename("time_point" = "time") %>%
  arrange(PG.ProteinGroups, k_site, time)

MCF7_sum$EG.ModifiedSequence <- gsub("\\+42", "Acetyl (K)", MCF7_sum$EG.ModifiedSequence)
MCF7_sum$EG.ModifiedSequence <- gsub("\\+57", "\\+C2+H3+N+O", MCF7_sum$EG.ModifiedSequence)

HCT_sum <- pep_HCT %>%
  filter(stoich_corr <= 1) %>%
  group_by(time_point, PG.ProteinGroups, PG.ProteinDescriptions, PG.Genes,
           PEP.StrippedSequence, EG.ModifiedPeptide, k_site) %>%
  summarise(stoich_mean = mean(stoich_corr, na.rm = TRUE),
            stoich_sd = sd(stoich_corr, na.rm = TRUE),
            stoich_n = n()) %>%
  ungroup() %>%
  rename("time_point" = "time") %>%
  arrange(PG.ProteinGroups, k_site, time)


############
HCT_wide_mean <- HCT_sum %>%
  select(-stoich_sd, -stoich_n) %>%
  spread(time, stoich_mean) %>%
  rename("0" = "Hr0_mean",
         "0.25" = "Hr0.25_mean",
         "1" = "Hr1_mean",
         "2" = "Hr2_mean",
         "4" = "Hr4_mean")
HCT_wide_sd <- HCT_sum %>%
  select(-stoich_mean, -stoich_n) %>%
  spread(time, stoich_sd) %>%
  rename("0" = "Hr0_sd",
         "0.25" = "Hr0.25_sd",
         "1" = "Hr1_sd",
         "2" = "Hr2_sd",
         "4" = "Hr4_sd")
HCT_wide_n <- HCT_sum %>%
  select(-stoich_sd, -stoich_mean) %>%
  spread(time, stoich_n) %>%
  rename("0" = "Hr0_n",
         "0.25" = "Hr0.25_n",
         "1" = "Hr1_n",
         "2" = "Hr2_n",
         "4" = "Hr4_n")
HCT_wide <- merge(HCT_wide_mean, HCT_wide_sd)
HCT_wide <- merge(HCT_wide, HCT_wide_n)
rm(HCT_wide_mean, HCT_wide_sd, HCT_wide_n)
HCT_wide$change_0.25 <- HCT_wide$Hr0.25_mean - HCT_wide$Hr0_mean
HCT_wide$change_1 <- HCT_wide$Hr1_mean - HCT_wide$Hr0_mean
HCT_wide$change_2 <- HCT_wide$Hr2_mean - HCT_wide$Hr0_mean
HCT_wide$change_4 <- HCT_wide$Hr4_mean - HCT_wide$Hr0_mean

HCT_wide <- HCT_wide %>% select(PG.ProteinGroups, PG.Genes, k_site, everything())

MCF7_wide_mean <- MCF7_sum %>%
  select(-stoich_sd, -stoich_n) %>%
  spread(time, stoich_mean) %>%
  rename("0" = "Hr0_mean",
         "1" = "Hr1_mean",
         "2" = "Hr2_mean",
         "4" = "Hr4_mean")
MCF7_wide_sd <- MCF7_sum %>%
  select(-stoich_mean, -stoich_n) %>%
  spread(time, stoich_sd) %>%
  rename("0" = "Hr0_sd",
         "1" = "Hr1_sd",
         "2" = "Hr2_sd",
         "4" = "Hr4_sd")
MCF7_wide_n <- MCF7_sum %>%
  select(-stoich_sd, -stoich_mean) %>%
  spread(time, stoich_n) %>%
  rename("0" = "Hr0_n",
         "1" = "Hr1_n",
         "2" = "Hr2_n",
         "4" = "Hr4_n")
MCF7_wide <- merge(MCF7_wide_mean, MCF7_wide_sd)
MCF7_wide <- merge(MCF7_wide, MCF7_wide_n)
rm(MCF7_wide_mean, MCF7_wide_sd, MCF7_wide_n)
MCF7_wide$short_change <- MCF7_wide$Hr1_mean - MCF7_wide$Hr0_mean
MCF7_wide$mid_change <- MCF7_wide$Hr2_mean - MCF7_wide$Hr0_mean
MCF7_wide$long_change <- MCF7_wide$Hr4_mean - MCF7_wide$Hr0_mean

MCF7_wide <- MCF7_wide %>% select(PG.ProteinGroups, PG.Genes, k_site, everything())


#### COMBINING THE TWO EXPERIMENTS ALL TOGETHER
MCF7_sum$experiment <- "MCF7"
HCT_sum$experiment <- "HCT116"

MCF7_sum <- MCF7_sum %>%
  rename( "EG.StrippedSequence" = "PEP.StrippedSequence",
          "EG.ModifiedSequence" = "EG.ModifiedPeptide")

combined <- rbind(MCF7_sum, HCT_sum)

# spread <- combined %>%
#   select(time, PG.ProteinGroups, PG.ProteinDescriptions, PG.Genes,
#          PEP.StrippedSequence, EG.ModifiedPeptide, k_site, stoich_mean, experiment) %>%
#   spread(key = experiment, value = stoich_mean)
# 
# ggplot((spread), aes(HCT116, MCF7)) +
#   geom_point(aes(color = as.factor(time)))

colnames(MCF7_wide)[7:21] <- paste(colnames(MCF7_wide)[7:21], "MCF7", sep = "_")
colnames(HCT_wide)[7:25] <- paste(colnames(HCT_wide)[7:25], "HCT", sep = "_")

MCF7_wide$experiment <- NULL
HCT_wide$experiment <- NULL

merge <- merge(MCF7_wide, HCT_wide, all = TRUE)

# ggplot(subset(HCT_sum, PG.Genes == "FIS1" & k_site == "K89"), aes(time, stoich_mean)) +
#   geom_errorbar(aes(ymin = stoich_mean-stoich_sd, ymax=stoich_mean+stoich_sd), width=.1) +
#   geom_line() +
#   geom_point( color = "red") +
#   theme_light(base_size = 10) +
#   scale_y_continuous(limits = c(0:1)) +
#   labs(title = "SSHCT116", y = "Mean Stoichiometry", x = "Time (h)")
# #ggsave("HCT_Fis1_K89.pdf", height = 3, width = 3)
# 
# ggplot(subset(MCF7_sum, PG.Genes == "FIS1" & k_site == "K89"), aes(time, stoich_mean)) +
#   geom_errorbar(aes(ymin=stoich_mean-stoich_sd, ymax=stoich_mean+stoich_sd), width=.1) +
#   geom_line() +
#   geom_point(color = "red") +
#   theme_light(base_size = 10) +
#   scale_y_continuous(limits = c(0:1)) +
#   labs(title = "SSMCF7", y = "Mean Stoichiometry", x = "Time (h)")
# #ggsave("MCF7_Fis1_K89.pdf", height = 3, width = 3)

ggplot(subset(combined, (PG.Genes == "FIS1" & k_site == "K89") | 
                (PG.Genes == "SRXN1" & k_site == "K116")), 
       aes(time, stoich_mean*100, group = experiment)) +
  geom_errorbar(aes(color = experiment, ymin = (stoich_mean*100 - stoich_sd*100), ymax = (stoich_mean*100 + stoich_sd*100))) +
  geom_line(aes(color = experiment), lwd = 1) +
  geom_point(aes(color = experiment), size = 1.5) +
  scale_color_brewer(palette = "Set1") +
  facet_wrap(~PG.Genes, nrow = 2) +
  theme_light(base_size = 10) +
  labs(y = "Average Stoichiometry (%)")
#ggsave("ExampleProteins_SerumStimulation_errorbars.pdf", height = 3, width = 3)

################################
#Merging the wide formats with stats
MCF7_wide_stats <- read.csv("./Serum_Stimulated_MCF7/Cleaned/20190402_SSMCF7_stoich_wide_stats_sigalltimepts.csv")
HCT_wide_stats <- read.csv("./Serum_Stimulated_HCT116/Cleaned/20190305_HCT116_stoich_wide_stats_sigalltimepts.csv")

MCF7_wide_stats$EG.ModifiedSequence <- gsub("\\+42", "Acetyl (K)", MCF7_wide_stats$EG.ModifiedSequence)
MCF7_wide_stats$EG.ModifiedSequence <- gsub("\\+57", "\\+C2+H3+N+O", MCF7_wide_stats$EG.ModifiedSequence)

MCF7_wide_stats <- MCF7_wide_stats %>% 
  rename( "EG.StrippedSequence" = "PEP.StrippedSequence",
          "EG.ModifiedSequence" = "EG.ModifiedPeptide")

# HCT_wide_stats <- HCT_wide_stats %>% 
#   select(-Hr0_sd, -Hr0.25_sd, -Hr1_sd, -Hr2_sd, -Hr4_sd,
#          -Hr0_n, -Hr0.25_n, -Hr1_n, -Hr2_n, -Hr4_n)
MCF7_wide_stats$unique_name <- NULL

colnames(MCF7_wide_stats)[7:27] <- paste(colnames(MCF7_wide_stats)[7:27], "MCF7", sep = "_")
colnames(HCT_wide_stats)[7:35] <- paste(colnames(HCT_wide_stats)[7:35], "HCT", sep = "_")

merge2 <- merge(MCF7_wide_stats, HCT_wide_stats, all = TRUE)

#write.table(merge2, "20190402_MCF7_HCT116_combined_with_stats_sigalltimepts.csv", sep = ",", row.names = FALSE)

sig_both <- merge2 %>% 
  subset((pval_time_0_1_MCF7 <= 0.05 | pval_time_0_2_MCF7 <= 0.05 | pval_time_0_4_MCF7 <= 0.05 |
            pval_time_1_2_MCF7 <= 0.05 | pval_time_1_4_MCF7 <= 0.05 | pval_time_2_4_MCF7 <= 0.05) &
           (pval_time_0_0.25_HCT <= 0.05 | pval_time_0_1_HCT <= 0.05 | pval_time_0_2_HCT <= 0.05 | pval_time_0_4_HCT <= 0.05 |
              pval_time_0.25_1_HCT <= 0.05 | pval_time_0.25_2_HCT <= 0.05 | pval_time_0.25_4_HCT <= 0.05 |
              pval_time_1_2_HCT <= 0.05 | pval_time_1_4_HCT <= 0.05 | pval_time_2_4_HCT <= 0.05))
#write.table(sig_both, "SS_HCT_MCF7_Ksite_Sig_BothExp.csv", sep = ",", row.names = FALSE)

changing <- sig_both %>% 
  subset((abs(short_change_MCF7) >= 0.05 | abs(med_change_MCF7) >= 0.05 | abs(long_change_MCF7) >= 0.05) &
           (abs(change_0.25_HCT) >= 0.05 | abs(change_1_HCT) >= 0.05 | 
              abs(change_2_HCT) >= 0.05 | abs(change_4_HCT) >= 0.05))


ggplot(subset(combined, (PG.Genes == "FIS1" & k_site == "K89") | (PG.Genes == "SRXN1" & k_site == "K116")), 
       aes(time, stoich_mean*100, group = experiment)) +
  geom_errorbar(aes(color = experiment, ymin = (stoich_mean*100 - stoich_sd*100), ymax = (stoich_mean*100 + stoich_sd*100))) +
  geom_line(aes(color = experiment), lwd = 1) +
  geom_point(aes(color = experiment), size = 1.5) +
  scale_color_brewer(palette = "Set1") +
  facet_wrap(~PG.Genes, nrow = 2) +
  theme_light(base_size = 10) +
  labs(y = "Average Stoichiometry (%)")
#ggsave("ExampleProteins_SerumStimulation_errorbars.pdf", height = 3, width = 3)

## VENN DIAGRAMS
v <- venneuler(c(A=450, B=1800, "A&B"=230))
plot(v)



