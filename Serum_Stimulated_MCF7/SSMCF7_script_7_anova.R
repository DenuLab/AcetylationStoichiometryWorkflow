library(dplyr)
library(tidyr)
library(readxl)
library(ggplot2)
library(purrr)
library(stringr)

ms1 <- read.csv("./Serum_Stimulated_MCF7/Cleaned/20181101_SS_tc_tidy_peptide_stoich_output_REPROCESS_1sthalf.csv")
ms2 <- read.csv("./Serum_Stimulated_MCF7/Cleaned/20181101_SS_tc_tidy_peptide_stoich_output_REPROCESS_2ndhalf.csv")
ms <- rbind(ms1, ms2)
rm(ms1, ms2);gc()

ms$EG.ModifiedSequence <- ms$EG.ModifiedSequence %>% 
  str_replace_all("(_)$", "") %>% 
  str_replace_all("^_", "")

ms <- ms %>% subset(stoich_corr <= 1)

mss <- split(ms, paste(ms$PG.ProteinGroups, ms$k_site, ms$EG.ModifiedSequence, sep = "_"))


# Looking for sites that change at any time point
tmpfn <- function(dat) {
  if(!any(dat$time_point == 0))
    return(NULL)
  if(length(unique(dat$time_point)) < 2)
    return(NULL)
  fit <- lm(stoich_corr ~ factor(time_point), dat)
  rsq <- summary(fit)$r.squared
  if(rsq < 0.99)
    pval <- anova(fit)[,"Pr(>F)"][1]
  else
    pval <- NA
  coefs <- t(coef(fit))
  data.frame(coefs,rsq, pval)
}
(msdiff <- mss %>%
    map(tmpfn) %>%
    compact() %>% 
    bind_rows(.id = "k_site") %>%
    rename(
      "control" = "X.Intercept.",
      "time1" = "factor.time_point.1",
      "time2" = "factor.time_point.2",
      "time4" = "factor.time_point.4") %>%
    mutate_if(is.numeric, function(x) round(x, 4)) %>%
    arrange(pval)) %>%
  filter(pval <= .1)

# Looking for significant changes in the first time_point
tmpfn <- function(dat) {
  dat <- dat %>% filter(time_point < 2)
  if(length(unique(dat$time_point)) < 2)
    return(NULL)
  fit <- lm(stoich_corr ~ factor(time_point), dat)
  rsq <- summary(fit)$r.squared
  if(rsq < 0.99)
    pval <- anova(fit)[,"Pr(>F)"][1]
  else
    pval <- NA
  coefs <- t(coef(fit))
  data.frame(coefs,rsq, pval)
}
(msdiff2 <- mss %>%
    map(tmpfn) %>%
    compact() %>% 
    bind_rows(.id = "k_site") %>%
    rename(
      "control" = "X.Intercept.",
      "time1" = "factor.time_point.1") %>%
    mutate_if(is.numeric, function(x) round(x, 4)) %>%
    arrange(pval)) %>%
  filter(pval <= .05)

# My attempt to change things
tmpfn <- function(dat) {
  dat <- dat %>% filter(time_point == 0 | time_point == 2)
  if(length(unique(dat$time_point)) < 2)
    return(NULL)
  fit <- lm(stoich_corr ~ factor(time_point), dat)
  rsq <- summary(fit)$r.squared
  if(rsq < 0.99)
    pval <- anova(fit)[,"Pr(>F)"][1]
  else
    pval <- NA
  coefs <- t(coef(fit))
  data.frame(coefs,rsq, pval)
}
(msdiff3 <- mss %>%
    map(tmpfn) %>%
    compact() %>% 
    bind_rows(.id = "k_site") %>%
    rename(
      "control" = "X.Intercept.",
      "time2" = "factor.time_point.2") %>%
    mutate_if(is.numeric, function(x) round(x, 4)) %>%
    arrange(pval)) %>%
  filter(pval <= .05)

# Sig 0 to 4
tmpfn <- function(dat) {
  dat <- dat %>% filter(time_point == 0 | time_point == 4)
  if(length(unique(dat$time_point)) < 2)
    return(NULL)
  fit <- lm(stoich_corr ~ factor(time_point), dat)
  rsq <- summary(fit)$r.squared
  if(rsq < 0.99)
    pval <- anova(fit)[,"Pr(>F)"][1]
  else
    pval <- NA
  coefs <- t(coef(fit))
  data.frame(coefs,rsq, pval)
}
(msdiff4 <- mss %>%
    map(tmpfn) %>%
    compact() %>% 
    bind_rows(.id = "k_site") %>%
    rename(
      "control" = "X.Intercept.",
      "time4" = "factor.time_point.4") %>%
    mutate_if(is.numeric, function(x) round(x, 4)) %>%
    arrange(pval)) %>%
  filter(pval <= .05)

# Sig 1 to 2
tmpfn <- function(dat) {
  dat <- dat %>% filter(time_point == 1 | time_point == 2)
  if(length(unique(dat$time_point)) < 2)
    return(NULL)
  fit <- lm(stoich_corr ~ factor(time_point), dat)
  rsq <- summary(fit)$r.squared
  if(rsq < 0.99)
    pval <- anova(fit)[,"Pr(>F)"][1]
  else
    pval <- NA
  coefs <- t(coef(fit))
  data.frame(coefs,rsq, pval)
}
(msdiff5 <- mss %>%
    map(tmpfn) %>%
    compact() %>% 
    bind_rows(.id = "k_site") %>%
    rename(
      "control" = "X.Intercept.",
      "time2" = "factor.time_point.2") %>%
    mutate_if(is.numeric, function(x) round(x, 4)) %>%
    arrange(pval)) %>%
  filter(pval <= .05)

# Sig 1 to 4
tmpfn <- function(dat) {
  dat <- dat %>% filter(time_point == 1 | time_point == 4)
  if(length(unique(dat$time_point)) < 2)
    return(NULL)
  fit <- lm(stoich_corr ~ factor(time_point), dat)
  rsq <- summary(fit)$r.squared
  if(rsq < 0.99)
    pval <- anova(fit)[,"Pr(>F)"][1]
  else
    pval <- NA
  coefs <- t(coef(fit))
  data.frame(coefs,rsq, pval)
}
(msdiff6 <- mss %>%
    map(tmpfn) %>%
    compact() %>% 
    bind_rows(.id = "k_site") %>%
    rename(
      "control" = "X.Intercept.",
      "time4" = "factor.time_point.4") %>%
    mutate_if(is.numeric, function(x) round(x, 4)) %>%
    arrange(pval)) %>%
  filter(pval <= .05)

# Sig 2 to 4
tmpfn <- function(dat) {
  dat <- dat %>% filter(time_point == 2 | time_point == 4)
  if(length(unique(dat$time_point)) < 2)
    return(NULL)
  fit <- lm(stoich_corr ~ factor(time_point), dat)
  rsq <- summary(fit)$r.squared
  if(rsq < 0.99)
    pval <- anova(fit)[,"Pr(>F)"][1]
  else
    pval <- NA
  coefs <- t(coef(fit))
  data.frame(coefs,rsq, pval)
}
(msdiff7 <- mss %>%
    map(tmpfn) %>%
    compact() %>% 
    bind_rows(.id = "k_site") %>%
    rename(
      "control" = "X.Intercept.",
      "time4" = "factor.time_point.4") %>%
    mutate_if(is.numeric, function(x) round(x, 4)) %>%
    arrange(pval)) %>%
  filter(pval <= .05)

#### 
names(msdiff)[1] <- "unique_name"
names(msdiff)[6] <- "rsq_all"
names(msdiff)[7] <- "pval_all"
msdiff <- msdiff %>% select(c(1,6:7))

names(msdiff2)[1] <- "unique_name"
names(msdiff2)[4] <- "rsq_time_0_1"
names(msdiff2)[5] <- "pval_time_0_1"
msdiff2 <- msdiff2 %>% select(c(1,4:5))

names(msdiff3)[1] <- "unique_name"
names(msdiff3)[4] <- "rsq_time_0_2"
names(msdiff3)[5] <- "pval_time_0_2"
msdiff3 <- msdiff3 %>% select(c(1,4:5))

names(msdiff4)[1] <- "unique_name"
names(msdiff4)[4] <- "rsq_time_0_4"
names(msdiff4)[5] <- "pval_time_0_4"
msdiff4 <- msdiff4 %>% select(c(1,4:5))

names(msdiff5)[1] <- "unique_name"
names(msdiff5)[4] <- "rsq_time_1_2"
names(msdiff5)[5] <- "pval_time_1_2"
msdiff5 <- msdiff5 %>% select(c(1,4:5))

names(msdiff6)[1] <- "unique_name"
names(msdiff6)[4] <- "rsq_time_1_4"
names(msdiff6)[5] <- "pval_time_1_4"
msdiff6 <- msdiff6 %>% select(c(1,4:5))

names(msdiff7)[1] <- "unique_name"
names(msdiff7)[4] <- "rsq_time_2_4"
names(msdiff7)[5] <- "pval_time_2_4"
msdiff7 <- msdiff7 %>% select(c(1,4:5))

ssmcf7_sig <- merge(msdiff, msdiff2, by = "unique_name", all.x = TRUE)
ssmcf7_sig <- merge(ssmcf7_sig, msdiff3, by = "unique_name", all.x = TRUE)
ssmcf7_sig <- merge(ssmcf7_sig, msdiff4, by = "unique_name", all.x = TRUE)
ssmcf7_sig <- merge(ssmcf7_sig, msdiff5, by = "unique_name", all.x = TRUE)
ssmcf7_sig <- merge(ssmcf7_sig, msdiff6, by = "unique_name", all.x = TRUE)
ssmcf7_sig <- merge(ssmcf7_sig, msdiff7, by = "unique_name", all.x = TRUE)

#write.table(ssmcf7_sig, "./Serum_Stimulated_MCF7/Cleaned/20190402_SSMCF7_sig_alltime_tomergewith_peptidedata.csv", sep = ",", row.names = FALSE)

# labeling the replicates
ms$rep <- NA
ms$rep[which(ms$sample_id == 1)] <- 1
ms$rep[which(ms$sample_id == 2)] <- 2
ms$rep[which(ms$sample_id == 3)] <- 3
ms$rep[which(ms$sample_id == 4)] <- 1
ms$rep[which(ms$sample_id == 5)] <- 2
ms$rep[which(ms$sample_id == 6)] <- 3
ms$rep[which(ms$sample_id == 7)] <- 1
ms$rep[which(ms$sample_id == 8)] <- 2
ms$rep[which(ms$sample_id == 9)] <- 3
ms$rep[which(ms$sample_id == 10)] <- 1
ms$rep[which(ms$sample_id == 11)] <- 2
ms$rep[which(ms$sample_id == 12)] <- 3

# Rearranging columns and deleting sample_id
ms <- ms %>% 
  select(time_point, rep, everything()) %>% 
  select(-sample_id) %>% 
  arrange(PG.ProteinGroups, k_site, time_point, rep)

# Cleaning modified sequence
ms$EG.ModifiedSequence <- ms$EG.ModifiedSequence %>% 
  str_replace_all("(_)$", "") %>% 
  str_replace_all("^_", "")

# Summarizing stoich data
ms_sum <- ms %>% 
  group_by(PG.ProteinGroups, PG.ProteinDescriptions, PG.Genes, 
           EG.StrippedSequence, EG.ModifiedSequence, k_site, time_point) %>% 
  summarise(stoich_mean = mean(stoich_corr, na.rm = TRUE),
            stoich_sd = sd(stoich_corr, na.rm = TRUE),
            stoich_n = n()) %>% 
  ungroup() %>% 
  rename("time" = "time_point") %>% 
  arrange(PG.ProteinGroups, k_site, time)

ms_sum$unique_name <- paste(ms_sum$PG.ProteinGroups, ms_sum$k_site, ms_sum$EG.ModifiedSequence, sep = "_")

ms_sum <- merge(ms_sum, ssmcf7_sig, all.x = TRUE, by = "unique_name")

ms_wide <- ms_sum %>% 
  select(-stoich_sd, -stoich_n) %>% 
  spread(time, stoich_mean)

# Counting time point observations
ms_wide$tp_count <- apply(ms_wide[16:19], 1, function(x){
  length(which(!is.na(x)))
})

ms_wide_mean <- ms_sum %>% 
  select(PG.ProteinGroups, PG.Genes, k_site, PG.ProteinDescriptions, time,
         EG.StrippedSequence, EG.ModifiedSequence, stoich_mean, pval_time_0_1,
         pval_time_0_2, pval_time_0_4, pval_time_1_2, pval_time_1_4, pval_time_2_4) %>% 
  filter(!is.na(k_site)) %>% 
  spread(time, stoich_mean) %>% 
  rename("Hr0_mean" = "0",
         "Hr1_mean" = "1",
         "Hr2_mean" = "2",
         "Hr4_mean" = "4")

ms_wide_sd <- ms_sum %>% 
  select(PG.ProteinGroups, PG.Genes, k_site, PG.ProteinDescriptions, time,
         EG.StrippedSequence, EG.ModifiedSequence, stoich_sd, pval_time_0_1,
         pval_time_0_2, pval_time_0_4, pval_time_1_2, pval_time_1_4, pval_time_2_4) %>% 
  filter(!is.na(k_site)) %>% 
  spread(time, stoich_sd) %>% 
  rename("Hr0_sd" = "0",
         "Hr1_sd" = "1",
         "Hr2_sd" = "2",
         "Hr4_sd" = "4")

ms_wide_n <- ms_sum %>% 
  select(PG.ProteinGroups, PG.Genes, k_site, PG.ProteinDescriptions, time,
         EG.StrippedSequence, EG.ModifiedSequence, stoich_n, pval_time_0_1,
         pval_time_0_2, pval_time_0_4, pval_time_1_2, pval_time_1_4, pval_time_2_4) %>% 
  filter(!is.na(k_site)) %>% 
  spread(time, stoich_n) %>% 
  rename("Hr0_n" = "0",
         "Hr1_n" = "1",
         "Hr2_n" = "2",
         "Hr4_n" = "4")

ms_wide <- merge(ms_wide_mean, ms_wide_sd, all = TRUE)
ms_wide <- merge(ms_wide, ms_wide_n, all = TRUE)

rm(ms_wide_mean, ms_wide_sd, ms_wide_n);gc()

ms_wide$short_change <- ms_wide$Hr1_mean - ms_wide$Hr0_mean
ms_wide$med_change <- ms_wide$Hr2_mean - ms_wide$Hr0_mean
ms_wide$long_change <- ms_wide$Hr4_mean - ms_wide$Hr0_mean

#write.table(ms_wide, "./Serum_Stimulated_MCF7/Cleaned/20190402_SSMCF7_stoich_wide_stats_sigalltimepts.csv", sep = ",", row.names = FALSE)

