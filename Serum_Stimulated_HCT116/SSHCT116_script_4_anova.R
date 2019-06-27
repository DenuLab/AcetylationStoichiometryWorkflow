options(stringsAsFactors = FALSE)
library(dplyr)
library(tidyr)
library(readxl)
library(ggplot2)
library(purrr)
library(stringr)

ms <- read.csv("./Serum_Stimulated_HCT116/Cleaned/20190220_SSHCT116_tc_tidy_peptide_stoich_output.csv")

ms$EG.ModifiedPeptide <- ms$EG.ModifiedPeptide %>% 
  str_replace_all("(_)$", "") %>% 
  str_replace_all("^_", "")

ms <- ms %>% subset(stoich_corr <= 1)

mss <- split(ms, paste(ms$PG.ProteinGroups, ms$k_site, ms$EG.ModifiedPeptide, sep = "_"))


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
      "time0.25" = "factor.time_point.0.25",
      "time1" = "factor.time_point.1",
      "time2" = "factor.time_point.2",
      "time4" = "factor.time_point.4") %>%
    mutate_if(is.numeric, function(x) round(x, 4)) %>%
    arrange(pval)) %>%
  filter(pval <= .1)

# Looking for significant changes between 0 and 0.25
tmpfn <- function(dat) {
  dat <- dat %>% filter(time_point < 1)
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
(msdiff_0.25 <- mss %>%
    map(tmpfn) %>%
    compact() %>% 
    bind_rows(.id = "k_site") %>%
    rename(
      "control" = "X.Intercept.",
      "time0.25" = "factor.time_point.0.25") %>%
    mutate_if(is.numeric, function(x) round(x, 4)) %>%
    arrange(pval)) %>%
  filter(pval <= .05)

# Looking for significant changes between 0 and 1
tmpfn <- function(dat) {
  dat <- dat %>% filter(time_point == 0 | time_point == 1)
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
(msdiff_1 <- mss %>%
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
(msdiff_2 <- mss %>%
    map(tmpfn) %>%
    compact() %>% 
    bind_rows(.id = "k_site") %>%
    rename(
      "control" = "X.Intercept.",
      "time2" = "factor.time_point.2") %>%
    mutate_if(is.numeric, function(x) round(x, 4)) %>%
    arrange(pval)) %>%
  filter(pval <= .05)

# My attempt to change things - between 0 and 4
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
(msdiff_4 <- mss %>%
    map(tmpfn) %>%
    compact() %>% 
    bind_rows(.id = "k_site") %>%
    rename(
      "control" = "X.Intercept.",
      "time4" = "factor.time_point.4") %>%
    mutate_if(is.numeric, function(x) round(x, 4)) %>%
    arrange(pval)) %>%
  filter(pval <= .05)

#Significance between 0.25 and 1
tmpfn <- function(dat) {
  dat <- dat %>% filter(time_point == 0.25 | time_point == 1)
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
(msdiff_5 <- mss %>%
    map(tmpfn) %>%
    compact() %>% 
    bind_rows(.id = "k_site") %>%
    rename(
      "control" = "X.Intercept.",
      "time1" = "factor.time_point.1") %>%
    mutate_if(is.numeric, function(x) round(x, 4)) %>%
    arrange(pval)) %>%
  filter(pval <= .05)

#Significance between 0.25 and 2
tmpfn <- function(dat) {
  dat <- dat %>% filter(time_point == 0.25 | time_point == 2)
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
(msdiff_6 <- mss %>%
    map(tmpfn) %>%
    compact() %>% 
    bind_rows(.id = "k_site") %>%
    rename(
      "control" = "X.Intercept.",
      "time2" = "factor.time_point.2") %>%
    mutate_if(is.numeric, function(x) round(x, 4)) %>%
    arrange(pval)) %>%
  filter(pval <= .05)

#Significance between 0.25 and 4
tmpfn <- function(dat) {
  dat <- dat %>% filter(time_point == 0.25 | time_point == 4)
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
(msdiff_7 <- mss %>%
    map(tmpfn) %>%
    compact() %>% 
    bind_rows(.id = "k_site") %>%
    rename(
      "control" = "X.Intercept.",
      "time4" = "factor.time_point.4") %>%
    mutate_if(is.numeric, function(x) round(x, 4)) %>%
    arrange(pval)) %>%
  filter(pval <= .05)

#Significance between 1 and 2
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
(msdiff_8 <- mss %>%
    map(tmpfn) %>%
    compact() %>% 
    bind_rows(.id = "k_site") %>%
    rename(
      "control" = "X.Intercept.",
      "time2" = "factor.time_point.2") %>%
    mutate_if(is.numeric, function(x) round(x, 4)) %>%
    arrange(pval)) %>%
  filter(pval <= .05)

#Significance between 1 and 4
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
(msdiff_9 <- mss %>%
    map(tmpfn) %>%
    compact() %>% 
    bind_rows(.id = "k_site") %>%
    rename(
      "control" = "X.Intercept.",
      "time2" = "factor.time_point.4") %>%
    mutate_if(is.numeric, function(x) round(x, 4)) %>%
    arrange(pval)) %>%
  filter(pval <= .05)

#Significance between 2 and 4
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
(msdiff_10 <- mss %>%
    map(tmpfn) %>%
    compact() %>% 
    bind_rows(.id = "k_site") %>%
    rename(
      "control" = "X.Intercept.",
      "time2" = "factor.time_point.4") %>%
    mutate_if(is.numeric, function(x) round(x, 4)) %>%
    arrange(pval)) %>%
  filter(pval <= .05)


#
#### 
# names(msdiff)[1] <- "unique_name"
# names(msdiff)[6] <- "rsq_all"
# names(msdiff)[7] <- "pval_all"
# msdiff <- msdiff %>% select(c(1,6:7))

names(msdiff_0.25)[1] <- "unique_name"
names(msdiff_0.25)[4] <- "rsq_time_0_0.25"
names(msdiff_0.25)[5] <- "pval_time_0_0.25"
msdiff_0.25 <- msdiff_0.25 %>% select(c(1,4:5))

names(msdiff_1)[1] <- "unique_name"
names(msdiff_1)[4] <- "rsq_time_0_1"
names(msdiff_1)[5] <- "pval_time_0_1"
msdiff_1 <- msdiff_1 %>% select(c(1,4:5))

names(msdiff_2)[1] <- "unique_name"
names(msdiff_2)[4] <- "rsq_time_0_2"
names(msdiff_2)[5] <- "pval_time_0_2"
msdiff_2 <- msdiff_2 %>% select(c(1,4:5))

names(msdiff_4)[1] <- "unique_name"
names(msdiff_4)[4] <- "rsq_time_0_4"
names(msdiff_4)[5] <- "pval_time_0_4"
msdiff_4 <- msdiff_4 %>% select(c(1,4:5))

names(msdiff_5)[1] <- "unique_name"
names(msdiff_5)[4] <- "rsq_time_0.25_1"
names(msdiff_5)[5] <- "pval_time_0.25_1"
msdiff_5 <- msdiff_5 %>% select(c(1,4:5))

names(msdiff_6)[1] <- "unique_name"
names(msdiff_6)[4] <- "rsq_time_0.25_2"
names(msdiff_6)[5] <- "pval_time_0.25_2"
msdiff_6 <- msdiff_6 %>% select(c(1,4:5))

names(msdiff_7)[1] <- "unique_name"
names(msdiff_7)[4] <- "rsq_time_0.25_4"
names(msdiff_7)[5] <- "pval_time_0.25_4"
msdiff_7 <- msdiff_7 %>% select(c(1,4:5))

names(msdiff_8)[1] <- "unique_name"
names(msdiff_8)[4] <- "rsq_time_1_2"
names(msdiff_8)[5] <- "pval_time_1_2"
msdiff_8 <- msdiff_8 %>% select(c(1,4:5))

names(msdiff_9)[1] <- "unique_name"
names(msdiff_9)[4] <- "rsq_time_1_4"
names(msdiff_9)[5] <- "pval_time_1_4"
msdiff_9 <- msdiff_9 %>% select(c(1,4:5))

names(msdiff_10)[1] <- "unique_name"
names(msdiff_10)[4] <- "rsq_time_2_4"
names(msdiff_10)[5] <- "pval_time_2_4"
msdiff_10 <- msdiff_10 %>% select(c(1,4:5))

sshct_sig <- merge(msdiff_0.25, msdiff_1, by = "unique_name", all = TRUE)
sshct_sig <- merge(sshct_sig, msdiff_2, by = "unique_name", all = TRUE)
sshct_sig <- merge(sshct_sig, msdiff_4, by = "unique_name", all = TRUE)
sshct_sig <- merge(sshct_sig, msdiff_5, by = "unique_name", all = TRUE)
sshct_sig <- merge(sshct_sig, msdiff_6, by = "unique_name", all = TRUE)
sshct_sig <- merge(sshct_sig, msdiff_7, by = "unique_name", all = TRUE)
sshct_sig <- merge(sshct_sig, msdiff_8, by = "unique_name", all = TRUE)
sshct_sig <- merge(sshct_sig, msdiff_9, by = "unique_name", all = TRUE)
sshct_sig <- merge(sshct_sig, msdiff_10, by = "unique_name", all = TRUE)

# write.table(sshct_sig, "./Serum_Stimulated_HCT116/Cleaned/20190318_sshct_sig_tomergewith_peptidedata_alltimecompare.csv", sep = ",", row.names = FALSE)

# labeling the replicates
ms$rep <- ms$sample_id
ms$sample_id <- NA
ms$sample_id[which(ms$rep == 1 & ms$time_point == 0)] <- 1
ms$sample_id[which(ms$rep == 2 & ms$time_point == 0)] <- 2
ms$sample_id[which(ms$rep == 3 & ms$time_point == 0)] <- 3
ms$sample_id[which(ms$rep == 1 & ms$time_point == 0.25)] <- 4
ms$sample_id[which(ms$rep == 2 & ms$time_point == 0.25)] <- 5
ms$sample_id[which(ms$rep == 3 & ms$time_point == 0.25)] <- 6
ms$sample_id[which(ms$rep == 1 & ms$time_point == 1)] <- 7
ms$sample_id[which(ms$rep == 2 & ms$time_point == 1)] <- 8
ms$sample_id[which(ms$rep == 3 & ms$time_point == 1)] <- 9
ms$sample_id[which(ms$rep == 1 & ms$time_point == 2)] <- 10
ms$sample_id[which(ms$rep == 2 & ms$time_point == 2)] <- 11
ms$sample_id[which(ms$rep == 3 & ms$time_point == 2)] <- 12
ms$sample_id[which(ms$rep == 1 & ms$time_point == 4)] <- 13
ms$sample_id[which(ms$rep == 2 & ms$time_point == 4)] <- 14
ms$sample_id[which(ms$rep == 3 & ms$time_point == 4)] <- 15

# Rearranging columns and deleting sample_id
ms <- ms %>% 
  select(time_point, rep, everything()) %>% 
  select(-sample_id) %>% 
  arrange(PG.ProteinGroups, k_site, time_point, rep)

# Cleaning modified sequence
ms$EG.ModifiedPeptide <- ms$EG.ModifiedPeptide %>% 
  str_replace_all("(_)$", "") %>% 
  str_replace_all("^_", "")

# Summarizing stoich data
ms_sum <- ms %>% 
  group_by(PG.ProteinGroups, PG.ProteinDescriptions, PG.Genes, 
           PEP.StrippedSequence, EG.ModifiedPeptide, k_site, time_point) %>% 
  summarise(stoich_mean = mean(stoich_corr, na.rm = TRUE),
            stoich_sd = sd(stoich_corr, na.rm = TRUE),
            stoich_n = n()) %>% 
  ungroup() %>% 
  rename("time" = "time_point") %>% 
  arrange(PG.ProteinGroups, k_site, time)

ms_sum$unique_name <- paste(ms_sum$PG.ProteinGroups, ms_sum$k_site, ms_sum$EG.ModifiedPeptide, sep = "_")

ms_sum <- merge(ms_sum, sshct_sig, all.x = TRUE, by = "unique_name")

ms_wide <- ms_sum %>% 
  select(-stoich_sd, -stoich_n) %>% 
  spread(time, stoich_mean)

# Counting time point observations
ms_wide$tp_count <- apply(ms_wide[16:19], 1, function(x){
  length(which(!is.na(x)))
})

ms_wide_mean <- ms_sum %>% 
  select(PG.ProteinGroups, PG.Genes, k_site, PG.ProteinDescriptions, time,
         PEP.StrippedSequence, EG.ModifiedPeptide, stoich_mean, pval_time_0_0.25, pval_time_0_1,
         pval_time_0_2, pval_time_0_4, pval_time_0.25_1, pval_time_0.25_2, pval_time_0.25_4,
         pval_time_1_2, pval_time_1_4, pval_time_2_4) %>% 
  filter(!is.na(k_site)) %>% 
  spread(time, stoich_mean) %>% 
  rename("Hr0_mean" = "0",
         "Hr0.25_mean" = "0.25",
         "Hr1_mean" = "1",
         "Hr2_mean" = "2",
         "Hr4_mean" = "4")

ms_wide_sd <- ms_sum %>% 
  select(PG.ProteinGroups, PG.Genes, k_site, PG.ProteinDescriptions, time,
         PEP.StrippedSequence, EG.ModifiedPeptide, stoich_sd, pval_time_0_0.25, pval_time_0_1,
         pval_time_0_2, pval_time_0_4, pval_time_0.25_1, pval_time_0.25_2, pval_time_0.25_4,
         pval_time_1_2, pval_time_1_4, pval_time_2_4) %>% 
  filter(!is.na(k_site)) %>% 
  spread(time, stoich_sd) %>% 
  rename("Hr0_sd" = "0",
         "Hr0.25_sd" = "0.25",
         "Hr1_sd" = "1",
         "Hr2_sd" = "2",
         "Hr4_sd" = "4")

ms_wide_n <- ms_sum %>% 
  select(PG.ProteinGroups, PG.Genes, k_site, PG.ProteinDescriptions, time,
         PEP.StrippedSequence, EG.ModifiedPeptide, stoich_n, pval_time_0_0.25, pval_time_0_1,
         pval_time_0_2, pval_time_0_4, pval_time_0.25_1, pval_time_0.25_2, pval_time_0.25_4,
         pval_time_1_2, pval_time_1_4, pval_time_2_4) %>% 
  filter(!is.na(k_site)) %>% 
  spread(time, stoich_n) %>% 
  rename("Hr0_n" = "0",
         "Hr0.25_n" = "0.25",
         "Hr1_n" = "1",
         "Hr2_n" = "2",
         "Hr4_n" = "4")

ms_wide <- merge(ms_wide_mean, ms_wide_sd, all = TRUE)
ms_wide <- merge(ms_wide, ms_wide_n, all = TRUE)

ms_wide$change_0.25 <- ms_wide$Hr0.25_mean - ms_wide$Hr0_mean
ms_wide$change_1 <- ms_wide$Hr1_mean - ms_wide$Hr0_mean
ms_wide$change_2 <- ms_wide$Hr2_mean - ms_wide$Hr0_mean
ms_wide$change_4 <- ms_wide$Hr4_mean - ms_wide$Hr0_mean

rm(ms_wide_mean, ms_wide_sd, ms_wide_n);gc()

#write.table(ms_wide, "./Serum_Stimulated_HCT116/Cleaned/20190305_HCT116_stoich_wide_stats_sigalltimepts.csv", sep = ",", row.names = FALSE)


