
# Libraries ---------------------------------------------------------------

options(stringsAsFactors = FALSE)
library(tidyverse)
library(ggridges)


# Importing data ----------------------------------------------------------

frg_sum <- read.csv("./HSC_StoichCurve/Cleaned/20181113_HSC_HPRP_tidy_frg_stoich_output_REPROCESSED.csv")
pep_sum <- read.csv("./HSC_StoichCurve/Cleaned/20181113_HSC_HPRP_tidy_peptide_stoich_output_REPROCESSED.csv")


# Formatting data ---------------------------------------------------------

# Filtering H_corr peaks
frg_sum <- frg_sum %>% 
  filter(H_corr_sum > 0)

pep_sum <- pep_sum %>% 
  filter(H_corr_sum > 0)

# Boxplot data - fragment ions
x <- frg_sum %>% 
  select(1:15, stoich_obs) %>% 
  mutate(corrected = FALSE) %>% 
  rename("stoich" = "stoich_obs")

y <- frg_sum %>% 
  select(1:15, stoich_corr) %>% 
  mutate(corrected = TRUE) %>% 
  rename("stoich" = "stoich_corr")
boxplot_frg <- rbind(x, y)
rm(x,y);gc()

# Boxplot data - peptides
x <- pep_sum %>% 
  select(1:9, stoich_obs) %>% 
  mutate(corrected = FALSE) %>% 
  rename("stoich" = "stoich_obs")

y <- pep_sum %>% 
  select(1:9, stoich_corr) %>% 
  mutate(corrected = TRUE) %>% 
  rename("stoich" = "stoich_corr")
boxplot_pep <- rbind(x, y)
rm(x,y);gc()


# Plots -------------------------------------------------------------------


# Corrected vs Un-corrected Boxplot
### Figure 2B
ggplot(boxplot_frg) + 
  geom_boxplot(aes(x = factor(stoich_input), y = stoich, fill = corrected), lwd = 0.2, outlier.color = NA) + 
  labs(title = "Fragment ion stoichiometry curve", 
       x = "Input stoichiometry", 
       y = "Corrected Stoichiometry") +
  scale_y_continuous(limits = c(0, 1.2)) +
  scale_fill_brewer(palette = "Paired") +
  theme_light(base_size = 14)
#ggsave(filename = "figures/Frg_stoich_curve_boxplot_square_20180731.pdf", height = 5, width = 6)

# Joy Log10(L/H) plot
### Figure 2D
ggplot(frg_sum) + 
  geom_density_ridges(aes(y = as.factor(stoich_input), x = log10(L_sum / H_corr_sum), 
                          fill = factor(stoich_input)), scale = 3) +
  labs(x = expression(Log[10]~Peak~Area-(Light/Heavy)),
       y = "Stoichiometry Input") +
  scale_x_continuous(limits = c(-4,4)) +
  theme_light(base_size = 14) +
  guides(fill = FALSE) +
  scale_fill_brewer(palette = "Spectral")
#ggsave(filename = "figures/Frg_stoich_curve_ridgeplot_20180731.pdf", height = 8, width = 4)

# Comparison of Q.values
ggplot(frg_sum) +
  geom_density(aes(x = log10(EG.Qvalue_L_rms), color = "L")) +
  geom_density(aes(x = log10(EG.Qvalue_H_rms), color = "H")) +
  facet_wrap(~stoich_input, scales = "free") +
  scale_color_brewer(palette = "Set1")

# Log2 Peak Areas - Fragment ions
ggplot(frg_sum) +
  geom_density(aes(x = log2(L_sum), color = "L")) +
  geom_density(aes(x = log2(H_sum), color = "H")) +
  geom_density(aes(x = log2(H_corr_sum), color = "Hc")) +
  scale_x_continuous(limits = c(0,30)) +
  scale_color_brewer(palette = "Set1") +
  facet_wrap(~stoich_input) +
  theme_light(base_size = 14)
#ggsave(filename = "figures/Frg_peak_area_distributions.pdf", height = 5, width = 6.5)

# Density plot of fragment stoichiometry
ggplot(frg_sum) +
  geom_density(aes(x = stoich_corr, color = "TRUE")) +
  geom_density(aes(x = stoich_obs, color = "FALSE")) +
  scale_color_brewer(palette = "Set1") +
  facet_wrap(~stoich_input, scales = "free") +
  theme_light(base_size = 12) +
  labs(#title = "Stoichiometry Distributions",
    x = "Stoichiometry",
    color = "Corrected") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#ggsave(filename = "figures/Frg_stoich_distribution.pdf", height = 5, width = 6.5)


# Peptide -----------------------------------------------------------------

# Corrected vs Un-corrected Boxplot
ggplot(boxplot_pep) +
  geom_boxplot(aes(x = factor(stoich_input), y = stoich, fill = corrected), outlier.color = NA) +
  labs(title = "Peptide stoichiometry curve",
       x = "Input stoichiometry",
       y = "Corrected Stoichiometry") +
  scale_y_continuous(limits = c(0, 1.2)) +
  scale_fill_brewer(palette = "Set1") +
  theme_light(base_size = 14)

# Joy Log10(L/H) plot
ggplot(pep_sum) +
  geom_density_ridges(aes(y = as.factor(stoich_input), x = log10(L_sum / H_corr_sum),
                          fill = factor(stoich_input)), scale = 3) +
  labs(x = expression(Log[10]~Peak~Area-(Light/Heavy)),
       y = "Stoichiometry Input") +
  scale_x_continuous(limits = c(-4,4)) +
  theme_light(base_size = 14) +
  guides(fill = FALSE) +
  scale_fill_brewer(palette = "Spectral")

# Comparison of Q.values
ggplot(pep_sum) +
  geom_density(aes(x = log10(EG.Qvalue_L_rms), color = "L")) +
  geom_density(aes(x = log10(EG.Qvalue_H_rms), color = "H")) +
  facet_wrap(~stoich_input, scales = "free")

# Log2 Peak Areas - Peptide
ggplot(pep_sum) +
  geom_density(aes(x = log2(L_sum), color = "L")) +
  geom_density(aes(x = log2(H_sum), color = "H")) +
  geom_density(aes(x = log2(H_corr_sum), color = "Hc")) +
  scale_x_continuous(limits = c(0,30)) +
  scale_color_brewer(palette = "Set1") +
  facet_wrap(~stoich_input) +
  theme_light(base_size = 14)

# Density plot of fragment stoichiometry
ggplot(pep_sum) +
  geom_density(aes(x = stoich_corr, color = "TRUE")) +
  geom_density(aes(x = stoich_obs, color = "FALSE")) +
  scale_color_brewer(palette = "Set1") +
  facet_wrap(~stoich_input, scales = "free") +
  labs(title = "Stoichiometry Distributions",
       x = "Stoichiometry",
       color = "Corrected")

