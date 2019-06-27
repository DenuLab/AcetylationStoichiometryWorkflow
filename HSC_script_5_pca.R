
# Libraries ---------------------------------------------------------------


options(stringsAsFactors = FALSE)
library(tidyverse)
library(svglite)
#install.packages("factoextra")
library(factoextra)


# Importing data ----------------------------------------------------------


frg_sum <- read.csv("./HSC_StoichCurve/Cleaned/20181113_HSC_HPRP_tidy_frg_stoich_output_REPROCESSED.csv")


# Formatting data ---------------------------------------------------------


# Filtering data
frg_sum <- frg_sum %>% 
  select(-sample_id) %>% 
  filter(H_corr_sum > 0, F.FrgZ == 1)

data_pca <- frg_sum %>% 
  select(stoich_input, PG.ProteinGroups, k_site, EG.ModifiedSequence, 
         F.FrgIon, F.FrgLossType, stoich_corr) %>% 
  unite(frg_id, PG.ProteinGroups, EG.ModifiedSequence, 
        k_site, F.FrgIon, F.FrgLossType, sep = "_") %>% 
  # gather(type, value, 3) %>% 
  spread(frg_id, stoich_corr)


# Removing columns with NAs
colname_na <- colnames(data_pca[2:ncol(data_pca)])[apply(data_pca[2:ncol(data_pca)], 2, anyNA)]
data_pca[colname_na] <- NULL


# PCA ---------------------------------------------------------------------

# Performing the PCA
pca <- prcomp(data_pca[2:ncol(data_pca)], center = T, scale. = T)


# # Generic plotting
plot(pca, type = "l")
summary(pca)


pca_plot <- data.frame(pca$x, stoich = data_pca$stoich_input)


# Plots -------------------------------------------------------------------


# PCA plot
### Figure 2C
ggplot(pca_plot, aes(PC1, PC2, color = factor(stoich))) +
  geom_point(size = 5) +
  geom_point(size = 5, shape = 1, color = "black") +
  theme_light(base_size = 12) +
  scale_color_brewer(palette = "Spectral") +
  guides(color = FALSE) +
  expand_limits(x = c(-12,10), y = c(-2,4)) +
  labs(#title = "Principle Component Analysis", 
    colour = "Stoich\nInput") +
  geom_text(aes(label = stoich), color = "black", hjust = 1.25, vjust = -1.75) +
  theme(text = element_text(size=15))
#ggsave(filename = "figures/pca_stoichiometry_curve2.pdf", height = 4, width = 4, useDingbats = FALSE)



# Generic plotting
pca_sum <- summary(pca)
pca_df <- as.data.frame(pca_sum$importance)
pca_df <- pca_df %>% 
  mutate(condition = rownames(.)) %>% 
  gather(temp, value, 1:(ncol(.)-1)) 
pca_df$number <- as.numeric(substring(pca_df$temp, 
                                      first = 3, 
                                      last = nchar(pca_df$temp)))

# Proportion of variance for PCA
ggplot(pca_df %>% filter(condition == "Proportion of Variance")) +
  geom_point(aes(x = reorder(temp, number), y = value, group = condition)) +
  geom_line(aes(x = reorder(temp, number), y = value, group = condition)) +
  ggthemes::theme_pander(base_size = 14) +
  labs(title = "Principle Components Analysis",
       subtitle = "Proportion of Variance",
       y = "Variance",
       x = "Principle Component") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#ggsave(filename = "figures/PCA_proportion_of_variance.pdf")

# Proportion of variance for PCA
ggplot(pca_df %>% filter(condition == "Cumulative Proportion")) +
  geom_point(aes(x = reorder(temp, number), y = value, group = condition)) +
  geom_line(aes(x = reorder(temp, number), y = value, group = condition)) +
  ggthemes::theme_pander(base_size = 14) +
  labs(title = "Principle Components Analysis",
       subtitle = "Cumulative Proportion of Variance",
       y = "Cumulative Variance",
       x = "Principle Component") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#ggsave(filename = "figures/PCA_cumulative_proportion.pdf")

