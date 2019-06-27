
# Libraries ---------------------------------------------------------------

options(stringsAsFactors = FALSE)
library(tidyverse)
library(dplyr)
library(ggridges)
library(broom)
library(purrr)
library(viridisLite)


# Importing data ----------------------------------------------------------

frg_sum <- read.csv("./HSC_StoichCurve/Cleaned/20181113_HSC_HPRP_tidy_frg_stoich_output_REPROCESSED.csv")

# Formatting data ---------------------------------------------------------

# Filtering H_corr peaks
frg_sum <- frg_sum %>% 
  select(-sample_id, k_count) %>% 
  filter(H_corr_sum > 0,
         F.FrgZ == 1)


# functions ---------------------------------------------------------------

# Function to count the number of time points
stoich_input_count <- function(x){
  length(unique(x$stoich_input))
}

# linear model function for the nested dataframe
lm_stoich <- function(x){
  lm(x$stoich_corr ~ x$stoich_input, 
     data = x,
     na.action = na.exclude)
}


# Nesting data ------------------------------------------------------------

data <- frg_sum %>%
  group_by(PG.ProteinGroups, PG.ProteinDescriptions, PG.Genes,
           EG.StrippedSequence, EG.ModifiedSequence, k_site) %>%
  nest() %>%
  mutate(n_frg = map_dbl(data, nrow),
         n_stoich_input = map_dbl(data, stoich_input_count)) %>%
  filter(n_stoich_input >= 9)


# Running the linear model
data <- data %>%
  mutate(lm_data = map(data, lm_stoich)) %>% 
  arrange(n_stoich_input)


# Summary Statistics ------------------------------------------------------


# Unnesting data
stoich_models <- data %>%
  mutate(tidy_data = map(lm_data, tidy),
         glance_data = map(lm_data, glance)) %>%
  unnest(glance_data) %>%
  unnest(tidy_data) %>%
  arrange(-n_stoich_input, PG.ProteinGroups)


# Formatting data ---------------------------------------------------------


# Extracting slope and intercept
stoich_tidy_glance <- stoich_models %>%
  select(1:8, adj.r.squared, p.value, term, estimate) %>%
  spread(term, estimate) %>%
  rename("intercept" = "(Intercept)",
         "slope" = "x$stoich_input") %>%
  arrange(-n_stoich_input, PG.ProteinGroups)

# Extracting the standard error
temp_df <- stoich_models %>%
  select(1:8, term, std.error) %>%
  spread(term, std.error) %>%
  rename("intercept_std.error" = "(Intercept)",
         "slope_std.error" = "x$stoich_input") %>%
  arrange(-n_stoich_input, PG.ProteinGroups)

# combining datasets
stoich_tidy_glance <- merge(stoich_tidy_glance, temp_df)
stoich_tidy_glance <- stoich_tidy_glance %>%
  arrange(-n_stoich_input, PG.ProteinGroups)
rm(temp_df);gc()

### data for high confidence scatter plots
data_facet <- stoich_models %>% 
  filter(n_stoich_input == 11) %>% 
  arrange(p.value) %>% 
  dplyr::slice(1:50)

data_facet <- frg_sum[which(frg_sum$EG.ModifiedSequence %in% data_facet$EG.ModifiedSequence),]


# model summary -----------------------------------------------------------


model_tidy <- data %>% 
  mutate(tidy_data = map(lm_data, tidy)) %>%
  unnest(tidy_data) %>%
  arrange(-n_stoich_input, PG.ProteinGroups)

model_glance <- data %>% 
  mutate(glance_data = map(lm_data, glance)) %>%
  unnest(glance_data) %>%
  select(-data, -lm_data) %>% 
  arrange(-n_stoich_input, PG.ProteinGroups)

model_augment <- data %>% 
  mutate(augment_data = map(lm_data, augment)) %>% 
  unnest(augment_data) %>% 
  arrange(-n_stoich_input, PG.ProteinGroups)


# Plots -------------------------------------------------------------------
stoich_tidy_glance2 <- stoich_tidy_glance %>% arrange(adj.r.squared)

# Linear regression analysis - fraction data
### Figure 2E
ggplot(data = data.frame(x = seq(0,100, length.out = 10), 
                         y = seq(0,1, length.out = 10)), 
       aes(x = x, y = y)) + 
  geom_point(color = "white") +
  geom_abline(data = stoich_tidy_glance2, lwd = 0.4,
              aes(slope = slope, intercept = intercept, color = adj.r.squared)) +
  scale_color_distiller(palette = "RdBu") +
  # scale_color_gradientn(colors = viridis(100, direction = -1, option = "D")) +
  labs(title = "Linear regression analysis",
       y = "Corrected Stoichiometry",
       x = "Input Stoichiometry",
       color = expression(Adj.R^2)) +
  theme_light(base_size = 14) +
  scale_y_continuous(limits = c(0,1))
#ggsave(filename = "figures/stoich_curve_linear_regression_rdbl4_again.pdf", width = 7.5, height = 5)



# Distribution of R2 value
### Figure 2F
ggplot(stoich_tidy_glance) +
  geom_density(aes(x = adj.r.squared)) +#, color = factor(n_stoich_input))) +
  theme_light(base_size = 14) +
  labs(title = "Global R2 distribution",
       x = expression(Adjusted~R^2),
       color = "No. of \ndata points") +
  annotate("text", x = 0.8, y = 9, label = paste("median adj.R2 is ",round(median(
    stoich_tidy_glance$adj.r.squared), digits = 3)))
#ggsave(filename = "figures/Adj_r_squared_density.pdf")

ggplot(stoich_tidy_glance) +
  geom_boxplot(aes(x = factor(n_stoich_input), y = adj.r.squared)) +
  theme_light(base_size = 14) +
  labs(y = expression(Adjusted~R~squared),
       x = "Stoichiometry input")
#ggsave(filename = "figures/Adj_r_squared_boxplot.pdf", height = 5, width = 7)

# fraction data scatter plot
ggplot(data_facet, aes(x = stoich_input, y = stoich_corr)) +
  geom_jitter(width = 3, size = 1) +
  geom_smooth(method = "lm", fullrange = TRUE) +
  facet_wrap(~PG.ProteinGroups) +
  #scale_y_continuous(limits = c(0:1)) +
  theme_light(base_size = 14) +
  labs(y = "Corrected Stoichiometry",
       x = "Input Stoichiometry")
#ggsave(filename = "figures/Top_proteins.pdf")

# Residual analysis - histogram
ggplot(model_augment) +
  geom_histogram(aes(x = .resid), bins = 75, color = "black", fill = "grey") +
  theme_light(base_size = 14) +
  labs(#title = "Residual analysis",
    x = "Residuals")
#ggsave(filename = "figures/stoich_curve_residual_analysis.pdf", height = 5, width = 7)

