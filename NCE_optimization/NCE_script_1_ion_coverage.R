
# Loading libraries -------------------------------------------------------


options(stringsAsFactors = FALSE)
library(tidyverse)
library(stringr)


# Importing data ----------------------------------------------------------


data <- read.csv("./NCE_optimization/20181214_MQ_NCE_msms.csv")


# Formatting data ---------------------------------------------------------

# Selecting top protein ID
data$Proteins <- unlist(lapply(data$Proteins, function(x){
  unlist(strsplit(x, split = ";"))[1]
}))

# Filtering
data$contaminant <- FALSE
data$contaminant[grep("CON_", data$Proteins)] <- TRUE
data$contaminant[grep("REV_", data$Proteins)] <- TRUE
data$contaminant[grep("Biognosys", data$Proteins)] <- TRUE

# Formatting the data
data <- data %>% 
  filter(contaminant == FALSE) %>% 
  select(Raw.file, Sequence, Modifications, Proteins, Matches, contaminant) %>% 
  mutate(NCE = str_sub(Raw.file, 26, 27),
         pep_length = str_length(Sequence),
         last_aa = str_sub(Sequence, -1, -1),
         acetyl = str_detect(Modifications, "Acetyl"),
         k_count = str_count(Sequence, "K"),
         y_total = str_count(Matches, "y"),
         b_total = str_count(Matches, "b"),
         y_ion = y_total / ((pep_length - 1) * 3),
         b_ion = b_total / ((pep_length - 1) * 3))

# Long format of data
data_long <- data %>% gather(ion_type, coverage, 14:15)


# Plots -------------------------------------------------------------------


# Acetyl peptide boxplot
## Figure 3B
ggplot(data_long %>% filter(acetyl == TRUE)) + 
  geom_boxplot(aes(x = NCE, y = coverage, fill = ion_type)) + 
  theme_light(base_size = 15) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "NCE Optimization",
       subtitle = "All acetyl peptides",
       y = "Ion Coverage\n(observed / total possible)",
       x = "NCE",
       fill = "Ion type") +
  scale_y_continuous(limits = c(0,1))
#ggsave(filename = "figures/NCE_opt_acetyl_peptides.png")

