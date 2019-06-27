
# Loading libraries -------------------------------------------------------

options(stringsAsFactors = FALSE)
library(tidyverse)
library(ggthemes)

# Importing data ----------------------------------------------------------


peptides <- read.csv("./NCE_optimization/MaxQuant_NCE_peptides.csv")

# Formatting data ---------------------------------------------------------


names(peptides)[54] <- 15
names(peptides)[55] <- 20
names(peptides)[56] <- 25
names(peptides)[57] <- 30
names(peptides)[58] <- 35
names(peptides)[59] <- 40
names(peptides)[60] <- 45
names(peptides)[61] <- 50

# Removing the contaminants and decoy observations
peptides$contaminants <- FALSE
peptides$contaminants[grep("CON_", peptides$Leading.razor.protein)] <- TRUE
peptides$contaminants[grep("REV_", peptides$Leading.razor.protein)] <- TRUE


# formatting the data
pep_short <- peptides %>% 
  filter(contaminants == FALSE) %>% 
  select(c(Sequence, Amino.acid.before, Last.amino.acid, `15`:`50`)) %>% 
  gather(NCE, count, 4:11) %>% 
  filter(!is.na(count))


# Plots -------------------------------------------------------------------


# Bar graph
ggplot(pep_short[which(pep_short$Last.amino.acid == "K" | 
                         pep_short$Last.amino.acid == "R" | 
                         pep_short$Last.amino.acid == "E"), ]) + 
  geom_bar(aes(x = NCE, fill = Last.amino.acid), position = "stack") + 
  labs(title = "NCE Optimization", y = "count", x = "NCE", fill = "C-terminal\namino acid") + 
  scale_fill_brewer(palette = "Set1") + 
  theme_light(base_size = 15)
#ggsave(filename = "figures/NCE_optimization_bar_graph.png")

# c-terminal amino acid bar graph
ggplot(pep_short %>% filter(NCE == 25)) +
  geom_bar(aes(x = Last.amino.acid)) +
  ggthemes::theme_pander(base_size = 15) +
  labs(x = "C-terminal amino acid",
       title = "Acetyl-peptidome cleavage by \nTrypsin & Glu-C (NCE-25)")
#ggsave(filename = "figures/c_term_aa_bargraph.png")

