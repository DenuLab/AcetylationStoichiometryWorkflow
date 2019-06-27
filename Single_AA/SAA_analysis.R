
# Libraries ---------------------------------------------------------------


options(stringsAsFactors = FALSE)
library(tidyverse)
library(ggthemes)
library(rio)


# Data Import -------------------------------------------------------------


data <- import("./Single_AA/Jing_data.csv")
data$compartment <- factor(data$compartment, 
                           levels = c("cytosol", "mitochondria", "histone", "Nuclear nonhistone"))

ack <- import("./Single_AA/acetyl_lys_M1.csv")
ack$compartment <- factor(ack$compartment, 
                          levels = c("cytosol", "mitochondria", "histone", "Nuclear nonhistone"))


# Plots -------------------------------------------------------------------


# Acetyl lysine plot
### Figure  4D
ggplot(data %>% filter(mod == "acetyl"), aes(x = compartment, y = mean)) +
  geom_bar(stat = "identity") + 
  geom_errorbar(aes(ymin = mean - std_dev, ymax = mean + std_dev), width = 0.1) +
  labs(title = "Acetyl",
       y = "Modified Lysine / Lysine") +
  theme_pander(base_size = 14)
# ggsave(filename = "figures/Acetyl_bar2.pdf", height = 2, width = 3.5)

# Me1 plot
ggplot(data %>% filter(mod == "Me1"), aes(x = compartment, y = mean)) +
  geom_bar(stat = "identity") + 
  geom_errorbar(aes(ymin = mean - std_dev, ymax = mean + std_dev), width = 0.1) +
  labs(title = "Monomethyl",
         y = "Modified Lysine / Lysine") +
  theme_pander(base_size = 14)
# ggsave(filename = "figures/Me1_bar.png")


# Me2 plot
ggplot(data %>% filter(mod == "Me2"), aes(x = compartment, y = mean)) +
  geom_bar(stat = "identity") + 
  geom_errorbar(aes(ymin = mean - std_dev, ymax = mean + std_dev), width = 0.1) +
  labs(title = "Dimethyl",
       y = "Modified Lysine / Lysine") +
  theme_pander(base_size = 14)
# ggsave(filename = "figures/Me2_bar.png")


# Me3 plot
ggplot(data %>% filter(mod == "Me3"), aes(x = compartment, y = mean)) +
  geom_bar(stat = "identity") + 
  geom_errorbar(aes(ymin = mean - std_dev, ymax = mean + std_dev), width = 0.1) +
  labs(title = "Trimethyl",
       y = "Modified Lysine / Lysine") +
  theme_pander(base_size = 14)
# ggsave(filename = "figures/Me3_bar.png")


