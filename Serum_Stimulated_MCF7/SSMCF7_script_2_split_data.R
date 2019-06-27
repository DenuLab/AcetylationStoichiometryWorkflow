options(stringsAsFactors = FALSE)
library(tidyverse)

### Data Import

#Spectronaut output
data <- read.csv("./Serum_Stimulated_MCF7/Raw/20181031_SSMCF7_stoich_Spectronaut10_reprocess_correctFASTA_Report.csv")

# Splitting the data into two halves to more quickly process through the data
data1 <- data %>% subset(time_point == 0 | time_point == 1)
data2 <- data %>% subset(time_point == 2 | time_point == 4)

# write.table(data1, "./Serum_Stimulated_MCF7/Raw/20181031_SSMCF7_stoich_Spectronaut10_reprocess_correctFASTA_1sthalf.csv", sep = ",", row.names = FALSE)
# write.table(data2, "./Serum_Stimulated_MCF7/Raw/20181031_SSMCF7_stoich_Spectronaut10_reprocess_correctFASTA_2ndhalf.csv", sep = ",", row.names = FALSE)
