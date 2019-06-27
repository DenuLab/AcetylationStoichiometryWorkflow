# Loading Libraries -------------------------------------------------------

options(stringsAsFactors = FALSE)
library(tidyverse)
library(stringr)

# Reading & Formatting Data -----------------------------------------------


# Reading in data
spec_lib <- read.csv("./Spectral_Libraries/20181029_MCF7_Library_Remake_correctFASTA_inflated.csv")

# Selecting and arranging columns
spec_lib <- spec_lib %>% 
  select(PrecursorCharge, IntModifiedPeptide, ModifiedPeptide, StrippedPeptide, LabeledPeptide,
         PrecursorMz, FragmentLossType, FragmentNumber, FragmentType, FragmentCharge, FragmentMz,
         UniProtIds, ProteinDescription, Genes) %>% 
  arrange(StrippedPeptide, FragmentMz)

spec_lib$k_count <- str_count(spec_lib$StrippedPeptide, pattern = "K")
spec_lib$acetyl_count <- str_count(spec_lib$ModifiedPeptide, pattern = "Acetyl")


# Generating the fragment Ion Sequence
spec_lib$fragment_ion_sequence <- substring(spec_lib$StrippedPeptide, 
                                            ifelse(spec_lib$FragmentType == "b", 1, (nchar(spec_lib$StrippedPeptide) - (spec_lib$FragmentNumber - 1))), 
                                            ifelse(spec_lib$FragmentType == "y", nchar(spec_lib$StrippedPeptide), spec_lib$FragmentNumber))

# Counting lysines on fragment ion sequence
spec_lib$frag_k_count <- str_count(spec_lib$fragment_ion_sequence, "K")

# Adding Isotope Label ----------------------------------------------------


# Creating a fragment ion dataframe
frag_1k <- spec_lib %>% 
  select(-LabeledPeptide, -PrecursorMz) %>% 
  distinct(PrecursorCharge, StrippedPeptide, ModifiedPeptide, fragment_ion_sequence, FragmentMz, .keep_all = TRUE) %>% 
  arrange(PrecursorCharge, StrippedPeptide, ModifiedPeptide,
          fragment_ion_sequence, FragmentLossType, FragmentMz) %>% 
  filter(acetyl_count > 0, frag_k_count == 1)

# Adding the heavy or light isotope label
frag_1k$isotope_label <- NA

ptm <- proc.time()

for(i in seq_along(frag_1k$isotope_label)){
  
  # Conditional statement
  if(is.na(frag_1k$isotope_label[i]) &
     frag_1k$StrippedPeptide[i] == frag_1k$StrippedPeptide[i+1] &
     frag_1k$ModifiedPeptide[i] == frag_1k$ModifiedPeptide[i+1] &
     frag_1k$PrecursorCharge[i] == frag_1k$PrecursorCharge[i+1] &
     frag_1k$fragment_ion_sequence[i] == frag_1k$fragment_ion_sequence[i+1] &
     frag_1k$FragmentLossType[i] == frag_1k$FragmentLossType[i+1]){
    
    # Do this if all conditions are met
    frag_1k$isotope_label[i] <- "L"
    frag_1k$isotope_label[i+1] <- "H"
  }
}

proc.time() - ptm
# Estimated: 60s

# Adding isotope label to spec_lib df
spec_lib <- merge(spec_lib, frag_1k, all = TRUE)


spec_lib$isotope_label[which(is.na(spec_lib$isotope_label))] <- "L"


# Renaming columns to match Spectronaut output
spec_lib_isotope <- spec_lib %>%
  rename("FG.Charge" = "PrecursorCharge",
         "EG.ModifiedPeptide" = "IntModifiedPeptide",
         "EG.StrippedSequence" = "StrippedPeptide",
         "F.FrgLossType" = "FragmentLossType",
         "F.FrgNum" = "FragmentNumber",
         "F.FrgType" = "FragmentType",
         "F.Charge" = "FragmentCharge",
         "FG.PrecMz" = "PrecursorMz",
         "F.FrgMz" = "FragmentMz")


# Writing tables ----------------------------------------------------------

#write.table(spec_lib_isotope, file = "./Spectral_Libraries/20181031_tidy_mcf7_spec_lib_spectronaut10_correctFASTA_isotope.csv", sep = ",", row.names = FALSE)
