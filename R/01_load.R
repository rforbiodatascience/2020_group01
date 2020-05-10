# Clear workspace
# ------------------------------------------------------------------------------
rm(list = ls())

# Load libraries
# ------------------------------------------------------------------------------
library("tidyverse")
library(readxl)

# Define functions
# ------------------------------------------------------------------------------
#file.choose()
source(file = "R/99_project_functions.R")

# Load data
# ------------------------------------------------------------------------------
# BARRACODA FILES
#-------------------------

path_ct26 <- "data/raw/barracoda_output_CT26.xlsx"

all_ct26_barracoda_raw <- path_ct26 %>% 
  # function to import all sheets
  excel_sheets() %>% 
  # give names to each sheet
  set_names() %>% 
  # apply read_excel to each sheet, and add the number to the colum sheet
  map_df(~ read_excel(path = path_ct26, sheet = .x), .id = "sheet") 

path_4t1 <- "data/raw/barracoda_output_4T1.xlsx"

all_4t1_barracoda_raw <- path_4t1 %>% 
  excel_sheets() %>% 
  set_names() %>% 
  map_df(~ read_excel(path = path_4t1, sheet = .x), .id = "sheet")

# MUPEXI FILES
#-------------------------
mupexi_ct26 <- read_xlsx(path = "data/raw/ct26_library_mupexi.xlsx") %>% 
  # remove extra columns from previous handling
  select(-identifier, -Mut_peptide.y, -Allele) %>% 
  # convert Mut_MHCrank_EL and Expression level to numeric so we can join both files
  mutate(Mut_MHCrank_EL = as.numeric(Mut_MHCrank_EL),
         Expression_Level = as.numeric(Expression_Level))

mupexi_4t1 <- read_xlsx(path = "data/raw/4T1_library_mupexi.xlsx") %>% 
  select(-identifier, -"...1")
#-------------------------

# PE_population_info and buffycoat_HLA_info file
#-------------------------
sample_info <- read_xlsx(path = "data/raw/sample_info.xlsx")
buffycoat.HLA_info <- read_xlsx(path = "data/raw/buffycoatHLA_info.xlsx")
#-------------------------

# Wrangle data
# ------------------------------------------------------------------------------
all_barracoda <- full_join(all_ct26_barracoda_raw, all_4t1_barracoda_raw)

# Wrong HLA annotation: H-2XX instead of H2-XX as mupexi. 
all_barracoda %>% group_by(HLA) %>%
  count()

all_barracoda <-  all_barracoda %>% 
  # Rename it so that columns can be merged
  mutate(HLA = str_replace(HLA, "^H-2", "H2-")) %>% 
  # Add identifier colum to merge w/ mupexi
  mutate(identifier = paste(HLA, Sequence, sep = "_"))

all_mupexi <- full_join(mupexi_4t1, mupexi_ct26) %>% 
  # identifier column to merge with barracoda - HLA_peptidename
  mutate(identifier = paste(HLA_allele, Mut_peptide, sep = "_"))


# Merged dataset
# ------------------------------------------------------------------------------
mupexi_barracoda <- left_join(all_barracoda, all_mupexi, by = "identifier")

my_data <- left_join(mupexi_barracoda, sample_info) %>% #bring PE_population info of each sample into barracoda_mupexi file 
  # add HLA_match column if HLA in any of the buffycoat_HLA info file
  mutate(HLA_match = case_when(HLA %in% buffycoat.HLA_info$HLA1 | 
                                 HLA %in% buffycoat.HLA_info$HLA2 | 
                                 HLA %in% buffycoat.HLA_info$HLA3~"yes",
                               TRUE ~ "no"))

# Write data
# ------------------------------------------------------------------------------
write_tsv(x = my_data,
          path = "data/01_my_data.tsv")
