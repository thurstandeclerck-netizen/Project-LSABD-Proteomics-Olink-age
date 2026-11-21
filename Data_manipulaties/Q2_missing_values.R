# Question 2: where are the missing values, and why are they missing?
# Input:  the three csv files in Data_Proteins/
# Output: numbers in the console + 2 figures in output/

library(tidyverse)

npx <- read_csv(here::here("Data_Proteins", "proteins_npx_raw.csv"))
metadata <- read_csv(here::here("Data_Proteins", "sample_metadata.csv"))
annotation <- read_csv(here::here("Data_Proteins", "protein_annotation.csv"))


# --- Part A: missing values in the protein data -----------------------------

# Total number of missing values
sum(is.na(npx))

# Number of missing values per protein, only the proteins with missing values
missing_per_protein <- colSums(is.na(npx))
missing_per_protein[missing_per_protein > 0]

# Is missing linked to the batch? Three examples:
# CD97 = old name, ADGRE5 = new name of the same protein
table(metadata$olink_batch_id, is.na(npx$CD97))
table(metadata$olink_batch_id, is.na(npx$ADGRE5))
table(metadata$olink_batch_id, is.na(npx$BID))

# Figure: number of missing values per protein
missing_table <- tibble(
  protein = names(missing_per_protein),
  n_missing = missing_per_protein
) |>
  filter(n_missing > 0)

ggplot(missing_table, aes(x = n_missing, y = protein)) +
  geom_col() +
  labs(
    title = "Missing values per protein",
    x = "Number of samples with a missing value", y = NULL
  )
ggsave(here::here("output", "q2_missing_per_protein.png"))


# --- Part B: values below the limit of detection (LOD) ----------------------

# Number of proteins below the LOD in more than half of the samples
sum(annotation$pct_below_lod > 50)

# Is it the same for young and old? (if yes, filtering is fair to both)
metadata |>
  group_by(subject.ageGroup) |>
  summarise(median_pct_below_lod = median(pct_below_lod))

# Figure: % of samples below the LOD for each protein
ggplot(annotation, aes(x = pct_below_lod)) +
  geom_histogram(binwidth = 5) +
  geom_vline(xintercept = 50, linetype = "dashed") +
  labs(
    title = "Most proteins are rarely below the LOD",
    x = "% of samples below the LOD", y = "Number of proteins"
  )
ggsave(here::here("output", "q2_pct_below_lod.png"))


# --- Part C: missing values in the metadata ---------------------------------

# Number of missing values per variable, only the variables with missing values
missing_per_variable <- colSums(is.na(metadata))
missing_per_variable[missing_per_variable > 0]

# CMV serology: only measured at some visits?
table(metadata$sample.visitName, is.na(metadata$cmv.igg_serology))

# Total protein: stopped being measured after a certain year?
table(metadata$sample.drawYear, is.na(metadata$chem.protein))

# BMI: missing as often for young as for old?
table(metadata$subject.ageGroup, is.na(metadata$am.bmi))
