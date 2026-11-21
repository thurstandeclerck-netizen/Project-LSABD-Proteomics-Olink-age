# Part 1, question 1: how many samples and features are there, and what type
# is each variable?
# Input:  the three csv files in Data_Proteins/ (they are not modified)
# Output: text in the console only, nothing is written to disk

library(readr)
library(dplyr)

# --- Read the data ----------------------------------------------------------
# here() finds the project folder by itself, so no setwd() is needed
npx <- read_csv(
  here::here("Data_Proteins", "proteins_npx_raw.csv"),
  show_col_types = FALSE
)
metadata <- read_csv(
  here::here("Data_Proteins", "sample_metadata.csv"),
  show_col_types = FALSE
)
annotation <- read_csv(
  here::here("Data_Proteins", "protein_annotation.csv"),
  show_col_types = FALSE
)

# --- 1. Number of samples ---------------------------------------------------
# One row in the NPX data = one blood sample
cat("Number of samples (rows in NPX):", nrow(npx), "\n")
cat("Number of unique sample IDs:", n_distinct(npx$sample_kit_id), "\n")

# Does every sample have metadata? Otherwise we can't link it to age
cat(
  "All samples have metadata:",
  all(npx$sample_kit_id %in% metadata$sample_kit_id), "\n"
)

# --- 2. Number of people ----------------------------------------------------
# Important: samples from the same person are not independent
cat("Number of people:", n_distinct(metadata$subject.subjectGuid), "\n")

samples_per_person <- metadata |>
  count(subject.subjectGuid, name = "n_samples")
cat("Samples per person:\n")
print(summary(samples_per_person$n_samples))

# --- 3. Number of protein features ------------------------------------------
# The first column is the sample ID, the rest are proteins
n_protein_columns <- ncol(npx) - 1
cat("\nNumber of protein columns:", n_protein_columns, "\n")

# Old + new gene name of the same protein share the same olink_id
cat("Number of unique olink_ids:", n_distinct(annotation$olink_id), "\n")

# Some proteins were measured in several panels (same assay name)
proteins_in_several_panels <- annotation |>
  count(assay, name = "n_panels") |>
  filter(n_panels > 1)
print(proteins_in_several_panels)

# A protein measured in several panels counts once: subtract the extra copies
extra_copies <- sum(proteins_in_several_panels$n_panels - 1)
unique_proteins <- n_distinct(annotation$olink_id) - extra_copies
cat("Number of unique proteins:", unique_proteins, "\n")

# --- 4. Number of metadata variables ----------------------------------------
cat("\nNumber of metadata columns (incl. sample ID):", ncol(metadata), "\n")

# These two columns look duplicated; if TRUE one of them can be dropped
cat(
  "subject_guid identical to subject.subjectGuid:",
  identical(metadata$subject_guid, metadata$subject.subjectGuid), "\n"
)

# --- 5. Type of each variable -----------------------------------------------
# R only knows 'character' and 'numeric'. Whether a variable is nominal,
# ordinal, continuous or discrete is your own decision; the number of unique
# values helps: few text values = categorical, many decimals = continuous.
variable_types <- tibble(
  variable = names(metadata),
  r_type = sapply(metadata, class),
  n_unique_values = sapply(metadata, n_distinct)
)
print(variable_types, n = Inf)
