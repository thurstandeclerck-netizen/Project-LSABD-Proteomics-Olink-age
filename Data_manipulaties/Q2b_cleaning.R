# Cleaning: apply the decisions from question 2 to the datasets.
# Input:  the three csv files in Data_Proteins/ (they are not changed)
# Output: npx_clean.csv, metadata_clean.csv, annotation_clean.csv in output/

library(tidyverse)

npx <- read_csv(here::here("Data_Proteins", "proteins_npx_raw.csv"))
metadata <- read_csv(here::here("Data_Proteins", "sample_metadata.csv"))
annotation <- read_csv(here::here("Data_Proteins", "protein_annotation.csv"))


# --- Step 1: merge the 9 renamed proteins -----------------------------------
# coalesce(new, old) takes the new value, and the old value where new is NA
npx_merged <- npx |>
  mutate(
    ADGRE5 = coalesce(ADGRE5, CD97),
    GOLM2 = coalesce(GOLM2, CASC4),
    DARS1 = coalesce(DARS1, DARS),
    CEP43 = coalesce(CEP43, FGFR1OP),
    CEP20 = coalesce(CEP20, FOPNL),
    HARS1 = coalesce(HARS1, HARS),
    IRAG2 = coalesce(IRAG2, LRMP),
    CCN4 = coalesce(CCN4, WISP1),
    CCN5 = coalesce(CCN5, WISP2)
  )


# --- Step 2: decide which proteins we keep ----------------------------------
annotation_clean <- annotation |>
  # the old names are missing in 82% of samples, so this removes them
  filter(pct_missing < 50) |>
  # remove proteins below the detection limit in more than half the samples
  filter(pct_below_lod <= 50) |>
  # CXCL8, IL6, TNF are in 4 panels: keep the Neurology version only,
  # because it has the fewest values below the detection limit
  filter(!(assay %in% c("CXCL8", "IL6", "TNF")) | panel == "Neurology") |>
  # names without the _OID part, e.g. IL6_OID20911 becomes IL6
  mutate(protein = assay)

nrow(annotation_clean)


# --- Step 3: keep only those proteins in the NPX data -----------------------
npx_clean <- npx_merged |>
  rename(CXCL8 = CXCL8_OID20997, IL6 = IL6_OID20911, TNF = TNF_OID20848) |>
  select(sample_kit_id, all_of(annotation_clean$protein))

dim(npx_clean)


# --- Step 4: clean the metadata ---------------------------------------------
metadata_clean <- metadata |>
  # remove variables that are mostly empty or not reliable
  select(
    -cmv.igg_serology, -cmv.igg_serology_interpretation,
    -infl.anti_ccp3, -infl.anti_ccp31,
    -infl.rf_iga_result, -infl.rf_iga_interpretation,
    -infl.rf_igm_result, -infl.rf_igm_interpretation,
    -chem.egfr_aa, -chem.egfr_non_aa, -chem.protein,
    -n_missing_npx, -subject_guid
  ) |>
  # a weight of 0 kg is impossible: make it missing
  mutate(am.weight = if_else(am.weight == 0, NA, am.weight))

# Adult height does not change: give each person the median of their heights.
# Then calculate BMI again: weight (kg) / height (m) squared
metadata_clean <- metadata_clean |>
  group_by(subject.subjectGuid) |>
  mutate(am.height = median(am.height, na.rm = TRUE)) |>
  ungroup() |>
  mutate(am.bmi = am.weight / (am.height / 100)^2)

dim(metadata_clean)


# --- Step 5: checks ---------------------------------------------------------
# Missing values in the proteins: should be much lower than the 10486 before
sum(is.na(npx_clean))

# Samples must be in the same order in both files
all(npx_clean$sample_kit_id == metadata_clean$sample_kit_id)

# Height (cm), weight (kg) and BMI should be in a logical range
summary(metadata_clean$am.height)
summary(metadata_clean$am.weight)
summary(metadata_clean$am.bmi)


# --- Step 6: save the clean datasets ----------------------------------------
write_csv(npx_clean, here::here("output", "npx_clean.csv"))
write_csv(metadata_clean, here::here("output", "metadata_clean.csv"))
write_csv(annotation_clean, here::here("output", "annotation_clean.csv"))
