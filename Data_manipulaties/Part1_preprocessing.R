# Part 1: explore, pre-process and quality-control the protein data
# Input:  the three csv files in Data_Proteins/ (they are never changed)
# Output: figures and the clean datasets in output/ (all start with "part1_")
#
# The script follows the guiding questions of Part 1 in order.
# Run it from top to bottom: each question uses the results of the previous one.

library(tidyverse)

npx <- read_csv(here::here("Data_Proteins", "proteins_npx_raw.csv"))
metadata <- read_csv(here::here("Data_Proteins", "sample_metadata.csv"))
annotation <- read_csv(here::here("Data_Proteins", "protein_annotation.csv"))


# =============================================================================
# QUESTION 1: How many samples and features are there,
#             and what type is each variable?
# =============================================================================

# Rows = samples, columns = variables
dim(npx) # 867 samples; 1 ID column + 1481 proteins
dim(metadata) # 867 samples; 76 variables

# Number of different people (one person gave several samples)
n_distinct(metadata$subject.subjectGuid)

# How many samples does each person have?
samples_per_person <- count(metadata, subject.subjectGuid)
summary(samples_per_person$n)

# Old + new gene name of the same protein share one olink_id
n_distinct(annotation$olink_id)

# Proteins that were measured in more than one panel
count(annotation, assay) |>
  filter(n > 1)

# Type of each variable: <chr> = text, <dbl> = number
glimpse(metadata)


# =============================================================================
# QUESTION 2: Where are the missing values, why are they missing,
#             and how do we handle them?
# =============================================================================

# --- 2a. Missing values in the protein data ---------------------------------
sum(is.na(npx))

# Missing values per protein, only the proteins with missing values
missing_per_protein <- colSums(is.na(npx))
missing_per_protein[missing_per_protein > 0]

# Missing values are linked to the batch:
# CD97 (old name of ADGRE5) is only measured in the first batch,
# BID is missing in the whole batch Q-04064
table(metadata$olink_batch_id, is.na(npx$CD97))
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
ggsave(here::here("output", "part1_q2_missing_per_protein.png"))

# --- 2b. Values below the limit of detection (LOD) --------------------------
# Number of proteins below the LOD in more than half of the samples
sum(annotation$pct_below_lod > 50)

# Same for young and old? If yes, filtering on LOD is fair to both groups
metadata |>
  group_by(subject.ageGroup) |>
  summarise(median_pct_below_lod = median(pct_below_lod))

ggplot(annotation, aes(x = pct_below_lod)) +
  geom_histogram(binwidth = 5) +
  geom_vline(xintercept = 50, linetype = "dashed") +
  labs(
    title = "Most proteins are rarely below the LOD",
    x = "% of samples below the LOD", y = "Number of proteins"
  )
ggsave(here::here("output", "part1_q2_pct_below_lod.png"))

# --- 2c. Missing values in the metadata -------------------------------------
missing_per_variable <- colSums(is.na(metadata))
missing_per_variable[missing_per_variable > 0]

# CMV serology was only measured at some visits (missing by design)
table(metadata$sample.visitName, is.na(metadata$cmv.igg_serology))

# Total protein was not measured anymore in 2022 (missing by design)
table(metadata$sample.drawYear, is.na(metadata$chem.protein))

# BMI is missing as often for young as for old (looks random)
table(metadata$subject.ageGroup, is.na(metadata$am.bmi))

# --- 2d. Handle the missing values: clean the data --------------------------

# Step 1: merge the 9 renamed proteins.
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

# Step 2: decide which proteins we keep
annotation_clean <- annotation |>
  # the old names are missing in 82% of samples, so this removes them
  filter(pct_missing < 50) |>
  # remove proteins below the LOD in more than half of the samples
  filter(pct_below_lod <= 50) |>
  # CXCL8, IL6, TNF are in 4 panels: keep only the Neurology version,
  # because it has the fewest values below the LOD
  filter(!(assay %in% c("CXCL8", "IL6", "TNF")) | panel == "Neurology") |>
  # names without the _OID part, e.g. IL6_OID20911 becomes IL6
  mutate(protein = assay)

nrow(annotation_clean) # number of proteins we keep

# Step 3: keep only those proteins in the protein data
npx_clean <- npx_merged |>
  rename(CXCL8 = CXCL8_OID20997, IL6 = IL6_OID20911, TNF = TNF_OID20848) |>
  select(sample_kit_id, all_of(annotation_clean$protein))

# Step 4: clean the metadata
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

# Adult height does not change: give each person the median of their heights,
# then calculate BMI again: weight (kg) / height (m) squared
metadata_clean <- metadata_clean |>
  group_by(subject.subjectGuid) |>
  mutate(am.height = median(am.height, na.rm = TRUE)) |>
  ungroup() |>
  mutate(am.bmi = am.weight / (am.height / 100)^2)

# Checks: fewer missing values, and samples in the same order in both tables
sum(is.na(npx_clean))
all(npx_clean$sample_kit_id == metadata_clean$sample_kit_id)


# =============================================================================
# QUESTION 3: Are there samples or proteins that behave unlike the rest?
#             Technical error or interesting biology?
# =============================================================================

# --- 3a. Samples: PCA -------------------------------------------------------
# PCA summarises all proteins into a few axes (PC1, PC2, ...).
# PCA can't handle NA, so we only use proteins without missing values
npx_values <- select(npx_clean, -sample_kit_id)
npx_complete <- npx_values[, colSums(is.na(npx_values)) == 0]
pca <- prcomp(npx_complete, scale. = TRUE)
summary(pca)$importance[2, 1:5] # share of all differences per axis

pca_scores <- bind_cols(metadata_clean, as_tibble(pca$x[, 1:2]))

# The samples group by batch: a technical effect
ggplot(pca_scores, aes(x = PC1, y = PC2, colour = olink_batch_id)) +
  geom_point(alpha = 0.6) +
  labs(title = "PCA before batch correction", colour = "Batch")
ggsave(here::here("output", "part1_q3_pca_before_correction.png"))

# --- 3b. Samples: extreme values --------------------------------------------
# scale() turns every protein into z-scores (mean 0, sd 1).
# |z| > 5 = more than 5 standard deviations away from the mean
z_scores <- scale(npx_values)
metadata_clean <- metadata_clean |>
  mutate(n_extreme = rowSums(abs(z_scores) > 5, na.rm = TRUE))

summary(metadata_clean$n_extreme)

# Two samples have more than 20 extreme values
metadata_clean |>
  filter(n_extreme > 20) |>
  select(sample_kit_id, subject.subjectGuid, sample.visitName, n_extreme)

# The other samples of these two people are normal,
# so the problem is in the sample (technical), not in the person
metadata_clean |>
  filter(subject.subjectGuid %in% c("BR1018", "BR1043")) |>
  select(subject.subjectGuid, sample.visitName, n_extreme)

# Remove the two bad samples
bad_samples <- c("KT00397", "KT02473")
npx_clean <- filter(npx_clean, !(sample_kit_id %in% bad_samples))
metadata_clean <- filter(metadata_clean, !(sample_kit_id %in% bad_samples))

# --- 3c. Proteins -----------------------------------------------------------
# Which proteins vary the most? apply(..., 2, sd) = sd of every column
sd_per_protein <- apply(select(npx_clean, -sample_kit_id), 2, sd, na.rm = TRUE)
head(sort(sd_per_protein, decreasing = TRUE), 10)

# Metadata and proteins in one table, to make figures
data <- left_join(metadata_clean, npx_clean, by = "sample_kit_id")

# IL10 differs between batches -> technical
ggplot(data, aes(x = olink_batch_id, y = IL10)) +
  geom_boxplot() +
  labs(title = "IL10 differs between batches (technical)", x = "Batch")
ggsave(here::here("output", "part1_q3_IL10_per_batch.png"))

# PAEP is high only in young women -> biology (hormones)
ggplot(data, aes(x = subject.biologicalSex, y = PAEP, fill = subject.ageGroup)) +
  geom_boxplot() +
  labs(title = "PAEP is high only in young women (biology)", x = NULL, fill = "Age group")
ggsave(here::here("output", "part1_q3_PAEP_sex_age.png"))

# PNLIPRP2 has two groups of people -> biology (most likely genetic)
ggplot(data, aes(x = PNLIPRP2)) +
  geom_histogram(binwidth = 0.5) +
  labs(title = "PNLIPRP2 has two groups (biology)", y = "Number of samples")
ggsave(here::here("output", "part1_q3_PNLIPRP2.png"))


# =============================================================================
# QUESTION 4: Which pre-processing steps are needed?
# =============================================================================
# Filtering: already done in question 2 (proteins) and 3 (samples).
# Standardisation: only needed for PCA, done by prcomp(..., scale. = TRUE).

# --- 4a. Log-transformation: not needed -------------------------------------
# NPX is already on a log2 scale: on that scale IL6 is symmetric,
# on the original scale (2^NPX) it is skewed
ggplot(npx_clean, aes(x = IL6)) +
  geom_histogram(bins = 40) +
  labs(title = "IL6 in NPX (log2 scale): symmetric")
ggsave(here::here("output", "part1_q4_IL6_log_scale.png"))

ggplot(npx_clean, aes(x = 2^IL6)) +
  geom_histogram(bins = 40) +
  labs(title = "IL6 on the original scale: skewed")
ggsave(here::here("output", "part1_q4_IL6_original_scale.png"))

# --- 4b. Normalisation: correct for batch -----------------------------------
# Long format: one row per sample x protein, with the batch of the sample
npx_long <- npx_clean |>
  pivot_longer(-sample_kit_id, names_to = "protein", values_to = "npx") |>
  left_join(select(metadata_clean, sample_kit_id, olink_batch_id),
    by = "sample_kit_id"
  )

# Per batch and per protein: subtract the median of that batch,
# so every batch gets median 0 and the batches become comparable
npx_long <- npx_long |>
  group_by(olink_batch_id, protein) |>
  mutate(npx_corrected = npx - median(npx, na.rm = TRUE)) |>
  ungroup()

# Back to wide format: one row per sample, one column per protein
npx_corrected <- npx_long |>
  select(sample_kit_id, protein, npx_corrected) |>
  pivot_wider(names_from = protein, values_from = npx_corrected)

all(npx_corrected$sample_kit_id == metadata_clean$sample_kit_id)


# =============================================================================
# QUESTION 5: Do samples cluster by biology or by a technical factor?
# =============================================================================

# The same PCA as in question 3, now on the batch-corrected data
npx_values_corrected <- select(npx_corrected, -sample_kit_id)
npx_complete_corrected <- npx_values_corrected[, colSums(is.na(npx_values_corrected)) == 0]
pca_corrected <- prcomp(npx_complete_corrected, scale. = TRUE)
summary(pca_corrected)$importance[2, 1:5]

pca_scores_corrected <- bind_cols(metadata_clean, as_tibble(pca_corrected$x[, 1:5]))

# R-squared of lm() = share of an axis explained by a factor (0 to 1)
summary(lm(PC1 ~ olink_batch_id, data = pca_scores_corrected))$r.squared # batch
summary(lm(PC2 ~ subject.ageGroup, data = pca_scores_corrected))$r.squared # age
summary(lm(PC4 ~ subject.biologicalSex, data = pca_scores_corrected))$r.squared # sex
summary(lm(PC2 ~ cmv_status, data = pca_scores_corrected))$r.squared # CMV
summary(lm(PC2 ~ sample.visitName, data = pca_scores_corrected))$r.squared # visit
summary(lm(PC2 ~ subject.subjectGuid, data = pca_scores_corrected))$r.squared # person

# Batch: the batches are now mixed (the correction worked)
ggplot(pca_scores_corrected, aes(x = PC1, y = PC2, colour = olink_batch_id)) +
  geom_point(alpha = 0.6) +
  labs(title = "PCA after batch correction: batch", colour = "Batch")
ggsave(here::here("output", "part1_q5_pca_batch.png"))

# Age group: young and old are separated along PC2
ggplot(pca_scores_corrected, aes(x = PC1, y = PC2, colour = subject.ageGroup)) +
  geom_point(alpha = 0.6) +
  labs(title = "PCA after batch correction: age group", colour = "Age group")
ggsave(here::here("output", "part1_q5_pca_age_group.png"))

# Sex: men and women are separated along PC4
ggplot(pca_scores_corrected, aes(x = PC2, y = PC4, colour = subject.biologicalSex)) +
  geom_point(alpha = 0.6) +
  labs(title = "PCA after batch correction: sex", colour = "Sex")
ggsave(here::here("output", "part1_q5_pca_sex.png"))


# =============================================================================
# QUESTION 6: Is the pre-processed dataset suitable for the next parts?
# =============================================================================
# Yes, for comparing young and old (also with sex and CMV), but:
# - age is only 2 groups (25-37 and 55-67): no gradual age effects
# - samples of the same person are not independent: use one sample per
#   person, or a model that accounts for repeated samples
# - for statistical tests: use npx_clean with batch in the model,
#   e.g. lm(protein ~ age group + batch); use npx_corrected for PCA and plots

# Save the final datasets for the next parts
write_csv(npx_clean, here::here("output", "part1_npx_clean.csv"))
write_csv(npx_corrected, here::here("output", "part1_npx_corrected.csv"))
write_csv(metadata_clean, here::here("output", "part1_metadata_clean.csv"))
