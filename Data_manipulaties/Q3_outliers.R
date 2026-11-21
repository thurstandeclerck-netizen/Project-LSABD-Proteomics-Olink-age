# Question 3: are there samples or proteins that behave unlike the rest?
# Is that a technical error or interesting biology?
# Input:  the clean files in output/ (made by Q2b_cleaning.R)
# Output: numbers in the console + 5 figures in output/

library(tidyverse)

npx <- read_csv(here::here("output", "npx_clean.csv"))
metadata <- read_csv(here::here("output", "metadata_clean.csv"))

# Put metadata and proteins in one table, matched on sample ID
data <- left_join(metadata, npx, by = "sample_kit_id")

# Only the protein values (without the sample ID column)
npx_values <- select(npx, -sample_kit_id)


# --- Part A: samples that behave differently --------------------------------

# A1. PCA summarises 1300 proteins into a few axes (PC1, PC2, ...).
# PCA can't handle NA, so we use only proteins without missing values
npx_complete <- npx_values[, colSums(is.na(npx_values)) == 0]
pca <- prcomp(npx_complete, scale. = TRUE)
summary(pca)$importance[2, 1:5] # share of variation per axis

# Add the PCA scores to the metadata, so we can colour by batch
pca_scores <- bind_cols(metadata, as_tibble(pca$x[, 1:2]))

ggplot(pca_scores, aes(x = PC1, y = PC2, colour = olink_batch_id)) +
  geom_point(alpha = 0.6) +
  labs(title = "PCA of all samples, coloured by batch", colour = "Batch")
ggsave(here::here("output", "q3_pca_batch.png"))

# A2. Count per sample how many proteins have an extreme value.
# scale() turns every protein into z-scores (mean 0, sd 1);
# |z| > 5 means: more than 5 standard deviations away from the mean
z_scores <- scale(npx_values)
metadata <- metadata |>
  mutate(n_extreme = rowSums(abs(z_scores) > 5, na.rm = TRUE))

summary(metadata$n_extreme)

ggplot(metadata, aes(x = n_extreme)) +
  geom_histogram(binwidth = 1) +
  labs(
    title = "Number of extreme protein values per sample",
    x = "Number of proteins with |z| > 5", y = "Number of samples"
  )
ggsave(here::here("output", "q3_extreme_values_per_sample.png"))

# Which samples have more than 20 extreme values?
metadata |>
  filter(n_extreme > 20) |>
  select(sample_kit_id, subject.subjectGuid, sample.visitName, n_extreme)

# Are the other samples of these people normal? (if yes: problem of 1 sample)
metadata |>
  filter(subject.subjectGuid %in% c("BR1018", "BR1043")) |>
  select(subject.subjectGuid, sample.visitName, n_extreme)


# --- Part B: proteins that behave differently -------------------------------

# B1. Which proteins vary the most? apply(..., 2, sd) = sd of every column
sd_per_protein <- apply(npx_values, 2, sd, na.rm = TRUE)
head(sort(sd_per_protein, decreasing = TRUE), 10)

# B2. IL10: very different between batches -> technical
data |>
  group_by(olink_batch_id) |>
  summarise(median_IL10 = median(IL10))

ggplot(data, aes(x = olink_batch_id, y = IL10)) +
  geom_boxplot() +
  labs(title = "IL10 differs strongly between batches", x = "Batch", y = "IL10 (NPX)")
ggsave(here::here("output", "q3_IL10_per_batch.png"))

# B3. GH1 (growth hormone) and PAEP: different between men and women -> biology
data |>
  group_by(subject.biologicalSex, subject.ageGroup) |>
  summarise(
    median_GH1 = median(GH1), median_PAEP = median(PAEP),
    .groups = "drop" # remove the grouping afterwards
  )

ggplot(data, aes(x = subject.biologicalSex, y = PAEP, fill = subject.ageGroup)) +
  geom_boxplot() +
  labs(
    title = "PAEP is high only in young women", x = NULL,
    y = "PAEP (NPX)", fill = "Age group"
  )
ggsave(here::here("output", "q3_PAEP_sex_age.png"))

# B4. PNLIPRP2: a small group of samples has a very low value
ggplot(data, aes(x = PNLIPRP2)) +
  geom_histogram(binwidth = 0.5) +
  labs(title = "PNLIPRP2 has two groups", x = "PNLIPRP2 (NPX)", y = "Number of samples")
ggsave(here::here("output", "q3_PNLIPRP2.png"))
