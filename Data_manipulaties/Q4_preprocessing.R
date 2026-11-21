# Question 4: which pre-processing steps does the data need?
# Input:  the clean files in output/ (made by Q2b_cleaning.R)
# Output: npx_corrected.csv + metadata_final.csv + 3 figures in output/

library(tidyverse)

npx <- read_csv(here::here("output", "npx_clean.csv"))
metadata <- read_csv(here::here("output", "metadata_clean.csv"))


# --- Step 1: filtering - remove the 2 bad samples from question 3 -----------
bad_samples <- c("KT00397", "KT02473")
npx <- filter(npx, !(sample_kit_id %in% bad_samples))
metadata <- filter(metadata, !(sample_kit_id %in% bad_samples))
nrow(npx)


# --- Step 2: log-transformation - is it needed? -----------------------------
# NPX is already on a log2 scale. Compare IL6 as NPX and on the original
# scale (2^NPX): the original scale is skewed, the NPX scale is not
ggplot(npx, aes(x = IL6)) +
  geom_histogram(bins = 40) +
  labs(title = "IL6 in NPX (log2 scale): symmetric", x = "IL6 (NPX)")
ggsave(here::here("output", "q4_IL6_log_scale.png"))

ggplot(npx, aes(x = 2^IL6)) +
  geom_histogram(bins = 40) +
  labs(title = "IL6 on the original scale: skewed", x = "IL6 (2^NPX)")
ggsave(here::here("output", "q4_IL6_original_scale.png"))


# --- Step 3: normalisation - correct for batch ------------------------------
# Long format: one row per sample x protein, with the batch of the sample
npx_long <- npx |>
  pivot_longer(-sample_kit_id, names_to = "protein", values_to = "npx") |>
  left_join(select(metadata, sample_kit_id, olink_batch_id), by = "sample_kit_id")

# Per batch and per protein: subtract the median of that batch.
# Afterwards every batch has median 0, so batches are comparable
npx_long <- npx_long |>
  group_by(olink_batch_id, protein) |>
  mutate(npx_corrected = npx - median(npx, na.rm = TRUE)) |>
  ungroup()

# Back to wide format: one row per sample, one column per protein
npx_corrected <- npx_long |>
  select(sample_kit_id, protein, npx_corrected) |>
  pivot_wider(names_from = protein, values_from = npx_corrected)

# Check: same samples in the same order as the metadata
all(npx_corrected$sample_kit_id == metadata$sample_kit_id)


# --- Step 4: check the correction with a PCA --------------------------------
# Standardisation: scale. = TRUE gives every protein mean 0 and sd 1,
# so every protein counts equally in the PCA
npx_values <- select(npx_corrected, -sample_kit_id)
npx_complete <- npx_values[, colSums(is.na(npx_values)) == 0]
pca <- prcomp(npx_complete, scale. = TRUE)
pca_scores <- bind_cols(metadata, as_tibble(pca$x[, 1:2]))

ggplot(pca_scores, aes(x = PC1, y = PC2, colour = olink_batch_id)) +
  geom_point(alpha = 0.6) +
  labs(title = "PCA after batch correction, coloured by batch", colour = "Batch")
ggsave(here::here("output", "q4_pca_after_correction.png"))


# --- Step 5: save -----------------------------------------------------------
write_csv(npx_corrected, here::here("output", "npx_corrected.csv"))
write_csv(metadata, here::here("output", "metadata_final.csv"))
