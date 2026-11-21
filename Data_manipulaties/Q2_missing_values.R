# Part 1, question 2: where are the missing values, why are they missing,
# and do we need to handle them?
# Input:  the three csv files in Data_Proteins/ (they are not modified)
# Output: text in the console + 3 figures in output/

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

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
dir.create(here::here("output"), showWarnings = FALSE)


# =============================================================================
# PART A: missing values in the protein data (NPX)
# =============================================================================

# --- A1. How much is missing in total? --------------------------------------
npx_values <- select(npx, -sample_kit_id)
cat(
  "Missing NPX values:", sum(is.na(npx_values)), "of",
  nrow(npx_values) * ncol(npx_values), "\n"
)
cat("Percentage:", round(mean(is.na(npx_values)) * 100, 2), "%\n")

# --- A2. Long format: one row per sample x protein --------------------------
# pivot_longer stacks the 1481 protein columns on top of each other, so we
# can count per batch and per protein with group_by()
batch_per_sample <- select(metadata, sample_kit_id, olink_batch_id)
npx_long <- npx |>
  left_join(batch_per_sample, by = "sample_kit_id") |>
  pivot_longer(
    cols = -c(sample_kit_id, olink_batch_id),
    names_to = "protein", values_to = "npx"
  )

# --- A3. Which proteins are missing, and in which batch? --------------------
missing_per_batch <- npx_long |>
  group_by(olink_batch_id, protein) |>
  summarise(pct_missing = mean(is.na(npx)) * 100, .groups = "drop")

# Only proteins with at least one missing value, otherwise the plot is unreadable
proteins_with_missing <- missing_per_batch |>
  filter(pct_missing > 0) |>
  distinct(protein)
cat(
  "\nNumber of proteins with at least one missing value:",
  nrow(proteins_with_missing), "\n"
)

missing_plot_data <- missing_per_batch |>
  filter(protein %in% proteins_with_missing$protein)

plot_npx_missing <- ggplot(
  missing_plot_data,
  aes(x = protein, y = olink_batch_id, fill = pct_missing)
) +
  geom_tile(colour = "white", linewidth = 0.5) +
  scale_fill_gradient(low = "#f1f5f9", high = "#1e40af", name = "% missing") +
  labs(
    title = "Missing NPX values occur in whole batches",
    subtitle = "Only proteins with at least one missing value",
    x = NULL, y = "Olink batch"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid = element_blank()
  )
ggsave(here::here("output", "q2_npx_missing_per_batch.png"),
  plot_npx_missing,
  width = 10, height = 4, dpi = 150
)

# --- A4. Old and new gene names: do they complement each other? -------------
# pct_missing > 50 are the old names; same_assay_as gives the new name
renamed <- annotation |>
  filter(pct_missing > 50) |>
  select(old_name = protein, new_name = same_assay_as)
print(renamed)

old_values <- as.matrix(npx[renamed$old_name])
new_values <- as.matrix(npx[renamed$new_name])
# If they complement each other perfectly, both numbers are 0
cat(
  "Samples where old AND new name are measured:",
  sum(!is.na(old_values) & !is.na(new_values)), "\n"
)
cat(
  "Samples where old AND new name are missing:",
  sum(is.na(old_values) & is.na(new_values)), "\n"
)

# In which batch was the old name used?
print(table(metadata$olink_batch_id, measured_under_old_name = !is.na(npx$CD97)))

# --- A5. Is the column n_missing_npx in the metadata correct? ---------------
own_count <- rowSums(is.na(npx_values))
cat(
  "\nn_missing_npx equal to our own count:",
  all(own_count == metadata$n_missing_npx), "\n"
)
print(table(own_count, n_missing_npx = metadata$n_missing_npx))


# =============================================================================
# PART B: below the limit of detection (LOD)
# =============================================================================
# Values below the LOD exist, but can't be told apart from background noise.

cat(
  "\nProteins with > 50% of samples below the LOD:",
  sum(annotation$pct_below_lod > 50), "\n"
)
print(annotation |>
  group_by(panel) |>
  summarise(n_proteins = n(), more_than_50pct_below_lod = sum(pct_below_lod > 50)))

# Important: does % below LOD differ between young and old? If not, an LOD
# filter does not bias the young vs old comparison
print(metadata |>
  group_by(subject.ageGroup) |>
  summarise(median_pct_below_lod = median(pct_below_lod)))

plot_lod <- ggplot(annotation, aes(x = pct_below_lod)) +
  geom_histogram(
    binwidth = 5, boundary = 0, fill = "#1e40af",
    colour = "white"
  ) +
  geom_vline(xintercept = 50, linetype = "dashed", colour = "grey30") +
  facet_wrap(~panel) +
  labs(
    title = "Most proteins are rarely below the limit of detection",
    subtitle = "Dashed line = 50%, the cut-off we use for filtering",
    x = "% of samples below the LOD", y = "Number of proteins"
  ) +
  theme_minimal(base_size = 11)
ggsave(here::here("output", "q2_pct_below_lod.png"),
  plot_lod,
  width = 8, height = 5, dpi = 150
)


# =============================================================================
# PART C: missing values in the metadata
# =============================================================================

metadata_missing <- tibble(
  variable = names(metadata),
  n_na = colSums(is.na(metadata))
) |>
  mutate(pct_na = round(n_na / nrow(metadata) * 100, 1)) |>
  filter(n_na > 0) |>
  arrange(desc(n_na))
print(metadata_missing, n = Inf)

plot_metadata_missing <- ggplot(
  metadata_missing,
  aes(x = pct_na, y = reorder(variable, pct_na))
) +
  geom_col(fill = "#1e40af", width = 0.7) +
  labs(
    title = "Missing values per metadata variable",
    subtitle = "Only variables with at least one NA",
    x = "% of samples with NA", y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.y = element_blank())
ggsave(here::here("output", "q2_metadata_missing.png"),
  plot_metadata_missing,
  width = 7, height = 9, dpi = 150
)

# --- C1. Are the lab values missing in the same samples? --------------------
cat("\nAre ALT (chemistry) and white blood cells (blood count) missing together?\n")
print(table(alt_na = is.na(metadata$chem.alt), wbc_na = is.na(metadata$bc.wbc)))

# --- C2. Is something missing 'by design' of the study? ---------------------
# CMV serology per visit: only measured at certain visits?
print(table(metadata$sample.visitName, cmv_na = is.na(metadata$cmv.igg_serology)))

# Total protein per year: did they stop measuring after a certain year?
print(table(metadata$sample.drawYear, protein_na = is.na(metadata$chem.protein)))

# eGFR exists in two versions, depending on race
print(table(metadata$subject.race, egfr_na = is.na(metadata$chem.egfr_aa)))

# --- C3. BMI: is missingness related to age? --------------------------------
print(table(metadata$subject.ageGroup, bmi_na = is.na(metadata$am.bmi)))

# Adult height does not change: is it known at another visit?
height_check <- metadata |>
  group_by(subject.subjectGuid) |>
  mutate(height_ever_measured = any(!is.na(am.height))) |>
  ungroup()
cat("Samples without height:", sum(is.na(height_check$am.height)), "\n")
cat(
  "...of which height is known at another visit:",
  sum(is.na(height_check$am.height) & height_check$height_ever_measured), "\n"
)
