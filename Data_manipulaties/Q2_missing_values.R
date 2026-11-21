# Part 1, vraag 2: waar zitten de ontbrekende waarden, waarom ontbreken ze,
# en moeten we er iets mee doen?
# Input:  de drie csv-bestanden in Data_Proteins/ (worden niet aangepast)
# Output: tekst in de console + 3 figuren in output/

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
annotatie <- read_csv(
  here::here("Data_Proteins", "protein_annotation.csv"),
  show_col_types = FALSE
)
dir.create(here::here("output"), showWarnings = FALSE)


# =============================================================================
# DEEL A: ontbrekende waarden in de eiwitdata (NPX)
# =============================================================================

# --- A1. Hoeveel ontbreekt er in totaal? -----------------------------------
npx_waarden <- select(npx, -sample_kit_id)
cat(
  "Ontbrekende NPX-waarden:", sum(is.na(npx_waarden)), "van",
  nrow(npx_waarden) * ncol(npx_waarden), "\n"
)
cat("Percentage:", round(mean(is.na(npx_waarden)) * 100, 2), "%\n")

# --- A2. Lang formaat: één rij per staal x eiwit ---------------------------
# pivot_longer zet de 1481 eiwitkolommen onder elkaar, zodat we per batch
# en per eiwit kunnen tellen met group_by()
batch_per_staal <- select(metadata, sample_kit_id, olink_batch_id)
npx_lang <- npx |>
  left_join(batch_per_staal, by = "sample_kit_id") |>
  pivot_longer(
    cols = -c(sample_kit_id, olink_batch_id),
    names_to = "eiwit", values_to = "npx"
  )

# --- A3. Welke eiwitten missen, en in welke batch? -------------------------
missing_per_batch <- npx_lang |>
  group_by(olink_batch_id, eiwit) |>
  summarise(pct_missing = mean(is.na(npx)) * 100, .groups = "drop")

# Alleen eiwitten die ergens iets missen, anders wordt de figuur onleesbaar
eiwitten_met_missing <- missing_per_batch |>
  filter(pct_missing > 0) |>
  distinct(eiwit)
cat(
  "\nAantal eiwitten met minstens één ontbrekende waarde:",
  nrow(eiwitten_met_missing), "\n"
)

missing_figuur_data <- missing_per_batch |>
  filter(eiwit %in% eiwitten_met_missing$eiwit)

figuur_npx_missing <- ggplot(
  missing_figuur_data,
  aes(x = eiwit, y = olink_batch_id, fill = pct_missing)
) +
  geom_tile(colour = "white", linewidth = 0.5) +
  scale_fill_gradient(low = "#f1f5f9", high = "#1e40af", name = "% missing") +
  labs(
    title = "Ontbrekende NPX-waarden zitten in vaste batches",
    subtitle = "Alleen de eiwitten met minstens één ontbrekende waarde",
    x = NULL, y = "Olink-batch"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid = element_blank()
  )
ggsave(here::here("output", "q2_npx_missing_per_batch.png"),
  figuur_npx_missing,
  width = 10, height = 4, dpi = 150
)

# --- A4. Oude en nieuwe gennamen: vullen ze elkaar aan? --------------------
# pct_missing > 50 zijn de oude namen; same_assay_as geeft de nieuwe naam
hernoemd <- annotatie |>
  filter(pct_missing > 50) |>
  select(oude_naam = protein, nieuwe_naam = same_assay_as)
print(hernoemd)

oud <- as.matrix(npx[hernoemd$oude_naam])
nieuw <- as.matrix(npx[hernoemd$nieuwe_naam])
# Als ze elkaar perfect aanvullen, zijn beide getallen 0
cat("Stalen waar oud EN nieuw gemeten is:", sum(!is.na(oud) & !is.na(nieuw)), "\n")
cat("Stalen waar oud EN nieuw ontbreekt:", sum(is.na(oud) & is.na(nieuw)), "\n")

# In welke batch is de oude naam gebruikt?
print(table(metadata$olink_batch_id, gemeten_onder_oude_naam = !is.na(npx$CD97)))

# --- A5. Klopt de kolom n_missing_npx uit de metadata? ---------------------
eigen_telling <- rowSums(is.na(npx_waarden))
cat(
  "\nn_missing_npx gelijk aan eigen telling:",
  all(eigen_telling == metadata$n_missing_npx), "\n"
)
print(table(eigen_telling, n_missing_npx = metadata$n_missing_npx))


# =============================================================================
# DEEL B: onder de detectielimiet (LOD)
# =============================================================================
# Waarden onder de LOD zijn er wel, maar niet te onderscheiden van ruis.

cat(
  "\nEiwitten met > 50% van de stalen onder de LOD:",
  sum(annotatie$pct_below_lod > 50), "\n"
)
print(annotatie |>
  group_by(panel) |>
  summarise(aantal = n(), boven_50pct_onder_lod = sum(pct_below_lod > 50)))

# Belangrijk: verschilt % onder LOD tussen jong en oud? Zo niet, dan maakt
# een LOD-filter de vergelijking jong vs oud niet scheef
print(metadata |>
  group_by(subject.ageGroup) |>
  summarise(mediaan_pct_onder_lod = median(pct_below_lod)))

figuur_lod <- ggplot(annotatie, aes(x = pct_below_lod)) +
  geom_histogram(
    binwidth = 5, boundary = 0, fill = "#1e40af",
    colour = "white"
  ) +
  geom_vline(xintercept = 50, linetype = "dashed", colour = "grey30") +
  facet_wrap(~panel) +
  labs(
    title = "De meeste eiwitten liggen zelden onder de detectielimiet",
    subtitle = "Stippellijn = 50%, de grens die we gebruiken om te filteren",
    x = "% stalen onder de LOD", y = "Aantal eiwitten"
  ) +
  theme_minimal(base_size = 11)
ggsave(here::here("output", "q2_pct_onder_lod.png"),
  figuur_lod,
  width = 8, height = 5, dpi = 150
)


# =============================================================================
# DEEL C: ontbrekende waarden in de metadata
# =============================================================================

metadata_missing <- tibble(
  variabele = names(metadata),
  aantal_na = colSums(is.na(metadata))
) |>
  mutate(pct_na = round(aantal_na / nrow(metadata) * 100, 1)) |>
  filter(aantal_na > 0) |>
  arrange(desc(aantal_na))
print(metadata_missing, n = Inf)

figuur_meta <- ggplot(
  metadata_missing,
  aes(x = pct_na, y = reorder(variabele, pct_na))
) +
  geom_col(fill = "#1e40af", width = 0.7) +
  labs(
    title = "Ontbrekende waarden per metadata-variabele",
    subtitle = "Alleen variabelen met minstens één NA",
    x = "% stalen met NA", y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.y = element_blank())
ggsave(here::here("output", "q2_metadata_missing.png"),
  figuur_meta,
  width = 7, height = 9, dpi = 150
)

# --- C1. Missen de labo-waarden in steeds dezelfde stalen? -----------------
cat("\nOntbreken ALT (chemie) en witte bloedcellen (bloedbeeld) in dezelfde stalen?\n")
print(table(alt_na = is.na(metadata$chem.alt), wbc_na = is.na(metadata$bc.wbc)))

# --- C2. Ontbreekt iets 'door het ontwerp' van de studie? ------------------
# CMV-serologie per bezoek: alleen bij bepaalde bezoeken gemeten?
print(table(metadata$sample.visitName, cmv_na = is.na(metadata$cmv.igg_serology)))

# Totaal eiwit per jaar: gestopt met meten na een bepaald jaar?
print(table(metadata$sample.drawYear, protein_na = is.na(metadata$chem.protein)))

# eGFR bestaat in twee versies, afhankelijk van ras
print(table(metadata$subject.race, egfr_na = is.na(metadata$chem.egfr_aa)))

# --- C3. BMI: hangt het ontbreken samen met leeftijd? ----------------------
print(table(metadata$subject.ageGroup, bmi_na = is.na(metadata$am.bmi)))

# Lengte verandert niet bij volwassenen: is ze gekend bij een ander bezoek?
lengte_check <- metadata |>
  group_by(subject.subjectGuid) |>
  mutate(lengte_ooit_gemeten = any(!is.na(am.height))) |>
  ungroup()
cat("Stalen zonder lengte:", sum(is.na(lengte_check$am.height)), "\n")
cat(
  "...waarvan lengte bij een ander bezoek wel gekend:",
  sum(is.na(lengte_check$am.height) & lengte_check$lengte_ooit_gemeten), "\n"
)
