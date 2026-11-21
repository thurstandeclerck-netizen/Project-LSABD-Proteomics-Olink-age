# Part 1, vraag 1: hoeveel stalen en features zijn er, en welk type heeft
# elke variabele?
# Input:  de drie csv-bestanden in Data_Proteins/ (worden niet aangepast)
# Output: alleen tekst in de console, er wordt niets weggeschreven

library(readr)
library(dplyr)

# --- Inlezen ----------------------------------------------------------------
# here() zoekt zelf de projectmap, dus geen setwd() nodig
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

# --- 1. Aantal stalen -------------------------------------------------------
# Eén rij in de NPX-data = één bloedstaal
cat("Aantal stalen (rijen NPX):", nrow(npx), "\n")
cat("Aantal unieke staal-ID's:", n_distinct(npx$sample_kit_id), "\n")

# Heeft elk staal ook metadata? Anders kunnen we het niet koppelen aan leeftijd
cat(
  "Alle stalen hebben metadata:",
  all(npx$sample_kit_id %in% metadata$sample_kit_id), "\n"
)

# --- 2. Aantal personen -----------------------------------------------------
# Belangrijk: stalen van dezelfde persoon zijn niet onafhankelijk
cat("Aantal personen:", n_distinct(metadata$subject.subjectGuid), "\n")

stalen_per_persoon <- metadata |>
  count(subject.subjectGuid, name = "aantal_stalen")
cat("Stalen per persoon:\n")
print(summary(stalen_per_persoon$aantal_stalen))

# --- 3. Aantal eiwit-features ----------------------------------------------
# De eerste kolom is de staal-ID, de rest zijn eiwitten
aantal_eiwitkolommen <- ncol(npx) - 1
cat("\nAantal eiwitkolommen:", aantal_eiwitkolommen, "\n")

# Oude + nieuwe gennaam voor hetzelfde eiwit delen dezelfde olink_id
cat("Aantal unieke olink_id's:", n_distinct(annotatie$olink_id), "\n")

# Sommige eiwitten zijn in meerdere panels gemeten (zelfde assay-naam)
eiwitten_in_meerdere_panels <- annotatie |>
  count(assay, name = "aantal_panels") |>
  filter(aantal_panels > 1)
print(eiwitten_in_meerdere_panels)

# Elk eiwit in meerdere panels telt maar één keer: trek de extra kopieën af
extra_kopieen <- sum(eiwitten_in_meerdere_panels$aantal_panels - 1)
unieke_eiwitten <- n_distinct(annotatie$olink_id) - extra_kopieen
cat("Aantal unieke eiwitten:", unieke_eiwitten, "\n")

# --- 4. Aantal metadata-variabelen -----------------------------------------
cat("\nAantal metadata-kolommen (incl. staal-ID):", ncol(metadata), "\n")

# Deze twee kolommen lijken dubbel; als TRUE mag er één weg
cat(
  "subject_guid gelijk aan subject.subjectGuid:",
  identical(metadata$subject_guid, metadata$subject.subjectGuid), "\n"
)

# --- 5. Type van elke variabele --------------------------------------------
# R kent alleen 'character' en 'numeric'. Of iets nominaal, ordinaal,
# continu of discreet is, beslis je zelf; het aantal unieke waarden helpt:
# weinig unieke tekstwaarden = categorisch, veel kommagetallen = continu.
variabele_types <- tibble(
  variabele = names(metadata),
  r_type = sapply(metadata, class),
  aantal_unieke_waarden = sapply(metadata, n_distinct)
)
print(variabele_types, n = Inf)
