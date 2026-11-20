# Verwijdert eiwitten van slechte kwaliteit uit de eiwit-annotatie.
# Input:  Data_Proteins/protein_annotation.csv (wordt NIET aangepast)
# Output: Data_proteins_gefilterd/protein_annotation_gefilterd.csv

library(dplyr)
library(readr)

# Drempels op één plek, zodat het team ze makkelijk kan aanpassen
max_pct_missing <- 50 # meer dan de helft van de stalen zonder meting
max_pct_below_lod <- 50 # meer dan de helft van de stalen onder detectielimiet

eiwitten_ruw <- read_csv(
  here::here("Data_Proteins", "protein_annotation.csv"),
  show_col_types = FALSE
)
cat("Start:", nrow(eiwitten_ruw), "eiwitten\n")

# --- Controles zonder treffers (alleen ter info) ---------------------------
# Ontbrekende waarden in de kolommen die we nodig hebben?
cat(
  "NA's in kernkolommen:",
  sum(is.na(select(
    eiwitten_ruw, protein, olink_id, panel,
    npx_median, pct_below_lod, pct_missing
  ))), "\n"
)
# Percentages moeten tussen 0 en 100 liggen
cat(
  "Percentages buiten 0-100:",
  sum(eiwitten_ruw$pct_below_lod < 0 | eiwitten_ruw$pct_below_lod > 100 |
    eiwitten_ruw$pct_missing < 0 | eiwitten_ruw$pct_missing > 100), "\n"
)

# --- Filter 1: te veel ontbrekende metingen --------------------------------
# Dit zijn precies de oude gennamen (bv. CD97 -> ADGRE5) van eiwitten die
# twee keer in de data staan: de oude naam mist in 82% van de stalen.
eiwitten_genoeg_metingen <- eiwitten_ruw |>
  filter(pct_missing <= max_pct_missing)
cat("Na filter 1 (missing):", nrow(eiwitten_genoeg_metingen), "\n")

# Controle: nu mag elke olink_id nog maar één keer voorkomen
cat(
  "Dubbele olink_id's over:",
  sum(duplicated(eiwitten_genoeg_metingen$olink_id)), "\n"
)

# --- Filter 2: te vaak onder de detectielimiet (LOD) -----------------------
# Onder de LOD is het signaal niet te onderscheiden van achtergrondruis
eiwitten_meestal_boven_lod <- eiwitten_genoeg_metingen |>
  filter(pct_below_lod <= max_pct_below_lod)
cat("Na filter 2 (LOD):", nrow(eiwitten_meestal_boven_lod), "\n")

# --- Filter 3: mediaan NPX onder mediaan LOD -------------------------------
# Voor het 'typische' staal ligt het signaal dan onder de ruisgrens, ook al
# zit pct_below_lod net onder 50%
eiwitten_signaal_ok <- eiwitten_meestal_boven_lod |>
  filter(npx_median >= lod_median)
cat("Na filter 3 (mediaan < LOD):", nrow(eiwitten_signaal_ok), "\n")

# --- Filter 4: zelfde eiwit gemeten in meerdere panels ---------------------
# CXCL8, IL6 en TNF zitten in alle 4 de panels. We houden per eiwit (assay)
# de meting met het minst waarden onder de LOD; bij gelijkstand de eerste.
eiwitten_gefilterd <- eiwitten_signaal_ok |>
  group_by(assay) |>
  slice_min(pct_below_lod, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(protein = assay) # naam zonder _OID-achtervoegsel, bv. IL6_OID20101
cat("Na filter 4 (dubbele panels):", nrow(eiwitten_gefilterd), "\n")

# Controle: elk eiwit (assay) en elke olink_id nog maar één keer
cat(
  "Dubbele assays over:", sum(duplicated(eiwitten_gefilterd$assay)),
  "| dubbele olink_id's over:", sum(duplicated(eiwitten_gefilterd$olink_id)),
  "\n"
)

# --- Resultaat bewaren in een aparte map, de ruwe data blijft onaangeroerd -
dir.create(here::here("Data_proteins_gefilterd"), showWarnings = FALSE)
write_csv(
  eiwitten_gefilterd,
  here::here("Data_proteins_gefilterd", "protein_annotation_gefilterd.csv")
)

# Overzicht: hoeveel eiwitten blijven er per panel over?
print(count(eiwitten_gefilterd, panel))
