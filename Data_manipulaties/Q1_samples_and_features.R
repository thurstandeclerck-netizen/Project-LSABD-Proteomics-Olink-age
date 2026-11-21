# Question 1: how many samples and features are there, and what type is
# each variable?
# Input:  the three csv files in Data_Proteins/
# Output: numbers in the console

library(tidyverse)

# Read the three datasets
npx <- read_csv(here::here("Data_Proteins", "proteins_npx_raw.csv"))
metadata <- read_csv(here::here("Data_Proteins", "sample_metadata.csv"))
annotation <- read_csv(here::here("Data_Proteins", "protein_annotation.csv"))

# Number of rows and columns (rows = samples, columns = variables)
dim(npx)
dim(metadata)
dim(annotation)

# Number of different people: one person gave several samples
n_distinct(metadata$subject.subjectGuid)

# How many samples does each person have?
samples_per_person <- count(metadata, subject.subjectGuid)
summary(samples_per_person$n)

# Some proteins appear twice (old + new gene name) but share one olink_id
n_distinct(annotation$olink_id)

# Some proteins were measured in more than one panel
count(annotation, assay) |>
  filter(n > 1)

# Type of each variable: <chr> = text, <dbl> = number
glimpse(metadata)
