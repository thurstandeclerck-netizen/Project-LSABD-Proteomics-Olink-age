# Question 5: do samples cluster by biology or by a technical factor?
# Input:  npx_corrected.csv and metadata_final.csv in output/ (from Q4)
# Output: numbers in the console + 3 figures in output/

library(tidyverse)

npx <- read_csv(here::here("output", "npx_corrected.csv"))
metadata <- read_csv(here::here("output", "metadata_final.csv"))


# --- Step 1: PCA on the batch-corrected data --------------------------------
npx_values <- select(npx, -sample_kit_id)
npx_complete <- npx_values[, colSums(is.na(npx_values)) == 0]
pca <- prcomp(npx_complete, scale. = TRUE)
summary(pca)$importance[2, 1:5] # share of variation per axis

# Add the first 5 PCA axes to the metadata
pca_scores <- bind_cols(metadata, as_tibble(pca$x[, 1:5]))


# --- Step 2: which factor explains an axis? ---------------------------------
# R-squared of lm() = share of an axis explained by a factor (0 to 1)
summary(lm(PC1 ~ olink_batch_id, data = pca_scores))$r.squared # technical
summary(lm(PC2 ~ subject.ageGroup, data = pca_scores))$r.squared # age
summary(lm(PC4 ~ subject.biologicalSex, data = pca_scores))$r.squared # sex
summary(lm(PC2 ~ cmv_status, data = pca_scores))$r.squared # CMV
summary(lm(PC2 ~ sample.visitName, data = pca_scores))$r.squared # visit
summary(lm(PC2 ~ subject.subjectGuid, data = pca_scores))$r.squared # person


# --- Step 3: figures --------------------------------------------------------
# Age group: young and old are separated along PC2
ggplot(pca_scores, aes(x = PC1, y = PC2, colour = subject.ageGroup)) +
  geom_point(alpha = 0.6) +
  labs(title = "PCA after batch correction: age group", colour = "Age group")
ggsave(here::here("output", "q5_pca_age_group.png"))

# Sex: men and women are separated along PC4
ggplot(pca_scores, aes(x = PC2, y = PC4, colour = subject.biologicalSex)) +
  geom_point(alpha = 0.6) +
  labs(title = "PCA after batch correction: sex", colour = "Sex")
ggsave(here::here("output", "q5_pca_sex.png"))

# Person: samples of the same person lie close together (5 people as example)
five_people <- c("BR1001", "BR1002", "BR1003", "BR2004", "BR2005")
ggplot(
  filter(pca_scores, subject.subjectGuid %in% five_people),
  aes(x = PC2, y = PC3, colour = subject.subjectGuid)
) +
  geom_point(size = 3) +
  labs(
    title = "Samples of the same person tend to lie close together",
    colour = "Person"
  )
ggsave(here::here("output", "q5_pca_person.png"))
