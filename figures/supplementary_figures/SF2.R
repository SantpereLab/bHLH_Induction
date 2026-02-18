
################ NON-CANONICAL EBOXES VS CANONICAL EBOXES DISTRIBUTION ############
####################################################################################
######## Supplemetary figure 2


### sumarizar
system("
  
  cd ~/proneural/chip-seq

  echo '
  # Load the input fixed noncanonical_eboxes RDS file
  args <- commandArgs(trailingOnly = TRUE)
  input_file <- args[1]

  library(dplyr) 

  study <- readRDS(input_file)
  study$ebox_center <- study$ebox_start + 3
  study$dist_summit <- study$summit - study$ebox_center

  # Create a new column for the 5bp window
  study <- study %>%
    mutate(window = cut(dist_summit, breaks = seq(-400, 400, by = 5), include.lowest = TRUE, labels = seq(-400, 395, by = 5)))

  # Count the number of eboxes in each window, grouped by central_dinucleotide
  summary <- study %>% group_by(window, central_dinucleotide) %>% summarise(count = n()) %>% ungroup()

  # Ensure all windows and central_dinucleotides are represented, even those without eboxes
  all_windows <- expand.grid(window = seq(-400, 395, by = 5), central_dinucleotide = unique(study$central_dinucleotide))

  summary <- summary %>%
    mutate(window = as.numeric(as.character(window))) %>%
    full_join(all_windows, by = c("window", "central_dinucleotide")) %>%
    mutate(count = ifelse(is.na(count), 0, count))

  saveRDS(summary, file = paste0("summary_", input_file))

  print(summary)
  ' > summarise_by_dinucleotide.R


  module load R
  
  cd human
  for file in neikes_mycmax_mcf7_treatment_1000nm_eboxes_canonical.RDS
  do
  sbatch --wrap="Rscript ../summarise_by_dinucleotide.R $file"
  done

  for file in neikes_mycmax_mcf7_treatment_1000nm_eboxes_noncanonical.RDS
  do
  sbatch --wrap="Rscript ../summarise_by_dinucleotide.R $file"
  done

  cd ../mouse

  for file in *canonical.RDS
  do
  sbatch --wrap="Rscript ../summarise_by_dinucleotide.R $file"
  done

  for file in *noncanonical.RDS
  do
  sbatch --wrap="Rscript ../summarise_by_dinucleotide.R $file"
  done

")



system("
  srun --partition=bigmem -n 1 --mpi=none --nodes=1 --pty bash -i
  module load R
  R
")


setwd("~/proneural/chip-seq")




joined <- data.frame()
for (study in studies)  {
  print(study)
  canonical_file <- c(
  list.files(path = "~/proneural/chip-seq/human", pattern = paste0("summary_", study, "_eboxes_canonical.RDS$"), full.names = TRUE),
  list.files(path = "~/proneural/chip-seq/mouse", pattern = paste0("summary_", study, "_eboxes_canonical.RDS$"), full.names = TRUE)
  )

  noncanonical_file <- c(
  list.files(path = "~/proneural/chip-seq/human", pattern = paste0("summary_", study, "_eboxes_noncanonical.RDS$"), full.names = TRUE),
  list.files(path = "~/proneural/chip-seq/mouse", pattern = paste0("summary_", study, "_eboxes_noncanonical.RDS$"), full.names = TRUE)
  )

  tryCatch({
  canon <- readRDS(canonical_file)
  noncanon <- readRDS( unique(noncanonical_file) )

  canon$ebox <- "canonical"
  noncanon$ebox <- "noncanonical"
  bound <- rbind(canon, noncanon)

  bound$study <- study

  joined <- rbind(joined, bound)
  }, error = function(e) {
  message("Error reading RDS files for study: ", study)
  })
}

library(dplyr)
joined_neikes <- joined %>% distinct()
joined<-rbind(joined, joined_neikes)

saveRDS(joined, "summary_canonical_noncanonical_by_dinucleotide.RDS")

library(dplyr)
setwd("~/proneural/chip-seq")
joined <- readRDS("summary_canonical_noncanonical_by_dinucleotide.RDS")

# Load the study mapping
study_mapping <- readRDS("study_mapping.RDS")
# Replace the studies column with the renamed studies
joined <- joined %>%
  left_join(study_mapping, by = c("study" = "study")) %>%
  mutate(study = renamed_study) %>%
  select(-renamed_study)


studies <- readRDS("filtered_studies.RDS")


joined <- joined %>% filter(study %in% studies)

library(dplyr)
library(ggridges)
library(ggplot2)

# Adjust the window size to 20bp and calculate the mean
joined <- joined %>%
  mutate(window = cut(window, breaks = seq(-400, 400, by = 20), include.lowest = TRUE, labels = seq(-400, 380, by = 20))) %>%
  group_by(study, ebox, central_dinucleotide, window) %>%
  summarise(count = mean(count, na.rm = TRUE), .groups = "drop")

# Normalize the count of motifs within each study with respect to the maximum
joined <- joined %>%
  group_by(study) %>%
  mutate(normalized_count = count / max(count, na.rm = TRUE)) %>%
  ungroup()

# Filter the data to include only distances less than 200
# Plot the ridge plot with canonical distribution on top
# Reorder the studies according to the specified order
orden <- readRDS("orden.RDS")
joined <- joined %>%
  mutate(study = factor(study, levels = unique(orden)))

joined <- joined %>% filter(study %in% studies)

joined$central_dinucleotide<-factor(joined$central_dinucleotide, levels = c("CAT-CAT","CAT-CAG","CAG-CAG","CAG-CAC","CAC-CAC","CAT-CAC","CAA-CAT","CAA-CAG","CAA-CAC","CAA-CAA") )

# Remove rows with NA central_dinucleotide
joined <- joined %>% filter(!is.na(central_dinucleotide))


pdf("~/Dropbox/induction/1_motifs_qc/canonical_vs_noncanonical_distribution_by_dinucleotide.pdf", height = 100, width = 30)
ggplot(joined, aes(x = as.numeric(as.character(window)), y = normalized_count, group = interaction(study, ebox, central_dinucleotide), fill = interaction(ebox, central_dinucleotide))) +
  geom_area(alpha = 0.7, position = "identity", color = "black") +
  facet_grid(study ~ central_dinucleotide, scales = "fixed") +
  scale_fill_manual(values = c(
    "canonical.CAT-CAT" = "black",  "canonical.CAT-CAG" = "firebrick2", 
    "canonical.CAG-CAG" = "blue",  
    "canonical.CAG-CAC" = "chartreuse1", 
    "canonical.CAC-CAC" = "goldenrod2", 
    "canonical.CAT-CAC" = "darkorange1", 
    "canonical.CAA-CAT" = "hotpink1", 
    "canonical.CAA-CAG" = "azure4", 
    "canonical.CAA-CAC" = "seagreen", 
    "canonical.CAA-CAA" = "yellow2",
    "noncanonical.CAT-CAT" = "grey", 
    "noncanonical.CAT-CAG" = "grey", "noncanonical.CAG-CAG" = "grey", "noncanonical.CAG-CAC" = "grey",   "noncanonical.CAC-CAC" = "grey", "noncanonical.CAT-CAC" = "grey", 
    "noncanonical.CAA-CAT" = "grey", 
    "noncanonical.CAA-CAG" = "grey",  "noncanonical.CAA-CAC" = "grey",  "noncanonical.CAA-CAA" = "grey"
  )) +
  theme_void() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text.y = element_text(angle = 0, size = 10),
    panel.spacing.y = unit(0.5, "lines"),
    legend.position = "bottom"
  ) +
  labs(
    title = "Canonical vs Non-Canonical E-boxes by Central Dinucleotide (Area)",
    x = "Distance from Summit (20bp bins)",
    y = "",
    fill = "E-box Type"
  )
dev.off()


