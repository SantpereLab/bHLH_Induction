######## CHROMATIN ACCESSIBILITY IN PEAKS ####################
############################################################
#### Main figure 2a



# hacer hetamaps que sean como ribbons. de la accesibilidad en ventanas pequeñitas. y luego encima de eso poner los cuadritos de los cuantiles

setwd('~/proneural/chip-seq')
# Read all files that start with "chromatin"


filtered_studies <- readRDS("filtered_studies.RDS")
study_mapping <- readRDS("study_mapping.RDS")
filtered_studies_renamed<-study_mapping %>% filter(renamed_study %in% filtered_studies) %>% pull(study)

chromatin_files <- list.files(pattern = "^chromatin_.*1000bp_around_summit\\.bed$")
chromatin_files<- paste0("chromatin_", filtered_studies_renamed, "_confidence_1000bp_around_summit.bed") 

# Initialize an empty list to store mean accessibilities
mean_accessibilities <- list()

# Loop through each file and process it
for (file in chromatin_files) {
  print(file)
  # Read the file
  data <- read.table(file, header = FALSE)
  
  # Sort by the max accessibility, assuming its in the 5th column
  data <- data[order(data$V5, decreasing = TRUE), ]
  
  # Divide into 271 groups. Because is the minumum number of peaks of a study
  groups <- split(data, cut(seq(nrow(data)), 271, labels = FALSE))
  
  # Calculate the mean accessibility for each group
  mean_accessibility <- sapply(groups, function(group) mean(group$V5))
  
  # Store the result in the list
  mean_accessibilities[[file]] <- mean_accessibility
}



# Combine the results into a data frame
mean_accessibilities_df <- do.call(cbind, mean_accessibilities)
colnames(mean_accessibilities_df) <- chromatin_files

# Print the result
print(mean_accessibilities_df)

# Convert the data frame to a matrix
heatmap_matrix <- as.matrix(mean_accessibilities_df)


# Apply log transformation and handle -Inf values
log_transformed_matrix <- log1p(heatmap_matrix)
log_transformed_matrix[is.infinite(log_transformed_matrix)] <- NA

# Replace NA values with the minimum value in the matrix
min_value <- min(log_transformed_matrix, na.rm = TRUE)
log_transformed_matrix[is.na(log_transformed_matrix)] <- min_value

# Print the log-transformed matrix

# Clean up column names
colnames(log_transformed_matrix) <- gsub("chromatin_|_confidence_1000bp_around_summit.bed", "", colnames(log_transformed_matrix))

orden<-readRDS("orden.RDS")


# Load the study mapping
study_mapping <- readRDS("study_mapping.RDS")

# Rename the columns of log_transformed_matrix based on study_mapping
colnames(log_transformed_matrix) <- study_mapping[match(colnames(log_transformed_matrix), study_mapping$study), "renamed_study"]


orden<-orden[which(orden %in% colnames(log_transformed_matrix))]
log_transformed_matrix <-log_transformed_matrix[,orden]

library(pheatmap)
library(viridis)
library(grid)



# Generate the heatmap using pheatmap without column names
pdf("~/Dropbox/induction/2_pre_induction_chromatin/rpkms_capped_in_peaks_quantiles_divided.pdf", height = 15, width = 9, useDingbats = FALSE)
  # Cap the matrix values at 5
  log_transformed_matrix_cap <- log_transformed_matrix
  log_transformed_matrix_cap[log_transformed_matrix_cap > 5] <- 5

  # Draw heatmap without borders
  p <- pheatmap(
    t(log_transformed_matrix_cap),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = viridis::mako(256),
    scale = "none",
    gaps_col = c(68, 136, 204),
    show_colnames = FALSE,
    border_color = NA,
    main = "Maximum accessibility in regions, measured as log (RPKMS +1)",
    silent = TRUE
  )

  # Overlay horizontal white separators robustly via gtable (no viewport assumptions)
  library(grid)
  library(gtable)

  g <- p$gtable
  mat_idx <- which(g$layout$name == "matrix")
  nr <- nrow(t(log_transformed_matrix_cap))

  if (length(mat_idx) == 1 && nr > 1) {
    ys <- 1 - (1:(nr - 1)) / nr
    segs <- lapply(ys, function(y)
      segmentsGrob(x0 = unit(0, "npc"), x1 = unit(1, "npc"),
                   y0 = unit(y, "npc"), y1 = unit(y, "npc"),
                   gp = gpar(col = "white", lwd = 0.5))
    )
    g <- gtable_add_grob(
      g, grobs = segs,
      t = g$layout$t[mat_idx], l = g$layout$l[mat_idx],
      b = g$layout$b[mat_idx], r = g$layout$r[mat_idx],
      z = Inf, name = paste0("row_sep_", seq_along(segs))
    )
  } else {
    # Fallback: draw borders if matrix grob not found
    p <- pheatmap(
      t(log_transformed_matrix_cap),
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      color = viridis::mako(256),
      scale = "none",
      gaps_col = c(68, 136, 204),
      show_colnames = FALSE,
      border_color = "white",
      main = "Maximum accessibility in regions, measured as log (RPKMS +1)",
      silent = TRUE
    )
    g <- p$gtable
  }

  grid.newpage()
  grid.draw(g)
dev.off()



################################################################
######## Main figure 2b ########################################
################################################################


setwd("~/proneural/chip-seq")

read_half <- function(studies, direction = c("upstream", "downstream")) {
  direction <- match.arg(direction)
  dfs <- lapply(studies, function(study) {
    paths <- file.path(
      c("~/proneural/chip-seq/human", "~/proneural/chip-seq/mouse"),
      paste0(study, "_confidence_", direction, "_half_nochrM_eboxes.RDS")
    )
    existing <- paths[file.exists(paths)]
    if (!length(existing)) {
      message("No files found for study: ", study)
      return(NULL)
    }
    df <- do.call(rbind, lapply(existing, readRDS))
    df$study <- study
    df
  })
  dfs <- Filter(Negate(is.null), dfs)
  if (length(dfs)) do.call(rbind, dfs) else data.frame()
}

upstream_joined_studies  <- read_half(studies, "upstream")
downstream_joined_studies <- read_half(studies, "downstream")

## convert region column. from SRX3231473.05_peak_1_DOWNHALF to peak_1
upstream_joined_studies$region <- sub(".*_(peak_\\d+)_UPHALF", "\\1", upstream_joined_studies$region)
downstream_joined_studies$region <- sub(".*_(peak_\\d+)_DOWNHALF", "\\1", downstream_joined_studies$region)

## convert ebox_seq to central_dinucleotide. for example, CATATG to CAT-CAT. 

# Map motif hexamers to central dinucleotide names
motif_map <- c(
  "CATATG" = "CAT-CAT",
  "CAGCTG" = "CAG-CAG",
  "CACGTG" = "CAC-CAC",
  "CAATTG" = "CAA-CAA",
  "CAGATG" = "CAT-CAG", "CATCTG" = "CAT-CAG",
  "CACATG" = "CAT-CAC", "CATGTG" = "CAT-CAC",
  "CAAATG" = "CAA-CAT", "CATTTG" = "CAA-CAT",
  "CAGGTG" = "CAG-CAC", "CACCTG" = "CAG-CAC",
  "CAACTG" = "CAA-CAG", "CAGTTG" = "CAA-CAG",
  "CAAGTG" = "CAA-CAC", "CACTTG" = "CAA-CAC"
)
upstream_joined_studies$central_dinucleotide <- motif_map[upstream_joined_studies$ebox_seq]
downstream_joined_studies$central_dinucleotide <- motif_map[downstream_joined_studies$ebox_seq]

flanking_joined_studies <- rbind(
  upstream_joined_studies ,
  downstream_joined_studies
)

#### count the number of ocurrences of each central_dinucleotide per peak. 
library(dplyr)
flanking_counts <- flanking_joined_studies %>%
  filter(!is.na(central_dinucleotide)) %>%
  group_by(study, region, central_dinucleotide) %>%
  summarise(count = n(), .groups = "drop")


flanking_counts$study<-study_mapping[match(flanking_counts$study, study_mapping$study), "renamed_study"]

### rename columns
flanking_counts<- flanking_counts %>%
  rename(peak_id = region, flanking_count = count)

### join flanking_counts and dinucleotide_counts by study, peak_id and central_dinucleotide
### dinucleotide_counts has all the fields. if a combination not present in flanking_counts, fill with 0
peak_and_flanking_counts <- dinucleotide_counts %>%
  left_join(flanking_counts, by = c("study", "peak_id", "central_dinucleotide")) %>%
  mutate(flanking_count = tidyr::replace_na(flanking_count, 0L)) %>%
  dplyr::select(study, peak_id, central_dinucleotide, count, flanking_count, peak_length, log_accessibility)

#### for each peak, and for each central_dinucleotide (ebox variant), create an index of motif enrichment over background
## for that, first add a pseudocount of +1 to the dinucleotide count and flanking count.
## then do the division, and then the log.
### log( (count+1) / (flanking count+1) )
peak_and_flanking_counts <- peak_and_flanking_counts %>%
  mutate(
    count_plus1 = count + 1,
    flanking_count_plus1 = flanking_count + 1,
    motif_enrichment = log(count_plus1 / flanking_count_plus1)
  )

saveRDS(peak_and_flanking_counts, "~/peak_and_flanking_counts.RDS")


peak_and_flanking_counts<- readRDS("~/peak_and_flanking_counts.RDS")


###########################
######## LM ###############
###########################

# Multivariate fits per study: log_accessibility ~ motif_enrichment for all motifs (one predictor per motif)
suppressPackageStartupMessages({
  if (!requireNamespace("broom", quietly = TRUE)) install.packages("broom")
  library(dplyr)
  library(tidyr)
  library(broom)
})


motif_levels <- c("CAT-CAT","CAT-CAG","CAG-CAG","CAG-CAC","CAC-CAC",
                  "CAT-CAC","CAA-CAT","CAA-CAG","CAA-CAC","CAA-CAA")


# Build peak-level wide design: one row per (study, peak_id) with columns motif_enrichment_<motif>
design_wide <- peak_and_flanking_counts %>%
  select(study, peak_id, central_dinucleotide, motif_enrichment, log_accessibility, peak_length) %>%
  distinct() %>%
  mutate(central_dinucleotide = factor(central_dinucleotide, levels = motif_levels)) %>%
  pivot_wider(
    names_from = central_dinucleotide,
    values_from = motif_enrichment,
    names_prefix = "motif_enrichment_"
  )





# Fit per-study multivariate LM and extract coefficients for each motif
count_fits <- dplyr::bind_rows(lapply(split(design_wide, design_wide$study), function(df_study) {
  df_study <- df_study %>% dplyr::filter(!is.na(log_accessibility))
  if (nrow(df_study) < 5) return(NULL)


  predictors_all <- paste0("motif_enrichment_", motif_levels)
  preds <- intersect(predictors_all, colnames(df_study))
  if (!length(preds)) return(NULL)


  # Drop zero-variance predictors and replace NA with 0
  preds <- preds[vapply(df_study[preds], function(x) stats::sd(x, na.rm = TRUE) > 0, logical(1))]
  if (!length(preds)) return(NULL)
  df_study[preds] <- lapply(df_study[preds], function(x) { x[is.na(x)] <- 0; x })


  # Backtick predictors with '-' so the formula parses correctly
  rhs <- paste(sprintf("`%s`", preds), collapse = " + ")
  fit <- tryCatch(lm(as.formula(paste("log_accessibility ~", rhs)),
                     data = df_study, na.action = na.exclude),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)


  td <- broom::tidy(fit)
  # Keep only motif terms (handle backticked names)
  keep_terms <- c(preds, sprintf("`%s`", preds))
  td <- td[td$term %in% keep_terms, , drop = FALSE]
  if (!nrow(td)) return(NULL)


  term_clean <- gsub("^`|`$", "", td$term)
  data.frame(
    study = unique(df_study$study),
    central_dinucleotide = sub("^motif_enrichment_", "", term_clean),
    slope = td$estimate,
    p_value = td$p.value,
    stringsAsFactors = FALSE
  )
}))



# P-values with BH correction across all comparisons (multivariate)
dinucleotide_pvals <- count_fits %>%
  dplyr::select(study, central_dinucleotide, slope, p_value) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(q_count = p.adjust(p_value, method = "BH"))


## combine motif count slopes and centrality slopes to do balloon plot. balloon size counts and color centrality

dinucleotide_pvals$slope<- -dinucleotide_pvals$slope






############################### HEATMAPS

# Two heatmaps: one for centrality slope, one for count slope (nmotifs)
# Overlay significance (q < 0.05) as an asterisk

library(ggplot2)
library(dplyr)
library(tidyr)
library(grid)
library(pheatmap)
library(RColorBrewer)

setwd("~/proneural/chip-seq")
df <- dinucleotide_pvals
orden <- readRDS("orden.RDS")
orden<- orden[which(orden %in% unique(df$study))]
motif_levels <- c("CAT-CAT","CAT-CAG","CAG-CAG","CAG-CAC","CAC-CAC","CAT-CAC","CAA-CAT","CAA-CAG","CAA-CAC","CAA-CAA")

# Prepare matrices for heatmaps
df <- df %>%
  mutate(
    study = factor(study, levels = orden),
    central_dinucleotide = factor(central_dinucleotide, levels = motif_levels)
  )


  # Count slope matrix (use dinucleotide_pvals columns: slope, q_count)
  count_mat <- df %>%
    select(study, central_dinucleotide, slope) %>%
    pivot_wider(names_from = central_dinucleotide, values_from = slope) %>%
    as.data.frame()
  rownames(count_mat) <- count_mat$study
  count_mat <- as.matrix(count_mat[,-1])
  count_mat <- count_mat[orden, motif_levels]

  # Significance matrix from q_count in dinucleotide_pvals
  sig_count <- df %>%
    mutate(sig = ifelse(!is.na(q_count) & q_count < 0.01, "*", "")) %>%
    select(study, central_dinucleotide, sig) %>%
    pivot_wider(names_from = central_dinucleotide, values_from = sig) %>%
    as.data.frame()
  rownames(sig_count) <- sig_count$study
  sig_count <- as.matrix(sig_count[,-1])
  sig_count <- sig_count[orden, motif_levels]

  # Color palette: diverging with white at 0
  make_diverging <- function(mat, n = 200,
                             col_min = "#2166AC",
                             col_mid = "white",
                             col_max = "#B2182B") {
    rng <- range(mat, na.rm = TRUE)
    min_v <- rng[1]; max_v <- rng[2]
    if (min_v == max_v) {
      return(list(colors = rep(col_mid, n), breaks = seq(min_v, max_v, length.out = n + 1)))
    }
    if (min_v >= 0 || max_v <= 0) {
      cols <- colorRampPalette(c(col_min, col_mid, col_max))(n)
      return(list(colors = cols, breaks = seq(min_v, max_v, length.out = n + 1)))
    }
    frac_neg <- abs(min_v) / (abs(min_v) + max_v)
    n_neg <- max(5, floor(n * frac_neg))
    n_pos <- n - n_neg
    neg_cols_full <- colorRampPalette(c(col_min, col_mid))(n_neg + 1)
    neg_cols <- neg_cols_full[1:n_neg]
    pos_cols_full <- colorRampPalette(c(col_mid, col_max))(n_pos + 1)
    pos_cols <- pos_cols_full[2:(n_pos + 1)]
    colors <- c(neg_cols, pos_cols)
    breaks <- c(seq(min_v, 0, length.out = n_neg + 1),
                seq(0, max_v, length.out = n_pos + 1)[-1])
    list(colors = colors, breaks = breaks)
  }



library(pheatmap )
# Count slope heatmap
count_spec <- make_diverging(count_mat, n = 200)
pdf("~/Dropbox/induction/2_pre_induction_chromatin/no_bins/multivariate/accessibility~motif_enrichment_0.01.pdf",
  height = 15, width = 9, useDingbats = FALSE)
pheatmap(
  count_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = count_spec$colors,
  breaks = count_spec$breaks,
  main = "Slope of Motif Count (per study/motif)",
  display_numbers = sig_count,
  number_color = "black",
  fontsize_number = 14,
  fontsize_row = 8,
  fontsize_col = 8,
  scale = "none"
)
dev.off()






############ DISTANCE BETWEEN EBOXES IN QUANTILES #######
#########################################################
##### Main figure 2c ####################################


setwd("~/proneural/chip-seq")
combined_data<-readRDS("combined_data.RDS")
studies<-readRDS("filtered_studies.RDS")

# Rename the studies based on a mapping file
study_mapping <- readRDS("study_mapping.RDS")

# Update the study names in the combined_data
combined_data <- combined_data %>%
  mutate(study = study_mapping[match(study, study_mapping$study), "renamed_study"])

combined_data<-combined_data %>% filter(study %in% studies)

library(dplyr)
combined_data$ebox_center<-combined_data$ebox_start+3
combined_data$dist_summit<-combined_data$ebox_center-combined_data$summit



library(dplyr)
combined_data$window_start<-combined_data$summit-as.integer(combined_data$peak_length/2)
combined_data$window_end<-combined_data$summit+as.integer(combined_data$peak_length/2)
combined_data$ebox_center<-combined_data$ebox_start+3


# Filter motifs that are within the window or keep rows with NA in ebox_center
filtered_data <- combined_data %>%
  filter(is.na(ebox_center) | (ebox_center >= window_start & ebox_center <= window_end))


# Convert quantiles into groups 1 to 4
combined_data_quantiles <- combined_data %>%
  mutate(quantile = case_when(
    quantile <= 25 ~ 1,
    quantile > 25 & quantile <= 50 ~ 2,
    quantile > 50 & quantile <= 75 ~ 3,
    quantile > 75 ~ 4
))

# Calculate distances between eboxes within each peak for each quantile
distances <- combined_data_quantiles %>%
  group_by(study, peak_id, quantile) %>%
  arrange(ebox_center) %>%
  mutate(distance_to_next = lead(ebox_center) - ebox_center) %>%
  filter(!is.na(distance_to_next)) %>%
  ungroup()

distances$distance_to_next <- distances$distance_to_next - 6

saveRDS(distances, file = "distance_summary_between_eboxes_quantiles.RDS")


distances<-readRDS("distance_summary_between_eboxes_quantiles.RDS")

studies<-readRDS("filtered_studies.RDS")
distances<-distances %>% filter(study %in% studies)


library(dplyr)
# Filter distances to include only those between 1 and 15
filtered_distances <- distances %>%
  filter(distance_to_next >= 1 & distance_to_next <= 15)

# Calculate the proportion of each distance with respect to the rest of distances within each study and quantile
distance_proportions <- filtered_distances %>%
  group_by(study, quantile) %>%
  mutate(total_count = n()) %>%
  ungroup() %>%
  group_by(study, quantile, distance_to_next) %>%
  summarise(count = n(), total_count = first(total_count), .groups = "drop") %>%
  mutate(proportion = count / total_count)

# Save the summary to an RDS file
saveRDS(distance_proportions, file = "distance_proportions_summary_quantiles.RDS")



# probar a hacer solo con las proporciones de la distancia de 5 y 6


setwd("~/proneural/chip-seq")
library(dplyr)
distance_proportions<-readRDS("distance_proportions_summary_quantiles.RDS")

# Rename the studies based on a mapping file
study_mapping <- readRDS("study_mapping.RDS")

# Update the study names in the distance_proportions
distance_proportions <- distance_proportions %>%
  mutate(study = study_mapping[match(study, study_mapping$study), "renamed_study"])

studies<-readRDS("filtered_studies.RDS")
distance_proportions<-distance_proportions %>% filter(study %in% studies)

library(dplyr)
distance_proportions_5<-distance_proportions %>% filter(distance_to_next==5)

library(ggplot2)

distance_proportions_5$quantile<-5-distance_proportions_5$quantile
# Convert quantile to a factor for better plotting
distance_proportions_5$quantile <- factor(distance_proportions_5$quantile, levels = c(1, 2, 3, 4)) 


# Reshape the data for a heatmap
heatmap_data <- distance_proportions_5 %>%
  select(study, quantile, proportion) %>%
  tidyr::pivot_wider(names_from = quantile, values_from = proportion, values_fill = 0)

# Convert to matrix for heatmap
heatmap_matrix <- as.matrix(heatmap_data[,-1])
rownames(heatmap_matrix) <- heatmap_data$study

orden<-readRDS("orden.RDS")
head(heatmap_matrix)


orden<-orden[which(orden %in% rownames(heatmap_matrix))]
heatmap_matrix <-heatmap_matrix[orden,rev(c(1,2,3,4))]

# Save the plot as a PDF
pdf("~/Dropbox/induction/2_pre_induction_chromatin/proportion_distance_5_across_quantiles.pdf", height = 10, width = 4, useDingbats = FALSE)

# Generate the heatmap
library(pheatmap)
pheatmap(
  heatmap_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = viridis::rocket(100),
  main = "Proportion of Distance 6 Across Quantiles",
  fontsize_row = 8,
  fontsize_col = 8,
  fontsize = 4
)
dev.off()

## the same plot but with proportion of 6bp


# Proportion of 6bp distance across quantiles
distance_proportions_6 <- distance_proportions %>% filter(distance_to_next == 6)

distance_proportions_6$quantile <- 5 - distance_proportions_6$quantile
distance_proportions_6$quantile <- factor(distance_proportions_6$quantile, levels = c(1, 2, 3, 4))

# Reshape the data for a heatmap
heatmap_data_6 <- distance_proportions_6 %>%
  select(study, quantile, proportion) %>%
  tidyr::pivot_wider(names_from = quantile, values_from = proportion, values_fill = 0)

# Convert to matrix for heatmap
heatmap_matrix_6 <- as.matrix(heatmap_data_6[,-1])
rownames(heatmap_matrix_6) <- heatmap_data_6$study

orden <- readRDS("orden.RDS")
orden <- orden[which(orden %in% rownames(heatmap_matrix_6))]
heatmap_matrix_6 <- heatmap_matrix_6[orden, rev(c(1,2,3,4))]

library(pheatmap)
pdf("~/Dropbox/induction/2_pre_induction_chromatin/proportion_distance_6_across_quantiles.pdf", height = 10, width = 4, useDingbats = FALSE)
pheatmap(
  heatmap_matrix_6,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = viridis::rocket(100),
  main = "Proportion of Distance 7 Across Quantiles",
  fontsize_row = 8,
  fontsize_col = 8,
  fontsize = 4
)
dev.off()



# Proportion of 9bp distance across quantiles
distance_proportions_9 <- distance_proportions %>% filter(distance_to_next == 9)

distance_proportions_9$quantile <- 5 - distance_proportions_9$quantile
distance_proportions_9$quantile <- factor(distance_proportions_9$quantile, levels = c(1, 2, 3, 4))

# Reshape the data for a heatmap
heatmap_data_9 <- distance_proportions_9 %>%
  select(study, quantile, proportion) %>%
  tidyr::pivot_wider(names_from = quantile, values_from = proportion, values_fill = 0)

# Convert to matrix for heatmap
heatmap_matrix_9 <- as.matrix(heatmap_data_9[,-1])
rownames(heatmap_matrix_9) <- heatmap_data_9$study

orden <- readRDS("orden.RDS")
orden <- orden[which(orden %in% rownames(heatmap_matrix_9))]
heatmap_matrix_9 <- heatmap_matrix_9[orden, rev(c(1,2,3,4))]

library(pheatmap)
pdf("~/Dropbox/induction/2_pre_induction_chromatin/proportion_distance_9_across_quantiles.pdf", height = 10, width = 4, useDingbats = FALSE)
pheatmap(
  heatmap_matrix_9,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = viridis::rocket(100),
  main = "Proportion of Distance 9 Across Quantiles",
  fontsize_row = 8,
  fontsize_col = 8,
  fontsize = 4
)
dev.off()


# Proportion of 10bp distance across quantiles
distance_proportions_10 <- distance_proportions %>% filter(distance_to_next == 10)

distance_proportions_10$quantile <- 5 - distance_proportions_10$quantile
distance_proportions_10$quantile <- factor(distance_proportions_10$quantile, levels = c(1, 2, 3, 4))

# Reshape the data for a heatmap
heatmap_data_10 <- distance_proportions_10 %>%
  select(study, quantile, proportion) %>%
  tidyr::pivot_wider(names_from = quantile, values_from = proportion, values_fill = 0)

# Convert to matrix for heatmap
heatmap_matrix_10 <- as.matrix(heatmap_data_10[,-1])
rownames(heatmap_matrix_10) <- heatmap_data_10$study

orden <- readRDS("orden.RDS")
orden <- orden[which(orden %in% rownames(heatmap_matrix_10))]
heatmap_matrix_10 <- heatmap_matrix_10[orden, rev(c(1,2,3,4))]

library(pheatmap)
pdf("~/Dropbox/induction/2_pre_induction_chromatin/proportion_distance_10_across_quantiles.pdf", height = 10, width = 4, useDingbats = FALSE)
pheatmap(
  heatmap_matrix_10,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = viridis::rocket(100),
  main = "Proportion of Distance 10 Across Quantiles",
  fontsize_row = 8,
  fontsize_col = 8,
  fontsize = 4
)
dev.off()











