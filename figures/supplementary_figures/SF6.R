#############################  Supplementary figure 6 a ############################
###################################################################################




library(dplyr)

# lee_ascl1_48h no because is reanalisis of wapinski
fc_studies<-c("aydin_ascl1_12h","aydin_neurog2_48h","barfeld_myc","casey_ascl1_24h","casey_ascl2_24h","casey_myod1_24h","conerly_myf5","conerly_myod","fong_neurod2_induction","lee_myod1_48h","li_tw2_gm_induction","li_tw2_induction","lin_myc_t24","manandhar_myod1","park_ascl1_dox","pereira_astrocytes_ngn2","smith_neurog2_induction","walz_myc_dox","wang_ascl1","wapinski_ascl1_48hr","zeid_mycn_6hr")

setwd("~/proneural")
i<-1
for (study in fc_studies) {
  print(study)
  filee<-paste0("intersect_foldchange_induction_",study,"_join.bed")
  tabl<-read.table(filee) %>% select(V4,V5)
  colnames(tabl)<-c("fc","intersect")
  tabl$intersect[which(tabl$intersect==".")]<-"no"
  tabl$intersect[which(tabl$intersect!="no")]<-"yes"
  tabl$study<-study
  
  if (i==1) { foldchanges<-tabl }
  else { foldchanges<-rbind(foldchanges,tabl) }
  i<-i+1
}

library(ggplot2)

foldchanges$log_fc<- log(foldchanges$fc)


# Rename the studies based on a mapping file
study<-unique(foldchanges$study)
study_renamed<-c("ASCL1 -- EB [Aydin et al. 2019]", "NEUROG2 -- EB [Aydin et al. 2019]", "MYC -- LNCaP [Barfeld et al. 2017]", 
                      "ASCL1 -- mESC [Casey et al. 2018]", "ASCL2 -- mESC [Casey et al. 2018]", "MYOD1 -- mESC [Casey et al. 2018]",
                      "MYF5 -- MEF [Conerly et al. 2016]", "MYOD1 -- MEF [Conerly et al. 2016]", "NEUROD2 -- P19 [Fong et al. 2012]",
                      "MYOD1 -- MEF [Lee et al. 2020]",
                      "TWIST2 -- Myoblast (GM) [Li et al. 2019]", "TWIST2 -- Myoblast (DM) [Li et al. 2019]",
                      "MYC -- P493-6 [Lin et al. 2012]", "MYOD1 -- Human fibroblasts [Manandhar et al. 2017]", 
                      "ASCL1 -- G523NS [Park et al. 2017]",
                      "NEUROG2 -- Astrocytes [Pereira et al. 2024]",
                      "NEUROG2 -- MRC5 [Smith et al. 2016]", "MYC -- U2OS [Walz et al. 2014]",
                      "AS2kbCL1 -- GI-MEN [Wang et al. 2023]", "ASCL1 -- MEF [Wapinski et al. 2013]",
                      "MYCN -- SHEP [Zeid et al. 2018]")

study_mapping <- data.frame(study = study, renamed_study = study_renamed)

# Update the study names in the foldchanges data frame
foldchanges <- foldchanges %>%
  mutate(study = study_mapping[match(study, study_mapping$study), "renamed_study"])




library(tidyr)

# Calculate delta of the median fold change for each study, excluding outliers
median_deltas <- foldchanges %>%
  group_by(study, intersect) %>%
  filter(log_fc >= -4 & log_fc <= 4) %>%
  summarize(median_fc = median(log_fc, na.rm = TRUE)) %>%
  spread(intersect, median_fc) %>%
  mutate(delta_median = yes - no)

# Print the median deltas
print(median_deltas)

# Create a barplot of the median deltas, ordered by delta_median
median_deltas <- median_deltas %>%
  arrange(delta_median)

jj <- ggplot(median_deltas, aes(x = reorder(study, delta_median), y = delta_median)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Median Fold Change Deltas by Study",
       x = "Study",
       y = "Delta Median Fold Change")







# Reorder the studies in foldchanges based on the delta_median order
foldchanges$study <- factor(foldchanges$study, levels = median_deltas$study)

# Create the reordered boxplot with reversed order and white background
gigi_reordered <- ggplot(foldchanges, aes(x = forcats::fct_rev(study), y = log_fc, fill = intersect) ) +
  geom_boxplot(outlier.size = 0.01, outlier.color = "grey20") + 
  ylim(-4, 4) +  
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 

# Save the reordered boxplot to a PDF file
pdf(file = "~/Dropbox/induction/5_remodeling/boxplots_foldchange_induction_bound_vs_not_bound.pdf", height = 10, width = 20, useDingbats = F)
  print(gigi_reordered)
dev.off()







##################### Supplementary fig. 6b ######################
##################################################################


setwd("~/proneural")
# Read all *_motif_archetype_median_summary.rds files, only those with intersect == "yes"
library(dplyr)

# List all RDS files matching the pattern in the current directory (or set your path)
rds_files <- list.files(pattern = "_intersect_yes_motif_archetype_median_summary\\.rds$")

# Read and filter for intersect == "yes"
motif_summaries <- lapply(rds_files, function(f) {
  df <- readRDS(f)
})

# Combine into a single data frame
motif_summaries_yes <- bind_rows(motif_summaries)


# Extract quantile and study name from the file column
motif_summaries_yes <- motif_summaries_yes %>%
  mutate(
    quantile = sub(".*_quantil_([0-9])_intersect_yes.*$", "\\1", file),
    study = sub("juancho_combined_vierstra_with_distance_peakid_([a-zA-Z0-9_]+)_quantil_\\d+_intersect_yes.*$", "\\1", file)
  )

motif_summaries_yes$quantile<-as.numeric(motif_summaries_yes$quantile)+1

motif_summaries_yes<-motif_summaries_yes %>% filter(quantile==4)

motif_summaries_yes$centrality<-500-motif_summaries_yes$median_abs_distance

motif_summaries_yes <- motif_summaries_yes %>%
  group_by(study) %>%
  mutate(centrality_scaled = (centrality - min(centrality, na.rm = TRUE)) / (max(centrality, na.rm = TRUE) - min(centrality, na.rm = TRUE))) %>%
  ungroup()

library(ggplot2)
library(tidyr)

# Prepare data for heatmap: rows = motif_name, columns = quantile, values = centrality_scaled
motif_heatmap_data <- motif_summaries_yes %>%
  select(study, motif_name, quantile, centrality_scaled)


mouse<-motif_summaries_yes %>% filter(study %in% c( "wapinski_ascl1_48hr","aydin_ascl1_12h", "aydin_neurog2_48h", "casey_ascl1_24h")  )

shared_archetypes<-unique(mouse$motif_name)

## some motif archetypes non annotated as E-boxes are similar to eboxes. so annotate them like that
# Rename motif_name values in motif_heatmap_data according to the mapping
motif_rename_map <- c(
  "HEN1" = "HEN1(CAGCTG)",
  "ZNF317" = "ZNF317(CAGCTG)",
  "ZNF331" = "ZNF331(CAGCTG)",
  "ZNF563" = "ZNF563(CAGCTG)",
  "SCRT1" = "SCRT1(CAGGTG)",
  "SNAI2" = "SNAI2(CAGGTG)",
  "MIES" = "MIES(CAGGTG)",
  "HIF" = "HIF(NACGTG)",
  "CREB3/XBP1" = "CREB3/XBP1(GACGTG)",
  "GMEB2/2" = "GMEB2/2(TACGTA)",
  "AHR" = "AHR(NGCGTG)",
  "ZNF134" = "ZNF134(CAGTTG)",
  "MYB/5" = "MYB/5(CAGTTG)",
  "OSR2" = "OSR2(CWGCTN)",
  "PRDM9" = "PRDM9(YAGYAN)"
)

motif_heatmap_data$motif_name <- ifelse(
  motif_heatmap_data$motif_name %in% names(motif_rename_map),
  motif_rename_map[motif_heatmap_data$motif_name],
  motif_heatmap_data$motif_name
)

## rename shared archetypes too
shared_archetypes<-ifelse(
  shared_archetypes %in% names(motif_rename_map),
  motif_rename_map[shared_archetypes],
  shared_archetypes
)




# I set the colors based on the central dinucleotide of the E-box
ebox_colors <- c(
  "HEN1(CAGCTG)" = "steelblue1", "ZNF317(CAGCTG)" = "royalblue1", "Ebox/CAGCTG" = "royalblue3", "ZNF331(CAGCTG)" = "royalblue4", 
  "ZNF563(CAGCTG)" = "navyblue", "Ebox/CAGATGG" = "red1", "Ebox/CATATG" = "black", "Ebox/CACCTG" = "darkolivegreen4", 
  "SCRT1(CAGGTG)" = "olivedrab1", "SNAI2(CAGGTG)" = "chartreuse1", "MIES(CAGGTG)" = "chartreuse3", "Ebox/CACGTG/1" = "goldenrod4", 
  "Ebox/CACGTG/2" = "goldenrod3", "HIF(NACGTG)" = "lightgoldenrod1", "CREB3/XBP1(GACGTG)" = "lemonchiffon3", 
  "GMEB2/2(TACGTA)" = "wheat2", "AHR(NGCGTG)" = "khaki2", "ZNF134(CAGTTG)" = "violetred2", "MYB/5(CAGTTG)" = "violet", 
  "OSR2(CWGCTN)" = "purple3", "PRDM9(YAGYAN)" = "plum4" )



# Remove unwanted studies from motif_heatmap_data
motif_heatmap_data <- motif_heatmap_data %>%
  filter(!(study %in% c("manandhar_myod1", "lee_ascl1_48h")))



# Rename the studies based on a mapping file
study<-unique(motif_heatmap_data$study)

print(study)


study_renamed<- c( "ASCL1 -- EB [Aydin et al. 2019]", "NEUROG2 -- EB [Aydin et al. 2019]", 
                      "ASCL1 -- mESC [Casey et al. 2018]", "ASCL2 -- mESC [Casey et al. 2018]", 
                      "MYOD1 -- mESC [Casey et al. 2018]", "MYOD1 -- MEF [Lee et al. 2020]",
                      "ASCL1 -- G523NS [Park et al. 2017]", "NEUROG2 -- Astrocytes [Pereira et al. 2024]",
                      "ASCL1 -- GI-MEN [Wang et al. 2023]", "ASCL1 -- MEF [Wapinski et al. 2013]" )

# rename study in motif_heatmap_data 
study_mapping <- data.frame(study = study, renamed_study = study_renamed)

motif_heatmap_data <- motif_heatmap_data %>%
  mutate(study = study_mapping[match(study, study_mapping$study), "renamed_study"])


studies <- unique(motif_heatmap_data$study)


# Get the top 50 archetype motifs with highest centrality in each study, then combine and get unique motif names.
## archetypes must be shared between mouse and human

top_motifs_all <- motif_heatmap_data %>% filter(motif_name %in% shared_archetypes) %>%
  filter(quantile == 4) %>%
  group_by(study) %>%
  arrange(desc(centrality_scaled)) %>%
  slice_head(n = 50) %>%
  ungroup() %>%
  pull(motif_name) %>%
  unique()

print(top_motifs_all)


motif_heatmap_data<- motif_heatmap_data %>%
  filter(motif_name %in% top_motifs_all)

# only 4th quantile
motif_heatmap_data <- motif_heatmap_data %>%
  filter(quantile == 4)




# Prepare data for heatmap: rows = studies, columns = motif archetypes, values = centrality_scaled
library(tidyr)
library(ggplot2)
library(ggtext)

# Spread the data to wide format: rows = study, columns = motif_name
heatmap_matrix <- motif_heatmap_data %>%
  select(study, motif_name, centrality_scaled) %>%
  tidyr::pivot_wider(names_from = motif_name, values_from = centrality_scaled)

# Convert to matrix for heatmap plotting (remove study column)
rownames_heatmap <- heatmap_matrix$study
heatmap_matrix <- as.matrix(heatmap_matrix[,-1])
rownames(heatmap_matrix) <- rownames_heatmap

# Prepare motif label colors
motif_names <- colnames(heatmap_matrix)
motif_labels <- sapply(motif_names, function(m) {
  col <- if (m %in% names(ebox_colors)) ebox_colors[[m]] else "#888888"
  paste0("<span style='color:", col, ";'>", m, "</span>")
})
names(motif_labels) <- motif_names

# Convert to long format for ggplot
heatmap_long <- as.data.frame(heatmap_matrix)
heatmap_long$study <- rownames(heatmap_matrix)
heatmap_long <- tidyr::pivot_longer(heatmap_long, -study, names_to = "motif_name", values_to = "centrality_scaled")



  # Plot heatmap with clustered columns (motifs)
  library(pheatmap)

  # Prepare matrix for clustering (motif_name as columns, study as rows)
  # Use the same heatmap_matrix as above
  # Replace NA with 0 or another value if you prefer
  mat_for_clust <- heatmap_matrix
  mat_for_clust[is.na(mat_for_clust)] <- 0

  # Cluster columns (motifs) and get the order
  col_clust <- hclust(dist(t(mat_for_clust)))
  motif_order_clust <- colnames(mat_for_clust)[col_clust$order]

  # Use ggtext::element_markdown to color x-axis labels via HTML <span> tags
  library(ggtext)

  motif_labels_clust <- sapply(motif_order_clust, function(m) {
    col <- if (m %in% names(ebox_colors)) ebox_colors[[m]] else "#888888"
    paste0("<span style='color:", col, ";'>", m, "</span>")
  })
  names(motif_labels_clust) <- motif_order_clust

  heatmap_long$motif_name <- factor(heatmap_long$motif_name, levels = motif_order_clust)







  # This will color the axis labels but not allow for rich text formatting.
  pdf(file = "~/Dropbox/induction/5_remodeling/archetypes around atac summits/archetypes_atac_intersect_centrality_foldchange_4thquantileonly_top50.pdf", width = 12, height = 4, useDingbats = FALSE)

  motif_axis_colors_vec <- sapply(levels(heatmap_long$motif_name), function(m) {
    if (m %in% names(ebox_colors)) ebox_colors[[m]] else "#888888"
  })

  gg2 <- ggplot(heatmap_long, aes(x = motif_name, y = study, fill = centrality_scaled)) +
    geom_tile() +
    scale_fill_viridis_c(option = "cividis", na.value = "#888888") +
    labs(title = "Motif Archetype Centrality (4th Quantile)", x = "Motif Archetype", y = "Study", fill = "Centrality (scaled)") +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 7),
      axis.text.x = element_text(size = 6, angle = 62, hjust = 1, vjust = 1, colour = motif_axis_colors_vec)
    ) +
    scale_x_discrete(drop = FALSE)

  print(gg2)
dev.off()



####################### THE SAME BUT WITH THE ATAC-SEQ PEAKS THAT DO NOT INTERSECT ##############
#################################################################################################
##################### Supplementary figure 6c



setwd("~/proneural")
# Read all *_motif_archetype_median_summary.rds files, only those with intersect == "yes"
library(dplyr)

# List all RDS files matching the pattern in the current directory (or set your path)
rds_files <- list.files(pattern = "_intersect_no_motif_archetype_median_summary\\.rds$")

# Read and filter for intersect == "yes"
motif_summaries <- lapply(rds_files, function(f) {
  df <- readRDS(f) } )



# Combine into a single data frame
motif_summaries_no <- bind_rows(motif_summaries)

# Extract quantile and study name from the file column
motif_summaries_no <- motif_summaries_no %>%
  mutate(
    quantile = sub(".*_quantil_([0-9])_intersect_no.*$", "\\1", file),
    study = sub("juancho_combined_vierstra_with_distance_peakid_([a-zA-Z0-9_]+)_quantil_\\d+_intersect_no.*$", "\\1", file)
  )

motif_summaries_no$quantile <- as.numeric(motif_summaries_no$quantile) + 1
motif_summaries_no$centrality <- 500 - motif_summaries_no$median_abs_distance

motif_summaries_no<-motif_summaries_no %>% filter(motif_name!=".")

motif_summaries_no <- motif_summaries_no %>% filter(quantile == 4)

motif_summaries_no <- motif_summaries_no %>%
  group_by(study) %>%
  mutate(centrality_scaled = (centrality - min(centrality, na.rm = TRUE)) / (max(centrality, na.rm = TRUE) - min(centrality, na.rm = TRUE))) %>%
  ungroup()

# Apply motif renaming
motif_summaries_no$motif_name <- ifelse(
  motif_summaries_no$motif_name %in% names(motif_rename_map),
  motif_rename_map[motif_summaries_no$motif_name],
  motif_summaries_no$motif_name
)

# Remove unwanted studies
motif_summaries_no <- motif_summaries_no %>%
  filter(!(study %in% c("manandhar_myod1", "lee_ascl1_48h")))

# Rename studies
study <- unique(motif_summaries_no$study)
study_mapping <- data.frame(study = study, renamed_study = study_renamed)
motif_summaries_no <- motif_summaries_no %>%
  mutate(study = study_mapping[match(study, study_mapping$study), "renamed_study"])

studies <- unique(motif_summaries_no$study)





# 2. Combined heatmap for 4th quantile, top 50 motifs per study
top_motifs_all_no <- motif_summaries_no %>% filter(motif_name %in% shared_archetypes) %>%
  filter(quantile == 4) %>%
  group_by(study) %>%
  arrange(desc(centrality_scaled)) %>%
  slice_head(n = 50) %>%
  ungroup() %>%
  pull(motif_name) %>%
  unique()

motif_summaries_no <- motif_summaries_no %>%
  filter(motif_name %in% top_motifs_all_no, quantile == 4)

heatmap_matrix_no <- motif_summaries_no %>%
  select(study, motif_name, centrality_scaled) %>%
  tidyr::pivot_wider(names_from = motif_name, values_from = centrality_scaled)

rownames_heatmap_no <- heatmap_matrix_no$study
heatmap_matrix_no <- as.matrix(heatmap_matrix_no[,-1])
rownames(heatmap_matrix_no) <- rownames_heatmap_no

motif_names_no <- colnames(heatmap_matrix_no)
motif_labels_no <- sapply(motif_names_no, function(m) {
  col <- if (m %in% names(ebox_colors)) ebox_colors[[m]] else "#888888"
  paste0("<span style='color:", col, ";'>", m, "</span>")
})
names(motif_labels_no) <- motif_names_no

heatmap_long_no <- as.data.frame(heatmap_matrix_no)
heatmap_long_no$study <- rownames(heatmap_matrix_no)
heatmap_long_no <- tidyr::pivot_longer(heatmap_long_no, -study, names_to = "motif_name", values_to = "centrality_scaled")

library(pheatmap)
mat_for_clust_no <- heatmap_matrix_no
mat_for_clust_no[is.na(mat_for_clust_no)] <- 0
col_clust_no <- hclust(dist(t(mat_for_clust_no)))
motif_order_clust_no <- colnames(mat_for_clust_no)[col_clust_no$order]
motif_labels_clust_no <- motif_labels_no[motif_order_clust_no]
heatmap_long_no$motif_name <- factor(heatmap_long_no$motif_name, levels = motif_order_clust_no)
pdf(file = "~/Dropbox/induction/5_remodeling/archetypes_atac_nointersect_centrality_foldchange_4thquantileonly_top50.pdf", width = 12, height = 4, useDingbats = FALSE)
motif_axis_colors_vec_no <- sapply(levels(heatmap_long_no$motif_name), function(m) {
  if (m %in% names(ebox_colors)) ebox_colors[[m]] else "#888888"
})
gg <- ggplot(heatmap_long_no, aes(x = motif_name, y = study, fill = centrality_scaled)) +
  geom_tile() +
  scale_fill_viridis_c(option = "cividis", na.value = "#888888") +
  labs(title = "Motif Archetype Centrality (4th Quantile, Not Bound)", x = "Motif Archetype", y = "Study", fill = "Centrality (scaled)") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 6, angle = 62, hjust = 1, vjust = 1, colour = motif_axis_colors_vec_no)
  ) +
  scale_x_discrete(drop = FALSE)
print(gg)
dev.off()




