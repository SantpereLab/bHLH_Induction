
setwd("~/proneural/chip-seq")
studies <- readRDS("filtered_studies.RDS")
study_mapping <- readRDS("study_mapping.RDS")
# Rename studies using the mapping
studies <- study_mapping[match(studies, study_mapping$renamed_study), "study"]

# Display the first few rows of the renamed studies
head(studies)

methylation_files <- list.files(pattern = "\\.methylation\\.bed$", recursive = TRUE)

# Initialize an empty list to store data frames
meth_data_list <- list()

# Loop through each file and process it
for (file in methylation_files) {
  print(file)
  # Read the file
  data <- read.table(file, header = FALSE)
  
  # Add a column for the study
  data$study <- gsub(".*\\/|\\.methylation\\.bed$", "", file)

  # Store the data frame in the list
  meth_data_list[[file]] <- data
}


# Combine all data frames into one
combined_meth_data <- do.call(rbind, meth_data_list)


# The peak IDs are incorrect. Extract the peak IDs from the chromosome and summit

combined_meth_data$summit <- as.numeric(combined_meth_data$V2)+as.numeric(combined_meth_data$V10)

combined_meth_data$peakid_new<-paste0(combined_meth_data$V1,"-",combined_meth_data$summit)

# Print the combined data frame
library(dplyr)
combined_meth_data<-combined_meth_data %>% select(V11,study,peakid_new)

colnames(combined_meth_data) <- c("meth", "study","peakid_new")

combined_meth_data$study <- gsub("_confidence\\.bed$", "", combined_meth_data$study)


studies_meth<-c("matsuda_neurod1_microglia", "pataskar_neurod1_mesc", "kim_neurod1_hek_treatment_pcag", "fong_neurod2_mef", "pereira_neurog2_astrocytes", "kim_twist1_hek_treatment_pcag", "kim_twist1_hek_treatment_pcdna", "mazzoni_olig2_mesc", "lin_mesp1_mesc_treatment_24h", "fong_myod1_mef", "casey_myod1_mesc", "lee_myod1_mef", "lee_myod1_mesc", "conerly_myod1_mef", "manandhar_myod1_fibroblasts", "conerly_myf5_mef", "casey_ascl1_mesc", "lee_ascl1_mesc", "lee_ascl1_mef", "wapinski_ascl1_mef", "vainorius_ascl1_mesc", "woods_ascl1_neuroblastoma_wt31_treatment", "woods_ascl1_neuroblastoma_wt51_treatment", "joung_ascl1_hesc", "casey_ascl2_mesc", "soufi_myc_bj_treatment", "croci_myc_mef", "depretis_myc_mef_treatment_2h", "weber_hey1_mesc_treatment", "weber_hey2_mesc_treatment")


combined_chromatin_data<-readRDS("combined_chromatin_data.RDS")

combined_chromatin_data<-combined_chromatin_data %>% filter(study %in% studies_meth)

combinee<-merge(combined_chromatin_data, combined_meth_data, by = c("peakid_new","study"), all.x=TRUE   )

combinee <- combinee %>% filter(study!="thomas_myc_hek_treatment")
saveRDS(combinee, "meth_accessibility.RDS")


################### PLOT METHYLATION DISTRIBUTION IN QUANTILES #####################
####################################################################################
###### Supplementary figure 5a

setwd("~/proneural/chip-seq")
combinee<-readRDS("meth_accessibility.RDS")

library(dplyr)
combinee <- combinee %>%
  mutate(quantile = case_when(
    quantile <= 25 ~ 4,
    quantile > 25 & quantile <= 50 ~ 3,
    quantile > 50 & quantile <= 75 ~ 2,
    quantile > 75 ~ 1 ))


library(ggridges)
library(viridis)
viridis_soft <- rev(viridis(100, option = "G")[15:95] )

library(gridExtra)
study_mapping <- readRDS("study_mapping.RDS")
# Rename studies using the mapping
combinee<- combinee %>% mutate(study = study_mapping[match(study, study_mapping$study), "renamed_study"])
# Ensure the studies are in the correct order

studies_meth<-unique(combinee$study)
orden<-readRDS("orden.RDS")
studies_meth <- orden[which(orden %in% studies_meth)]


studies_meth <- gsub("NEUROD1 -- HEK295 \\[Kim et al\\. 2023\\]", "NEUROD1 -- HEK293 [Kim et al. 2023]", studies_meth)
studies_meth <- gsub("TWIST1 -- HEK293 \\[Kim et al\\. 2023\\]", "TWIST1 -- HEK293 [Kim et al. 2023] (1)", studies_meth)
studies_meth<- gsub("TWIST1 -- HEK294 \\[Kim et al\\. 2023\\]", "TWIST1 -- HEK293 [Kim et al. 2023] (2)", studies_meth)

combinee$study<-gsub("NEUROD1 -- HEK295 \\[Kim et al\\. 2023\\]", "NEUROD1 -- HEK293 [Kim et al. 2023]", combinee$study)
combinee$study<-gsub("TWIST1 -- HEK293 \\[Kim et al\\. 2023\\]", "TWIST1 -- HEK293 [Kim et al. 2023] (1)", combinee$study)
combinee$study<-gsub("TWIST1 -- HEK294 \\[Kim et al\\. 2023\\]", "TWIST1 -- HEK293 [Kim et al. 2023] (2)", combinee$study)


library(ggplot2)
# Create a list to store individual plots
plot_list <- list()

for (estudio in studies_meth) {
  
  rafa <- combinee %>% filter(study == estudio) %>% select(peakid_new, meth, quantile) %>% distinct()

  rafa <- rafa %>% select(meth, quantile) %>% filter(meth != ".")
  rafa$meth <- as.numeric(rafa$meth)
  
  if (max(rafa$meth) < 2) { rafa$meth <- rafa$meth * 100 }
   
  rafa$quantile <- factor(rafa$quantile, levels = c("4", "3", "2", "1"))

  gigi <- ggplot(rafa, aes(x = meth, y = quantile, fill = ..x..)) +
  geom_density_ridges_gradient(alpha = 0.7, scale = 1.5) +
  scale_fill_gradientn(colors = viridis_soft) +
  scale_x_continuous(limits = c(0, 100)) +
  theme_minimal() + labs(
    title = estudio,
    x = "Mean % mCpG",
    y = "Quantile"
  ) + theme(
    legend.position = "none",
    axis.text.y = element_text(size = 12),
    panel.border = element_rect(color = "white", fill = NA, size = 1)
  ) + scale_x_continuous(limits = c(0, 100))

  # Add the plot to the list
  plot_list[[estudio]] <- gigi
}

# aqui lo he dejado. queda esto y hacer bien las distribuciones de los motivos metilacion
saveRDS(plot_list, "plot_list.RDS")

plot_list<-readRDS("plot_list.RDS")


# Arrange all plots in a grid
pdf("~/Dropbox/induction/4_methylation/distribution_methylation_in_quantiles.pdf", useDingbats = F, height = 45, width = 15)
grid.arrange(grobs = plot_list, ncol = 3) # Adjust ncol for desired layout
dev.off()






######### MOTIFS IN METHYLATION + ACCESSIBILITY ##########
##########################################################
######## Supplementary figure 5b



library(dplyr)
setwd("~/proneural/chip-seq")

# eboxes + methylation + accessibility. 400bp around summit
methylationm_accessibility_peaks<-readRDS("meth_accessibility.RDS") %>% distinct()

## canonical eboxes 400bp around summits
canonical_eboxes<-readRDS("combined_canonical_eboxes.RDS")
canonical_eboxes$peakid_new <- paste0(canonical_eboxes$chrom, "-", canonical_eboxes$summit)

# assign methylation and accessibility to eboxes
canonical_eboxes %>% left_join(methylationm_accessibility_peaks, by = c("peakid_new", "study")) -> canonical_eboxes_meth_atac

canonical_eboxes_meth_atac<- canonical_eboxes_meth_atac %>% filter(meth != ".")

canonical_eboxes_meth_atac$meth<-as.numeric(canonical_eboxes_meth_atac$meth) 


## obtain max methylation for each study
max_meth <- canonical_eboxes_meth_atac %>%
  group_by(study) %>%
  summarise(max_meth = max(meth, na.rm = TRUE), .groups = "drop")

max_meth %>% print(n=30)


# if the maximum is 1, multiply by 100
canonical_eboxes_meth_atac <- canonical_eboxes_meth_atac %>%
  left_join(max_meth, by = "study") %>%
  mutate(meth = ifelse(max_meth == 1, meth * 100, meth)) %>%
  select(-max_meth)

# convert methylation into 0-10% 10-20% .. etc
canonical_eboxes_meth_atac$meth_window <- cut(canonical_eboxes_meth_atac$meth, 
                                          breaks = seq(0, 100, by = 10), 
                                          include.lowest = TRUE, 
                                          labels = paste0(seq(0, 90, by = 10), "-", seq(10, 100, by = 10)))


canonical_eboxes_meth_atac$ebox_center<-canonical_eboxes_meth_atac$ebox_start+3
canonical_eboxes_meth_atac$dist_summit<-canonical_eboxes_meth_atac$ebox_center-canonical_eboxes_meth_atac$summit

head(canonical_eboxes_meth_atac)

saveRDS(canonical_eboxes_meth_atac, "canonical_eboxes_meth_atac.RDS")



##### plot central dinucleotides in 4th quantile of methylation.
################################################################

# both centrality and num

# load natural length of peaks
jj<-readRDS("jj.RDS")

head(canonical_eboxes_meth_atac)

canonical_eboxes_meth_atac_natural_length<-merge(canonical_eboxes_meth_atac, jj, by = c("peakid_new", "study"), all.x = TRUE)

saveRDS(canonical_eboxes_meth_atac_natural_length, "canonical_eboxes_meth_atac_natural_length.RDS")




# load
canonical_eboxes_meth_atac_natural_length<-readRDS("canonical_eboxes_meth_atac_natural_length.RDS")



## convert 1 to 100 quantiles to 1 to 4
canonical_eboxes_meth_atac_natural_length$quantile <- case_when(
  canonical_eboxes_meth_atac_natural_length$quantile <= 25 ~ 4,
  canonical_eboxes_meth_atac_natural_length$quantile > 25 & canonical_eboxes_meth_atac_natural_length$quantile <= 50 ~ 3,
  canonical_eboxes_meth_atac_natural_length$quantile > 50 & canonical_eboxes_meth_atac_natural_length$quantile <= 75 ~ 2,
  canonical_eboxes_meth_atac_natural_length$quantile > 75 ~ 1
)

canonical_eboxes_meth_atac_natural_length_quantile4<- canonical_eboxes_meth_atac_natural_length %>%
  filter(quantile == 4) %>% distinct()

# change meth windows. do from 0-25%, 25-50%, 50-75%, 75-100% . from meth values
canonical_eboxes_meth_atac_natural_length_quantile4$meth_window<- cut(canonical_eboxes_meth_atac_natural_length_quantile4$meth, 
                                          breaks = seq(0, 100, by = 25), 
                                          include.lowest = TRUE, 
                                          labels = paste0(seq(0, 75, by = 25), "-", seq(25, 100, by = 25)))


# for each central dinucleotide and study and meth window, get the mean distance to the summit
centralities<-canonical_eboxes_meth_atac_natural_length_quantile4 %>% group_by(study, central_dinucleotide, meth_window) %>%
  summarise(centrality = 400-mean(abs(dist_summit), na.rm = TRUE) ) 
# scale centrality to 0-1. with respect to the maximum and minimum of each study
centralities<- centralities %>%
  group_by(study) %>%
  mutate(centrality = (centrality - min(centrality, na.rm = TRUE)) / (max(centrality, na.rm = TRUE) - min(centrality, na.rm = TRUE))) %>%
  ungroup()

# for each central dinucleotide and study and meth window, get the mean number of motifs per peak. taking the natural length of peaks

#first get the total number of peaks per study and meth window
npeaks_meth_window<- canonical_eboxes_meth_atac_natural_length_quantile4 %>%
  select(study, meth_window, peakid_new) %>% distinct() %>% group_by(study,meth_window) %>% summarise(npeaks=n())

npeaks_meth_window %>% arrange(npeaks)

# get the eboxes that are within the natural peak length. within the windows
canonical_eboxes_meth_atac_natural_length_quantile4 %>% filter(ebox_center >= window_start & ebox_center <= window_end ) -> canonical_eboxes_meth_atac_natural_length_quantile4_within_natural_length

# get the number of motifs (each different central_dinucleotide) per meth/window/study
nmotifs_meth_windows <- canonical_eboxes_meth_atac_natural_length_quantile4_within_natural_length %>%
  group_by(study, central_dinucleotide, meth_window) %>%
  summarise(nmotifs = n(), .groups = "drop") %>%
  left_join(npeaks_meth_window, by = c("study", "meth_window")) %>%
  mutate(nmotifs_per_peak_corrected = nmotifs / npeaks) %>%
  select(-npeaks)

# scale between 0 and 1 with respect to the maximum and minimum of each study
nmotifs_meth_windows <- nmotifs_meth_windows %>%
  group_by(study) %>%
  mutate(nmotifs_per_peak_corrected_scaled = (nmotifs_per_peak_corrected - min(nmotifs_per_peak_corrected, na.rm = TRUE)) / 
         (max(nmotifs_per_peak_corrected, na.rm = TRUE) - min(nmotifs_per_peak_corrected, na.rm = TRUE))) %>%
  ungroup()

# join the table of the centrality and the table of the number of motifs per peak
nmotifs_centrality_meth_windows <- nmotifs_meth_windows %>%
  left_join(centralities, by = c("study", "central_dinucleotide", "meth_window"))

nmotifs_centrality_meth_windows <- nmotifs_centrality_meth_windows %>% select(study, central_dinucleotide, meth_window, nmotifs_per_peak_corrected_scaled, centrality)
# complete. all study, central, central_dinucleotide, meth_window combinationns mus be filled. if not present, put 0 in the other fields
nmotifs_centrality_meth_windows <- nmotifs_centrality_meth_windows %>%
  tidyr::complete(study, central_dinucleotide, meth_window, fill = list(nmotifs_per_peak_corrected_scaled = 0, centrality = 0)) %>%
  arrange(study, central_dinucleotide, meth_window)

## do balloon plot. size of baloon representing number of motifs. and color representing centrality.
## do a balloon plot for each study, in a loop. in each "row", central dinucleotides, and in the "columns" meth_windows

# aqui lo he dejado. tengo que ponerlo bonito

saveRDS(nmotifs_centrality_meth_windows, "nmotifs_centrality_meth_windows.RDS")

setwd("~/proneural/chip-seq")
nmotifs_centrality_meth_windows <- readRDS("nmotifs_centrality_meth_windows.RDS")


library(ggplot2)
library(viridis)
library(dplyr)
library(tidyr)


# rename studies
study_mapping <- readRDS("study_mapping.RDS")
nmotifs_centrality_meth_windows <- nmotifs_centrality_meth_windows %>%
  mutate(study = study_mapping[match(study, study_mapping$study), "renamed_study"])

# correct HEK
nmotifs_centrality_meth_windows$study <- gsub("NEUROD1 -- HEK295 \\[Kim et al\\. 2023\\]", "NEUROD1 -- HEK293 [Kim et al. 2023]", nmotifs_centrality_meth_windows$study)
nmotifs_centrality_meth_windows$study <- gsub("TWIST1 -- HEK293 \\[Kim et al\\. 2023\\]", "TWIST1 -- HEK293 [Kim et al. 2023] (1)", nmotifs_centrality_meth_windows$study)
nmotifs_centrality_meth_windows$study <- gsub("TWIST1 -- HEK294 \\[Kim et al\\. 2023\\]", "TWIST1 -- HEK293 [Kim et al. 2023] (2)", nmotifs_centrality_meth_windows$study)
# Ensure the studies are in the correct order
orden<-readRDS("orden.RDS")
orden<-orden[which(orden %in% unique(nmotifs_centrality_meth_windows$study) )]

nmotifs_centrality_meth_windows$central_dinucleotide <- factor(nmotifs_centrality_meth_windows$central_dinucleotide, 
                                levels = c("CAT-CAT", "CAT-CAG", "CAG-CAG", "CAG-CAC", "CAC-CAC", 
                                           "CAT-CAC", "CAA-CAT", "CAA-CAG", "CAA-CAC", "CAA-CAA"))


# Reverse the order of meth_window factor levels for plotting
nmotifs_centrality_meth_windows$meth_window <- factor(
  nmotifs_centrality_meth_windows$meth_window,
  levels = rev(levels(nmotifs_centrality_meth_windows$meth_window))
)



library(cowplot)

# Create a list to store plots
plot_list <- list()

for(estudio in orden) {
  study_data <- nmotifs_centrality_meth_windows %>% dplyr::filter(study == estudio)
  
  p <- ggplot(study_data, aes(x = central_dinucleotide, y = meth_window)) +
    geom_point(aes(size = nmotifs_per_peak_corrected_scaled, color = centrality), alpha = 0.7) +
    scale_size_continuous(range = c(1, 10), name = "Normalized # E-boxes") +
    scale_color_viridis_c(option = "cividis", name = "Normalized Distance") +
    theme_minimal(base_size = 7) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
      axis.text.y = element_text(size = 6),
      axis.title.x = element_text(size = 7),
      axis.title.y = element_text(size = 7),
      plot.title = element_text(size = 8),
      legend.position = "none"
    ) +
    labs(
      title = estudio,
      x = "Central Dinucleotide",
      y = "Methylation Window"
    )
  plot_list[[estudio]] <- p
}

# Arrange all plots in a grid using cowplot
ncol_grid <- 3
nrow_grid <- ceiling(length(plot_list) / ncol_grid)

pdf("~/Dropbox/induction/4_methylation/centrality_nmotifs_meth_windows_4th_quantil.pdf", useDingbats = F, height = 2.2 * nrow_grid, width = 4.5 * ncol_grid)
  plot_grid(plotlist = plot_list, ncol = ncol_grid)
dev.off()


