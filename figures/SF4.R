################# CENTRAL DINUCLEOTIDES IN PRE INDUCTION CHROMATIN ################
###################################################################################
############## Supplementary figure 4a

setwd("~/proneural/chip-seq")

canonical_eboxes<-readRDS("combined_canonical_eboxes.RDS")


# Initialize an empty list to store data frames
chromatin_data_list <- list()

options(scipen = 999) # Disable scientific notation for large numbers
# Loop through each file and process it
for (file in chromatin_files) {
  print(file)
  # Read the file
  data <- read.table(file, header = FALSE)
  
  # Add a column for the study
  data$study <- gsub("chromatin_|_confidence_1000bp_around_summit.bed", "", file)
  
  # Store the data frame in the list
  chromatin_data_list[[file]] <- data
}

# Combine all data frames into one
combined_chromatin_data <- do.call(rbind, chromatin_data_list)

# The peak IDs are incorrect. Extract the peak IDs from the chromosome and summit
combined_chromatin_data$summit<-as.numeric(combined_chromatin_data$V2+500)
combined_chromatin_data$peakid_new<-paste0(combined_chromatin_data$V1,"-",combined_chromatin_data$summit)

# Print the combined data frame
library(dplyr)
combined_chromatin_data<-combined_chromatin_data %>% select(V5,study,peakid_new)

colnames(combined_chromatin_data) <- c("accessibility", "study","peakid_new")

# Order each study by accessibility and create 100 quantiles
combined_chromatin_data <- combined_chromatin_data %>%
  group_by(study) %>%
  arrange(desc(as.numeric(accessibility))) %>%
  mutate(quantile = ntile(accessibility, 100)) %>%
  ungroup()

saveRDS(combined_chromatin_data,"combined_chromatin_data.RDS")
saveRDS(combined_chromatin_data,"~/proneural/combined_chromatin_data.RDS")


canonical_eboxes$peakid_new<-paste0(canonical_eboxes$chrom,"-",canonical_eboxes$summit)

# Combine combined_chromatin_data and canonical_eboxes by "peak_id" and "study" columns
combined_data <- merge(combined_chromatin_data, canonical_eboxes, by = c("peakid_new", "study"), all.x=TRUE)

saveRDS(combined_data,"combined_data.RDS")



####### PLOT DISTRIBUTION IN 4 QUANTILES
setwd("~/proneural/chip-seq")
combined_data<-readRDS("combined_data.RDS")

# Rename the studies based on a mapping file
study_mapping <- readRDS("study_mapping.RDS")

# Update the study names in the combined_data
combined_data <- combined_data %>%
  mutate(study = study_mapping[match(study, study_mapping$study), "renamed_study"])


studies<-readRDS("filtered_studies.RDS")
library(dplyr)
combined_data<-combined_data %>% filter(study %in% studies)



library(dplyr)
combined_data$ebox_center<-combined_data$ebox_start+3

combined_data$dist_summit<-combined_data$ebox_center-combined_data$summit

 # Convert quantiles into groups 1 to 4
 combined_data_quantiles <- combined_data %>%
   mutate(quantile = case_when(
     quantile <= 25 ~ 4,
     quantile > 25 & quantile <= 50 ~ 3,
     quantile > 50 & quantile <= 75 ~ 2,
     quantile > 75 ~ 1
   ))




library(ggplot2)
library(ggridges)
library(tidyr)


tabla <- combined_data %>%
  group_by(study, central_dinucleotide, dist_summit, quantile) %>%
  summarise(ocurrencias = n(), .groups = "drop")

tabla <- tabla %>%
  group_by(study, central_dinucleotide) %>% # Group by necessary columns
  tidyr::complete(dist_summit,quantile, fill = list(ocurrencias = 0)) %>% # Fill missing values
  ungroup() # Ungroup after completion



tabla<-tabla %>% mutate(central_dinucleotide = factor(central_dinucleotide, levels = c("CAT-CAT","CAT-CAG","CAG-CAG","CAG-CAC","CAC-CAC","CAT-CAC","CAA-CAT","CAA-CAG","CAA-CAC","CAA-CAA") )  )

#tabla<-tabla %>% filter(quantile %in% c(1, 4))

tabla <- tabla %>%
  mutate(bin = floor(dist_summit / 20) * 20) %>% # Crear bins de 5 bp
  group_by(study, central_dinucleotide,quantile, bin) %>%
  summarise(ocurrencias = mean(ocurrencias, na.rm = TRUE), .groups = "drop") %>%
  rename(dist_summit = bin) 


tabla <- tabla %>% filter(abs(dist_summit)<=395)

maxii<-tabla %>% group_by(study,central_dinucleotide) %>% summarise(maximo=max(ocurrencias))

tabla<-merge(maxii, tabla, by= c("study","central_dinucleotide"))

tabla$ocurrencias_escaladas<-tabla$ocurrencias/tabla$maximo



# Remove rows with NA central_dinucleotide
tabla <- tabla %>% filter(!is.na(central_dinucleotide) )

# Update study names based on the specified changes
tabla <- tabla %>%
  mutate(study = case_when(
    study == "NEUROD1 -- HEK295 [Kim et al. 2023]" ~ "NEUROD1 -- HEK293 [Kim et al. 2023]",
    study == "TWIST1 -- HEK293 [Kim et al. 2023]" ~ "TWIST1 -- HEK293 (1) [Kim et al. 2023]",
    study == "TWIST1 -- HEK294 [Kim et al. 2023]" ~ "TWIST1 -- HEK293 (2) [Kim et al. 2023]",
    TRUE ~ study
  ))






# Add a column for the maximum occurrences within each study and central dinucleotide
tabla_scaled <- tabla %>%
  group_by(study, central_dinucleotide) %>%
  mutate(max_ocurrencias = max(ocurrencias, na.rm = TRUE)) %>%
  ungroup()

# Scale occurrences using the new max_ocurrencias column
tabla_scaled <- tabla_scaled %>%
  mutate(ocurrencias_scaled = ifelse(max_ocurrencias > 0, ocurrencias / max_ocurrencias, 0))

### the same plot but 400bp around summit
tabla_scaled_400bp <- tabla_scaled %>%
  filter(abs(dist_summit) <= 200) # Filter for distances within 200 bp


pdf("~/Dropbox/induction/2_pre_induction_chromatin/distribution_central_dinucleotides_1st_vs_4th_quantiles_800bpwindow_motif_scaled.pdf", height = 100, width = 30, useDingbats = FALSE)
tabla_scaled$quantile <- as.factor(tabla_scaled$quantile)
gg <- ggplot(tabla_scaled %>% filter(quantile %in% c("1", "4")), 
  aes(x = dist_summit, y = ocurrencias_scaled, fill = quantile, color = quantile)) +
  geom_area(alpha = 0.5, position = "identity") +
  facet_grid(study ~ central_dinucleotide, scales = "free_y") +
  scale_fill_manual(values = c("1" = "blue", "4" = scales::alpha("azure", 0.1))) +
  scale_color_manual(values = c("1" = "black", "4" = "black")) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text.y = element_text(angle = 0, size = 10),
    panel.spacing.y = unit(1, "lines"),
    legend.position = "bottom"
  ) +
  labs(
    title = "Distribution of Central Dinucleotides Across Studies (Scaled by Study and Motif)",
    x = "Distance to summit (binned by 10 bp)",
    y = "Scaled occurrences",
    fill = "Quantile",
    color = "Quantile"
  )
print(gg)
dev.off()










#####################################################################
########### Supplementary figure 4b ###############
#####################################################################


estudios<-unique(dinucleotide_counts$study)

# List all studies in the "human" directory
human_studies <- list.files(path = "~/proneural/chip-seq/human", pattern = "_confidence.bed$", full.names = FALSE)
# Extract unique study names by removing the "_confidence.bed" suffix
human_studies <- unique(gsub("_confidence.bed$", "", human_studies))
print(human_studies)

mouse_studies <- list.files(path = "~/proneural/chip-seq/mouse", pattern = "_confidence.bed$", full.names = FALSE)
# Extract unique study names by removing the "_confidence.bed" suffix
mouse_studies <- unique(gsub("_confidence.bed$", "", mouse_studies))
print(mouse_studies)

estudios_human<-estudios[which(estudios %in% human_studies)]
estudios_mouse<-estudios[which(estudios %in% mouse_studies)]



tss_mm10<-GenomicFeatures::makeTxDbFromGFF("~/n/genomes/gencode.vM24.annotation.filtered.gtf", format="gtf")

tss_hg38<-GenomicFeatures::makeTxDbFromGFF("~/n/genomes/gencode.v47.annotation.filtered.gtf", format="gtf")




dist_tss<-function (peaks, tss ) {
  require(dplyr)
  require(GenomicRanges)
  require(ChIPseeker)
  colnames(peaks)<-c("chr","summit","peak_id")
  
  peaks$start<-as.numeric(peaks$summit)-1
  peaks$end<-as.numeric(peaks$summit)
  peaks$peak_id_2<-paste0(peaks$chr,":",peaks$start,"-",peaks$end)
  ranges<-makeGRangesFromDataFrame(peaks,start.field = "start", end.field = "end", seqnames.field = "chr" )
  options(ChIPseeker.ignore_1st_exon = TRUE)
  options(ChIPseeker.ignore_downstream = TRUE)
  anno<-annotatePeak(ranges, TxDb=tss)
  anno<-as.data.frame( anno@anno )
  ano<-anno %>% select(seqnames,start,end,distanceToTSS)
  colnames(ano)[1]<-"chr"
  ano$peak_id_2<-paste0(ano$chr,":",anno$start,"-",anno$end)
  ano<-ano %>% select(peak_id_2,distanceToTSS)
  ano$distanceToTSS<-abs(ano$distanceToTSS)
  return(ano)
}



################### sacar en cada quantil cuantos picos hay de promotor y de enhancer #######################


combined_data <- readRDS("merged_data.RDS")


peaks_human<- combined_data %>% filter(study %in% estudios_human) %>% select(peakid_new,study) %>% distinct()

peaks_human$chrom <- gsub("-.*", "", peaks_human$peakid_new)
peaks_human$summit <- as.numeric(gsub(".*-", "", peaks_human$peakid_new))
peaks_human$peak_id<-peaks_human$peakid_new



i<-1
for (studio in estudios_human) {
  print(studio)
  studioo<-peaks_human %>% filter(study==studio) 
  studio_summits<-studioo %>% select(chrom,summit, peak_id)
  jj<-dist_tss(studio_summits,tss_hg38)
  studioo$distance_tss<-jj$distanceToTSS
  
  if (i==1) {
    bind_filtered_tss_human<-studioo
  }
  else {
    bind_filtered_tss_human<-rbind(bind_filtered_tss_human,studioo)
  }
  i<-i+1
  gc()
}

saveRDS(bind_filtered_tss_human, file = "bind_filtered_tss_human.RDS")



# Process mouse studies
peaks_mouse <- combined_data %>% filter(study %in% estudios_mouse) %>% select(peakid_new, study) %>% distinct()

peaks_mouse$chrom <- gsub("-.*", "", peaks_mouse$peakid_new)
peaks_mouse$summit <- as.numeric(gsub(".*-", "", peaks_mouse$peakid_new))
peaks_mouse$peak_id <- peaks_mouse$peakid_new

i <- 1
for (studio in estudios_mouse) {
  print(studio)
  studioo <- peaks_mouse %>% filter(study == studio)
  studio_summits <- studioo %>% select(chrom, summit, peak_id)
  jj <- dist_tss(studio_summits, tss_mm10)
  studioo$distance_tss <- jj$distanceToTSS
  
  if (i == 1) {
    bind_filtered_tss_mouse <- studioo
  } else {
    bind_filtered_tss_mouse <- rbind(bind_filtered_tss_mouse, studioo)
  }
  i <- i + 1
  gc()
}

saveRDS(bind_filtered_tss_mouse, file = "bind_filtered_tss_mouse.RDS")


bind_filtered_tss_mouse<-readRDS("bind_filtered_tss_mouse.RDS")
bind_filtered_tss_human<-readRDS("bind_filtered_tss_human.RDS")


combined_tss_data <- rbind(bind_filtered_tss_human, bind_filtered_tss_mouse)
# Merge the combined data with the TSS distance data
combined_data <- combined_data %>%
  left_join(combined_tss_data %>% select(peakid_new, study,distance_tss), by = c("peakid_new", "study"))

combined_data[which(is.na(combined_data$distance_tss)),] %>% dim()

# Create a new column to classify peaks based on their distance to TSS
combined_data <- combined_data %>%
  mutate(region_type = case_when(
    abs(distance_tss) <= 2000 ~ "Promoter",
    abs(distance_tss) > 2000 & abs(distance_tss) <= 4000 ~ "Proximal Enhancer",
    abs(distance_tss) > 4000 ~ "Distal Enhancer",
    TRUE ~ NA_character_
  ))




# Count the number of peaks in each region type for each study and quantile
region_counts <- combined_data %>% select(study, quantile, region_type, peakid_new) %>% distinct() %>%
  group_by(study, quantile, region_type) %>%
  summarise(count = n(), .groups = "drop")

# Do a stacked barplot for each study. Each bar is a quantile, and the colors are the region types
library(ggplot2)
library(dplyr)
region_counts <- region_counts %>%
  mutate(region_type = factor(region_type, levels = c("Promoter", "Proximal Enhancer", "Distal Enhancer")))


# Update the study names in the region_counts
region_counts <- region_counts %>%
  mutate(study = study_mapping[match(study, study_mapping$study), "renamed_study"])

filtered_studies <- readRDS("filtered_studies.RDS")

# Filter region_counts to include only studies in filtered_studies
region_counts <- region_counts %>% filter(study %in% filtered_studies)



# Prepare data for heatmap: proportion of promoters per study and quantile
region_counts_heatmap <- region_counts %>%
  group_by(study, quantile) %>%
  mutate(total = sum(count)) %>%
  ungroup() %>%
  filter(region_type == "Promoter") %>%
  mutate(proportion_promoter = count / total) %>%
  select(study, quantile, proportion_promoter)

region_counts_heatmap$quantile<-5-region_counts_heatmap$quantile
# Reshape for heatmap
heatmap_data <- region_counts_heatmap %>%
  tidyr::pivot_wider(names_from = quantile, values_from = proportion_promoter, values_fill = 0)

heatmap_matrix <- as.matrix(heatmap_data[,-1])
rownames(heatmap_matrix) <- heatmap_data$study

orden <- readRDS("orden.RDS")
orden <- orden[orden %in% rownames(heatmap_matrix)]
heatmap_matrix <- heatmap_matrix[orden, c(4,3,2,1)]

library(pheatmap)
library(viridis)
pdf("~/Dropbox/induction/2_pre_induction_chromatin/promoter_proportion_heatmap.pdf", height = 10, width = 4, useDingbats = FALSE)
pheatmap(
  heatmap_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = viridis::rocket(100),
  main = "Proportion of Promoters Across Quantiles",
  fontsize_row = 8,
  fontsize_col = 8
)
dev.off()






###########################################
######## Supplementary figure 4c ##########
###########################################

# Filter for Twist studies
twist_data <- distance_proportions %>% filter(grepl("WIST", study))

# Convert quantile to a factor for better plotting
twist_data$quantile<-5-twist_data$quantile
twist_data$quantile <- factor(twist_data$quantile, levels = c(1, 2, 3, 4))



# Generate stacked barplots with numbers inside the bars
pdf("~/Dropbox/induction/2_pre_induction_chromatin/twist_studies_distance_proportions_in_quantiles_stacked.pdf", height = 15, width = 4, useDingbats = FALSE)
ggplot(twist_data, aes(x = quantile, y = proportion, fill = factor(distance_to_next))) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = distance_to_next, size = proportion), 
            position = position_stack(vjust = 0.5), 
            color = "salmon") +
  facet_wrap(~study, ncol = 1, scales = "free_y") +
  scale_fill_viridis_d(option = "mako", name = "Distance", direction = -1) +
  scale_size_continuous(range = c(2, 6), guide = "none") + # Adjust size range for text
  theme_minimal() +
  labs(
    title = "Proportion of Spacings Across Quantiles",
    x = "Quantile",
    y = "Proportion"
  ) +
  theme(
    strip.text = element_text(size = 10),
    legend.position = "bottom"
)
dev.off()









