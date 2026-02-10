######### MEDIAN DISTANCE FROM SUMMIT + motif rate ########
######################################################
## Main figure 1a



# Reorder the studies according to the specified order
orden <- c(
  "NEUROD1 -- Microglia [Matsuda et al. 2016]", "NEUROD1 -- mESC [Pataskar et al. 2016]", "NEUROD1 -- HEK295 [Kim et al. 2023]",
  "NEUROD2 -- P19 [Fong et al. 2012]", "NEUROD2 -- MEF [Fong et al. 2012]", "NEUROD2 -- P19 [Fong et al. 2015]",
  "NEUROG2 -- EB [Aydin et al. 2019]", "NEUROG2 -- EB [Velasco et al. 2016]",
  "NEUROG2 -- MRC5 [Smith et al. 2016]", "NEUROG2 -- mESC [Vainorius et al. 2023]",
  "NEUROG2 -- Astrocytes [Pereira et al. 2024]", 
  "OLIG2 -- mESC [Mazzoni et al. 2011]",
  "TWIST2 -- Myoblast (DM) [Li et al. 2019]", "TWIST2 -- Myoblast (GM) [Li et al. 2019]",
  "TWIST1 -- EB [Fan et al. 2012]", "TWIST1 -- MCF7 [Yu et al. 2023]", "TWIST1 -- BMSC [Lee et al. 2023]",
  "TWIST1 -- HEK293 [Kim et al. 2023]", "TWIST1 -- HEK294 [Kim et al. 2023]",
  "ATOH1 -- mESC [Costa et al. 2022]", "ATOH8 -- A549 [Liu et al. 2023]",    
  "MESP1 -- mESC [Lin et al. 2022]", "MSGN1 -- EB [Chalamalasetty et al. 2014]",
  "MYOD1 -- P19 [Fong et al. 2012]", "MYOD1 -- MEF [Fong et al. 2012]",
  "MYOD1 -- mESC [Casey et al. 2018]", "MYOD1 -- MEF [Lee et al. 2020]", "MYOD1 -- mESC [Lee et al. 2020]",
  "MYOD1 -- MEF [Conerly et al. 2016]", "MYOD1 -- Human fibroblasts [Manandhar et al. 2017]", "MYF5 -- MEF [Conerly et al. 2016]",
  "ASCL1 -- mESC [Casey et al. 2018]", "ASCL1 -- EB [Aydin et al. 2019]", "ASCL1 -- mESC [Lee et al. 2020]", "ASCL1 -- MEF [Lee et al. 2020]",
  "ASCL1 -- MEF [Wapinski et al. 2013]", "ASCL1 -- NPC [Wapinski et al. 2013]", "ASCL1 -- G523NS [Park et al. 2017]",
  "ASCL1 -- mESC [Vainorius et al. 2023]", "ASCL1 -- SH-SY5Y [Woods et al. 2023] (1)",
  "ASCL1 -- SH-SY5Y [Woods et al. 2023] (2)", "ASCL1 -- GI-MEN [Wang et al. 2023]",
  "ASCL1 -- hESC [Joung et al. 2023]", "ASCL1 -- NS5 [Raposo et al. 2015]","ASCL2 -- mESC [Casey et al. 2018]",
  "ASCL2 -- CD4+ T cells [Liu et al. 2023]", "TCF4 -- SU-DHL-2 [Jain et al. 2019]", "TCF4 -- TMD8 [Jain et al. 2019]",
  "MYC -- U2OS [Lorenzin et al. 2016]", "MYC -- BJ [Soufi et al. 2012]",
  "MYC -- P493-6 [Lin et al. 2012]", "MYC -- U2OS [See et al. 2022]",
  "MYC -- U2OS [Walz et al. 2014]", 
  "MYC -- LNCaP [Barfeld et al. 2017]",
  "MYC -- MEF [Croci et al. 2017]", "MYC -- MEF [Sabo et al. 2014]", "MYC -- P493-6 [Sabo et al. 2014] (1)",
  "MYC -- P493-6 [Sabo et al. 2014] (2)", "MYC -- MEF [De Pretis et al. 2017]", 
  "MYC -- MCF10A [Muthalagu et al. 2014]", "MYC -- ATCC CCL-221 [Tameire et al. 2018]",
  "MYC -- IMEC [Jaenicke et al. 2016]", "MYC -- H2171 [Liang et al. 2018]"," MYC -- HEK293 [Thomas et al. 2015]",
  "MYC -- U2OS [Jung et al. 2017]", "MYCN -- SHEP [Zeid et al. 2018]",
  "MYCN -- LAN5 [Upton et al. 2020]","MYCN -- NB1643 [Upton et al. 2020]","MYCN -- COGN415 [Upton et al. 2020]",
  "MYCN -- NGP [Bosse et al. 2017]", "MYCN -- Kelly [Bosse et al. 2017]", "MYCN -- NB1643 [Bosse et al. 2017]",
  "MYCN -- SHEP [Herold et al. 2019]", "MYC-MAX -- MCF-7 [Neikes et al. 2023]",
  "MITF -- 501mel [Louphrasitthiphol et al. 2020]", "HEY1 -- CM7/1 [Weber et al. 2016]",
  "HEY1 -- Cardiomiocytes [Weber et al. 2016]","HEY2 -- CM7/1 [Weber et al. 2016]","HEY2 -- Cardiomiocytes [Weber et al. 2016]",
  "HES1 -- mESC [Shang et al. 2016]","HES1 -- SW1353 [Sugita et al. 2015]",
  "SIM2 -- mESC [Letourneau et al. 2014] (1)","SIM2 -- mESC [Letourneau et al. 2014] (2)","SIM2 -- mESC [Letourneau et al. 2014] (3)",
  "SIM2 -- mESC [Letourneau et al. 2014] (4)","TFEB -- HEK293 [Gambardella et al. 2020]",
  "TFE3 -- HK-2 [Tai-Nagara et al. 2020]","BHLHE40 -- Jurkat [Emming et al. 2020] (1)",
  "BHLHE40 -- Jurkat [Emming et al. 2020] (2)"
)


saveRDS(orden,"orden.RDS")




canonical_eboxes_filtered<-readRDS("canonical_eboxes_filtered.RDS")

library(dplyr)
median_distances<-canonical_eboxes_filtered %>% select(central_dinucleotide,study,dist_summit) %>%
  group_by(study, central_dinucleotide) %>% filter(abs(dist_summit)<=200) %>% 
  summarise(median_dist = median(abs(dist_summit), na.rm = TRUE))


setwd("~/proneural/chip-seq")
canonical_eboxes<-readRDS("combined_canonical_eboxes.RDS")

canonical_eboxes <- canonical_eboxes %>%
  left_join(study_mapping, by = c("study" = "study")) %>%
  mutate(study = renamed_study) %>%
  select(-renamed_study)



library(dplyr)
canonical_eboxes$window_start<-canonical_eboxes$summit-as.integer(canonical_eboxes$peak_length/2)
canonical_eboxes$window_end<-canonical_eboxes$summit+as.integer(canonical_eboxes$peak_length/2)
canonical_eboxes$ebox_center<-canonical_eboxes$ebox_start+3


# Filter motifs that are within the window or keep rows with NA in ebox_center
filtered_data <- canonical_eboxes %>%
  filter(is.na(ebox_center) | (ebox_center >= window_start & ebox_center <= window_end))


filtered_data<-filtered_data %>% select(chrom, window_start, window_end, summit, ebox_center, central_dinucleotide, study, peak_id) 
colnames(filtered_data)[c(2,3)]<-c("peak_start","peak_end")

saveRDS(filtered_data, "combined_canonical_eboxes_filtered.RDS")


canonical_eboxes_filtered_peak<-readRDS("combined_canonical_eboxes_filtered.RDS")

total_eboxes_in_peaks <- canonical_eboxes_filtered_peak %>%
  group_by(study,central_dinucleotide) %>%
  summarise(num_eboxes = n(), .groups = "drop")


# Join total_eboxes_in_peaks and median_distances tables by study and central_dinucleotide
joined_data <- total_eboxes_in_peaks %>%
  inner_join(median_distances, by = c("study", "central_dinucleotide"))


# Normalize both the distance to the summit and num_eboxes with respect to the maximum and minimum in each study
joined_data <- joined_data %>%
  group_by(study) %>%
  filter(!all(is.na(median_dist)) & !all(is.na(num_eboxes))) %>%
  mutate(
    normalized_distance = scales::rescale(median_dist, to = c(0, 1), na.rm = TRUE),
    normalized_num_eboxes = scales::rescale(num_eboxes, to = c(0, 1), na.rm = TRUE)
  ) %>% ungroup()

joined_data$centrality<-1-joined_data$normalized_distance

orden<-readRDS("orden.RDS")
joined_data <- joined_data %>%
  mutate(study = factor(study, levels = rev(unique(orden))) )

joined_data<-joined_data %>% mutate(central_dinucleotide = factor(central_dinucleotide,
levels = c("CAT-CAT","CAT-CAG","CAG-CAG","CAG-CAC","CAC-CAC","CAT-CAC","CAA-CAT","CAA-CAG","CAA-CAC","CAA-CAA") )  )

joined_data$centrality<-1-joined_data$normalized_distance






library(ggplot2)
# Create the balloon plot


pdf("~/Dropbox/induction/1_motifs_qc/balloon_motifs_centrality_and_rate.pdf", height = 30, width = 9, useDingbats = F)
ggplot(joined_data, aes(x = central_dinucleotide, y = study)) +
  geom_point(aes(size = normalized_num_eboxes, color = centrality), alpha = 0.7) +
  scale_size_continuous(range = c(1, 10), name = "Normalized # E-boxes") +
  scale_color_viridis_c(option = "cividis", name = "Normalized Distance") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8)
  ) +
  labs(
    title = "Balloon Plot of E-box Metrics",
    x = "Central Dinucleotide",
    y = "Study"
  )
dev.off()



################## FLANKING NUCLEOTIDE #############################################
#################################################################################### 
### Main figure 1b


setwd("~/proneural/chip-seq")
library(dplyr)

# read all flanking nucleotide files. they are in the human and mouse subdirectories ending in flanking.RDS
flanking_files <- c(
  list.files(path = "~/proneural/chip-seq/human", pattern = "flanking.RDS$", full.names = TRUE),
  list.files(path = "~/proneural/chip-seq/mouse", pattern = "flanking.RDS$", full.names = TRUE)
)
flanking_data <- data.frame()
## read files
for (file in flanking_files) {
  print(file)flanking_data
  data <- readRDS(file)
  flanking_data <- rbind(flanking_data, data)
}

saveRDS(flanking_data, file = "combined_flanking_nucleotides.RDS")


setwd("~/proneural/chip-seq")
library(dplyr)

motifs<-readRDS("combined_flanking_nucleotides.RDS")




study_mapping <- readRDS("study_mapping.RDS")

filtered_studies <- readRDS("filtered_studies.RDS")

filtered_studies<-study_mapping$study[which(study_mapping$renamed_study %in% filtered_studies)]

library(dplyr)
motifs<-motifs %>% filter(study %in% filtered_studies)

# duplicate each entry of the motifs. but in one entry put one half site. and in the other entry put the other half site.



# Duplicate each entry of the motifs, but in one entry put one half site (before the "-")
flanking_data_half1 <- motifs %>%
  mutate(half_site = sub("-.*", "", flanking)) %>%
  select(chrom, start, end, summit, peak_length, peak_id, ebox_start, ebox_seq, real_start, real_end, half_site, study)


flanking_data_half2<- motifs %>%
  mutate(half_site = sub(".*-", "", flanking)) %>%
  select(chrom, start, end, summit, peak_length, peak_id, ebox_start, ebox_seq, real_start, real_end, half_site, study)

flanking_data_halfsite_separated<- rbind(flanking_data_half1, flanking_data_half2)
# Save the combined dataframe to an RDS file
saveRDS(flanking_data_halfsite_separated, file = "flanking_data_halfsite_separated.RDS")


# obtain the number of falf sites per real peak size. (within each study)
half_sites_within_peaks<- flanking_data_halfsite_separated %>% filter(!is.na(half_site)) %>% filter(ebox_start>real_start & (ebox_start+10)<real_end) %>%
  group_by(study, half_site) %>%
  summarise(num_half_sites = n(), .groups = "drop") 





# centrality in the resized 800bp peaks. centrality around the summit
flanking_data_halfsite_separated$ebox_center<-flanking_data_halfsite_separated$ebox_start+5
flanking_data_halfsite_separated$dist_summit<-flanking_data_halfsite_separated$summit-flanking_data_halfsite_separated$ebox_center

half_sites_centrality <- flanking_data_halfsite_separated %>% group_by(study,half_site) %>% summarise(median_dist = median(abs(dist_summit), na.rm = TRUE), .groups = "drop")

half_sites_centrality$centrality<-400-half_sites_centrality$median_dist

half_sites_centrality <- half_sites_centrality %>%
  group_by(study) %>%
  filter(!all(is.na(centrality))) %>%
  mutate(
    normalized_centrality = scales::rescale(centrality, to = c(0, 1), na.rm = TRUE)
  ) %>% ungroup()

# join half_sites_within_peaks and half_sites_centrality by study and half_site. if a field is empty, put 0
# all combinations of study and half_site should be present, even if the count is 0
half_sites_summary <- half_sites_within_peaks %>%
  right_join(half_sites_centrality, by = c("study", "half_site")) %>%
  mutate(num_half_sites = ifelse(is.na(num_half_sites), 0, num_half_sites)) %>%
  select(study, half_site, num_half_sites, normalized_centrality)

# rescale num_half_sites with respect to the maximum and minimum in each study
half_sites_summary <- half_sites_summary %>%
  group_by(study) %>%
  filter(!all(is.na(num_half_sites))) %>%
  mutate(
    normalized_num_half_sites = scales::rescale(num_half_sites, to = c(0, 1), na.rm = TRUE)
  ) %>% ungroup()

# do balloon plot showing centrality with color and number of half sites with size

# reorder variable of flanking dinucleotide. first the ones with CAT, then CAC, then CAG, then CAA (as the las 3 nucleotides)
ordered_halfsites<-c(  unique(half_sites_summary$half_site)[grep("CAT", unique(half_sites_summary$half_site))],
                       unique(half_sites_summary$half_site)[grep("CAG", unique(half_sites_summary$half_site))],
                       unique(half_sites_summary$half_site)[grep("CAC", unique(half_sites_summary$half_site))],
                       unique(half_sites_summary$half_site)[grep("CAA", unique(half_sites_summary$half_site))] )

half_sites_summary$half_site <- factor(half_sites_summary$half_site, levels = ordered_halfsites)

# rename studies
half_sites_summary <- half_sites_summary %>%
  left_join(study_mapping, by = c("study" = "study")) %>%
  mutate(study = renamed_study) %>%
  select(-renamed_study)

# Reorder the studies according to the specified order
orden <- readRDS("orden.RDS")


half_sites_summary <- half_sites_summary %>%
  mutate(study = factor(study, levels = rev(unique(orden))) )

half_sites_summary<-half_sites_summary %>% filter (!is.na(half_site) & !is.na(study))




library(ggplot2)
pdf("~/Dropbox/induction/1_motifs_qc/balloon_flanking_nucleotides_centrality_and_rate.pdf", height = 30, width = 10, useDingbats = F)
ggplot(half_sites_summary, aes(x = half_site, y = study)) +
  geom_point(aes(size = normalized_num_half_sites, color = normalized_centrality), alpha = 0.7) +
  scale_size_continuous(range = c(1, 10), name = "Normalized # Half Sites") +
  scale_color_viridis_c(option = "cividis", name = "Normalized Centrality") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8)
  ) +
  labs(
    title = "Balloon Plot of Flanking Nucleotides Metrics",
    x = "Flanking Nucleotide Half Site",
    y = "Study"
  )
dev.off()






####### NUMBER OF MOTIFS PER PEAK STACKED ########
##################################################
## Main fig 1 c



setwd("~/proneural/chip-seq")
studies<-readRDS("filtered_studies.RDS")

library(dplyr)

npeaks_motifs_join <- data.frame()

for (study in studies) {
  print(study)
  study_files <- c(
    list.files(path = "~/proneural/chip-seq/human", pattern = paste0(study, "_eboxes.RDS$"), full.names = TRUE),
    list.files(path = "~/proneural/chip-seq/mouse", pattern = paste0(study, "_eboxes.RDS$"), full.names = TRUE)
  )
  for (file in study_files) {
    data <- readRDS(file)
    
    data <- data %>%
      mutate(window_start = summit - peak_length / 2,
         window_end = summit + peak_length / 2 )

    data$ebox_center<-data$ebox_start+3
    data <- data %>% filter(ebox_center >= window_start & ebox_center <= window_end)
    

    data <- data %>% filter(grepl("CA..TG", ebox_seq)) %>% distinct()
  
    peak_counts <- data %>%
      group_by(peak_id) %>%
      summarise(num_eboxes = n())
    
        data_original <- readRDS(file)
 
    missing_peaks <- setdiff(data_original$peak_id, peak_counts$peak_id)
    
    if (length(missing_peaks) > 0) {
      missing_data <- data.frame(peak_id = missing_peaks, num_eboxes = 0)
      peak_counts <- rbind(peak_counts, missing_data)
    }

    npeaks_motifs <- peak_counts %>%
      group_by(num_eboxes) %>%
      summarise(n = n())
    
    npeaks_motifs$study <- study

    npeaks_motifs_join <- rbind(npeaks_motifs_join,npeaks_motifs)
  }
}

saveRDS(npeaks_motifs_join, file = "motif_counts_per_peak.RDS")


setwd("~/proneural/chip-seq")
npeaks_motifs_join <- readRDS("motif_counts_per_peak.RDS")
# Plot the number of motifs per peak

# Replace the studies column with the renamed studies
study_mapping <- data.frame(study = studies, renamed_study = studies_renamed)

npeaks_motifs_join <- npeaks_motifs_join %>%
  left_join(study_mapping, by = c("study" = "study")) %>%
  mutate(study = renamed_study) %>%
  select(-renamed_study)


studies_filtered<-readRDS("filtered_studies.RDS")


npeaks_motifs_join<-npeaks_motifs_join %>% filter(study %in% studies_filtered)

library(ggplot2)
# Cap the number of motifs at 4 or more
npeaks_motifs_join$num_eboxes <- ifelse(npeaks_motifs_join$num_eboxes >= 4, ">3", as.character(npeaks_motifs_join$num_eboxes))

# Convert num_eboxes to a factor to maintain the order
npeaks_motifs_join$num_eboxes <- factor(npeaks_motifs_join$num_eboxes, levels = c("0", "1", "2", "3", ">3") )

# Calculate the proportion of peaks with different numbers of motifs
npeaks_motifs_join <- npeaks_motifs_join %>%
  group_by(study) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

# Aggregate rows with num_eboxes == ">3" by study
npeaks_motifs_join <- npeaks_motifs_join %>%
  group_by(study, num_eboxes) %>%
  summarise(n = sum(n), proportion = sum(proportion), .groups = "drop")


orden<-readRDS("orden.RDS")
# Reorder the studies according to the specified order
npeaks_motifs_join <- npeaks_motifs_join %>%
  mutate(study = factor(study, levels = rev(orden) ))




  # Plot the stacked barplot without grey lines inside the bars
  pdf("~/Dropbox/induction/1_motifs_qc/canonical_eboxes_per_peak.pdf", height = 28, width = 10, useDingbats = F)
  ggplot(npeaks_motifs_join, aes(x = study, y = proportion, fill = num_eboxes)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_brewer(palette = "Blues") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 14), axis.text.y = element_text(hjust = 1, size = 14)) +
    labs(title = "Proportion of Peaks with Different Numbers of Motifs",
         x = "Study",
         y = "Proportion",
         fill = "Number of Motifs") +
    coord_flip()
  dev.off()


###### SPACING BETWEEN MOTIFS ###################
################################################
##### Main figure 1d


setwd("~/proneural/chip-seq")
canon<-readRDS("canonical_eboxes_filtered.RDS")

canon$peak_start<-canon$summit-as.integer(canon$peak_length/2)
canon$peak_end<-canon$summit+as.integer(canon$peak_length/2)
canon$ebox_center<-canon$ebox_start+3
library(dplyr)
canonical_eboxes_filtered<-canon %>% filter(ebox_center >= peak_start & ebox_center <= peak_end)

canonical_eboxes_filtered<-canonical_eboxes_filtered %>% distinct()

# Calculate distances between eboxes within each peak
distances <- canonical_eboxes_filtered %>%
  group_by(study, peak_id) %>%
  arrange(ebox_center) %>%
  mutate(distance_to_next = lead(ebox_center) - ebox_center) %>%
  filter(!is.na(distance_to_next) ) %>%
  ungroup()

distances$distance_to_next<-distances$distance_to_next-6


# Save the summary to an RDS file
saveRDS(distances, file = "distance_summary_between_eboxes.RDS")



setwd("~/proneural/chip-seq")
distances<-readRDS("distance_summary_between_eboxes.RDS")

distances<-distances %>% filter(study %in% studies_filtered)
# Print the summary
print(distance_summary)


library(dplyr)
# Filter distances to include only those between 1 and 15
filtered_distances <- distances %>%
  filter(distance_to_next >= 1 & distance_to_next <= 15)

# Calculate the proportion of distances within each study
distance_proportions <- filtered_distances %>%
  group_by(study, distance_to_next) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

# Save the summary to an RDS file
saveRDS(distance_proportions, file = "distance_proportions_summary.RDS")

# Print the summary
print(distance_proportions)

# Load the distance proportions summary
distance_proportions <- readRDS("distance_proportions_summary.RDS")
distance_proportions<-distance_proportions %>% filter(study %in% studies_filtered)


# Reorder the studies according to the specified order
distance_proportions <- distance_proportions %>%
  mutate(study = factor(study, levels = rev(orden)))

# Plot the heatmap
library(ggplot2)
library(viridis)

pdf("~/Dropbox/induction/1_motifs_qc/distance_between_eboxes_heatmap.pdf", height = 10, width = 7, useDingbats = F)
  ggplot(distance_proportions, aes(x = distance_to_next, y = study, fill = proportion)) +
    geom_tile(color = "white") +
    scale_fill_viridis_c(option = "rocket", na.value = "grey50") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Heatmap of Distances Between E-boxes",
        x = "Distance to Next E-box",
        y = "Study",
        fill = "Proportion")
dev.off()



#####################################################################
########### Main figure 1e ##########################################
#####################################################################



# hacerlo para todos los picos en general, sin dividir por quantiles
region_counts_all <- region_counts %>% 
  group_by(study, region_type) %>%
  summarise(count = sum(count), .groups = "drop")

pdf("~/Dropbox/induction/1_motifs_qc/promoters_enhancers_distribution.pdf", height = 6, width = 17, useDingbats = FALSE)
ggplot(region_counts_all %>% mutate(region_type = factor(region_type, levels = c("Distal Enhancer", "Proximal Enhancer", "Promoter"))), 
  aes(x = study, y = count, fill = region_type)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c("Promoter" = "brown1", "Proximal Enhancer" = "goldenrod1", "Distal Enhancer" = "moccasin")) +
  theme_minimal() +
  labs(
    title = "Distribution of Region Types Across Studies",
    x = "Study",
    y = "Proportion of Peaks",
    fill = "Region Type"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.y = element_text(angle = 0)
  )
dev.off()