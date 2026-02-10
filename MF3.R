
#######################################################
#################### Main figure 3b ###################
#######################################################


########### un plot que sumarice la media y la dispersion.


combineee <- combineee %>% select(peakid_new, meth, quantile,study) %>% distinct()


combineee <- combineee %>%
  mutate(quantile = case_when(
    quantile <= 25 ~ 4,
    quantile > 25 & quantile <= 50 ~ 3,
    quantile > 50 & quantile <= 75 ~ 2,
    quantile > 75 ~ 1 ))

# Adjust meth for studies where some values are greater than 1
combine_corrected <- combineee %>%
  filter(meth != ".") %>% # Remove rows where meth is "."
  mutate(meth = as.numeric(meth)) 

combine_corrected2<-data.frame()
for (estudio in unique(combine_corrected$study) ) {
  rafa <- combine_corrected %>% filter(study == estudio) 
  if (max(rafa$meth) < 2) { rafa$meth <- rafa$meth * 100 }
  combine_corrected2<-rbind(combine_corrected2, rafa)
}

summary_stats <- combine_corrected2 %>%
  mutate(meth = as.numeric(meth)) %>% # Ensure 'meth' is numeric
  group_by(study, quantile) %>%
  summarise(
    mean_meth = mean(meth, na.rm = TRUE),
    sd_meth = sd(meth, na.rm = TRUE),
    .groups = "drop" # Avoid grouped output warning
  )

study_mapping <- readRDS("study_mapping.RDS")
summary_stats$study <- study_mapping[match(summary_stats$study, study_mapping$study), "renamed_study"]

orden<-readRDS("orden.RDS")
orden<-orden[which(orden %in% summary_stats$study)]
summary_stats$study <- factor(summary_stats$study, levels = rev(orden) )
summary_stats <- summary_stats %>% filter(!is.na(study))

viridis_soft <- rev(viridis(100, option = "G")[1:100] )
library(ggplot2)




### the same but plot a heatmap of the median methylation

combineee <- combinee %>% select(peakid_new, meth, quantile,study) %>% distinct()



# Adjust meth for studies where some values are greater than 1
combine_corrected <- combineee %>%
  filter(meth != ".") %>% # Remove rows where meth is "."
  mutate(meth = as.numeric(meth))
combine_corrected2<-data.frame()
for (estudio in unique(combine_corrected$study) ) {
  rafa <- combine_corrected %>% filter(study == estudio) 
  if (max(rafa$meth) < 2) { rafa$meth <- rafa$meth * 100 }
  combine_corrected2<-rbind(combine_corrected2, rafa)
}



summary_stats <- combine_corrected2 %>%
  mutate(meth = as.numeric(meth)) %>% # Ensure 'meth' is numeric
  group_by(study, quantile) %>%
  summarise(
    median_meth = median(meth, na.rm = TRUE),
    .groups = "drop" # Avoid grouped output warning
  )


orden<-readRDS("orden.RDS")
orden<-orden[which(orden %in% summary_stats$study)]
summary_stats$study <- factor(summary_stats$study, levels = rev(orden) )
summary_stats <- summary_stats %>% filter(!is.na(study))
viridis_soft <- rev(viridis(100, option = "G")[1:100] )
library(ggplot2)



pdf("~/Dropbox/induction/4_methylation/summary_median_meth_quantiles_heatmap.pdf",
    useDingbats = F, height = 7, width = 7)
# Create a heatmap
ggplot(summary_stats, aes(x = quantile, y = study, fill = median_meth)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = viridis_soft, name = "Median % mCpG") +
  theme_minimal() +
  labs(
    title = "Heatmap of Median Methylation Levels",
    x = "Quantile",
    y = "Study"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )
dev.off()








#################### Main fig. 3c ####################
######################################################


########## plot only myc and CAT-CAC and CAC-CAC in quantile 4. 

canonical_eboxes_meth_atac<-readRDS("canonical_eboxes_meth_atac.RDS")

canonical_eboxes_meth_atac_myc <- canonical_eboxes_meth_atac[
  grepl("myc|hey1|hey2", canonical_eboxes_meth_atac$study, ignore.case = TRUE),
]

# Convert quantiles to 1-4 (1 = most accessible, 4 = least accessible)
canonical_eboxes_meth_atac_myc$quantile <- dplyr::case_when(
  canonical_eboxes_meth_atac_myc$quantile <= 25 ~ 4,
  canonical_eboxes_meth_atac_myc$quantile > 25 & canonical_eboxes_meth_atac_myc$quantile <= 50 ~ 3,
  canonical_eboxes_meth_atac_myc$quantile > 50 & canonical_eboxes_meth_atac_myc$quantile <= 75 ~ 2,
  canonical_eboxes_meth_atac_myc$quantile > 75 ~ 1
)

canonical_eboxes_meth_atac_4th_quantile<-canonical_eboxes_meth_atac_myc %>% filter(quantile==4)

canonical_eboxes_meth_atac_4th_quantile$meth_window <- factor(canonical_eboxes_meth_atac_4th_quantile$meth_window, levels = c(
  "0-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-90", "90-100"
))

# get number of peaks (number of unique peakid_new) in each meth window for each study
npeaks_meth_window<- canonical_eboxes_meth_atac_4th_quantile %>%
  select(study, meth_window, peakid_new) %>% distinct() %>% group_by(study,meth_window) %>% summarise(npeaks=n())
 
tabla<-canonical_eboxes_meth_atac_4th_quantile %>% select(central_dinucleotide,study,dist_summit,meth_window)


tabla <- tabla %>%
  group_by(study, central_dinucleotide, dist_summit,meth_window) %>%
  summarise(ocurrencias = n() )


library(ggplot2)
library(ggridges)
library(tidyr)


tabla <- tabla %>%
  mutate(bin = floor(dist_summit / 20) * 20) %>% # Crear bins de 5 bp
  group_by(study, central_dinucleotide,meth_window, bin) %>%
  summarise(ocurrencias = sum(ocurrencias, na.rm = TRUE), .groups = "drop") %>%
  rename(dist_summit = bin) 

tabla$dist_summit<- factor(tabla$dist_summit, levels = seq(-400, 400, by = 20)) # Convertir a factor para orden correcto

tabla <- tabla %>%
  group_by(study, central_dinucleotide) %>% # Group by necessary columns
  tidyr::complete(dist_summit,meth_window, fill = list(ocurrencias = 0)) %>% # Fill missing values
  ungroup() # Ungroup after completion



# correct occurences by the number of peaks in the meth windows
tabla <- tabla %>%
  left_join(npeaks_meth_window, by = c("study", "meth_window")) %>%
  mutate(ocurrencias = ocurrencias / npeaks) %>%
  select(-npeaks)


tabla$dist_summit <- as.numeric(as.character(tabla$dist_summit)) # Convert dist_summit to numeric for filtering
tabla <- tabla %>% filter(abs(dist_summit)<=200)

maxii<-tabla %>% group_by(study) %>% summarise(maximo=max(ocurrencias))

tabla<-merge(maxii, tabla, by="study")

tabla$ocurrencias_escaladas<-tabla$ocurrencias/tabla$maximo



# Remove rows with NA central_dinucleotide
tabla <- tabla %>% filter(!is.na(central_dinucleotide) )


study_mapping <- readRDS("study_mapping.RDS")

tabla <- tabla %>%
  mutate(study = study_mapping[match(study, study_mapping$study), "renamed_study"])




orden<-readRDS("orden.RDS")


tabla <- tabla %>%
  mutate(study = factor(study, levels = orden))


tabla<-tabla %>% filter(central_dinucleotide %in% c("CAT-CAC","CAC-CAC")  )



# Define colors for the central dinucleotides
dinucleotide_colors <- c(
   "CAC-CAC" = "goldenrod2", "CAT-CAC" = "darkorange1"
)

library(patchwork)

# Create a list to store plots for each study
plot_list <- list()

for (estudio in unique(tabla$study)) {
  tabla_estudio <- tabla %>% filter(study == estudio)

  gg <- ggplot(tabla_estudio, aes(x = dist_summit, y = ocurrencias_escaladas, group = meth_window, fill = central_dinucleotide)) +
    geom_area(alpha = 0.7, color = "black") +
    facet_grid(meth_window ~ central_dinucleotide, scales = "fixed") +
    scale_fill_manual(values = dinucleotide_colors) +
    theme_void() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      strip.text.y = element_text(angle = 0, size = 10),
      panel.spacing.y = unit(-4, "lines"),
      legend.position = "bottom"
    ) +
    labs(
      title = estudio,
      x = "Distance to summit (binned by 20 bp)",
      y = "Frequency",
      fill = "Dinucleótido central"
    ) + ylim(0, 1)
  plot_list[[estudio]] <- gg
}

# Combine all plots side by side
combined_plot <- wrap_plots(plot_list, nrow = 1, guides = "collect")

pdf("~/Dropbox/induction/4_methylation/caccac_catcac_myc_4th_quantile_meth_windows.pdf", height = 13, width = 5 * length(plot_list), useDingbats = F)
  print(combined_plot)
dev.off()





