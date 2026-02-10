##################################################################################
############# Main figure 4a #####################################################
##################################################################################






setwd("~/proneural/mnase")

mesc_dyads_chemical <- read.table("mm10.GSM2183909_unique.map_95pc.bed")

colnames(mesc_dyads_chemical) <- c("chr", "start", "end", "id")
mesc_dyads_chemical$start <- as.numeric(mesc_dyads_chemical$start)
mesc_dyads_chemical$end <- as.numeric(mesc_dyads_chemical$end)
colnames(mesc_dyads_chemical)[1] <- "chrom"
mesc_dyads_chemical <- mesc_dyads_chemical %>% select(chrom, start, end, id)



mesc_dyads_chemical_extended <- mesc_dyads_chemical

mesc_dyads_chemical_extended$start <- mesc_dyads_chemical_extended$start - 200
mesc_dyads_chemical_extended$end <- mesc_dyads_chemical_extended$end + 200


mesc_dyads_extended<-mesc_dyads_chemical_extended

mesc_dyads_extended<-mesc_dyads_extended[,-4]

mesc_dyads_extended$summit<-mesc_dyads_extended$start + 200


library(dplyr)
setwd("~/proneural/chip-seq")
canonical_eboxes<-readRDS("combined_canonical_eboxes.RDS")

# Rename the studies based on a mapping file
study_mapping <- readRDS("study_mapping.RDS")

# Update the study names in the combined_data
canonical_eboxes <- canonical_eboxes %>%
  mutate(study = study_mapping[match(study, study_mapping$study), "renamed_study"] )


studies<-readRDS("filtered_studies.RDS")

mesc_studies<- studies[ grep("mESC|EB|CM7",studies) ]



canonical_eboxes_mesc<- canonical_eboxes %>% filter( study %in% mesc_studies )  


summits_mesc<-canonical_eboxes_mesc %>% select(chrom, summit,study) %>% distinct()
summits_mesc$start<-summits_mesc$summit
summits_mesc$end<-summits_mesc$summit +1
## reorder columns
summits_mesc<- summits_mesc %>% select(chrom, start, end, summit, study)


library(valr)
i <- 1
for (estudio in mesc_studies) {
  print(estudio)
  gc()
  estudio_summits <- summits_mesc %>% filter(study == estudio )
  colnames(mesc_dyads_extended)[1] <- "chrom"
  dyads <- mesc_dyads_extended %>% select(chrom, start, end, summit) %>% distinct()
  interr <- bed_intersect(dyads, estudio_summits)
  interr <- interr %>% select(summit.x, summit.y)
  dist_dyad <- interr$summit.x - interr$summit.y
  dff <- data.frame("dist_dyad" = dist_dyad, "study" = estudio)
  if (i == 1) { data <- dff }
  if (i != 1) { data <- rbind(dff, data) }
  i <- i + 1
}

library(dplyr)

summ_distance <- data %>% group_by(study, dist_dyad) %>%
  summarise(nsummits = n())

library(tidyheatmaps)
summ_distance <- tibble(
  study = summ_distance$study,
  dist_dyad = summ_distance$dist_dyad,
  nsummits = summ_distance$nsummits
)

library(tidyr)
summ_distance <- summ_distance %>% complete(study, dist_dyad)
summ_distance$nsummits[is.na(summ_distance$nsummits)] <- 0

library(scales)
summ_distance <- summ_distance %>% group_by(study)

summ_distance_200bp_around_dyad <- summ_distance %>%
  filter(dist_dyad >= -100 & dist_dyad <= 100)

summ_distance <- summ_distance %>%
  mutate(bin_10bp = 10 * ifelse(dist_dyad >= 0, floor(dist_dyad/10 + 0.5), ceiling(dist_dyad/10 - 0.5))) %>%
  group_by(study, bin_10bp) %>%
  summarise(mean_nsummits = mean(nsummits))

summ_distance_200bp_around_dyad<- summ_distance_200bp_around_dyad %>%
mutate (bin_10bp = 10 * ifelse(dist_dyad >= 0, floor(dist_dyad/10 + 0.5), ceiling(dist_dyad/10 - 0.5))) %>%
group_by(study, bin_10bp) %>%
summarise(mean_nsummits = mean(nsummits))


summ_distance <- summ_distance %>%
  group_by(study) %>%
  mutate(mean_nsummits_scaled = rescale(mean_nsummits, to = c(0, 1)))

summ_distance_200bp_around_dyad <- summ_distance_200bp_around_dyad %>%
  group_by(study) %>%
  mutate(mean_nsummits_scaled = rescale(mean_nsummits, to = c(0, 1)))









######################## oct4 and sox2





oct4<-read.table("~/proneural/mnase/oct4_SRX236476.05.bed") %>% select(V1,V2,V10)

sox2<-read.table("~/proneural/mnase/sox2_SRX236477.05.bed") %>% select(V1,V2,V10)


oct4$summit<- oct4$V2 + oct4$V10

oct4$start<- oct4$summit
oct4$end<- oct4$summit +1
oct4<- oct4 %>% select(V1, start, end, summit)
colnames(oct4)[1]<-"chrom"

sox2$summit<- sox2$V2 + sox2$V10
sox2$start<- sox2$summit
sox2$end<- sox2$summit +1
sox2<- sox2 %>% select(V1, start, end, summit)
colnames(sox2)[1]<-"chrom"

#### plot the distribution of oct4 and sox2 summits around dyads


# For SOX2
colnames(mesc_dyads_extended)[1] <- "chrom"
dyads <- mesc_dyads_extended %>% select(chrom, start, end, summit) %>% distinct()
inter_sox2 <- bed_intersect(dyads, sox2)
inter_sox2 <- inter_sox2 %>% select(summit.x, summit.y)
dist_dyad_sox2 <- inter_sox2$summit.x - inter_sox2$summit.y
data_sox2 <- data.frame("dist_dyad" = dist_dyad_sox2, "tf" = "SOX2")

# For OCT4
inter_oct4 <- bed_intersect(dyads, oct4)
inter_oct4 <- inter_oct4 %>% select(summit.x, summit.y)
dist_dyad_oct4 <- inter_oct4$summit.x - inter_oct4$summit.y
data_oct4 <- data.frame("dist_dyad" = dist_dyad_oct4, "tf" = "OCT4")

# Combine
data_tf <- rbind(data_sox2, data_oct4)

summ_distance_tf <- data_tf %>%
  group_by(tf, dist_dyad) %>%
  summarise(nsummits = n(), .groups = "drop") %>%
  tidyr::complete(tf, dist_dyad, fill = list(nsummits = 0)) %>%
  mutate(bin_10bp = 10 * ifelse(dist_dyad >= 0, floor(dist_dyad/10 + 0.5), ceiling(dist_dyad/10 - 0.5))) %>%
  group_by(tf, bin_10bp) %>%
  summarise(mean_nsummits = mean(nsummits), .groups = "drop_last") %>%
  group_by(tf) %>%
  mutate(mean_nsummits_scaled = scales::rescale(mean_nsummits, to = c(0, 1))) %>%
  ungroup()





summ_distance_200bp_around_dyad

summ_distance_tf_200bp

join<- rbind(summ_distance_200bp_around_dyad,
             summ_distance_tf_200bp %>% mutate(study=tf) %>% select(-tf) )

orden<-readRDS("~/proneural/chip-seq/orden.RDS")

orden<- orden[ orden %in% unique(join$study) ]

orden<- c(orden, c("OCT4"), c("SOX2"))

join$study <- factor(join$study, rev(orden)  )




library(valr)
i <- 1
for (estudio in mesc_studies) {
  print(estudio)
  gc()
  estudio_summits <- summits_mesc %>% filter(study == estudio )
  colnames(mesc_dyads_extended)[1] <- "chrom"
  dyads <- mesc_dyads_extended %>% select(chrom, start, end, summit) %>% distinct()
  interr <- bed_intersect(dyads, estudio_summits)
  interr <- interr %>% select(summit.x, summit.y)
  dist_dyad <- interr$summit.x - interr$summit.y
  dff <- data.frame("dist_dyad" = dist_dyad, "study" = estudio)
  if (i == 1) { data <- dff }
  if (i != 1) { data <- rbind(dff, data) }
  i <- i + 1
}

library(dplyr)

summ_distance <- data %>% group_by(study, dist_dyad) %>%
  summarise(nsummits = n())

library(tidyheatmaps)
summ_distance <- tibble(
  study = summ_distance$study,
  dist_dyad = summ_distance$dist_dyad,
  nsummits = summ_distance$nsummits
)

library(tidyr)
summ_distance <- summ_distance %>% complete(study, dist_dyad)
summ_distance$nsummits[is.na(summ_distance$nsummits)] <- 0

library(scales)
summ_distance <- summ_distance %>% group_by(study)

summ_distance_200bp_around_dyad <- summ_distance %>%
  filter(dist_dyad >= -110 & dist_dyad <= 110)


summ_distance_200bp_around_dyad <- summ_distance_200bp_around_dyad %>%
  group_by(study) %>%
  mutate(nsummits_scaled = rescale(nsummits, to = c(0, 1)))


summ_distance_tf <- data_tf %>% filter(dist_dyad >= -110 & dist_dyad <= 110)  %>%
  group_by(tf, dist_dyad) %>%
  summarise(nsummits = n(), .groups = "drop") %>%
  tidyr::complete(tf, dist_dyad, fill = list(nsummits = 0)) %>%
  group_by(tf) %>%
  mutate(nsummits_scaled = scales::rescale(nsummits, to = c(0, 1))) %>%
  ungroup()

colnames(summ_distance_tf)[1]<-"study"

join_1bp<- rbind(summ_distance_200bp_around_dyad,summ_distance_tf)


saveRDS(join_1bp, "~/proneural/mnase/chipseq_summits_distributions_around_dyads_200bpwindow_joined_1bp.RDS")




join_1bp<-readRDS("~/proneural/mnase/chipseq_summits_distributions_around_dyads_200bpwindow_joined_1bp.RDS")

# Try different ways of smoothing nsummits_scaled values and plotting

library(ggplot2)
library(dplyr)
library(zoo)      # For rollmean
library(signal)   # For sgolayfilt (Savitzky-Golay)
library(tidyr)






# 2. Gaussian kernel smoothing (window = 11, sd = 3)
gaussian_smooth <- function(x, window = 11, sd = 3) {
  if (length(x) < window) return(rep(NA, length(x)))
  kern <- dnorm(seq(-(window-1)/2, (window-1)/2), sd = sd)
  kern <- kern / sum(kern)
  stats::filter(x, kern, sides = 2)
}
join_1bp_smooth_gauss <- join_1bp %>%
  group_by(study) %>%
  arrange(dist_dyad) %>%
  mutate(nsummits_gauss = gaussian_smooth(nsummits_scaled, window = 11, sd = 3)) %>%
  ungroup()


join_1bp_smooth_gauss$study <- factor(join_1bp_smooth_gauss$study, rev(orden)  )

pdf("~/Dropbox/induction/3_mnase/chemical/chipseq_summits_distributions_around_dyads_200bpwindow_joined_gauss.pdf", height = 4, width = 7, useDingbats = F)
ggplot(join_1bp_smooth_gauss, aes(x = dist_dyad, y = study, fill = nsummits_gauss)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(title = "Occurrences of TF summits around the dyad (Gaussian smoothed)",
       y = "TF") +
  scale_x_continuous(
    breaks = c(-100, -70, 0, 70, 100),
    labels = c("-100", "-70", "dyad", "70", "100"),
    expand = c(0, 0)
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8),
    axis.ticks.x = element_line(size = 0.5)
  )
dev.off()














