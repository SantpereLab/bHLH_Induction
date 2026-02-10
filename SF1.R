

######## TEST CANONICAL EBOX ENRICHMENT ########################
#################################################################
### Supplemental fig. 1a




# Function to filter CA..TG eboxes
filter_eboxes <- function(file) {
    data <- readRDS(file)
    filtered_data <- data[ grep("CA..TG", data$ebox_seq), ]
    saveRDS(filtered_data, file = gsub("\\_eboxes.RDS$", "_eboxes_canonical.RDS", file))
}

filter_noncanonical <- function(file) {
    data <- readRDS(file)
    filtered_data <- data[ -grep("CA..TG", data$ebox_seq), ]
    saveRDS(filtered_data, file = gsub("\\_eboxes.RDS$", "_eboxes_noncanonical.RDS", file))
}


setwd("~/proneural/chip-seq/human")

# List all RDS files in the human and mouse directories
human_files <- list.files(pattern = "\\_eboxes.RDS$", full.names = TRUE)

# Apply the filter function to each file using a loop
for (file in human_files) {
    print(file)
    filter_eboxes(file)
    filter_noncanonical(file)
}

setwd("~/proneural/chip-seq/mouse")
mouse_files <- list.files(pattern = "\\_eboxes.RDS$", full.names = TRUE)

# Apply the filter function to each file using a loop
for (file in mouse_files) {
    print(file)
    filter_
setwd("~/proneural/chip-seq/human")

# List all RDS files in the human and mouse directories
human_files <- list.files(pattern = "\\_eboxes.RDS$", full.names = TRUE)
eboxes(file)
    filter_noncanonical(file)
}





setwd("~/proneural/chip-seq")

# coger, para cada induccion, el timepoint con mas picos


studies_renamed <- c(
    "NEUROD1 -- Microglia [Matsuda et al. 2016]", "NEUROD1 -- mESC [Pataskar et al. 2016]", "NEUROD1 -- HEK295 [Kim et al. 2023]",
    "NEUROD2 -- P19 [Fong et al. 2012]", "NEUROD2 -- MEF [Fong et al. 2012]", "NEUROD2 -- P19 [Fong et al. 2015]",
    "NEUROG2 -- EB [Aydin et al. 2019]", "NEUROG2 -- EB [Velasco et al. 2016]",
    "NEUROG2 -- MRC5 [Smith et al. 2016]", "NEUROG2 -- mESC [Vainorius et al. 2023]",
    "NEUROG2 -- Astrocytes [Pereira et al. 2024]", 
    "TWIST2 -- Myoblast (DM) [Li et al. 2019]", "TWIST2 -- Myoblast (GM) [Li et al. 2019]",
    "TWIST1 -- EB [Fan et al. 2012]", "TWIST1 -- MCF7 [Yu et al. 2023]", "TWIST1 -- BMSC [Lee et al. 2023]",
    "TWIST1 -- HEK293 [Kim et al. 2023]", "TWIST1 -- HEK294 [Kim et al. 2023]",
    "OLIG2 -- mESC [Mazzoni et al. 2011]", "ATOH1 -- mESC [Costa et al. 2022]", "ATOH8 -- A549 [Liu et al. 2023]",    
    "MESP1 -- mESC [Lin et al. 2022]","MSGN1 -- EB [Chalamalasetty et al. 2014]",
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
    "MITF -- 501mel [Louphrasitthiphol et al. 2020]","HEY1 -- Cardiomiocytes [Weber et al. 2016]",
    "HEY1 -- CM7/1 [Weber et al. 2016]",
    "HEY2 -- Cardiomiocytes [Weber et al. 2016]",
    "HEY2 -- CM7/1 [Weber et al. 2016]",
    "HES1 -- mESC [Shang et al. 2016]","HES1 -- SW1353 [Sugita et al. 2015]",
    "SIM2 -- mESC [Letourneau et al. 2014] (1)","SIM2 -- mESC [Letourneau et al. 2014] (2)","SIM2 -- mESC [Letourneau et al. 2014] (3)",
    "SIM2 -- mESC [Letourneau et al. 2014] (4)","TFEB -- HEK293 [Gambardella et al. 2020]",
    "TFE3 -- HK-2 [Tai-Nagara et al. 2020]","BHLHE40 -- Jurkat [Emming et al. 2020] (1)",
    "BHLHE40 -- Jurkat [Emming et al. 2020] (2)"
)






studies <- c(
    "matsuda_neurod1_microglia", "pataskar_neurod1_mesc", "kim_neurod1_hek_treatment_pcag",
    "fong_neurod2_p19", "fong_neurod2_mef", "fong_neurod2_p19_treatment",
    "aydin_neurog2_eb_treatment_48h", "velasco_neurog2_eb_treatment_12h",
    "smith_neurog2_mrc5_treatment_1dpt", "vainorius_neurog2_mesc",
    "pereira_neurog2_astrocytes", 
    "li_twist2_myoblast_dm_treatment",
    "li_twist2_myoblast_gm_treatment", "fan_twist1_eb_treatment",
    "yu_twist1_mcf7_treatment", "lee_twist1_bmsc_treatment",
    "kim_twist1_hek_treatment_pcag", "kim_twist1_hek_treatment_pcdna",
    "mazzoni_olig2_mesc", "costa_atoh1_mesc", "liu_atoh8_lungcarcinoma",
    "lin_mesp1_mesc_treatment_24h", "chalamalasetty_msgn1_eb",
    "fong_myod1_p19", "fong_myod1_mef",
    "casey_myod1_mesc", "lee_myod1_mef", "lee_myod1_mesc",
    "conerly_myod1_mef", "manandhar_myod1_fibroblasts", "conerly_myf5_mef",
    "casey_ascl1_mesc", "aydin_ascl1_eb_treatment_12h", "lee_ascl1_mesc" ,"lee_ascl1_mef",
    "wapinski_ascl1_mef", "wapinski_ascl1_npc", "park_ascl1_glioblastoma_treatment",
    "vainorius_ascl1_mesc", "woods_ascl1_neuroblastoma_wt31_treatment",
    "woods_ascl1_neuroblastoma_wt51_treatment", "wang_ascl1_gi-men",
    "joung_ascl1_hesc", "raposo_ascl1_ns5", "casey_ascl2_mesc",
    "liu_ascl2_cd4_treatment", "jain_tcf4_sudhl", "jain_tcf4_tmd8",
    "lorenzin_myc_u2os_treatment", "soufi_myc_bj_treatment",
    "lin_myc_p493-6_treatment_24h", "see_myc_u2os_highmyc",
    "walz_myc_u2os_treatment",
    "barfeld_myc_lncap_treatment",
    "croci_myc_mef", "sabo_myc_mef_treatment", "sabo_myc_lymphoblastoma_highmyc",
    "sabo_myc_lymphoblastoma_treatment_24h", "depretis_myc_mef_treatment_2h",
    "muthalagu_myc_mcf10a", "tameire_myc_colorectaladenocarcinoma",
    "jaenicke_myc_imec", "liang_myc_h2171", "thomas_myc_hek_treatment",
    "jung_myc_u2os_treatment", "zeid_mycn_shep_treatment_6h",
    "upton_mycn_lan5", "upton_mycn_nb1643", "upton_mycn_cogn415",
    "bosse_mycn_ngp", "bosse_mycn_kelly", "bosse_mycn_nb1643",
    "herold_mycn_shep_treatment", "neikes_mycmax_mcf7_treatment_1000nm",
    "louphrasitthiphol_mitf_melanoma_treatment", "weber_hey1_cardiomyocytes_treatment",
    "weber_hey1_mesc_treatment", "weber_hey2_cardiomyocytes_treatment",
    "weber_hey2_mesc_treatment", "shang_hes1_mesc", "sugita_hes1_chrondrocytes",
    "letourneau_sim2_mesc_a6", "letourneau_sim2_mesc_b8", "letourneau_sim2_mesc_c4",
    "letourneau_sim2_mesc_eb3", "gambardella_tfeb_hek_treatment_18h",
    "tai-nagara_tfe3_hk-2", "emming_bhlhe40_jurkat_treatment_clone1",
    "emming_bhlhe40_jurkat_treatment_clone2"
)


canonical_eboxes <- data.frame()

for (study in studies) {
    print(study)
    study_files <- c(
      list.files(path = "~/proneural/chip-seq/human", pattern = paste0("^", study, ".*_eboxes_canonical.RDS$"), full.names = TRUE),
      list.files(path = "~/proneural/chip-seq/mouse", pattern = paste0("^",study, ".*_eboxes_canonical.RDS$"), full.names = TRUE)
    )
    for (file in study_files) {
        data <- readRDS(file)
        canonical_eboxes <- rbind(canonical_eboxes, data)
    }
}


# Save the combined dataframe to an RDS file
saveRDS(canonical_eboxes, file = "combined_canonical_eboxes.RDS")


setwd("~/proneural/chip-seq")

canonical_eboxes<-readRDS("combined_canonical_eboxes.RDS")



library(dplyr)
library(tidyr)
# Replace the studies column with the renamed studies
study_mapping <- data.frame(study = studies, renamed_study = studies_renamed)

saveRDS(study_mapping, "study_mapping.RDS")

canonical_eboxes <- canonical_eboxes %>%
  left_join(study_mapping, by = c("study" = "study")) %>%
  mutate(study = renamed_study) %>%
  select(-renamed_study)


library(dplyr)
canonical_eboxes$dist_summit<-as.numeric(canonical_eboxes$ebox_start)+3-as.numeric(canonical_eboxes$summit)
bind_close<-canonical_eboxes %>% filter(dist_summit<200) %>% filter(dist_summit>-200)

ruskis<-c()
studies<-c()
for (paper in studies_renamed ) {
  print(paper)
  summ_study<-bind_close %>% filter(study==paper)
  ruski<-ks.test(summ_study$dist_summit,alternative = "less", "punif", min = -199, max = 199)
  ruskis<-c(ruskis,ruski$p.value)
  studies<-c(studies,paper)
}

cor_table<-data.frame("ruskis"=ruskis,"studies"=studies)
cor_table<-cor_table %>% arrange(ruskis)


library(dplyr)
############### plotar 

summ_distance<-canonical_eboxes %>% group_by(study,dist_summit) %>%
  summarise(nmotifs=n())

library(tidyheatmaps)
summ_distance <- tibble(
  study = summ_distance$study,
  dist_summit = summ_distance$dist_summit,
  nmotifs = summ_distance$nmotifs
)



library(tidyr)
summ_distance <- summ_distance %>% filter(dist_summit >= -395 & dist_summit <= 395)
summ_distance<-summ_distance %>% complete(study, dist_summit)

summ_distance <- summ_distance %>% mutate(dist_group = floor(dist_summit / 5) * 5) %>% group_by(study, dist_group) %>% summarise(nmotifs = mean(nmotifs, na.rm = TRUE)) %>% rename(dist_summit = dist_group)

summ_distance$nmotifs[is.na(summ_distance$nmotifs)]<-0


library(scales)
summ_distance<-summ_distance %>% group_by(study) %>% mutate(nmotifs_scaled = rescale(nmotifs, to = c(0, 1)))

gc()

library(RColorBrewer)
summ_distance <- summ_distance %>% arrange(factor(study, levels = unique(cor_table$studies)))

library(ggplot2)


# Remove rows with missing values in the 'study' column
summ_distance <- summ_distance %>% filter(!is.na(study))

# Ensure summ_distance is properly formatted
if (!all(c("study", "dist_summit", "nmotifs_scaled") %in% colnames(summ_distance))) {
  stop("summ_distance must contain 'study', 'dist_summit', and 'nmotifs_scaled' columns.")
}

# Generate the heatmap
heatm <- tidyheatmap(
  df = summ_distance,
  rows = study,
  columns = dist_summit,
  values = nmotifs_scaled,
  scale = "none",
  colors = brewer.pal(5, "BuPu")
)

# Check if the heatmap object is valid before saving
if (!is.null(heatm)) {
  pdf("~/Dropbox/induction/1_motifs_qc/all_studies_centrality.pdf", height = 10, width = 6)
  print(heatm)
  dev.off()
} else {
  stop("Heatmap generation failed. Please check the input data.")
}





# Select the studies that appear to have the E-box centrally enriched. A ojo.

central_ebox_studies<-cor_table$studies[1:74]

saveRDS(central_ebox_studies,"filtered_studies.RDS")

central_ebox_studies<-readRDS("filtered_studies.RDS")


canonical_eboxes_filtered<-canonical_eboxes %>% filter(study %in% central_ebox_studies)

saveRDS(canonical_eboxes_filtered,"canonical_eboxes_filtered.RDS")






##### CO-BINDINDING MOTIF ARCHETYPES #########
##############################################
#### Supplementary figure 1c


system("

  module load BEDTools

  cd ~/proneural/chip-seq/mouse

  for file in matsuda_neurod1_microglia_confidence.bed pataskar_neurod1_mesc_confidence.bed fong_neurod2_p19_confidence.bed fong_neurod2_p19_control_confidence.bed fong_neurod2_p19_treatment_confidence.bed fong_neurod2_p19_treatment_confidence.bed fong_neurod2_mef_confidence.bed aydin_neurog2_eb_treatment_48h_confidence.bed velasco_neurog2_eb_treatment_12h_confidence.bed vainorius_neurog2_mesc_confidence.bed pereira_neurog2_astrocytes_confidence.bed pereira_neurog2_astrocytes_yy+_confidence.bed pereira_neurog2_astrocytes_yy-_confidence.bed pereira_neurog2_astrocytes_yy-_confidence.bed pereira_neurog2_astrocytes_yy+_confidence.bed li_twist2_myoblast_dm_treatment_confidence.bed li_twist2_myoblast_gm_treatment_confidence.bed fan_twist1_eb_treatment_confidence.bed mazzoni_olig2_mesc_confidence.bed lin_mesp1_mesc_treatment_24h_confidence.bed chalamalasetty_msgn1_eb_confidence.bed hershbach_myod1_mef_confidence.bed fong_myod1_p19_confidence.bed fong_myod1_mef_confidence.bed casey_myod1_mesc_confidence.bed lee_myod1_mef_confidence.bed lee_myod1_mesc_confidence.bed conerly_myod1_mef_confidence.bed conerly_myf5_mef_confidence.bed casey_ascl1_mesc_confidence.bed aydin_ascl1_eb_treatment_12h_confidence.bed lee_ascl1_mesc_confidence.bed lee_ascl1_mef_confidence.bed wapinski_ascl1_mef_confidence.bed wapinski_ascl1_npc_confidence.bed vainorius_ascl1_mesc_confidence.bed raposo_ascl1_ns5_confidence.bed casey_ascl2_mesc_confidence.bed liu_ascl2_cd4_treatment_confidence.bed kress_myc_liver_treatment_confidence.bed croci_myc_liver_treatment_confidence.bed croci_myc_mef_confidence.bed depretis_myc_mef_treatment_2h_confidence.bed weber_hey1_cardiomyocytes_treatment_confidence.bed weber_hey1_mesc_treatment_confidence.bed weber_hey2_cardiomyocytes_treatment_confidence.bed weber_hey2_mesc_treatment_confidence.bed
  do
    awk '{OFS="\t"; summit=$2+$10; start=summit-200; end=summit+200; if (start<0) start=0; print $1, start, end, $4}' $file > ${file%.bed}_400bp_around_summit.bed
    sbatch --partition=bigmem -w node03 --mem=300000 --wrap="bedtools intersect -a ${file%.bed}_400bp_around_summit.bed -b ~/n/genomes/mm10.archetype_motifs.v1.0.bed -wao | cut -f1,2,3,4,6,7,8 > ${file%.bed}_archetypes.bed"
  done


  cd ~/proneural/chip-seq/human
  for file in kim_neurod1_hek_treatment_pcag_confidence.bed smith_neurog2_mrc5_treatment_1dpt_confidence.bed kim_twist1_hek_treatment_pcag_confidence.bed kim_twist1_hek_treatment_pcdna_confidence.bed manandhar_myod1_fibroblasts_confidence.bed park_ascl1_glioblastoma_treatment_confidence.bed woods_ascl1_neuroblastoma_wt31_treatment_confidence.bed woods_ascl1_neuroblastoma_wt51_treatment_confidence.bed wang_ascl1_gi-men_confidence.bed joung_ascl1_hesc_confidence.bed jain_tcf4_sudhl_confidence.bed jain_tcf4_tmd8_confidence.bed lorenzin_myc_u2os_treatment_confidence.bed soufi_myc_bj_treatment_confidence.bed lin_myc_p493-6_treatment_24h_confidence.bed walz_myc_u2os_treatment_confidence.bed barfeld_myc_lncap_treatment_confidence.bed sabo_myc_lymphoblastoma_highmyc_confidence.bed sabo_myc_lymphoblastoma_treatment_24h_confidence.bed muthalagu_myc_mcf10a_confidence.bed tameire_myc_colorectaladenocarcinoma_confidence.bed jaenicke_myc_imec_confidence.bed liang_myc_h2171_confidence.bed jung_myc_u2os_treatment_confidence.bed zeid_mycn_shep_treatment_6h_confidence.bed upton_mycn_lan5_confidence.bed upton_mycn_nb1643_confidence.bed upton_mycn_cogn415_confidence.bed bosse_mycn_ngp_confidence.bed bosse_mycn_kelly_confidence.bed bosse_mycn_nb1643_confidence.bed herold_mycn_shep_treatment_confidence.bed neikes_mycmax_mcf7_treatment_1000nm_confidence.bed louphrasitthiphol_mitf_melanoma_treatment_confidence.bed emming_bhlhe40_jurkat_treatment_clone1_confidence.bed emming_bhlhe40_jurkat_treatment_clone2_confidence.bed
  do
    awk '{OFS="\t"; summit=$2+$10; start=summit-200; end=summit+200; if (start<0) start=0; print $1, start, end, $4}' $file > ${file%.bed}_400bp_around_summit.bed
    sbatch --partition=bigmem -w node03 --mem=300000 --wrap="bedtools intersect -a ${file%.bed}_400bp_around_summit.bed -b ~/n/genomes/hg38.archetype_motifs.v1.0.bed -wao | cut -f1,2,3,4,6,7,8 > ${file%.bed}_archetypes.bed"
  done

")


system("

  cd ~/proneural/chip-seq
  echo '
    args <- commandArgs(trailingOnly = TRUE)
    archetypes <- read.table(args[1])
    archetypes$summit <- archetypes$V2 + 200
    archetypes$motif_center <- as.integer((archetypes$V5+archetypes$V6)/2 )
    archetypes$dist_summit <- abs(archetypes$motif_center-archetypes$summit)
    library(dplyr)
    archetypes_summary <- archetypes %>%
      group_by(V7) %>%
      summarise(median_distance = median( abs(dist_summit) , na.rm = TRUE) )
    head(archetypes_summary)
    saveRDS(archetypes_summary, file = paste0("summary_", gsub(".bed", "", args[1]), ".RDS") )
  ' > summarize_archetypes.R


  module load R
  cd human
  for file in *archetypes.bed
  do
    sbatch --wrap="Rscript ../summarize_archetypes.R $file"
  done

  cd ../mouse
  for file in *archetypes.bed
  do
    sbatch --partition=normal --wrap="Rscript ../summarize_archetypes.R $file"
  done
")




setwd("~/proneural/chip-seq")

# List all RDS files for human and mouse separately
human_archetype_files <- list.files(path = "~/proneural/chip-seq/human", pattern = "summary.*archetypes\\.RDS$", full.names = TRUE)
mouse_archetype_files <- list.files(path = "~/proneural/chip-seq/mouse", pattern = "summary.*archetypes\\.RDS$", full.names = TRUE)

# Initialize empty dataframes for human and mouse
combined_human_archetypes <- data.frame()
combined_mouse_archetypes <- data.frame()

# Loop through each human file and combine them
for (file in human_archetype_files) {
  data <- readRDS(file)
  data$study <- gsub("summary_|_confidence_archetypes.RDS", "", basename(file))
  combined_human_archetypes <- rbind(combined_human_archetypes, data)
}

# Loop through each mouse file and combine them
for (file in mouse_archetype_files) {
  data <- readRDS(file)
  data$study <- gsub("summary_|_confidence_archetypes.RDS", "", basename(file))
  combined_mouse_archetypes <- rbind(combined_mouse_archetypes, data)
}


# Get unique archetypes for mouse. Because in human there are extra archetypes
archetypes <- unique(combined_mouse_archetypes$V7)[1:233]

combined_archetypes<-rbind(combined_human_archetypes,combined_mouse_archetypes)

# Ensure all archetypes are represented in each study
all_archetypes <- expand.grid(study = unique(combined_archetypes$study), V7 = archetypes)

# Merge with the combined_archetypes dataframe
combined_archetypes <- merge(all_archetypes, combined_archetypes, by = c("study", "V7"), all.x = TRUE)

# Replace NA values in median_distance with NA
combined_archetypes$median_distance[is.na(combined_archetypes$median_distance)] <- NA

# Save the combined dataframe to an RDS file
saveRDS(combined_archetypes, file = "combined_archetypes.RDS")



setwd("~/proneural/chip-seq")

combined_archetypes<-readRDS("combined_archetypes.RDS")
library(dplyr)


# Replace the studies column with the renamed studies
combined_archetypes <- combined_archetypes %>%
  left_join(study_mapping, by = c("study" = "study")) %>%
  mutate(study = renamed_study) %>%
  select(-renamed_study)


combined_archetypes<-combined_archetypes %>% filter(study %in% studies_filtered)



# Identify the top 10 most centrally enriched motifs for each study
top_motifs <- combined_archetypes %>%
  group_by(study) %>%
  top_n(-10, wt = median_distance) %>%
  ungroup()

# Filter the heatmap matrix to include only the top motifs
top_motifs_list <- unique(top_motifs$V7)
combined_archetypes <- combined_archetypes %>% filter(V7 %in% top_motifs_list)


combined_archetypes <- combined_archetypes %>% group_by(study) %>%
mutate(normalized_distance = scales::rescale( median_distance, to = c(0, 1) ) ) %>% ungroup()

combined_archetypes$scaled_centrality<-1-combined_archetypes$normalized_distance

orden<-readRDS("orden.RDS")

combined_archetypes <- combined_archetypes %>%
  mutate(study = factor(study, levels = rev(orden) ))

# Plot the heatmap
library(ggplot2)
library(viridis)

library(pheatmap)

library(tidyr)
# Prepare the data for clustering


heatmap_matrix <- combined_archetypes %>%
  select(study, V7, scaled_centrality) %>%
  pivot_wider(names_from = V7, values_from = scaled_centrality)

# Convert the heatmap_matrix to a matrix
heatmap_matrix <- as.matrix(heatmap_matrix[,-1])  # Remove the study column

# Set the row names to the study names
rownames(heatmap_matrix) <- combined_archetypes$study[!duplicated(combined_archetypes$study)]

heatmap_matrix <- heatmap_matrix[rownames(heatmap_matrix) %in% orden, ]


orden<-orden[which(orden %in% rownames(heatmap_matrix))]
heatmap_matrix <-heatmap_matrix[orden,]

# Manually annotate all the motifs that are canonical or non-canonical E-boxes,
# looking at the motif clusters at the vierstra webpage: https://resources.altius.org/~jvierstra/projects/motif-clustering/releases/v1.0/cluster_viz.html


## motifs that contain a canonical ebox:
"HEN1","ZNF317","Ebox/CAGCTG","ZNF331","ZNF563","Ebox/CAGATGG","Ebox/CATATG", 
"Ebox/CACCTG","SCRT1","SNAI2","MIES","Ebox/CACGTG/1","Ebox/CACGTG/2"
## non canonical eboxes:
"HIF","CREB3/XBP1","GMEB2/2","AHR","ZNF134","MYB/5","OSR2","PRDM9"


colnames(heatmap_matrix) <- sapply(colnames(heatmap_matrix), function(x) {
  switch(x,
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
    "PRDM9" = "PRDM9(YAGYAN)",
    x
  )
})

# pongo los colores en base al central dinucleotide de la Ebox
ebox_colors <- c(
  "HEN1(CAGCTG)" = "steelblue1", "ZNF317(CAGCTG)" = "royalblue1", "Ebox/CAGCTG" = "royalblue3", "ZNF331(CAGCTG)" = "royalblue4", 
  "ZNF563(CAGCTG)" = "navyblue", "Ebox/CAGATGG" = "red1", "Ebox/CATATG" = "black", "Ebox/CACCTG" = "darkolivegreen4", 
  "SCRT1(CAGGTG)" = "olivedrab1", "SNAI2(CAGGTG)" = "chartreuse1", "MIES(CAGGTG)" = "chartreuse3", "Ebox/CACGTG/1" = "goldenrod4", 
  "Ebox/CACGTG/2" = "goldenrod3", "HIF(NACGTG)" = "lightgoldenrod1", "CREB3/XBP1(GACGTG)" = "lemonchiffon3", 
  "GMEB2/2(TACGTA)" = "wheat2", "AHR(NGCGTG)" = "khaki2", "ZNF134(CAGTTG)" = "violetred2", "MYB/5(CAGTTG)" = "violet", 
  "OSR2(CWGCTN)" = "purple3", "PRDM9(YAGYAN)" = "plum4" )

# Assign grey color to archetypes not in the list
all_archetypes <- unique(combined_archetypes$V7)
missing_archetypes <- setdiff(all_archetypes, names(ebox_colors))
ebox_colors[missing_archetypes] <- "grey"

# Create a color annotation for the columns
annotation_col <- data.frame(archetype = colnames(heatmap_matrix))
rownames(annotation_col) <- colnames(heatmap_matrix)
annotation_colors <- list(archetype = ebox_colors)



# Ensure heatmap_matrix is numeric
heatmap_matrix <- as.matrix(heatmap_matrix)

rowww<-rownames(heatmap_matrix)
heatmap_matrix <- apply(heatmap_matrix, 2, as.numeric)
rownames(heatmap_matrix) <- rowww

heatmap_matrix <- heatmap_matrix[orden, ]


# Identify the most centrally enriched archetype for each study
most_central_archetypes <- combined_archetypes %>%
  group_by(study) %>%
  filter(scaled_centrality == max(scaled_centrality, na.rm = TRUE)) %>%
  select(study, V7, scaled_centrality) %>%
  arrange(desc(scaled_centrality))

# Save the results to a file
saveRDS(most_central_archetypes, file = "most_central_archetypes.RDS")

# Print the results
print(most_central_archetypes)

# Generate the clustered heatmap with color annotation and row names
pdf("~/Dropbox/induction/1_motifs_qc/archetypes_centrality.pdf", height = 10, width = 15, useDingbats = F)
  pheatmap(heatmap_matrix, 
    cluster_rows = FALSE, cluster_cols = T, treeheight_col = 0, color = viridis::viridis(100, option = "cividis"), 
    main = "Heatmap of Scaled Centrality for Archetypes", fontsize_row = 8, 
    fontsize_col = 8, show_colnames = TRUE, show_rownames = TRUE, annotation_col = annotation_col, annotation_colors = annotation_colors, 
    annotation_legend = TRUE, annotation_names_col = TRUE, annotation_position = "bottom", labels_row = rownames(heatmap_matrix))
dev.off()


##########  OVERLAP PRE VS POST INDUCTION PEAKS ###########
###################################################### 
######## Supplementary fig. 1c


setwd("~/proneural/chip-seq/human")

### read all the bed files for each study
all_peaks<-data.frame()
for (i in 1:nrow(df)){
  study<-df$studies[i]
  study_renamed<-df$studies_renamed[i]
  print(study)
  bed_file <- Sys.glob(paste0(study, "_confidence.bed"))
  if (length(bed_file) == 0) next
  peaks <- read.table(bed_file[1], header = FALSE)
  peaks$study<-study_renamed
  all_peaks<-rbind(all_peaks,peaks)
}

setwd("~/proneural/chip-seq/mouse")

### read all the bed files for each study
for (i in 1:nrow(df)){
  study<-df$studies[i]
  study_renamed<-df$studies_renamed[i]
  print(study)
  bed_file <- Sys.glob(paste0(study, "_confidence.bed"))
  if (length(bed_file) == 0) next
  peaks <- read.table(bed_file[1], header = FALSE)
  peaks$study<-study_renamed
  all_peaks<-rbind(all_peaks,peaks)
}


all_peaks<-all_peaks %>% distinct()



all_peaks$condition<- ifelse(grepl("control", all_peaks$study), "control", "treatment")
all_peaks$study<-gsub(" control", "", all_peaks$study)
all_peaks$study<-gsub(" treatment", "", all_peaks$study)



library(GenomicRanges)
library(dplyr)

# Prepare peak GRanges for overlap analysis
all_peaks$start <- all_peaks$V2
all_peaks$end <- all_peaks$V3
all_peaks$chr <- all_peaks$V1

# List of studies
studies <- unique(all_peaks$study)

overlap_results <- list()

for (study in studies) {
  control_peaks <- all_peaks %>% filter(study == !!study, condition == "control")
  treatment_peaks <- all_peaks %>% filter(study == !!study, condition == "treatment")
  
  gr_control <- GRanges(seqnames = control_peaks$chr, ranges = IRanges(start = control_peaks$start, end = control_peaks$end))
  gr_treatment <- GRanges(seqnames = treatment_peaks$chr, ranges = IRanges(start = treatment_peaks$start, end = treatment_peaks$end))
  
  # Find overlaps
  overlap_idx_control <- findOverlaps(gr_control, gr_treatment)
  overlap_control <- unique(queryHits(overlap_idx_control))
  overlap_treatment <- unique(subjectHits(overlap_idx_control))
  
  control_peaks$status <- ifelse(seq_len(nrow(control_peaks)) %in% overlap_control, "control+treatment", "only_control")
  treatment_peaks$status <- ifelse(seq_len(nrow(treatment_peaks)) %in% overlap_treatment, "control+treatment", "only_treatment")
  
  overlap_results[[study]] <- bind_rows(control_peaks, treatment_peaks)
}

overlap_peaks <- bind_rows(overlap_results)

# Summarize overlap proportions per study
plot_overlap <- overlap_peaks %>%
  group_by(study, status) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(study) %>%
  mutate(prop = count / sum(count))


plot_overlap<-plot_overlap %>% filter(study!="TWIST1 -- BMSC [Lee et al. 2023]")

setwd("~/proneural/chip-seq")
orden <- readRDS("orden.RDS")
orden<-orden[orden %in% unique(plot_overlap$study)]
plot_overlap$study <- factor(plot_overlap$study, levels = orden)


plot_overlap$status <- factor(plot_overlap$status, levels = c("only_control", "control+treatment","only_treatment") ) 

# Plot stacked barplot of overlap proportions
pdf("~/Dropbox/induction/6_control_vs_treatment/overlap_stacked_barplot.pdf", width=10, height=6, useDingbats=FALSE)
ggplot(plot_overlap, aes(x = study, y = prop, fill = status)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("control+treatment" = "grey70", "only_control" = "dodgerblue", "only_treatment" = "tomato")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Peak Overlap Proportions per Study", x = "Study", y = "Proportion", fill = "Peak Status")
dev.off()




######## PROMOTER VS ENAHNCERS Pre induction vs post induction chip-seq peaks #########
#######################################################################################
################################# Supplementary figure 1d


studies<-c(
"fong_neurod2_p19_control","fong_neurod2_p19_treatment",
"li_twist2_myoblast_dm_control","li_twist2_myoblast_dm_treatment",
"li_twist2_myoblast_gm_control","li_twist2_myoblast_gm_treatment",
"fan_twist1_eb_control","fan_twist1_eb_treatment",
"park_ascl1_glioblastoma_control","park_ascl1_glioblastoma_treatment",
"woods_ascl1_neuroblastoma_wt51_control","woods_ascl1_neuroblastoma_wt51_treatment",
"liu_ascl2_cd4_control","liu_ascl2_cd4_treatment",
"lorenzin_myc_u2os_control","lorenzin_myc_u2os_treatment",
"soufi_myc_bj_control","soufi_myc_bj_treatment",
"lin_myc_p493-6_control","lin_myc_p493-6_treatment_24h",
"walz_myc_u2os_control","walz_myc_u2os_treatment",
"barfeld_myc_lncap_control","barfeld_myc_lncap_treatment",
"croci_myc_liver_control","croci_myc_liver_treatment",
"sabo_myc_lymphoblastoma_lowmyc","sabo_myc_lymphoblastoma_highmyc",
"sabo_myc_lymphoblastoma_control","sabo_myc_lymphoblastoma_treatment_24h",
"depretis_myc_mef_control","depretis_myc_mef_treatment_2h",
"jung_myc_u2os_control","jung_myc_u2os_treatment",
"zeid_mycn_shep_control","zeid_mycn_shep_treatment_6h",
"herold_mycn_shep_control","herold_mycn_shep_treatment",
"weber_hey1_cardiomyocytes_control","weber_hey1_cardiomyocytes_treatment",
"weber_hey1_mesc_control","weber_hey1_mesc_treatment",
"weber_hey2_cardiomyocytes_control","weber_hey2_cardiomyocytes_treatment",
"weber_hey2_mesc_control","weber_hey2_mesc_treatment",
"emming_bhlhe40_jurkat_control_clone1","emming_bhlhe40_jurkat_treatment_clone1",
"emming_bhlhe40_jurkat_control_clone2","emming_bhlhe40_jurkat_treatment_clone2")

studies_renamed<-c(
"NEUROD2 -- P19 [Fong et al. 2012] control", "NEUROD2 -- P19 [Fong et al. 2012] treatment",
"TWIST2 -- Myoblast (DM) [Li et al. 2019] control", "TWIST2 -- Myoblast (DM) [Li et al. 2019] treatment",
"TWIST2 -- Myoblast (GM) [Li et al. 2019] control", "TWIST2 -- Myoblast (GM) [Li et al. 2019] treatment",
"TWIST1 -- EB [Fan et al. 2012] control", "TWIST1 -- EB [Fan et al. 2012] treatment",
"ASCL1 -- G523NS [Park et al. 2017] control", "ASCL1 -- G523NS [Park et al. 2017] treatment",
"ASCL1 -- SH-SY5Y [Woods et al. 2023] (2) control", "ASCL1 -- SH-SY5Y [Woods et al. 2023] (2) treatment",
"ASCL2 -- CD4+ T cells [Liu et al. 2023] control", "ASCL2 -- CD4+ T cells [Liu et al. 2023] treatment",
"MYC -- U2OS [Lorenzin et al. 2016] control", "MYC -- U2OS [Lorenzin et al. 2016] treatment",
"MYC -- BJ [Soufi et al. 2012] control", "MYC -- BJ [Soufi et al. 2012] treatment",
"MYC -- P493-6 [Lin et al. 2012] control", "MYC -- P493-6 [Lin et al. 2012] treatment",
"MYC -- U2OS [Walz et al. 2014] control", "MYC -- U2OS [Walz et al. 2014] treatment",
"MYC -- LNCaP [Barfeld et al. 2017] control", "MYC -- LNCaP [Barfeld et al. 2017] treatment",
"MYC -- MEF [Croci et al. 2017] control", "MYC -- MEF [Croci et al. 2017] treatment",
"MYC -- P493-6 [Sabo et al. 2014] (1) control", "MYC -- P493-6 [Sabo et al. 2014] (1) treatment",
"MYC -- P493-6 [Sabo et al. 2014] (2) control", "MYC -- P493-6 [Sabo et al. 2014] (2) treatment",
"MYC -- MEF [De Pretis et al. 2017] control", "MYC -- MEF [De Pretis et al. 2017] treatment",
"MYC -- U2OS [Jung et al. 2017] control", "MYC -- U2OS [Jung et al. 2017] treatment",
"MYCN -- SHEP [Zeid et al. 2018] control", "MYCN -- SHEP [Zeid et al. 2018] treatment",
"MYCN -- SHEP [Herold et al. 2019] control", "MYCN -- SHEP [Herold et al. 2019] treatment",
"HEY1 -- Cardiomiocytes [Weber et al. 2016] control", "HEY1 -- Cardiomiocytes [Weber et al. 2016] treatment",
"HEY1 -- CM7/1 [Weber et al. 2016] control", "HEY1 -- CM7/1 [Weber et al. 2016] treatment",
"HEY2 -- Cardiomiocytes [Weber et al. 2016] control", "HEY2 -- Cardiomiocytes [Weber et al. 2016] treatment",
"HEY2 -- CM7/1 [Weber et al. 2016] control", "HEY2 -- CM7/1 [Weber et al. 2016] treatment",
"BHLHE40 -- Jurkat [Emming et al. 2020] (1) control", "BHLHE40 -- Jurkat [Emming et al. 2020] (1) treatment",
"BHLHE40 -- Jurkat [Emming et al. 2020] (2) control", "BHLHE40 -- Jurkat [Emming et al. 2020] (2) treatment"
)


df<-data.frame( "studies_renamed"=studies_renamed,"studies"=studies)



## aqui lo he dejado

setwd("~/proneural/chip-seq/human")

### read all the bed files for each study
all_peaks<-data.frame()
for (i in 1:nrow(df)){
  study<-df$studies[i]
  study_renamed<-df$studies_renamed[i]
  print(study)
  bed_file <- Sys.glob(paste0(study, "_confidence.bed"))
  if (length(bed_file) == 0) next
  peaks <- read.table(bed_file[1], header = FALSE)
  peaks$study<-study_renamed
  all_peaks<-rbind(all_peaks,peaks)
}



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



all_peaks$summit<-all_peaks$V2+all_peaks$V10


all_peaks$peakid_new<-paste0(all_peaks$V1,"-",all_peaks$summit)

studies<- unique(all_peaks$study)

for(study in studies){
  print(study)
  idx<-which(all_peaks$study==study)
  all_peaks_dist_study<-dist_tss(all_peaks[idx,c("V1","summit","peakid_new")], tss_hg38)
  all_peaks_dist_study$study<-study
  if (study==studies[1]){
    all_peaks_dist<-all_peaks_dist_study
  } else {
    all_peaks_dist<-rbind(all_peaks_dist, all_peaks_dist_study)
  }
}

saveRDS(all_peaks_dist, file="all_peaks_dist_hg38.rds")



### the same for mouse

tss_mm10<-GenomicFeatures::makeTxDbFromGFF("~/n/genomes/gencode.vM24.annotation.filtered.gtf", format="gtf")

setwd("~/proneural/chip-seq/mouse")

### read all the bed files for each study
all_peaks<-data.frame()
for (i in 1:nrow(df)){
  study<-df$studies[i]
  study_renamed<-df$studies_renamed[i]
  print(study)
  bed_file <- Sys.glob(paste0(study, "_confidence.bed"))
  if (length(bed_file) == 0) next
  peaks <- read.table(bed_file[1], header = FALSE)
  peaks$study<-study_renamed
  all_peaks<-rbind(all_peaks,peaks)
}
all_peaks$summit<-all_peaks$V2+all_peaks$V10
all_peaks$peakid_new<-paste0(all_peaks$V1,"-",all_peaks$summit)
studies<- unique(all_peaks$study)
for(study in studies){
  print(study)
  idx<-which(all_peaks$study==study)
  all_peaks_dist_study<-dist_tss(all_peaks[idx,c("V1","summit","peakid_new")], tss_mm10)
  all_peaks_dist_study$study<-study
  if (study==studies[1]){
    all_peaks_dist<-all_peaks_dist_study
  } else {
    all_peaks_dist<-rbind(all_peaks_dist, all_peaks_dist_study)
  }
}
saveRDS(all_peaks_dist, file="all_peaks_dist_mm10.rds")


peaks_dist_hg38<-readRDS("~/proneural/chip-seq/human/all_peaks_dist_hg38.rds")
peaks_dist_mm10<-readRDS("~/proneural/chip-seq/mouse/all_peaks_dist_mm10.rds")

peaks_dist<-rbind(peaks_dist_hg38, peaks_dist_mm10)



peaks_dist <- peaks_dist %>%
  mutate(region_type = case_when(
    abs(distanceToTSS) <= 2000 ~ "Promoter",
    abs(distanceToTSS) > 2000 & abs(distanceToTSS) <= 4000 ~ "Proximal Enhancer",
    abs(distanceToTSS) > 4000 ~ "Distal Enhancer",
    TRUE ~ NA_character_
))


peaks_dist$condition<- ifelse(grepl("control", peaks_dist$study), "control", "treatment")
peaks_dist$study<-gsub(" control", "", peaks_dist$study)
peaks_dist$study<-gsub(" treatment", "", peaks_dist$study)




library(ggplot2)
library(dplyr)

# Summarize counts per study, condition, and region type
plot_data <- peaks_dist %>%
  group_by(study, condition, region_type) %>%
  summarise(count = n(), .groups = "drop")


plot_data<-plot_data %>%filter(study!="TWIST1 -- BMSC [Lee et al. 2023]")

orden<-readRDS("~/proneural/chip-seq/orden.RDS")

orden<-orden[orden %in% unique(plot_data$study)]

plot_data$study <- factor(plot_data$study, levels = orden)

plot_data$region_type <- factor(plot_data$region_type, levels = c("Distal Enhancer", "Proximal Enhancer", "Promoter"))

# Controls
pdf("~/Dropbox/induction/6_control_vs_treatment/promoter_enhancer_stacked_controls.pdf", width=10, height=6, useDingbats=FALSE)
ggplot(plot_data %>% filter(condition == "control") %>%
     group_by(study) %>%
     mutate(prop = count / sum(count)),
     aes(x = study, y = prop, fill = region_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Promoter" = "brown1", "Proximal Enhancer" = "goldenrod1", "Distal Enhancer" = "moccasin")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Region Type Distribution (Control)", x = "Study", y = "Proportion", fill = "Region Type")
dev.off()




# Treatments
pdf("~/Dropbox/induction/6_control_vs_treatment/promoter_enhancer_stacked_treatments.pdf", width=10, height=6, useDingbats=FALSE)
ggplot(plot_data %>% filter(condition == "treatment") %>%
     group_by(study) %>%
     mutate(prop = count / sum(count)),
     aes(x = study, y = prop, fill = region_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Promoter" = "brown1", "Proximal Enhancer" = "goldenrod1", "Distal Enhancer" = "moccasin")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Region Type Distribution (Treatment)", x = "Study", y = "Proportion", fill = "Region Type")
dev.off()












