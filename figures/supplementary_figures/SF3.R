####### NUMBER OF MOTIFS IN BACKGROUND SEQUENCES (flanking regions of the peaks) ##################
########################################################################################
##### Supplementary figure 3 a



setwd("~/proneural/chip-seq")

upstream_joined_studies<-data.frame()
for (study in studies) {
  paths <- c(
    file.path("~/proneural/chip-seq/human", paste0(study, "_confidence_upstream_half_nochrM_eboxes.RDS")),
    file.path("~/proneural/chip-seq/mouse", paste0(study, "_confidence_upstream_half_nochrM_eboxes.RDS"))
  )
  existing <- paths[file.exists(paths)]
  if (length(existing) == 0) {
    print(paste("No files found for study:", study))
    next
  }

  study_files <- if (length(existing) == 1) {
    readRDS(existing)
  } else {
    do.call(rbind, lapply(existing, readRDS))
  }
  # Add column with study name
  study_files$study <- study
   if (nrow(upstream_joined_studies) == 0) {
      upstream_joined_studies <- study_files
    next
   }
   else {
      upstream_joined_studies <- rbind(upstream_joined_studies, study_files)
   }
}


upstream_nmotifs_per_peak <- upstream_joined_studies %>%
  group_by(study, region) %>%
  summarise(num_eboxes = n()) %>%
  ungroup()




downstream_joined_studies<-data.frame()

for (study in studies) {
  paths <- c(
    file.path("~/proneural/chip-seq/human", paste0(study, "_confidence_downstream_half_nochrM_eboxes.RDS")),
    file.path("~/proneural/chip-seq/mouse", paste0(study, "_confidence_downstream_half_nochrM_eboxes.RDS"))
  )
  existing <- paths[file.exists(paths)]
  if (length(existing) == 0) {
    print(paste("No files found for study:", study))
    next
  }

  study_files <- if (length(existing) == 1) {
    readRDS(existing)
  } else {
    do.call(rbind, lapply(existing, readRDS))
  }
  # Add column with study name
  study_files$study <- study
   if (nrow(downstream_joined_studies) == 0) {
      downstream_joined_studies <- study_files
    next
   }
   else {
      downstream_joined_studies <- rbind(downstream_joined_studies, study_files)
   }
}

downstream_nmotifs_per_peak <- downstream_joined_studies %>%
  group_by(study, region) %>%
  summarise(num_eboxes = n()) %>%
  ungroup()


background_nmotifs_per_peak<- rbind(upstream_nmotifs_per_peak, downstream_nmotifs_per_peak)


background_nmotifs_per_peak$peak_id<- gsub("_UPHALF|_DOWNHALF", "", background_nmotifs_per_peak$region)

## count number of motifs per background peak (joining the two flanking regions)
background_nmotifs_summary <- background_nmotifs_per_peak %>% group_by(study, peak_id) %>% summarise(total_num_eboxes = sum(num_eboxes)) %>% ungroup()


setwd("~/proneural/chip-seq")

##### each study shoould have total_num_eboxes info for all peaks in that study. if a peak not present, put 0 as the count.
### complete missing peaks with 0 counts.

for (study in studies) {
  paths <- c(
    file.path("~/proneural/chip-seq/human", paste0(study, "_confidence.bed")),
    file.path("~/proneural/chip-seq/mouse", paste0(study, "_confidence.bed"))
  )
  existing <- paths[file.exists(paths)]
  if (length(existing) == 0) {
    message("No files found for study: ", study)
    next
  }
  study_files <- if (length(existing) == 1) {
    read.table(existing, header = FALSE, stringsAsFactors = FALSE)
  } else {
    do.call(rbind, lapply(existing, function(x) read.table(x, header = FALSE, stringsAsFactors = FALSE)))
  }
  # 4th column is peak_id
  peak_ids <- study_files[[4]]
  peaks_df <- data.frame(study = study, peak_id = peak_ids, stringsAsFactors = FALSE)
  if (!exists("study_peak_ids")) {
    study_peak_ids <- peaks_df
  } else {
    study_peak_ids <- rbind(study_peak_ids, peaks_df)
  }
}



study_peak_ids<- study_peak_ids %>% filter(study %in% unique(background_nmotifs_summary$study))

### join with background_nmotifs_summary to get all peaks, filling missing with 0
background_nmotifs_complete <- study_peak_ids %>%
  left_join(background_nmotifs_summary, by = c("study", "peak_id")) %>%
  mutate(total_num_eboxes = ifelse(is.na(total_num_eboxes), 0, total_num_eboxes))

### summarize number of peaks with different counts of motifs
background_nmotifs_counts <- background_nmotifs_complete %>%
  group_by(study, total_num_eboxes) %>%
  summarise(n = n()) %>%
  ungroup()




  library(dplyr)
  library(ggplot2)

  # Load mappings and filtered (renamed) studies
  study_mapping <- readRDS("study_mapping.RDS")
  filtered_renamed <- readRDS("filtered_studies.RDS")   # contains renamed study names
  orden <- readRDS("orden.RDS")                         # desired plotting order (renamed)

  # Attach renamed study, keep only those in filtered set
  background_nmotifs_counts <- background_nmotifs_counts %>%
    left_join(study_mapping, by = "study") %>%
    mutate(study = renamed_study) %>%
    select(-renamed_study) %>%
    filter(study %in% filtered_renamed)

  # Collapse counts ≥4
  background_nmotifs_counts$total_num_eboxes <- ifelse(
    background_nmotifs_counts$total_num_eboxes >= 4,
    ">3",
    as.character(background_nmotifs_counts$total_num_eboxes)
  )

  background_nmotifs_counts$total_num_eboxes <- factor(
    background_nmotifs_counts$total_num_eboxes,
    levels = c("0","1","2","3",">3")
  )

  # Compute proportions per study
  background_nmotifs_counts <- background_nmotifs_counts %>%
    group_by(study, total_num_eboxes) %>%
    summarise(n = sum(n), .groups = "drop") %>%
    group_by(study) %>%
    mutate(proportion = n / sum(n)) %>%
    ungroup()

  # Reorder studies
  background_nmotifs_counts <- background_nmotifs_counts %>%
    mutate(study = factor(study, levels = rev(orden)))

  # Plot
  pdf("~/Dropbox/induction/1_motifs_qc/flanking_background_motifs_per_peak_stacked.pdf",
      height = 28, width = 10, useDingbats = FALSE)
  ggplot(background_nmotifs_counts, aes(x = study, y = proportion, fill = total_num_eboxes)) +
    geom_bar(stat = "identity") +
    scale_fill_brewer(palette = "Blues", name = "# motifs") +
    coord_flip() +
    theme_classic() +
    theme(axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14)) +
    labs(title = "Proportion of Background Peak Halves by # Canonical E-boxes",
         x = "Study",
         y = "Proportion")
  dev.off()







#################### SPACING BETWEEN EBOXES IN BACKGROUND SEQUENCES ##################
######################################################################################
######### Supplementary figure 3b

####################### create 100 different sets of background sequences ###############


system("
######## create background sequences

# human
cd ~/proneural/chip-seq/human

module load Miniconda3

echo "
  input_bed=$1
  iteration=$2
  echo $file
  biasaway k -f ${input_bed%.bed}.fa -n 1 > iteration${iteration}_${input_bed%.bed}_background.fa
  biasaway g -f ${input_bed%.bed}.fa -b iteration${iteration}_${input_bed%.bed}_background.fa -r iteration${iteration}_${input_bed%.bed}_biasaway --nfold 1 --length --seed 42 > iteration${iteration}_${input_bed%.bed}_bg.fa
" > create_background.sh



for file in neikes_mycmax_mcf7_treatment_1000nm_confidence.bed barfeld_myc_lncap_treatment_confidence.bed bosse_mycn_kelly_confidence.bed park_ascl1_glioblastoma_treatment_confidence.bed bosse_mycn_nb1643_confidence.bed bosse_mycn_ngp_confidence.bed sabo_myc_lymphoblastoma_highmyc_confidence.bed emming_bhlhe40_jurkat_control_clone1_confidence.bed emming_bhlhe40_jurkat_treatment_clone1_confidence.bed sabo_myc_lymphoblastoma_treatment_24h_confidence.bed emming_bhlhe40_jurkat_treatment_clone2_confidence.bed sabo_myc_mef_treatment_confidence.bed herold_mycn_shep_treatment_confidence.bed see_myc_u2os_highmyc_confidence.bed jaenicke_myc_imec_confidence.bed jain_tcf4_sudhl_confidence.bed jain_tcf4_tmd8_confidence.bed smith_neurog2_mrc5_treatment_1dpt_confidence.bed joung_ascl1_hesc_confidence.bed jung_myc_u2os_treatment_confidence.bed soufi_myc_bj_treatment_confidence.bed kim_neurod1_hek_treatment_pcag_confidence.bed kim_twist1_hek_treatment_pcag_confidence.bed kim_twist1_hek_treatment_pcdna_confidence.bed tameire_myc_colorectaladenocarcinoma_confidence.bed lee_twist1_bmsc_treatment_confidence.bed thomas_myc_hek_treatment_confidence.bed liang_myc_h2171_confidence.bed upton_mycn_cogn415_confidence.bed upton_mycn_lan5_confidence.bed upton_mycn_nb1643_confidence.bed lin_myc_p493-6_treatment_24h_confidence.bed walz_myc_u2os_treatment_confidence.bed lorenzin_myc_u2os_control_confidence.bed wang_ascl1_gi-men_confidence.bed lorenzin_myc_u2os_treatment_confidence.bed wapinski_ascl1_mesc_confidence.bed louphrasitthiphol_mitf_melanoma_control_confidence.bed woods_ascl1_neuroblastoma_wt31_treatment_confidence.bed louphrasitthiphol_mitf_melanoma_treatment_confidence.bed manandhar_myod1_fibroblasts_confidence.bed woods_ascl1_neuroblastoma_wt51_treatment_confidence.bed muthalagu_myc_mcf10a_confidence.bed zeid_mycn_shep_treatment_6h_confidence.bed
do
  for i in {1..100}
  do
    sbatch --wrap="bash create_background.sh $file $i"
  done
done




for file in neikes_mycmax_mcf7_treatment_1000nm_confidence.bed barfeld_myc_lncap_treatment_confidence.bed bosse_mycn_kelly_confidence.bed park_ascl1_glioblastoma_treatment_confidence.bed bosse_mycn_nb1643_confidence.bed bosse_mycn_ngp_confidence.bed sabo_myc_lymphoblastoma_highmyc_confidence.bed emming_bhlhe40_jurkat_control_clone1_confidence.bed emming_bhlhe40_jurkat_treatment_clone1_confidence.bed sabo_myc_lymphoblastoma_treatment_24h_confidence.bed emming_bhlhe40_jurkat_treatment_clone2_confidence.bed sabo_myc_mef_treatment_confidence.bed herold_mycn_shep_treatment_confidence.bed see_myc_u2os_highmyc_confidence.bed jaenicke_myc_imec_confidence.bed jain_tcf4_sudhl_confidence.bed jain_tcf4_tmd8_confidence.bed smith_neurog2_mrc5_treatment_1dpt_confidence.bed joung_ascl1_hesc_confidence.bed jung_myc_u2os_treatment_confidence.bed soufi_myc_bj_treatment_confidence.bed kim_neurod1_hek_treatment_pcag_confidence.bed kim_twist1_hek_treatment_pcag_confidence.bed kim_twist1_hek_treatment_pcdna_confidence.bed tameire_myc_colorectaladenocarcinoma_confidence.bed lee_twist1_bmsc_treatment_confidence.bed thomas_myc_hek_treatment_confidence.bed liang_myc_h2171_confidence.bed upton_mycn_cogn415_confidence.bed upton_mycn_lan5_confidence.bed upton_mycn_nb1643_confidence.bed lin_myc_p493-6_treatment_24h_confidence.bed walz_myc_u2os_treatment_confidence.bed lorenzin_myc_u2os_control_confidence.bed wang_ascl1_gi-men_confidence.bed lorenzin_myc_u2os_treatment_confidence.bed wapinski_ascl1_mesc_confidence.bed louphrasitthiphol_mitf_melanoma_control_confidence.bed woods_ascl1_neuroblastoma_wt31_treatment_confidence.bed louphrasitthiphol_mitf_melanoma_treatment_confidence.bed manandhar_myod1_fibroblasts_confidence.bed woods_ascl1_neuroblastoma_wt51_treatment_confidence.bed muthalagu_myc_mcf10a_confidence.bed zeid_mycn_shep_treatment_6h_confidence.bed
do
  for i in {1..100}
  do
    outfile="iteration${i}_${file%.bed}_bg.fa"
    # Run only if the output file does not exist or is empty
    if [[ ! -s "$outfile" ]]; then
      rm iteration${i}_${file%.bed}_bg.fa
      rm iteration${i}_${file%.bed}_background.fa
      rm -r iteration${i}_${file%.bed}_biasaway
      sbatch --wrap="bash create_background.sh $file $i"
    fi
  done
done



# mouse
cd ~/proneural/chip-seq/mouse


for file in *confidence.bed
do
  for i in {1..100}
  do
    sbatch --partition=bigmem --wrap="bash ../human/create_background.sh $file $i"
  done
done


for file in liu_ascl2_cd4_treatment_confidence.bed matsuda_neurod1_microglia_confidence.bed mazzoni_olig2_mesc_confidence.bed croci_myc_mef_confidence.bed pataskar_neurod1_mesc_confidence.bed pereira_neurog2_astrocytes_confidence.bed depretis_myc_mef_treatment_2h_confidence.bed raposo_ascl1_ns5_confidence.bed vainorius_ascl1_mesc_confidence.bed vainorius_neurog2_mesc_confidence.bed fan_twist1_eb_treatment_confidence.bed velasco_neurog2_eb_treatment_12h_confidence.bed fong_myod1_mef_confidence.bed wapinski_ascl1_mef_confidence.bed fong_myod1_p19_confidence.bed wapinski_ascl1_npc_confidence.bed fong_neurod2_mef_confidence.bed fong_neurod2_p19_confidence.bed weber_hey1_cardiomyocytes_treatment_confidence.bed fong_neurod2_p19_treatment_confidence.bed weber_hey1_mesc_treatment_confidence.bed hershbach_myod1_mef_confidence.bed weber_hey2_cardiomyocytes_treatment_confidence.bed lee_ascl1_mef_confidence.bed weber_hey2_mesc_treatment_confidence.bed lee_ascl1_mesc_confidence.bed lee_myod1_mef_confidence.bed yu_twist1_mcf7_treatment_confidence.bed lee_myod1_mesc_confidence.bed
do
  for i in {1..100}
  do
    outfile="iteration${i}_${file%.bed}_bg.fa"
    if [[ ! -s "$outfile" ]]; then
      sbatch --partition=bigmem --wrap="bash ../human/create_background.sh $file $i"
    fi
  done
done

for file in liu_ascl2_cd4_treatment_confidence.bed matsuda_neurod1_microglia_confidence.bed mazzoni_olig2_mesc_confidence.bed croci_myc_mef_confidence.bed pataskar_neurod1_mesc_confidence.bed pereira_neurog2_astrocytes_confidence.bed depretis_myc_mef_treatment_2h_confidence.bed raposo_ascl1_ns5_confidence.bed vainorius_ascl1_mesc_confidence.bed vainorius_neurog2_mesc_confidence.bed fan_twist1_eb_treatment_confidence.bed velasco_neurog2_eb_treatment_12h_confidence.bed fong_myod1_mef_confidence.bed wapinski_ascl1_mef_confidence.bed fong_myod1_p19_confidence.bed wapinski_ascl1_npc_confidence.bed fong_neurod2_mef_confidence.bed fong_neurod2_p19_confidence.bed weber_hey1_cardiomyocytes_treatment_confidence.bed fong_neurod2_p19_treatment_confidence.bed weber_hey1_mesc_treatment_confidence.bed hershbach_myod1_mef_confidence.bed weber_hey2_cardiomyocytes_treatment_confidence.bed lee_ascl1_mef_confidence.bed weber_hey2_mesc_treatment_confidence.bed lee_ascl1_mesc_confidence.bed lee_myod1_mef_confidence.bed yu_twist1_mcf7_treatment_confidence.bed lee_myod1_mesc_confidence.bed
do
  for i in {1..100}
  do
    outfile="iteration${i}_${file%.bed}_bg.fa"
    # Run only if the output file does not exist or is empty
    if [[ ! -s "$outfile" ]]; then
      rm iteration${i}_${file%.bed}_bg.fa
      rm iteration${i}_${file%.bed}_background.fa
      rm -r iteration${i}_${file%.bed}_biasaway
      sbatch --wrap="bash create_background.sh $file $i"
    fi
  done
done




for i in {1..100}
do
  du -h iteration${i}_park_ascl1_glioblastoma_control_confidence_bg.fa
done



")



library(stringr)
library(dplyr)

system("
  #concatenate fastas of all iterations
  srun --partition=bigmem -n 1 --mpi=none --nodes=1 --pty bash -i
  
  cd ~/proneural/chip-seq/human

  for fasta in barfeld_myc_lncap_control_confidence neikes_mycmax_mcf7_treatment_1000nm_confidence barfeld_myc_lncap_treatment_confidence park_ascl1_glioblastoma_control_confidence bosse_mycn_kelly_confidence park_ascl1_glioblastoma_treatment_confidence bosse_mycn_nb1643_confidence sabo_myc_lymphoblastoma_control_confidence bosse_mycn_ngp_confidence sabo_myc_lymphoblastoma_highmyc_confidence emming_bhlhe40_jurkat_control_clone1_confidence sabo_myc_lymphoblastoma_lowmyc_confidence emming_bhlhe40_jurkat_control_clone2_confidence sabo_myc_lymphoblastoma_treatment_1h_confidence emming_bhlhe40_jurkat_treatment_clone1_confidence sabo_myc_lymphoblastoma_treatment_24h_confidence emming_bhlhe40_jurkat_treatment_clone2_confidence sabo_myc_mef_control_confidence herold_mycn_shep_control_confidence sabo_myc_mef_treatment_confidence herold_mycn_shep_treatment_confidence see_myc_u2os_highmyc_confidence jaenicke_myc_imec_confidence see_myc_u2os_lowmyc_confidence jain_tcf4_sudhl_confidence smith_neurog2_mrc5_treatment_05dpt_confidence jain_tcf4_tmd8_confidence smith_neurog2_mrc5_treatment_1dpt_confidence joung_ascl1_hesc_confidence smith_neurog2_mrc5_treatment_2dpt_confidence jung_myc_u2os_control_confidence soufi_myc_bj_control_confidence jung_myc_u2os_treatment_confidence soufi_myc_bj_treatment_confidence kim_neurod1_hek_treatment_pcag_confidence sugita_hes1_chrondrocytes_confidence kim_twist1_hek_treatment_pcag_confidence tai-nagara_tfe3_hk-2_confidence kim_twist1_hek_treatment_pcdna_confidence tameire_myc_colorectaladenocarcinoma_confidence lee_twist1_bmsc_control_confidence thomas_myc_hek_control_confidence lee_twist1_bmsc_treatment_confidence thomas_myc_hek_treatment_confidence liang_myc_h2171_confidence upton_mycn_cogn415_confidence lin_myc_p493-6_control_confidence upton_mycn_lan5_confidence lin_myc_p493-6_treatment_1h_confidence upton_mycn_nb1643_confidence lin_myc_p493-6_treatment_24h_confidence walz_myc_u2os_control_confidence liu_atoh8_lungcarcinoma_confidence walz_myc_u2os_treatment_confidence lorenzin_myc_u2os_control_confidence wang_ascl1_gi-men_confidence lorenzin_myc_u2os_treatment_confidence wapinski_ascl1_mesc_confidence louphrasitthiphol_mitf_melanoma_control_confidence woods_ascl1_neuroblastoma_wt31_treatment_confidence louphrasitthiphol_mitf_melanoma_treatment_confidence woods_ascl1_neuroblastoma_wt51_control_confidence manandhar_myod1_fibroblasts_confidence woods_ascl1_neuroblastoma_wt51_treatment_confidence muthalagu_myc_mcf10a_confidence zeid_mycn_shep_control_confidence neikes_mycmax_mcf7_treatment_0001nm_confidence zeid_mycn_shep_treatment_2h_confidence neikes_mycmax_mcf7_treatment_0010nm_confidence zeid_mycn_shep_treatment_6h_confidence
  do
    echo $fasta
    cat iteration{1..100}_${fasta}_bg.fa > all_iterations_${fasta}_bg.fa
  done

  cd ~/proneural/chip-seq/mouse
  for fasta in GSE211864_hershbach_myod1_mef_confidence letourneau_sim2_mesc_a6_confidence aydin_ascl1_eb_treatment_12h_confidence letourneau_sim2_mesc_b8_confidence aydin_ascl1_eb_treatment_48h_confidence letourneau_sim2_mesc_c4_confidence aydin_neurog2_eb_treatment_12h_confidence letourneau_sim2_mesc_eb3_confidence aydin_neurog2_eb_treatment_48h_confidence li_twist2_myoblast_dm_control_confidence casey_ascl1_mesc_confidence li_twist2_myoblast_dm_treatment_confidence casey_ascl2_mesc_confidence li_twist2_myoblast_gm_control_confidence casey_myod1_mesc_confidence li_twist2_myoblast_gm_treatment_confidence chalamalasetty_msgn1_eb_confidence lin_mesp1_mesc_treatment_12h_confidence conerly_myf5_mef_confidence lin_mesp1_mesc_treatment_24h_confidence conerly_myod1_mef_confidence liu_ascl2_cd4_control_confidence costa_atoh1_mesc_confidence liu_ascl2_cd4_treatment_confidence croci_myc_liver_control_confidence matsuda_neurod1_microglia_confidence croci_myc_liver_treatment_confidence mazzoni_olig2_mesc_confidence croci_myc_mef_confidence pataskar_neurod1_mesc_confidence depretis_myc_mef_control_confidence pereira_neurog2_astrocytes_confidence depretis_myc_mef_treatment_10min_confidence pereira_neurog2_astrocytes_yy+_confidence depretis_myc_mef_treatment_20min_confidence pereira_neurog2_astrocytes_yy-_confidence depretis_myc_mef_treatment_2h_confidence raposo_ascl1_ns5_confidence depretis_myc_mef_treatment_30min_confidence shang_hes1_mesc_confidence depretis_myc_mef_treatment_4h_confidence vainorius_ascl1_mesc_confidence fan_twist1_eb_control_confidence vainorius_neurog2_mesc_confidence fan_twist1_eb_treatment_confidence velasco_neurog2_eb_treatment_12h_confidence fong_myod1_mef_confidence wapinski_ascl1_mef_confidence fong_myod1_p19_confidence wapinski_ascl1_npc_confidence fong_neurod2_mef_confidence weber_hey1_cardiomyocytes_control_confidence fong_neurod2_p19_confidence weber_hey1_cardiomyocytes_treatment_confidence fong_neurod2_p19_control_confidence weber_hey1_mesc_control_confidence fong_neurod2_p19_treatment_confidence weber_hey1_mesc_treatment_confidence hershbach_myod1_mef_confidence weber_hey2_cardiomyocytes_control_confidence kress_myc_liver_control_confidence weber_hey2_cardiomyocytes_treatment_confidence kress_myc_liver_treatment_confidence weber_hey2_mesc_control_confidence lee_ascl1_mef_confidence weber_hey2_mesc_treatment_confidence lee_ascl1_mesc_confidence yu_twist1_mcf7_control_confidence lee_myod1_mef_confidence yu_twist1_mcf7_treatment_confidence lee_myod1_mesc_confidence
  do
    echo $fasta
    cat iteration{1..100}_${fasta}_bg.fa > all_iterations_${fasta}_bg.fa
  done

")

system("
  srun --partition=bigmem -n 1 --mpi=none --nodes=1 --pty bash -i
  module load R
  R
")

# Function to find E-boxes and calculate spacing
find_ebox_spacing <- function(fasta_file) {
  # Read the FASTA file
  fasta <- readLines(fasta_file)
  
  # Identify sequence headers
  header_idx <- which(grepl("^>", fasta))
  n_seqs <- length(header_idx)
  spacing_counts <- integer(0)
  
  # For each sequence, find E-box positions and spacings
  for (i in seq_len(n_seqs)) {
    start <- header_idx[i] + 1
    end <- if (i < n_seqs) header_idx[i + 1] - 1 else length(fasta)
    seq <- paste(fasta[start:end], collapse = "")
    ebox_positions <- str_locate_all(seq, "CA..TG")[[1]][, 1]
    if (length(ebox_positions) > 1) {
      spacings <- diff(ebox_positions) - 6
      filtered_spacings <- spacings[spacings >= 1 & spacings <= 15]
      spacing_counts <- c(spacing_counts, filtered_spacings)
    }
  }
  
  # Summarize the counts of each spacing
  if (length(spacing_counts) > 0) {
    spacing_summary <- as.data.frame(table(spacing_counts))
    colnames(spacing_summary) <- c("spacing", "count")
  } else {
    spacing_summary <- data.frame(spacing = integer(0), count = integer(0))
  }
  
  # Add the file name to the summary
  spacing_summary$file <- basename(fasta_file)
  
  return(spacing_summary)
}



# Directories containing *_bg.fa files
setwd("~/proneural/chip-seq")
fasta_files <- c(
  list.files(path = "mouse", pattern = "^all_iterations_.*_bg\\.fa$", full.names = TRUE),
  list.files(path = "human", pattern = "^all_iterations_.*_bg\\.fa$", full.names = TRUE)
)

fasta_files <- gsub("^(human/|mouse/)", "", fasta_files)

# Remove "all_iterations_" from the beginning and "_confidence_bg.fa" from the end
fasta_files <- gsub("^all_iterations_", "", fasta_files)
fasta_files <- gsub("_confidence_bg\\.fa$", "", fasta_files)


studies_filtered <- readRDS("filtered_studies.RDS")
study_mapping <- readRDS("study_mapping.RDS")
library(dplyr)
studies_filtered_renamed<- study_mapping %>% filter(renamed_study %in% studies_filtered) %>% pull(study)



library(stringr)
#spacing_results<-data.frame()



# Process each file in a loop and combine results
for (study_name in studies_filtered_renamed[17:length(studies_filtered_renamed)]) {
  # Construct the file path
  fasta_file <- c(
    file.path("mouse", paste0("all_iterations_",study_name, "_confidence_bg.fa")),
    file.path("human", paste0("all_iterations_",study_name, "_confidence_bg.fa"))
  )
  fasta_file <- fasta_file[file.exists(fasta_file)]

  print(paste("Processing file:", fasta_file))
  spacing_result <- find_ebox_spacing(fasta_file)
  if (nrow(spacing_result) > 0) {
    spacing_result$file <- basename(fasta_file)
    # Summarize the counts of each spacing for this file
    spacing_summary <- spacing_result %>%
      group_by(spacing) %>%
      summarise(count = sum(count), .groups = "drop") %>%
      mutate(proportion = count / sum(count) ) %>%
      mutate(study = gsub("_confidence_bg\\.fa", "", basename(fasta_file)))
    # Ensure all expected columns exist
    required_cols <- c("spacing", "count", "proportion", "study")
    for (col in required_cols) {
      if (!col %in% colnames(spacing_summary)) {
        spacing_summary[[col]] <- NA
      }
    }
    spacing_summary <- spacing_summary[, required_cols]
    spacing_results <- rbind(spacing_results, spacing_summary)
  } else {
    print(paste("Skipping file due to no E-boxes found:", fasta_file))
  }

  saveRDS(spacing_results, "ebox_spacing_background.RDS")

}




# Save the results to an RDS file
saveRDS(spacing_results, "ebox_spacing_background.RDS")




spacing_summary<-readRDS("ebox_spacing_background.RDS")

spacing_summary$study <- gsub("^all_iterations_", "", spacing_summary$study)

# Replace the studies column with the renamed studies
spacing_summary <- spacing_summary %>%
  left_join(study_mapping, by = c("study" = "study")) %>%
  mutate(study = renamed_study) %>%
  select(-renamed_study)

# Filter studies to include only the filtered ones
spacing_summary <- spacing_summary %>%
  filter(study %in% studies_filtered)

orden <- readRDS("orden.RDS")
# Reorder the studies according to the specified order
spacing_summary <- spacing_summary %>%
  mutate(study = factor(study, levels = rev(orden)))

library(ggplot2)
# Plot the heatmap
pdf("~/Dropbox/induction/1_motifs_qc/distance_between_eboxes_background_heatmap.pdf", height = 12, width = 7, useDingbats = F)
ggplot(spacing_summary, aes(x = spacing, y = study, fill = proportion)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "rocket", na.value = "grey50") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Heatmap of Distances Between E-boxes in Background Sequences",
       x = "Distance to Next E-box",
       y = "Study",
       fill = "Proportion")
dev.off()








