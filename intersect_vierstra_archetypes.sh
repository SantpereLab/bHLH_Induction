
########################## VIERSTRA SUMMITS ATAC ##########################
###########################################################################




# Obtain the regions witht he fold change bu with the summits of the peaks in the 4th column.
for file in aydin_ascl1_12h_join.bed aydin_neurog2_48h_join.bed casey_ascl1_24h_join.bed casey_ascl2_24h_join.bed casey_myod1_24h_join.bed lee_ascl1_48h_join.bed lee_myod1_48h_join.bed wapinski_ascl1_48hr_join.bed park_ascl1_dox_join.bed pereira_astrocytes_ngn2_join.bed wang_ascl1_join.bed manandhar_myod1_join.bed
do
  echo $file
  paste bedmap_pre*$file bedmap_post*$file | awk '{if ($11 == 0) print $1,$2,$3,$10,"inf"; else print $1,$2,$3,$10,$22/$11}' | tr " " "\t" > foldchange_induction_summits_${file}
done

for file in wapinski_ascl1_48hr_join.bed
do
  echo $file
  paste bedmap_pre*$file bedmap_post*$file | awk '{if ($11 == 0) print $1,$2,$3,$10,"inf"; else print $1,$2,$3,$10,$22/$11}' | tr " " "\t" > foldchange_induction_summits_${file}
done



## Only the ones that are ATAC-seq

# Intersect with the chip seq peaks using 
cd ~/proneural
bedtools intersect -wa -a foldchange_induction_summits_aydin_ascl1_12h_join.bed -b chip-seq/mouse/aydin_ascl1_eb_treatment_12h_confidence.bed > intersect_foldchange_induction_summits_aydin_ascl1_12h_join.bed
bedtools intersect -wa -a foldchange_induction_summits_aydin_neurog2_48h_join.bed -b chip-seq/mouse/aydin_neurog2_eb_treatment_48h_confidence.bed > intersect_foldchange_induction_summits_aydin_neurog2_48h_join.bed
bedtools intersect -wa -a foldchange_induction_summits_casey_ascl1_24h_join.bed -b chip-seq/mouse/casey_ascl1_mesc_confidence.bed > intersect_foldchange_induction_summits_casey_ascl1_24h_join.bed
bedtools intersect -wa -a foldchange_induction_summits_casey_ascl2_24h_join.bed -b chip-seq/mouse/casey_ascl2_mesc_confidence.bed > intersect_foldchange_induction_summits_casey_ascl2_24h_join.bed
bedtools intersect -wa -a foldchange_induction_summits_casey_myod1_24h_join.bed -b chip-seq/mouse/casey_myod1_mesc_confidence.bed > intersect_foldchange_induction_summits_casey_myod1_24h_join.bed
bedtools intersect -wa -a foldchange_induction_summits_lee_myod1_48h_join.bed -b chip-seq/mouse/lee_myod1_mef_confidence.bed > intersect_foldchange_induction_summits_lee_myod1_48h_join.bed
bedtools intersect -wa -a foldchange_induction_summits_lee_ascl1_48h_join.bed -b chip-seq/mouse/lee_ascl1_mef_confidence.bed > intersect_foldchange_induction_summits_lee_ascl1_48h_join.bed
bedtools intersect -wa -a foldchange_induction_summits_wapinski_ascl1_48hr_join.bed -b chip-seq/mouse/wapinski_ascl1_mef_confidence.bed > intersect_foldchange_induction_summits_wapinski_ascl1_48hr_join.bed
bedtools intersect -wa -a foldchange_induction_summits_park_ascl1_dox_join.bed -b chip-seq/human/park_ascl1_glioblastoma_treatment_confidence.bed > intersect_foldchange_induction_summits_park_ascl1_dox_join.bed
bedtools intersect -wa -a foldchange_induction_summits_pereira_astrocytes_ngn2_join.bed -b chip-seq/mouse/pereira_neurog2_astrocytes_confidence.bed > intersect_foldchange_induction_summits_pereira_astrocytes_ngn2_join.bed
bedtools intersect -wa -a foldchange_induction_summits_wang_ascl1_join.bed -b chip-seq/human/wang_ascl1_gi-men_confidence.bed > intersect_foldchange_induction_summits_wang_ascl1_join.bed


# The ones that do not intersect
bedtools intersect -v -a foldchange_induction_summits_aydin_ascl1_12h_join.bed -b chip-seq/mouse/aydin_ascl1_eb_treatment_12h_confidence.bed > no_intersect_foldchange_induction_summits_aydin_ascl1_12h_join.bed
bedtools intersect -v -a foldchange_induction_summits_aydin_neurog2_48h_join.bed -b chip-seq/mouse/aydin_neurog2_eb_treatment_48h_confidence.bed > no_intersect_foldchange_induction_summits_aydin_neurog2_48h_join.bed
bedtools intersect -v -a foldchange_induction_summits_casey_ascl1_24h_join.bed -b chip-seq/mouse/casey_ascl1_mesc_confidence.bed > no_intersect_foldchange_induction_summits_casey_ascl1_24h_join.bed
bedtools intersect -v -a foldchange_induction_summits_casey_ascl2_24h_join.bed -b chip-seq/mouse/casey_ascl2_mesc_confidence.bed > no_intersect_foldchange_induction_summits_casey_ascl2_24h_join.bed
bedtools intersect -v -a foldchange_induction_summits_casey_myod1_24h_join.bed -b chip-seq/mouse/casey_myod1_mesc_confidence.bed > no_intersect_foldchange_induction_summits_casey_myod1_24h_join.bed
bedtools intersect -v -a foldchange_induction_summits_lee_myod1_48h_join.bed -b chip-seq/mouse/lee_myod1_mef_confidence.bed > no_intersect_foldchange_induction_summits_lee_myod1_48h_join.bed
bedtools intersect -v -a foldchange_induction_summits_lee_ascl1_48h_join.bed -b chip-seq/mouse/lee_ascl1_mef_confidence.bed > no_intersect_foldchange_induction_summits_lee_ascl1_48h_join.bed
bedtools intersect -v -a foldchange_induction_summits_wapinski_ascl1_48hr_join.bed -b chip-seq/mouse/wapinski_ascl1_mef_confidence.bed > no_intersect_foldchange_induction_summits_wapinski_ascl1_48hr_join.bed
bedtools intersect -v -a foldchange_induction_summits_park_ascl1_dox_join.bed -b chip-seq/human/park_ascl1_glioblastoma_treatment_confidence.bed > no_intersect_foldchange_induction_summits_park_ascl1_dox_join.bed
bedtools intersect -v -a foldchange_induction_summits_pereira_astrocytes_ngn2_join.bed -b chip-seq/mouse/pereira_neurog2_astrocytes_confidence.bed > no_intersect_foldchange_induction_summits_pereira_astrocytes_ngn2_join.bed
bedtools intersect -v -a foldchange_induction_summits_wang_ascl1_join.bed -b chip-seq/human/wang_ascl1_gi-men_confidence.bed > no_intersect_foldchange_induction_summits_wang_ascl1_join.bed

# Resize intersect and no_intersect files to 1000bp around the summit
for file in intersect_foldchange_induction_summits_*wap*_join.bed no_intersect_foldchange_induction_summits_*wap*_join.bed
do
  awk 'BEGIN{OFS="\t"} {start=$2+$4-500; end=start+1000; print $1, start, end, $5, $6}' $file | sed 's/\t$//'  > resized_$file
done


cut -f1,2,3,4 ~/n/genomes/mm10.archetype_motifs.v1.0.bed > reduced_mm10.archetype_motifs.v1.0.bed

cut -f1,2,3,4 ~/n/genomes/hg38.archetype_motifs.v1.0.bed > reduced_hg38.archetype_motifs.v1.0.bed



# Mouse
for file in resized*inter*summits*aydin* resized*inter*summits*casey* resized*inter*summits*lee* resized*inter*summits*pereira* resized*inter*summits*wapinski*
do
  outfile="vierstra.$file"
  if [ ! -f "$outfile" ] || [ ! -s "$outfile" ]; then
    sbatch --partition=bigmem --mem=400000 --wrap="bedtools intersect -wao -b reduced_mm10.archetype_motifs.v1.0.bed -a $file > $outfile"
  fi
done

for file in resized*inter*summits*wapinski*
do
  outfile="vierstra.$file"
  if [ ! -f "$outfile" ] || [ ! -s "$outfile" ]; then
    sbatch --partition=bigmem --mem=400000 --wrap="bedtools intersect -wao -b reduced_mm10.archetype_motifs.v1.0.bed -a $file > $outfile"
  fi
done


# Human
for file in resized*inter*summits*park* resized*inter*summits*wang* 
do
  outfile="vierstra.$file"
  if [ ! -f "$outfile" ] || [ ! -s "$outfile" ]; then
    sbatch --partition=bigmem --mem=400000 --wrap="bedtools intersect -wao -b reduced_hg38.archetype_motifs.v1.0.bed -a $file > $outfile"
  fi
done


for file in vierstra.resized_intersect_*wap*.bed 
do
  sbatch --wrap="cut -f1,2,3,4,6,7,8 $file > reduced.$file"
done

for file in vierstra.resized_no_intersect_*wap*.bed 
do
  sbatch --wrap="cut -f1,2,3,4,6,7,8 $file > reduced.$file"
done


# Combine intersect and no_intersect files for each study and add a column indicating intersect status
for study in aydin_ascl1_12h aydin_neurog2_48h casey_ascl1_24h casey_ascl2_24h casey_myod1_24h lee_ascl1_48h lee_myod1_48h wapinski_ascl1_48hr park_ascl1_dox pereira_astrocytes_ngn2 wang_ascl1
do
  echo " $study"
  
  # Add intersect status to intersect files
  sbatch --partition=bigmem --wrap="yes 'yes' | head -n \$(wc -l < reduced.vierstra.resized_intersect_foldchange_induction_summits_${study}_join.bed) > temp_yes_${study}.txt && paste reduced.vierstra.resized_intersect_foldchange_induction_summits_${study}_join.bed temp_yes_${study}.txt > temp_intersect_${study}.bed && rm temp_yes_${study}.txt"
  
  # Add intersect status to no_intersect files
  sbatch --partition=bigmem --wrap="yes 'no' | head -n \$(wc -l < reduced.vierstra.resized_no_intersect_foldchange_induction_summits_${study}_join.bed) > temp_no_${study}.txt && paste reduced.vierstra.resized_no_intersect_foldchange_induction_summits_${study}_join.bed temp_no_${study}.txt > temp_no_intersect_${study}.bed && rm temp_no_${study}.txt"
done


# Combine intersect and no_intersect files for each study
for study in aydin_ascl1_12h aydin_neurog2_48h casey_ascl1_24h casey_ascl2_24h casey_myod1_24h lee_ascl1_48h lee_myod1_48h wapinski_ascl1_48hr park_ascl1_dox pereira_astrocytes_ngn2 wang_ascl1
do
  echo $study
  intersect_file="temp_intersect_${study}.bed"
  no_intersect_file="temp_no_intersect_${study}.bed"
  combined_file="combined_vierstra_${study}.bed"
  sbatch --partition=bigmem --wrap="cat $intersect_file $no_intersect_file > $combined_file"
done



srun --partition=bigmem -n 1 --mpi=none --nodes=1 --pty bash -i

# Add a column that computes the distance between the center of the motifs and the center of the regions
for study in aydin_ascl1_12h aydin_neurog2_48h casey_ascl1_24h casey_ascl2_24h casey_myod1_24h lee_ascl1_48h lee_myod1_48h wapinski_ascl1_48hr park_ascl1_dox pereira_astrocytes_ngn2 wang_ascl1
do
  combined_file="combined_vierstra_${study}.bed"
  output_file="combined_vierstra_with_distance_${study}.bed"
  
  echo $study
  
  sbatch --partition=bigmem --wrap="awk 'BEGIN{OFS=\"\t\"} {
    region_center = int((\$2 + \$3) / 2);
    motif_center = int((\$5 + \$6) / 2);
    distance = region_center - motif_center;
    print \$0, distance
  }' $combined_file > $output_file"
done








# Add a column with peak ID to each combined_vierstra_with_distance_ file
for study in aydin_ascl1_12h aydin_neurog2_48h casey_ascl1_24h casey_ascl2_24h casey_myod1_24h lee_ascl1_48h lee_myod1_48h wapinski_ascl1_48hr park_ascl1_dox pereira_astrocytes_ngn2 wang_ascl1
do
  combined_file="combined_vierstra_with_distance_${study}.bed"
  output_file="combined_vierstra_with_distance_peakid_${study}.bed"
  sbatch --partition=bigmem --wrap="awk 'BEGIN{OFS=\"\t\"} {peak_id = \$1 \":\" \$2 \"-\" \$3; print \$0, peak_id}' $combined_file > $output_file"
done


# Extract the first 4 columns and the peakid, then make unique and sort by column 4
for file in combined_vierstra_with_distance_peakid_*wap*.bed
do
  echo $file
  sbatch --wrap="cut -f1-3,4,10 $file | sort -k4,4 | uniq > unique_${file}"
done

##

# Split the unique_combined_vierstra_with_distance_peakid_ files into 4 subfiles of the same size, taking only column 4

for file in unique_combined_vierstra_with_distance_peakid_*wap*.bed
do
  echo $file
  total_lines=$(wc -l < $file)
  lines_per_file=$(( (total_lines + 3) / 4 ))
  split -l $lines_per_file -d -a 1 <(cut -f5 $file) ${file}_peakids_quantil_
done


for study in aydin_ascl1_12h aydin_neurog2_48h casey_ascl1_24h casey_ascl2_24h casey_myod1_24h lee_ascl1_48h lee_myod1_48h wapinski_ascl1_48hr park_ascl1_dox pereira_astrocytes_ngn2 wang_ascl1
do
  echo $study
  for quantil in {0..3}
  do
  sbatch --wrap="grep -f unique_combined_vierstra_with_distance_peakid_${study}.bed_peakids_quantil_${quantil} combined_vierstra_with_distance_peakid_${study}.bed > juancho_combined_vierstra_with_distance_peakid_${study}_quantil_${quantil}.bed"
  done
done

for study in aydin_ascl1_12h aydin_neurog2_48h casey_ascl1_24h casey_ascl2_24h casey_myod1_24h lee_ascl1_48h lee_myod1_48h wapinski_ascl1_48hr park_ascl1_dox pereira_astrocytes_ngn2 wang_ascl1
do
  for quantil in {0..3}
  do
    file="juancho_combined_vierstra_with_distance_peakid_${study}_quantil_${quantil}.bed"
    echo $file
    awk -F'\t' '$8=="yes"' "$file" > "${file%.bed}_intersect_yes.bed"
    awk -F'\t' '$8=="no"' "$file" > "${file%.bed}_intersect_no.bed"
  done
done





  echo '
  #!/usr/bin/env Rscript

  # Usage: Rscript summarize_motif_archetypes.R <input_file>
  # Example: Rscript summarize_motif_archetypes.R file_intersect_yes.bed

  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1) {
    stop("Usage: Rscript summarize_motif_archetypes.R <input_file>")
  }
  input_file <- args[1]
  output_rds <- sub("\\.bed$", "_motif_archetype_median_summary.rds", input_file)

  library(dplyr)

  # Read the input BED file
  df <- read.table(input_file, header = FALSE, stringsAsFactors = FALSE)
  # Columns: chrom, start, end, fc, motif_start, motif_end, motif_name, intersect, distance_summit, peakid (optional)
  if (ncol(df) == 10) {
    colnames(df) <- c("chrom", "start", "end", "fc", "motif_start", "motif_end", "motif_name", "intersect", "distance_summit", "peakid")
  } else {
    colnames(df) <- c("chrom", "start", "end", "fc", "motif_start", "motif_end", "motif_name", "intersect", "distance_summit")
  }
  df$abs_distance_summit <- abs(as.numeric(df$distance_summit))
  summary <- df %>%
    group_by(motif_name, intersect) %>%
    summarise(
      count = n(),
      median_abs_distance = median(abs_distance_summit, na.rm = TRUE),
      .groups = "drop"
    )
  summary$file <- basename(input_file)

  saveRDS(summary, file = output_rds)
  ' > summarize_motif_archetypes.R


module load R
for file in *wap*_intersect_yes.bed
do
  echo $file
  sbatch --wrap="Rscript summarize_motif_archetypes.R $file"
done

for file in *wap*_intersect_no.bed
do
  echo $file
  sbatch --wrap="Rscript summarize_motif_archetypes.R $file"
done

