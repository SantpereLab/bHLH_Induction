cd ~/proneural/chip-seq


echo '
find_noncanonical_mouse <- function(peaks, peaks_file, motif) {
    require(dplyr)
    require(valr)
    refindall <- function(seq, pattern) {
        m <- gregexpr(pattern, seq, perl = TRUE)[[1]]
        if (m[1] == -1) integer(0) else as.integer(m)
    }

    
    system(paste0("bedtools getfasta -fi ~/genomes/mm10.fa -bed ", peaks_file, " | grep -v chr > ", peaks_file, ".fa")) # nolint
    seqs <- read.table(paste0(peaks_file, ".fa"))
    seqs$V1 <- toupper(seqs$V1)
    
    ebox <- data.frame()
    for (i in seq_along(seqs$V1)) {
        ebox_chr <- peaks$chr[i]
        matches <- refindall(seqs$V1[i], motif)
        if (length(matches) == 0) next
        
        for (j in seq_along(matches)) {
            ebox <- rbind(ebox, data.frame(
                chrom = ebox_chr,
                start = peaks$start[i] + matches[j],
                end = peaks$start[i] + matches[j] + 6,
                ebox_seq = substr(seqs$V1[i], matches[j], matches[j] + 5)
            ))
        }
    }
    
    colnames(peaks)[1] <- "chrom"
    ebox <- ebox %>% mutate(start = as.numeric(start), end = as.numeric(end))
    inter <- valr::bed_intersect(peaks, ebox) %>%
        select(chrom, start = start.x, end = end.x, summit = summit.x, peak_length = peak_length.x, peak_id = peak_id.x, ebox_start = start.y, ebox_seq = ebox_seq.y)
    
    return(inter)
}

find_noncanonical_mouse_extended <- function(peaks, peakfile, motif) {
    require(dplyr)
    
    peaks$summit<-peaks$V10+peaks$V2
    peaks<-peaks %>% select(V1,V2,V3,summit)
    colnames(peaks)<-c("chr","start","end","summit")

    peaks <- peaks %>% filter(chr != "chrM")
    npeaks <- nrow(peaks)
    peaks <- peaks %>% mutate(
        peak_id = paste0("peak_", row_number()),
        peak_length = end - start,
        start = summit - 400,
        end = summit + 400
    )
    
    motifs <- do.call(rbind, lapply(seq(1, npeaks, by = 2000), function(i) {
        peaks1 <- peaks[i:min(i + 1999, npeaks), ]
        options(scipen = 999)
        write.table(peaks1, file = paste0("rep_", ceiling(i / 2000), "_", peakfile), sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)
        find_noncanonical_mouse(peaks = peaks1, peaks_file = paste0("rep_", ceiling(i / 2000), "_", peakfile), motif)
    }))



    motifs$central_dinucleotide <- substr(motifs$ebox_seq, 3, 4)
    #pairear los dinucleos
    paired_dinucleotides <- c("GA/TC", "GG/CC", "AA/TT", "CA/TG", "AC/GT", "AG/CT", "GC", "CG", "TA", "AT")
    motifs$central_dinucleotide <- sapply(motifs$central_dinucleotide, function(dinucleotide) {
        matched <- grep(dinucleotide, paired_dinucleotides, value = TRUE)
        if (length(matched) > 0) return(matched[1])
        return(dinucleotide)
    })

    #cambiar la notacion de los dinucleotides por la buena
    motifs$central_dinucleotide[motifs$central_dinucleotide == "GC"] <- "CAG-CAG"
    motifs$central_dinucleotide[motifs$central_dinucleotide == "CG"] <- "CAC-CAC"
    motifs$central_dinucleotide[motifs$central_dinucleotide == "GA/TC"] <- "CAT-CAG"
    motifs$central_dinucleotide[motifs$central_dinucleotide == "CA/TG"] <- "CAT-CAC"
    motifs$central_dinucleotide[motifs$central_dinucleotide == "TA"] <- "CAT-CAT"
    motifs$central_dinucleotide[motifs$central_dinucleotide == "AA/TT"] <- "CAA-CAT"
    motifs$central_dinucleotide[motifs$central_dinucleotide == "AC/GT"] <- "CAA-CAG"
    motifs$central_dinucleotide[motifs$central_dinucleotide == "AG/CT"] <- "CAA-CAC"
    motifs$central_dinucleotide[motifs$central_dinucleotide == "AT"] <- "CAA-CAA"
    motifs$central_dinucleotide[motifs$central_dinucleotide == "GG/CC"] <- "CAG-CAC"
    
    return(motifs)
}


args <- commandArgs(trailingOnly = TRUE)
file <- args[1]

motifs <- c(".A..TG", "C...TG", "CA...G", "CA..T.")
peaks <- read.table(file, header = FALSE)
all_eboxes <- data.frame()
for (motif in motifs) {
    eboxes <- find_noncanonical_mouse_extended(peaks, file, motif)
    all_eboxes <- rbind(all_eboxes, eboxes)
}

estudio<-gsub("_confidence.bed", "", file)
all_eboxes$study<-estudio
saveRDS(all_eboxes,paste0(estudio,"_eboxes.RDS") )
' > find_eboxes_mouse.R



echo '
find_noncanonical_human <- function(peaks, peaks_file, motif) {
    require(dplyr)
    require(valr)
    refindall <- function(seq, pattern) {
        m <- gregexpr(pattern, seq, perl = TRUE)[[1]]
        if (m[1] == -1) integer(0) else as.integer(m)
    }



    system(paste0("bedtools getfasta -fi ~/genomes/hg38.fa -bed ", peaks_file, " | grep -v chr > ", peaks_file, ".fa")) # nolint
    seqs <- read.table(paste0(peaks_file, ".fa"))
    seqs$V1 <- toupper(seqs$V1)
    ebox <- data.frame()
    for (i in seq_along(seqs$V1)) {
        ebox_chr <- peaks$chr[i]
        matches <- refindall(seqs$V1[i], motif)
        if (length(matches) == 0) next
        
        for (j in seq_along(matches)) {
            ebox <- rbind(ebox, data.frame(
                chrom = ebox_chr,
                start = peaks$start[i] + matches[j],
                end = peaks$start[i] + matches[j] + 6,
                ebox_seq = substr(seqs$V1[i], matches[j], matches[j] + 5)
            ))
        }
    }
    colnames(peaks)[1] <- "chrom"
    ebox <- ebox %>% mutate(start = as.numeric(start), end = as.numeric(end))
    print("Intersecting peaks with ebox")
    inter <- valr::bed_intersect(peaks, ebox) %>%
        select(chrom, start = start.x, end = end.x, summit = summit.x, peak_length = peak_length.x, peak_id = peak_id.x, ebox_start = start.y, ebox_seq = ebox_seq.y)
    
    print("Returning intersected data")
    return(inter)
}

find_noncanonical_human_extended <- function(peaks, peakfile, motif) {
    require(dplyr)
    
    peaks$summit<-peaks$V10+peaks$V2
    peaks<-peaks %>% select(V1,V2,V3,summit)
    colnames(peaks)<-c("chr","start","end","summit")

    peaks <- peaks %>% filter(chr != "chrM")
    npeaks <- nrow(peaks)
    peaks <- peaks %>% mutate(
        peak_id = paste0("peak_", row_number()),
        peak_length = end - start,
        start = summit - 400,
        end = summit + 400
    )
    
    motifs <- do.call(rbind, lapply(seq(1, npeaks, by = 2000), function(i) {
        peaks1 <- peaks[i:min(i + 1999, npeaks), ]
        options(scipen = 999)
        write.table(peaks1, file = paste0("rep_", ceiling(i / 2000), "_", peakfile), sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)
        find_noncanonical_human(peaks = peaks1, peaks_file = paste0("rep_", ceiling(i / 2000), "_", peakfile), motif)
    }))

    motifs$central_dinucleotide <- substr(motifs$ebox_seq, 3, 4)
    #pairear los dinucleos
    paired_dinucleotides <- c("GA/TC", "GG/CC", "AA/TT", "CA/TG", "AC/GT", "AG/CT", "GC", "CG", "TA", "AT")
    motifs$central_dinucleotide <- sapply(motifs$central_dinucleotide, function(dinucleotide) {
        matched <- grep(dinucleotide, paired_dinucleotides, value = TRUE)
        if (length(matched) > 0) return(matched[1])
        return(dinucleotide)
    })

    #cambiar la notacion de los dinucleotides por la buena
    motifs$central_dinucleotide[motifs$central_dinucleotide == "GC"] <- "CAG-CAG"
    motifs$central_dinucleotide[motifs$central_dinucleotide == "CG"] <- "CAC-CAC"
    motifs$central_dinucleotide[motifs$central_dinucleotide == "GA/TC"] <- "CAT-CAG"
    motifs$central_dinucleotide[motifs$central_dinucleotide == "CA/TG"] <- "CAT-CAC"
    motifs$central_dinucleotide[motifs$central_dinucleotide == "TA"] <- "CAT-CAT"
    motifs$central_dinucleotide[motifs$central_dinucleotide == "AA/TT"] <- "CAA-CAT"
    motifs$central_dinucleotide[motifs$central_dinucleotide == "AC/GT"] <- "CAA-CAG"
    motifs$central_dinucleotide[motifs$central_dinucleotide == "AG/CT"] <- "CAA-CAC"
    motifs$central_dinucleotide[motifs$central_dinucleotide == "AT"] <- "CAA-CAA"
    motifs$central_dinucleotide[motifs$central_dinucleotide == "GG/CC"] <- "CAG-CAC"
    
    return(motifs)
}

args <- commandArgs(trailingOnly = TRUE)
file <- args[1]

motifs <- c(".A..TG", "C...TG", "CA...G", "CA..T.")
peaks <- read.table(file, header = FALSE)
all_eboxes <- data.frame()
for (motif in motifs) {
    eboxes <- find_noncanonical_human_extended(peaks, file, motif)
    all_eboxes <- rbind(all_eboxes, eboxes)
}

estudio<-gsub("_confidence.bed", "", file)
all_eboxes$study<-estudio
saveRDS(all_eboxes,paste0(estudio,"_eboxes.RDS") )
' > find_eboxes_human.R



module load R
module load BEDTools

cd ~/proneural/chip-seq

cd human

for file in kim_neurod1_hek_treatment_pcag_confidence.bed smith_neurog2_mrc5_treatment_1dpt_confidence.bed kim_twist1_hek_treatment_pcag_confidence.bed kim_twist1_hek_treatment_pcdna_confidence.bed manandhar_myod1_fibroblasts_confidence.bed park_ascl1_glioblastoma_treatment_confidence.bed woods_ascl1_neuroblastoma_wt31_treatment_confidence.bed woods_ascl1_neuroblastoma_wt51_treatment_confidence.bed wang_ascl1_gi-men_confidence.bed joung_ascl1_hesc_confidence.bed jain_tcf4_sudhl_confidence.bed jain_tcf4_tmd8_confidence.bed lorenzin_myc_u2os_treatment_confidence.bed soufi_myc_bj_treatment_confidence.bed lin_myc_p493-6_treatment_24h_confidence.bed walz_myc_u2os_treatment_confidence.bed barfeld_myc_lncap_treatment_confidence.bed sabo_myc_lymphoblastoma_highmyc_confidence.bed sabo_myc_lymphoblastoma_treatment_24h_confidence.bed muthalagu_myc_mcf10a_confidence.bed tameire_myc_colorectaladenocarcinoma_confidence.bed jaenicke_myc_imec_confidence.bed liang_myc_h2171_confidence.bed jung_myc_u2os_treatment_confidence.bed zeid_mycn_shep_treatment_6h_confidence.bed upton_mycn_lan5_confidence.bed upton_mycn_nb1643_confidence.bed upton_mycn_cogn415_confidence.bed bosse_mycn_ngp_confidence.bed bosse_mycn_kelly_confidence.bed bosse_mycn_nb1643_confidence.bed herold_mycn_shep_treatment_confidence.bed neikes_mycmax_mcf7_treatment_1000nm_confidence.bed louphrasitthiphol_mitf_melanoma_treatment_confidence.bed emming_bhlhe40_jurkat_treatment_clone1_confidence.bed emming_bhlhe40_jurkat_treatment_clone2_confidence.bed
do
    sbatch --partition=normal --wrap="Rscript ../find_eboxes_human.R $file"
done


rm rep_*bed
rm rep_*fa

cd ../mouse

for file in lee_ascl1_mesc_confidence.bed matsuda_neurod1_microglia_confidence.bed pataskar_neurod1_mesc_confidence.bed fong_neurod2_p19_confidence.bed fong_neurod2_p19_control_confidence.bed fong_neurod2_p19_treatment_confidence.bed fong_neurod2_p19_treatment_confidence.bed fong_neurod2_mef_confidence.bed aydin_neurog2_eb_treatment_48h_confidence.bed velasco_neurog2_eb_treatment_12h_confidence.bed vainorius_neurog2_mesc_confidence.bed pereira_neurog2_astrocytes_confidence.bed pereira_neurog2_astrocytes_yy+_confidence.bed pereira_neurog2_astrocytes_yy-_confidence.bed pereira_neurog2_astrocytes_yy-_confidence.bed pereira_neurog2_astrocytes_yy+_confidence.bed li_twist2_myoblast_dm_treatment_confidence.bed li_twist2_myoblast_gm_treatment_confidence.bed fan_twist1_eb_treatment_confidence.bed mazzoni_olig2_mesc_confidence.bed lin_mesp1_mesc_treatment_24h_confidence.bed chalamalasetty_msgn1_eb_confidence.bed hershbach_myod1_mef_confidence.bed fong_myod1_p19_confidence.bed fong_myod1_mef_confidence.bed casey_myod1_mesc_confidence.bed lee_myod1_mef_confidence.bed lee_myod1_mesc_confidence.bed conerly_myod1_mef_confidence.bed conerly_myf5_mef_confidence.bed casey_ascl1_mesc_confidence.bed aydin_ascl1_eb_treatment_12h_confidence.bed lee_ascl1_mesc_confidence.bed lee_ascl1_mef_confidence.bed wapinski_ascl1_mef_confidence.bed wapinski_ascl1_npc_confidence.bed vainorius_ascl1_mesc_confidence.bed raposo_ascl1_ns5_confidence.bed casey_ascl2_mesc_confidence.bed liu_ascl2_cd4_treatment_confidence.bed kress_myc_liver_treatment_confidence.bed croci_myc_liver_treatment_confidence.bed croci_myc_mef_confidence.bed depretis_myc_mef_treatment_2h_confidence.bed weber_hey1_cardiomyocytes_treatment_confidence.bed weber_hey1_mesc_treatment_confidence.bed weber_hey2_cardiomyocytes_treatment_confidence.bed weber_hey2_mesc_treatment_confidence.bed
do
    sbatch --partition=normal --wrap="Rscript ../find_eboxes_mouse.R $file"
done


rm rep_*bed
rm rep_*fa



make_half_flank_peaks() {
    in="$1"
    base="${in%.bed}"
    up_file="${base}_upstream_half.bed"
    down_file="${base}_downstream_half.bed"
    : > "$up_file"
    : > "$down_file"
    awk -v up="$up_file" -v down="$down_file" 'BEGIN{OFS="\t"}{
        len=$3-$2
        h=int(len/2)
        up_start=$2-h
        if(up_start<0) up_start=0
        up_end=$2
        down_start=$3
        down_end=$3+h
        if(NF>=4){
            up_name=$4"_UPHALF"
            down_name=$4"_DOWNHALF"
        } else {
            up_name="UPHALF"
            down_name="DOWNHALF"
        }
        up_line=$1 OFS up_start OFS up_end OFS up_name
        down_line=$1 OFS down_start OFS down_end OFS down_name
        for(i=5;i<=NF;i++){
            up_line=up_line OFS $i
            down_line=down_line OFS $i
        }
        print up_line >> up
        print down_line >> down
    }' "$in"
}



cd ~/proneural/chip-seq/mouse
for file in lee_ascl1_mesc_confidence.bed matsuda_neurod1_microglia_confidence.bed pataskar_neurod1_mesc_confidence.bed fong_neurod2_p19_confidence.bed fong_neurod2_p19_control_confidence.bed fong_neurod2_p19_treatment_confidence.bed fong_neurod2_p19_treatment_confidence.bed fong_neurod2_mef_confidence.bed aydin_neurog2_eb_treatment_48h_confidence.bed velasco_neurog2_eb_treatment_12h_confidence.bed vainorius_neurog2_mesc_confidence.bed pereira_neurog2_astrocytes_confidence.bed pereira_neurog2_astrocytes_yy+_confidence.bed pereira_neurog2_astrocytes_yy-_confidence.bed pereira_neurog2_astrocytes_yy-_confidence.bed pereira_neurog2_astrocytes_yy+_confidence.bed li_twist2_myoblast_dm_treatment_confidence.bed li_twist2_myoblast_gm_treatment_confidence.bed fan_twist1_eb_treatment_confidence.bed mazzoni_olig2_mesc_confidence.bed lin_mesp1_mesc_treatment_24h_confidence.bed chalamalasetty_msgn1_eb_confidence.bed hershbach_myod1_mef_confidence.bed fong_myod1_p19_confidence.bed fong_myod1_mef_confidence.bed casey_myod1_mesc_confidence.bed lee_myod1_mef_confidence.bed lee_myod1_mesc_confidence.bed conerly_myod1_mef_confidence.bed conerly_myf5_mef_confidence.bed casey_ascl1_mesc_confidence.bed aydin_ascl1_eb_treatment_12h_confidence.bed lee_ascl1_mesc_confidence.bed lee_ascl1_mef_confidence.bed wapinski_ascl1_mef_confidence.bed wapinski_ascl1_npc_confidence.bed vainorius_ascl1_mesc_confidence.bed raposo_ascl1_ns5_confidence.bed casey_ascl2_mesc_confidence.bed liu_ascl2_cd4_treatment_confidence.bed kress_myc_liver_treatment_confidence.bed croci_myc_liver_treatment_confidence.bed croci_myc_mef_confidence.bed depretis_myc_mef_treatment_2h_confidence.bed weber_hey1_cardiomyocytes_treatment_confidence.bed weber_hey1_mesc_treatment_confidence.bed weber_hey2_cardiomyocytes_treatment_confidence.bed weber_hey2_mesc_treatment_confidence.bed
do
    echo $file
    make_half_flank_peaks "$file"
done

head li_twist2_myoblast_dm_treatment_confidence_upstream_half.bed
cd ~/proneural/chip-seq/human
for file in kim_neurod1_hek_treatment_pcag_confidence.bed smith_neurog2_mrc5_treatment_1dpt_confidence.bed kim_twist1_hek_treatment_pcag_confidence.bed kim_twist1_hek_treatment_pcdna_confidence.bed manandhar_myod1_fibroblasts_confidence.bed park_ascl1_glioblastoma_treatment_confidence.bed woods_ascl1_neuroblastoma_wt31_treatment_confidence.bed woods_ascl1_neuroblastoma_wt51_treatment_confidence.bed wang_ascl1_gi-men_confidence.bed joung_ascl1_hesc_confidence.bed jain_tcf4_sudhl_confidence.bed jain_tcf4_tmd8_confidence.bed lorenzin_myc_u2os_treatment_confidence.bed soufi_myc_bj_treatment_confidence.bed lin_myc_p493-6_treatment_24h_confidence.bed walz_myc_u2os_treatment_confidence.bed barfeld_myc_lncap_treatment_confidence.bed sabo_myc_lymphoblastoma_highmyc_confidence.bed sabo_myc_lymphoblastoma_treatment_24h_confidence.bed muthalagu_myc_mcf10a_confidence.bed tameire_myc_colorectaladenocarcinoma_confidence.bed jaenicke_myc_imec_confidence.bed liang_myc_h2171_confidence.bed jung_myc_u2os_treatment_confidence.bed zeid_mycn_shep_treatment_6h_confidence.bed upton_mycn_lan5_confidence.bed upton_mycn_nb1643_confidence.bed upton_mycn_cogn415_confidence.bed bosse_mycn_ngp_confidence.bed bosse_mycn_kelly_confidence.bed bosse_mycn_nb1643_confidence.bed herold_mycn_shep_treatment_confidence.bed neikes_mycmax_mcf7_treatment_1000nm_confidence.bed louphrasitthiphol_mitf_melanoma_treatment_confidence.bed emming_bhlhe40_jurkat_treatment_clone1_confidence.bed emming_bhlhe40_jurkat_treatment_clone2_confidence.bed
do
    echo $file
    make_half_flank_peaks "$file"
done


### takee each half-flank peak file, and generate a clean version without chrM
cd ~/proneural/chip-seq/human
for file in kim_neurod1_hek_treatment_pcag_confidence.bed smith_neurog2_mrc5_treatment_1dpt_confidence.bed kim_twist1_hek_treatment_pcag_confidence.bed kim_twist1_hek_treatment_pcdna_confidence.bed manandhar_myod1_fibroblasts_confidence.bed park_ascl1_glioblastoma_treatment_confidence.bed woods_ascl1_neuroblastoma_wt31_treatment_confidence.bed woods_ascl1_neuroblastoma_wt51_treatment_confidence.bed wang_ascl1_gi-men_confidence.bed joung_ascl1_hesc_confidence.bed jain_tcf4_sudhl_confidence.bed jain_tcf4_tmd8_confidence.bed lorenzin_myc_u2os_treatment_confidence.bed soufi_myc_bj_treatment_confidence.bed lin_myc_p493-6_treatment_24h_confidence.bed walz_myc_u2os_treatment_confidence.bed barfeld_myc_lncap_treatment_confidence.bed sabo_myc_lymphoblastoma_highmyc_confidence.bed sabo_myc_lymphoblastoma_treatment_24h_confidence.bed muthalagu_myc_mcf10a_confidence.bed tameire_myc_colorectaladenocarcinoma_confidence.bed jaenicke_myc_imec_confidence.bed liang_myc_h2171_confidence.bed jung_myc_u2os_treatment_confidence.bed zeid_mycn_shep_treatment_6h_confidence.bed upton_mycn_lan5_confidence.bed upton_mycn_nb1643_confidence.bed upton_mycn_cogn415_confidence.bed bosse_mycn_ngp_confidence.bed bosse_mycn_kelly_confidence.bed bosse_mycn_nb1643_confidence.bed herold_mycn_shep_treatment_confidence.bed neikes_mycmax_mcf7_treatment_1000nm_confidence.bed louphrasitthiphol_mitf_melanoma_treatment_confidence.bed emming_bhlhe40_jurkat_treatment_clone1_confidence.bed emming_bhlhe40_jurkat_treatment_clone2_confidence.bed
do
    up_file="${file%.bed}_upstream_half.bed"
    down_file="${file%.bed}_downstream_half.bed"
    clean_up_file="${up_file%.bed}_nochrM.bed"
    clean_down_file="${down_file%.bed}_nochrM.bed"
    grep -v "chrM" "$up_file" > "$clean_up_file"
    grep -v "chrM" "$down_file" > "$clean_down_file"
done


cd ~/proneural/chip-seq/mouse
for file in lee_ascl1_mesc_confidence.bed matsuda_neurod1_microglia_confidence.bed pataskar_neurod1_mesc_confidence.bed fong_neurod2_p19_confidence.bed fong_neurod2_p19_control_confidence.bed fong_neurod2_p19_treatment_confidence.bed fong_neurod2_p19_treatment_confidence.bed fong_neurod2_mef_confidence.bed aydin_neurog2_eb_treatment_48h_confidence.bed velasco_neurog2_eb_treatment_12h_confidence.bed vainorius_neurog2_mesc_confidence.bed pereira_neurog2_astrocytes_confidence.bed pereira_neurog2_astrocytes_yy+_confidence.bed pereira_neurog2_astrocytes_yy-_confidence.bed pereira_neurog2_astrocytes_yy-_confidence.bed pereira_neurog2_astrocytes_yy+_confidence.bed li_twist2_myoblast_dm_treatment_confidence.bed li_twist2_myoblast_gm_treatment_confidence.bed fan_twist1_eb_treatment_confidence.bed mazzoni_olig2_mesc_confidence.bed lin_mesp1_mesc_treatment_24h_confidence.bed chalamalasetty_msgn1_eb_confidence.bed hershbach_myod1_mef_confidence.bed fong_myod1_p19_confidence.bed fong_myod1_mef_confidence.bed casey_myod1_mesc_confidence.bed lee_myod1_mef_confidence.bed lee_myod1_mesc_confidence.bed conerly_myod1_mef_confidence.bed conerly_myf5_mef_confidence.bed casey_ascl1_mesc_confidence.bed aydin_ascl1_eb_treatment_12h_confidence.bed lee_ascl1_mesc_confidence.bed lee_ascl1_mef_confidence.bed wapinski_ascl1_mef_confidence.bed wapinski_ascl1_npc_confidence.bed vainorius_ascl1_mesc_confidence.bed raposo_ascl1_ns5_confidence.bed casey_ascl2_mesc_confidence.bed liu_ascl2_cd4_treatment_confidence.bed kress_myc_liver_treatment_confidence.bed croci_myc_liver_treatment_confidence.bed croci_myc_mef_confidence.bed depretis_myc_mef_treatment_2h_confidence.bed weber_hey1_cardiomyocytes_treatment_confidence.bed weber_hey1_mesc_treatment_confidence.bed weber_hey2_cardiomyocytes_treatment_confidence.bed weber_hey2_mesc_treatment_confidence.bed
do
    up_file="${file%.bed}_upstream_half.bed"
    down_file="${file%.bed}_downstream_half.bed"
    clean_up_file="${up_file%.bed}_nochrM.bed"
    clean_down_file="${down_file%.bed}_nochrM.bed"
    grep -v "chrM" "$up_file" > "$clean_up_file"
    grep -v "chrM" "$down_file" > "$clean_down_file"
done






echo "
find_background_eboxes_human <- function(peaks_file, motif = \"CA..TG\") {
    refindall <- function(seq, pattern) {
        m <- gregexpr(pattern, seq, perl = TRUE)[[1]]
        if (m[1] == -1) integer(0) else as.integer(m)
    }

    fasta_out <- paste0(peaks_file, \".fa.tab\")
    cmd <- paste(\"bedtools getfasta -fi ~/genomes/hg38.fa -bed\", shQuote(peaks_file), \"-tab >\", shQuote(fasta_out))
    system(cmd)

    tab <- read.table(fasta_out, sep = \"\\t\", stringsAsFactors = FALSE, quote = \"\", comment.char = \"\")
    unlink(fasta_out)

    if (nrow(tab) == 0) {
        return(data.frame())
    }

    seqs <- toupper(tab[[2]])
    peaks <- read.table(peaks_file, header = FALSE, stringsAsFactors = FALSE)
    n <- nrow(peaks)

    res <- vector(\"list\", n)
    k <- 0L
    for (i in seq_len(n)) {
        s <- seqs[i]
        pos <- refindall(s, motif)
        if (length(pos) == 0) next

        chrom <- peaks[i, 1]
        start0 <- as.integer(peaks[i, 2])
        name <- if (ncol(peaks) >= 4) peaks[i, 4] else paste0(\"region_\", i)

        rows <- lapply(pos, function(p) {
            gstart <- start0 + p
            data.frame(
                chrom = chrom,
                start = gstart,
                end = gstart + 6L,
                region = name,
                ebox_start = gstart,
                ebox_seq = substr(s, p, p + 5L),
                stringsAsFactors = FALSE
            )
        })
        for (r in rows) {
            k <- k + 1L
            res[[k]] <- r
        }
    }
    if (k == 0) return(data.frame())
    do.call(rbind, res[seq_len(k)])
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop(\"Usage: Rscript find_background_eboxes_human.R <bed_file>\")
}
infile <- args[1]
cat(\"Scanning E-boxes (CANNTG) in:\", infile, \"\\n\")
eboxes <- find_background_eboxes_human(infile, motif = \"CA..TG\")
base <- sub(\"\\\\.[^.]*$\", \"\", basename(infile))
outfile <- paste0(base, \"_eboxes.RDS\")
saveRDS(eboxes, outfile)
cat(\"Saved:\", outfile, \"with\", nrow(eboxes), \"matches\\n\")
" > find_background_eboxes_human.R









echo "
find_background_eboxes_mouse <- function(peaks_file, motif = \"CA..TG\") {
    refindall <- function(seq, pattern) {
        m <- gregexpr(pattern, seq, perl = TRUE)[[1]]
        if (m[1] == -1) integer(0) else as.integer(m)
    }


    fasta_out <- paste0(peaks_file, \".fa.tab\")
    cmd <- paste(\"bedtools getfasta -fi ~/genomes/mm10.fa -bed\", shQuote(peaks_file), \"-tab >\", shQuote(fasta_out))
    system(cmd)

    tab <- read.table(fasta_out, sep = \"\\t\", stringsAsFactors = FALSE, quote = \"\", comment.char = \"\")
    unlink(fasta_out)

    if (nrow(tab) == 0) {
        return(data.frame())
    }

    seqs <- toupper(tab[[2]])
    peaks <- read.table(peaks_file, header = FALSE, stringsAsFactors = FALSE)
    
    n <- nrow(peaks)

    res <- vector(\"list\", n)
    k <- 0L
    for (i in seq_len(n)) {
        s <- seqs[i]
        pos <- refindall(s, motif)
        if (length(pos) == 0) next

        chrom <- peaks[i, 1]
        start0 <- as.integer(peaks[i, 2])
        name <- if (ncol(peaks) >= 4) peaks[i, 4] else paste0(\"region_\", i)

        rows <- lapply(pos, function(p) {
            gstart <- start0 + p
            data.frame(
                chrom = chrom,
                start = gstart,
                end = gstart + 6L,
                region = name,
                ebox_start = gstart,
                ebox_seq = substr(s, p, p + 5L),
                stringsAsFactors = FALSE
            )
        })
        for (r in rows) {
            k <- k + 1L
            res[[k]] <- r
        }
    }
    if (k == 0) return(data.frame())
    do.call(rbind, res[seq_len(k)])
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop(\"Usage: Rscript find_background_eboxes_mouse.R <bed_file>\")
}
infile <- args[1]
cat(\"Scanning E-boxes (CANNTG) in:\", infile, \"\\n\")
eboxes <- find_background_eboxes_mouse(infile, motif = \"CA..TG\")
base <- sub(\"\\\\.[^.]*$\", \"\", basename(infile))
outfile <- paste0(base, \"_eboxes.RDS\")
saveRDS(eboxes, outfile)
cat(\"Saved:\", outfile, \"with\", nrow(eboxes), \"matches\\n\")
" > find_background_eboxes_mouse.R







##### now find e-boxes in half-flank peaks
# Find E-boxes in half-flank peaks (human) using explicit filenames (no globs)
cd ~/proneural/chip-seq/human

human_files=(
    kim_neurod1_hek_treatment_pcag_confidence
    smith_neurog2_mrc5_treatment_1dpt_confidence
    kim_twist1_hek_treatment_pcag_confidence
    kim_twist1_hek_treatment_pcdna_confidence
    manandhar_myod1_fibroblasts_confidence
    park_ascl1_glioblastoma_treatment_confidence
    woods_ascl1_neuroblastoma_wt31_treatment_confidence
    woods_ascl1_neuroblastoma_wt51_treatment_confidence
    wang_ascl1_gi-men_confidence
    joung_ascl1_hesc_confidence
    jain_tcf4_sudhl_confidence
    jain_tcf4_tmd8_confidence
    lorenzin_myc_u2os_treatment_confidence
    soufi_myc_bj_treatment_confidence
    lin_myc_p493-6_treatment_24h_confidence
    walz_myc_u2os_treatment_confidence
    barfeld_myc_lncap_treatment_confidence
    sabo_myc_lymphoblastoma_highmyc_confidence
    sabo_myc_lymphoblastoma_treatment_24h_confidence
    muthalagu_myc_mcf10a_confidence
    tameire_myc_colorectaladenocarcinoma_confidence
    jaenicke_myc_imec_confidence
    liang_myc_h2171_confidence
    jung_myc_u2os_treatment_confidence
    zeid_mycn_shep_treatment_6h_confidence
    upton_mycn_lan5_confidence
    upton_mycn_nb1643_confidence
    upton_mycn_cogn415_confidence
    bosse_mycn_ngp_confidence
    bosse_mycn_kelly_confidence
    bosse_mycn_nb1643_confidence
    herold_mycn_shep_treatment_confidence
    neikes_mycmax_mcf7_treatment_1000nm_confidence
    louphrasitthiphol_mitf_melanoma_treatment_confidence
    emming_bhlhe40_jurkat_treatment_clone1_confidence
    emming_bhlhe40_jurkat_treatment_clone2_confidence
)


for base in "${human_files[@]}"; do
    upstream="${base}_upstream_half_nochrM.bed"
    downstream="${base}_downstream_half_nochrM.bed"

    up_rds="${upstream%.bed}_eboxes.RDS"
    down_rds="${downstream%.bed}_eboxes.RDS"

    if [ -f "$upstream" ] && [ ! -f "$up_rds" ]; then
        sbatch --partition=normal --output="${upstream%.bed}_eboxes.log" --wrap="Rscript ../find_background_eboxes_human.R $upstream"
    fi
    if [ -f "$downstream" ] && [ ! -f "$down_rds" ]; then
        sbatch --partition=normal --output="${downstream%.bed}_eboxes.log" --wrap="Rscript ../find_background_eboxes_human.R $downstream"
    fi
done


cat kim_neurod1_hek_treatment_pcag_confidence_upstream_half_nochrM_eboxes.log
du -h kim_neurod1_hek_treatment_pcag_confidence_upstream_half_nochrM_eboxes.RDS

du -h fong_neurod2_p19_confidence_downstream_half_eboxes_nochrM.RDS

a

# Cleanup temporary chunk files (human)
rm -f rep_*bed rep_*fa



# Find E-boxes in half-flank peaks (mouse) using explicit filenames
cd ~/proneural/chip-seq/mouse

mouse_files=(
    lee_ascl1_mesc_confidence
    matsuda_neurod1_microglia_confidence
    pataskar_neurod1_mesc_confidence
    fong_neurod2_p19_confidence
    fong_neurod2_p19_control_confidence
    fong_neurod2_p19_treatment_confidence
    fong_neurod2_mef_confidence
    aydin_neurog2_eb_treatment_48h_confidence
    velasco_neurog2_eb_treatment_12h_confidence
    vainorius_neurog2_mesc_confidence
    pereira_neurog2_astrocytes_confidence
    pereira_neurog2_astrocytes_yy+_confidence
    pereira_neurog2_astrocytes_yy-_confidence
    li_twist2_myoblast_dm_treatment_confidence
    li_twist2_myoblast_gm_treatment_confidence
    fan_twist1_eb_treatment_confidence
    mazzoni_olig2_mesc_confidence
    lin_mesp1_mesc_treatment_24h_confidence
    chalamalasetty_msgn1_eb_confidence
    hershbach_myod1_mef_confidence
    fong_myod1_p19_confidence
    fong_myod1_mef_confidence
    casey_myod1_mesc_confidence
    lee_myod1_mef_confidence
    lee_myod1_mesc_confidence
    conerly_myod1_mef_confidence
    conerly_myf5_mef_confidence
    casey_ascl1_mesc_confidence
    aydin_ascl1_eb_treatment_12h_confidence
    lee_ascl1_mef_confidence
    wapinski_ascl1_mef_confidence
    wapinski_ascl1_npc_confidence
    vainorius_ascl1_mesc_confidence
    raposo_ascl1_ns5_confidence
    casey_ascl2_mesc_confidence
    liu_ascl2_cd4_treatment_confidence
    kress_myc_liver_treatment_confidence
    croci_myc_liver_treatment_confidence
    croci_myc_mef_confidence
    depretis_myc_mef_treatment_2h_confidence
    weber_hey1_cardiomyocytes_treatment_confidence
    weber_hey1_mesc_treatment_confidence
    weber_hey2_cardiomyocytes_treatment_confidence
    weber_hey2_mesc_treatment_confidence
)



for base in "${mouse_files[@]}"; do
    upstream="${base}_upstream_half_nochrM.bed"
    downstream="${base}_downstream_half_nochrM.bed"

    up_rds="${upstream%.bed}_eboxes.RDS"
    down_rds="${downstream%.bed}_eboxes.RDS"

    if [ -f "$upstream" ] && [ ! -f "$up_rds" ]; then
        sbatch --partition=normal --output="${upstream%.bed}_eboxes.log" --wrap="Rscript ../find_background_eboxes_mouse.R $upstream"
    fi
    if [ -f "$downstream" ] && [ ! -f "$down_rds" ]; then
        sbatch --partition=normal --output="${downstream%.bed}_eboxes.log" --wrap="Rscript ../find_background_eboxes_mouse.R $downstream"
    fi
done



cat fong_neurod2_p19_confidence_downstream_half_nochrM_eboxes.log

cat fong_neurod2_mef_confidence_downstream_half_nochrM_eboxes.log

cat depretis_myc_mef_treatment_2h_confidence_downstream_half_nochrM_eboxes.log

du -h fong_neurod2_p19_confidence_downstream_half_nochrM_eboxes.RDS

# Cleanup temporary chunk files (mouse)
rm -f rep_*bed rep_*fa




