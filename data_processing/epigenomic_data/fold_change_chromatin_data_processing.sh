############## CHROMATIN PEAK CALLING #########
###############################################


###### ATAC-seq

# Mouse
for file in aydin_0h_atac_1.sam.processed.bam aydin_12h_ascl1_atac.sam.processed.bam aydin_neurog2_12h_atac.sam.processed.bam \
casey_esc_uninduced_ascl1_atac.fastq.gz.sam.processed.bam casey_24h_ascl1_atac.fastq.gz.sam.processed.bam \
casey_esc_uninduced_ascl2_atac.fastq.gz.sam.processed.bam casey_24h_ascl2_atac.fastq.gz.sam.processed.bam \
casey_esc_uninduced_myod1_atac.fastq.gz.sam.processed.bam casey_24h_myod1_atac.fastq.gz.sam.processed.bam \
lee_mef_atac.sam.processed.bam lee_ascl1_48h_atac.sam.processed.bam lee_myod1_48h_atac.sam.processed.bam \
wapinski_MEF.sam.processed.bam wapinski_MEF+Ascl1_48hr_B1.sam.processed.bam
do
  sbatch --partition=bigmem --wrap="macs2 callpeak -f BAM -t $file -g mm --outdir . -n $file"
done



# faltan
cd ~/proneural/mouse/pairedd
for file in pereira_astrocytes_GSM6368722_ATACseq_GFP.sam.rmdup.bam pereira_astrocytes_GSM6368725_ATACseq_Ngn2.sam.rmdup.bam
do
  sbatch --partition=bigmem --mem=400000 --wrap="macs2 callpeak -f BAM -t $file -g mm --outdir . -n $file"
done


cd ~/proneural/mouse/single
for file in liu_GSM3659988_GFP-GM_K27ac_1.fastq.gz.sam.processed.bam li_Tw2-GM_K27ac_1.fastq.gz.sam.processed.bam
do
  sbatch --partition=bigmem --mem=400000 --wrap="macs2 callpeak --broad -f BAM -t $file -g mm --outdir . -n $file"
done

cd ~/proneural/mouse/single

## no me sale el peak calling del raposo, lo he probado tanto con mis bams como con bams descargados de ellos

for file in raposo_ns5control_dnase.bam raposo_ns5induced_dnase.bam
do
  sbatch --partition=bigmem --mem=600000 --wrap="macs2 callpeak -f BAM -t $file -g mm --outdir . -n $file"
done

cd ~/proneural/mouse/paired
for file in conerly_GSM1954034_H4Ac_MyoDs*bam conerly_GSM1954035_H4Ac_Myf5*bam mef_conerly_h4ac.sam.processed.bam
do
  sbatch --partition=bigmem --mem=400000 --wrap="macs2 callpeak --broad -f BAM -t $file -g mm --outdir . -n $file"
done

# Human
cd ~/proneural
for file in park_atac_control.fastq.gz.sam.processed.bam park_atac_dox.fastq.sam.processed.bam
do
  sbatch --partition=bigmem --mem=400000 --wrap="macs2 callpeak -f BAM -t $file -g hs --outdir . -n $file"
done

cd ~/proneural/human/paired
for file in wang_ATAC_GIMEN_ASCL1_Neg.sam.processed.bam wang_ATAC_GIMEN_ASCL1_Pos.sam.processed.bam
do
  sbatch --partition=bigmem --mem=400000 --wrap="macs2 callpeak -f BAM -t $file -g hs --outdir . -n $file"
done

##### PTM

# Mouse
for file in li_GFP-DM_K27ac_1.fastq.gz.sam.processed.bam li_Tw2-DM_K27ac_1.fastq.gz.sam.processed.bam \
fong_p19_control_ach4.fastq.gz.sam.processed.bam fong_p19_neurod2_induction_ach4.fastq.gz.sam.processed.bam \
do
  sbatch --partition=bigmem --wrap="macs2 callpeak --broad -f BAM -t $file -g mm --outdir . -n $file"
done

  
# Human
for file in smith_IP1_H3K27AC_MRC5_1DPT.fastq.gz.sam.processed.bam smith_IP1_H3K27AC_N_1DPT.fastq.gz.sam.processed.bam \
lin_P493-6_T0_H3K27AC.fastq.gz.sam.processed.bam lin_P493-6_T24_H3K27AC.fastq.gz.sam.processed.bam \
zeid_SHEP_0HR_H3K27AC.fastq.gz.sam.processed.bam zeid_SHEP_6HR_H3K27AC.fastq.gz.sam.processed.bam \
walz_U2OS_AcH4_+Dox_ChIPseq.fastq.gz.sam.processed.bam walz_U2OS_AcH4_-Dox_ChIPseq.fastq.gz.sam.processed.bam
do
  sbatch --partition=bigmem --wrap="macs2 callpeak --broad -f BAM -t $file -g hs --outdir . -n $file"
done

cd ~/proneural/human/single
for file in barfeld_H3K27ac_R1881.fastq.gz.sam.processed.bam barfeld_H3K27ac_R1881_Dox.fastq.gz.sam.processed.bam
do
  sbatch --partition=bigmem --mem=400000 --wrap="macs2 callpeak --broad -f BAM -t $file -g hs --outdir . -n $file"
done

for file in *manandhar_GSM2449222_Fb_dnase*bam *manandhar_GSM2449225_Fd_dnase*bam
do
  sbatch --partition=bigmem --mem=400000 --wrap="macs2 callpeak -f BAM -t $file -g hs --outdir . -n $file"
done

cd ~/proneural
mv */*/*Peak .



############ FOLD CHANGE COMPUTE #################
##################################################


### first concatenate pre-induction and post-induction peaks


module load BEDTools
concatenate_beds () {
  bedtools intersect -v -a $2 -b $1 > ~/proneural/missing.bed
  cat $1 ~/proneural/missing.bed | sort -k1,1 -k2,2n | grep -v JH | grep -v GL | grep -v random | grep -v chrM | grep -v KI | grep -v alt | uniq  > ~/proneural/$3.bed
}

concatenate_beds aydin_0h_atac_1.sam.processed.bam_peaks.narrowPeak aydin_12h_ascl1_atac.sam.processed.bam_peaks.narrowPeak "aydin_ascl1_12h_join"
concatenate_beds aydin_0h_atac_1.sam.processed.bam_peaks.narrowPeak aydin_48h_neurog2_atac.sam.processed.bam_peaks.narrowPeak "aydin_neurog2_48h_join"
concatenate_beds casey_esc_uninduced_ascl1_atac.fastq.gz.sam.processed.bam_peaks.narrowPeak casey_24h_ascl1_atac.fastq.gz.sam.processed.bam_peaks.narrowPeak "casey_ascl1_24h_join"
concatenate_beds casey_esc_uninduced_ascl2_atac.fastq.gz.sam.processed.bam_peaks.narrowPeak casey_24h_ascl2_atac.fastq.gz.sam.processed.bam_peaks.narrowPeak "casey_ascl2_24h_join"
concatenate_beds casey_esc_uninduced_myod1_atac.fastq.gz.sam.processed.bam_peaks.narrowPeak casey_24h_myod1_atac.fastq.gz.sam.processed.bam_peaks.narrowPeak "casey_myod1_24h_join"
concatenate_beds lee_mef_atac.sam.processed.bam_peaks.narrowPeak lee_myod1_48h_atac.sam.processed.bam_peaks.narrowPeak "lee_myod1_48h_join"
concatenate_beds lee_mef_atac.sam.processed.bam_peaks.narrowPeak lee_ascl1_48h_atac.sam.processed.bam_peaks.narrowPeak "lee_ascl1_48h_join"
concatenate_beds wapinski_MEF.sam.processed.bam_peaks.narrowPeak wapinski_MEF+Ascl1_48hr_B1.sam.processed.bam_peaks.narrowPeak "wapinski_ascl1_48hr_join"
concatenate_beds park_atac_control.fastq.gz.sam.processed.bam_peaks.narrowPeak park_atac_dox.fastq.sam.processed.bam_peaks.narrowPeak "park_ascl1_dox_join"
concatenate_beds pereira_astrocytes_GSM6368722_ATACseq_GFP.sam.rmdup.bam_peaks.narrowPeak pereira_astrocytes_GSM6368725_ATACseq_Ngn2.sam.rmdup.bam_peaks.narrowPeak "pereira_astrocytes_ngn2_join"
concatenate_beds wang_ATAC_GIMEN_ASCL1_Neg.sam.processed.bam_peaks.narrowPeak wang_ATAC_GIMEN_ASCL1_Pos.sam.processed.bam_peaks.narrowPeak "wang_ascl1_join"
concatenate_beds manandhar_GSM2449222_Fb_dnase.fastq.gz.sam.rmdup.sam.processed.bam_peaks.narrowPeak manandhar_GSM2449225_Fd_dnase.fastq.gz.sam.processed.bam_peaks.narrowPeak "manandhar_myod1_join"
concatenate_beds li_GFP-DM_K27ac_1.fastq.gz.sam.processed.bam_peaks.broadPeak li_Tw2-DM_K27ac_1.fastq.gz.sam.processed.bam_peaks.broadPeak "li_tw2_induction_join"
concatenate_beds fong_p19_control_ach4.fastq.gz.sam.processed.bam_peaks.broadPeak fong_p19_neurod2_induction_ach4.fastq.gz.sam.processed.bam_peaks.broadPeak "fong_neurod2_induction_join"
concatenate_beds smith_IP1_H3K27AC_MRC5_1DPT.fastq.gz.sam.processed.bam_peaks.broadPeak smith_IP1_H3K27AC_N_1DPT.fastq.gz.sam.processed.bam_peaks.broadPeak "smith_neurog2_induction_join"
concatenate_beds lin_P493-6_T0_H3K27AC.fastq.gz.sam.processed.bam_peaks.broadPeak lin_P493-6_T24_H3K27AC.fastq.gz.sam.processed.bam_peaks.broadPeak "lin_myc_t24_join"
concatenate_beds zeid_SHEP_0HR_H3K27AC.fastq.gz.sam.processed.bam_peaks.broadPeak zeid_SHEP_6HR_H3K27AC.fastq.gz.sam.processed.bam_peaks.broadPeak "zeid_mycn_6hr_join"
concatenate_beds walz_U2OS_AcH4_-Dox_ChIPseq.fastq.gz.sam.processed.bam_peaks.broadPeak walz_U2OS_AcH4_+Dox_ChIPseq.fastq.gz.sam.processed.bam_peaks.broadPeak "walz_myc_dox_join"
concatenate_beds liu_GSM3659988_GFP-GM_K27ac_1.fastq.gz.sam.processed.bam_peaks.broadPeak li_Tw2-GM_K27ac_1.fastq.gz.sam.processed.bam_peaks.broadPeak "li_tw2_gm_induction_join"
concatenate_beds mef_conerly_h4ac.sam.processed.bam_peaks.broadPeak conerly_GSM1954035_H4Ac_Myf5.sam.processed.bam_peaks.broadPeak "conerly_myf5_join"
concatenate_beds mef_conerly_h4ac.sam.processed.bam_peaks.broadPeak conerly_GSM1954034_H4Ac_MyoDs.sam.processed.bam_peaks.broadPeak "conerly_myod_join"
concatenate_beds barfeld_H3K27ac_R1881.fastq.gz.sam.processed.bam_peaks.broadPeak barfeld_H3K27ac_R1881_Dox.fastq.gz.sam.processed.bam_peaks.broadPeak "barfeld_myc_join"


concatenate_beds raposo_ns5control_dnase*bam_peaks.narrowPeak raposo_ns5induced_dnase*bam_peaks.narrowPeak "raposo_ascl1induced_join"




srun --partition=normal --mpi=none --nodes=1 --pty bash -i

cd ~/proneural


## resize the peaks to 1000bp centered at the summit. summit position within the peak is in the 10th column for narrowPeak.
## so summit position is col2 + col10
## and to resize, the new col2 is summit - 500 and new col3 is summit + 500
resize_narrowpeak () {
  awk -v OFS="\t" \'{summit=$2+$10; start=summit-500; if (start<0) start=0; end=summit+500; print $1,start,end,$4,$5,$6,$6,$7,$8,$9,$10}\' $1 > ~/proneural/${1%.bed}_1000bp.bed
}





bedmap_chromatin ( ) {
  bedtools map -o max -c 4 -a $1 -b /projects_ng/xabier/proneural/$2 > ~/proneural/bedmap_pre_${1}
  bedtools map -o max -c 4 -a $1 -b /projects_ng/xabier/proneural/$3 > ~/proneural/bedmap_post_${1}
}



module load BEDTools



bedmap_chromatin aydin_ascl1_12h_join.bed clean.aydin_0h_atac_2.sam.processed.bam.bdg clean.aydin_12h_ascl1_atac.sam.processed.bam.bdg
bedmap_chromatin aydin_neurog2_48h_join.bed clean.aydin_0h_atac_2.sam.processed.bam.bdg clean.aydin_48h_neurog2_atac.sam.processed.bam.bdg
bedmap_chromatin casey_ascl1_24h_join.bed clean.casey_esc_uninduced_ascl1_atac.fastq.gz.sam.processed.bam.bdg clean.casey_24h_ascl1_atac.fastq.gz.sam.processed.bam.bdg
bedmap_chromatin casey_ascl2_24h_join.bed clean.casey_esc_uninduced_ascl2_atac.fastq.gz.sam.processed.bam.bdg clean.casey_24h_ascl2_atac.fastq.gz.sam.processed.bam.bdg
bedmap_chromatin casey_myod1_24h_join.bed clean.casey_esc_uninduced_myod1_atac.fastq.gz.sam.processed.bam.bdg clean.casey_24h_myod1_atac.fastq.gz.sam.processed.bam.bdg
bedmap_chromatin lee_myod1_48h_join.bed clean.lee_mef_atac.sam.processed.bam.bdg clean.lee_myod1_48h_atac.sam.processed.bam.bdg
bedmap_chromatin lee_ascl1_48h_join.bed clean.lee_mef_atac.sam.processed.bam.bdg clean.lee_ascl1_48h_atac.sam.processed.bam.bdg
bedmap_chromatin wapinski_ascl1_48hr_join.bed clean.wapinski_MEF.sam.processed.bam.bdg clean.wapinski_MEF+Ascl1_48hr_B2.sam.processed.bam.bdg
bedmap_chromatin park_ascl1_dox_join.bed clean.park_atac_control.fastq.gz.sam.processed.bam.bdg clean.park_atac_dox.fastq.sam.processed.bam.bdg
bedmap_chromatin li_tw2_induction_join.bed clean.li_GFP-DM_K27ac_1.fastq.gz.sam.processed.bam.bdg clean.li_Tw2-DM_K27ac_1.fastq.gz.sam.processed.bam.bdg
bedmap_chromatin smith_neurog2_induction_join.bed clean.smith_IP1_H3K27AC_MRC5_1DPT.fastq.gz.sam.processed.bam.bdg clean.smith_IP1_H3K27AC_N_1DPT.fastq.gz.sam.processed.bam.bdg
bedmap_chromatin lin_myc_t24_join.bed clean.lin_P493-6_T0_H3K27AC.fastq.gz.sam.processed.bam.bdg clean.lin_P493-6_T24_H3K27AC.fastq.gz.sam.processed.bam.bdg
bedmap_chromatin zeid_mycn_6hr_join.bed clean.zeid_SHEP_0HR_H3K27AC.fastq.gz.sam.processed.bam.bdg clean.zeid_SHEP_6HR_H3K27AC.fastq.gz.sam.processed.bam.bdg
bedmap_chromatin fong_neurod2_induction_join.bed clean.fong_p19_control_ach4.fastq.gz.sam.processed.bam.bdg clean.fong_p19_neurod2_induction_ach4.fastq.gz.sam.processed.bam.bdg
bedmap_chromatin walz_myc_dox_join.bed clean.walz_U2OS_AcH4_-Dox_ChIPseq.fastq.gz.sam.processed.bam.bdg clean.walz_U2OS_AcH4_+Dox_ChIPseq.fastq.gz.sam.processed.bam.bdg
bedmap_chromatin barfeld_myc_join.bed clean.barfeld_H3K27ac_R1881.fastq.gz.sam.processed.bam.bdg clean.barfeld_H3K27ac_R1881_Dox.fastq.gz.sam.processed.bam.bdg
bedmap_chromatin pereira_astrocytes_ngn2_join.bed clean.pereira_astrocytes_GSM6368722_ATACseq_GFP.sam.rmdup.bam.bdg clean.pereira_astrocytes_GSM6368725_ATACseq_Ngn2.sam.rmdup.bam.bdg
bedmap_chromatin wang_ascl1_join.bed clean.wang_ATAC_GIMEN_ASCL1_Neg.sam.processed.bam.bdg clean.wang_ATAC_GIMEN_ASCL1_Pos.sam.processed.bam.bdg
bedmap_chromatin manandhar_myod1_join.bed clean.manandhar_GSM2449222_Fb_dnase*bam.bdg clean.*manandhar_GSM2449225_Fd_dnase*bam.bdg
bedmap_chromatin li_tw2_gm_induction_join.bed clean.liu_GSM3659988_GFP-GM_K27ac_1*bam.bdg clean.li_Tw2-GM_K27ac_1.fastq.gz.sam.processed.bam.bdg
bedmap_chromatin conerly_myf5_join.bed clean.mef_conerly_h4ac.sam.processed.bam.bdg clean.conerly_GSM1954035_H4Ac_Myf5*bam.bdg
bedmap_chromatin conerly_myod_join.bed clean.mef_conerly_h4ac.sam.processed.bam.bdg clean.conerly_GSM1954034_H4Ac_MyoDs*bam.bdg



raposo_ascl1induced_join.bed clean.raposo_ns5induced_dnase.bam.bdg clean.raposo_ns5control_dnase.bam.bdg



# Compute fold change
## Fold change between the maximum coverage pre induction within the peaks and the maximum coverage post induction



### broadPeak and narrowPeak files have different numbers of columns

for file in aydin_ascl1_12h_join.bed aydin_neurog2_48h_join.bed casey_ascl1_24h_join.bed casey_ascl2_24h_join.bed casey_myod1_24h_join.bed lee_ascl1_48h_join.bed lee_myod1_48h_join.bed wapinski_ascl1_48hr_join.bed park_ascl1_dox_join.bed pereira_astrocytes_ngn2_join.bed wang_ascl1_join.bed manandhar_myod1_join.bed
do
  echo $file
  paste bedmap_pre*$file bedmap_post*$file | awk '{if ($11 == 0) print $1,$2,$3,"inf"; else print $1,$2,$3,$22/$11}' | tr " " "\t" > foldchange_induction_${file}
done

for file in li_tw2_induction_join.bed smith_neurog2_induction_join.bed lin_myc_t24_join.bed zeid_mycn_6hr_join.bed fong_neurod2_induction_join.bed walz_myc_dox_join.bed barfeld_myc_join.bed li_tw2_gm_induction_join.bed conerly_myf5_join.bed conerly_myod_join.bed
do
  echo $file
  paste bedmap_pre*$file bedmap_post*$file | awk '{if ($10 == 0) print $1,$2,$3,"inf"; else print $1,$2,$3,$20/$10}' | tr " " "\t" > foldchange_induction_${file}
done

for file in wapinski_ascl1_48hr_join.bed
do
  echo $file
  paste bedmap_pre*$file bedmap_post*$file | awk '{if ($11 == 0) print $1,$2,$3,"inf"; else print $1,$2,$3,$22/$11}' | tr " " "\t" > foldchange_induction_${file}
done



########## Intersect with the chip seq peaks
############################################

cd ~/proneural
module load BEDTools
bedtools intersect -wao -a foldchange_induction_aydin_ascl1_12h_join.bed -b chip-seq/mouse/aydin_ascl1_eb_treatment_12h_confidence.bed > intersect_foldchange_induction_aydin_ascl1_12h_join.bed
bedtools intersect -wao -a foldchange_induction_aydin_neurog2_48h_join.bed -b chip-seq/mouse/aydin_neurog2_eb_treatment_48h_confidence.bed > intersect_foldchange_induction_aydin_neurog2_48h_join.bed
bedtools intersect -wao -a foldchange_induction_barfeld_myc_join.bed -b chip-seq/human/barfeld_myc_lncap_treatment_confidence.bed > intersect_foldchange_induction_barfeld_myc_join.bed
bedtools intersect -wao -a foldchange_induction_casey_ascl1_24h_join.bed -b chip-seq/mouse/casey_ascl1_mesc_confidence.bed > intersect_foldchange_induction_casey_ascl1_24h_join.bed
bedtools intersect -wao -a foldchange_induction_casey_ascl2_24h_join.bed -b chip-seq/mouse/casey_ascl2_mesc_confidence.bed > intersect_foldchange_induction_casey_ascl2_24h_join.bed
bedtools intersect -wao -a foldchange_induction_casey_myod1_24h_join.bed -b chip-seq/mouse/casey_myod1_mesc_confidence.bed > intersect_foldchange_induction_casey_myod1_24h_join.bed
bedtools intersect -wao -a foldchange_induction_conerly_myf5_join.bed -b chip-seq/mouse/conerly_myf5_mef_confidence.bed > intersect_foldchange_induction_conerly_myf5_join.bed
bedtools intersect -wao -a foldchange_induction_conerly_myod_join.bed -b chip-seq/mouse/conerly_myod1_mef_confidence.bed > intersect_foldchange_induction_conerly_myod_join.bed
bedtools intersect -wao -a foldchange_induction_fong_neurod2_induction_join.bed -b chip-seq/mouse/fong_neurod2_p19_confidence.bed > intersect_foldchange_induction_fong_neurod2_induction_join.bed
bedtools intersect -wao -a foldchange_induction_lee_myod1_48h_join.bed -b chip-seq/mouse/lee_myod1_mef_confidence.bed > intersect_foldchange_induction_lee_myod1_48h_join.bed
bedtools intersect -wao -a foldchange_induction_lee_ascl1_48h_join.bed -b chip-seq/mouse/lee_ascl1_mef_confidence.bed > intersect_foldchange_induction_lee_ascl1_48h_join.bed
bedtools intersect -wao -a foldchange_induction_li_tw2_gm_induction_join.bed -b chip-seq/mouse/li_twist2_myoblast_gm_treatment_confidence.bed > intersect_foldchange_induction_li_tw2_gm_induction_join.bed
bedtools intersect -wao -a foldchange_induction_li_tw2_induction_join.bed -b chip-seq/mouse/li_twist2_myoblast_dm_treatment_confidence.bed > intersect_foldchange_induction_li_tw2_induction_join.bed
bedtools intersect -wao -a foldchange_induction_lin_myc_t24_join.bed -b chip-seq/human/lin_myc_p493-6_treatment_24h_confidence.bed > intersect_foldchange_induction_lin_myc_t24_join.bed
bedtools intersect -wao -a foldchange_induction_manandhar_myod1_join.bed -b chip-seq/human/manandhar_myod1_fibroblasts_confidence.bed > intersect_foldchange_induction_manandhar_myod1_join.bed
bedtools intersect -wao -a foldchange_induction_park_ascl1_dox_join.bed -b chip-seq/human/park_ascl1_glioblastoma_treatment_confidence.bed > intersect_foldchange_induction_park_ascl1_dox_join.bed
bedtools intersect -wao -a foldchange_induction_pereira_astrocytes_ngn2_join.bed -b chip-seq/mouse/pereira_neurog2_astrocytes_confidence.bed > intersect_foldchange_induction_pereira_astrocytes_ngn2_join.bed
bedtools intersect -wao -a foldchange_induction_smith_neurog2_induction_join.bed -b chip-seq/human/smith_neurog2_mrc5_treatment_1dpt_confidence.bed  > intersect_foldchange_induction_smith_neurog2_induction_join.bed
bedtools intersect -wao -a foldchange_induction_walz_myc_dox_join.bed -b chip-seq/human/walz_myc_u2os_treatment_confidence.bed > intersect_foldchange_induction_walz_myc_dox_join.bed
bedtools intersect -wao -a foldchange_induction_wang_ascl1_join.bed -b chip-seq/human/wang_ascl1_gi-men_confidence.bed > intersect_foldchange_induction_wang_ascl1_join.bed
bedtools intersect -wao -a foldchange_induction_wapinski_ascl1_48hr_join.bed -b chip-seq/mouse/wapinski_ascl1_mef_confidence.bed > intersect_foldchange_induction_wapinski_ascl1_48hr_join.bed
bedtools intersect -wao -a foldchange_induction_zeid_mycn_6hr_join.bed -b chip-seq/human/zeid_mycn_shep_treatment_6h_confidence.bed > intersect_foldchange_induction_zeid_mycn_6hr_join.bed


