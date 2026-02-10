
#####################################################################
########### co-dependence foldchange / initial accessibility ########
#####################################################################


### use the 1000bp regions centered on the summits of the chip-seq peaks to compute the accessibility foldchange
setwd("~/proneural")






  srun --partition=normal --mpi=none --nodes=1 --pty bash -i

  cd ~/proneural/chip-seq

  bedmap_chromatin ( ) {
    grep -v chrM $1 | bedtools map -o max -c 4 -a - -b /projects_ng/xabier/proneural/$2 > ~/proneural/bedmap_pre_${1}
    grep -v chrM $1 | bedtools map -o max -c 4 -a - -b /projects_ng/xabier/proneural/$3 > ~/proneural/bedmap_post_${1}
  }

  ## los que las bdg ya estan en el projects  
  module load BEDTools

  bedmap_chromatin aydin_ascl1_eb_treatment_12h_confidence_1000bp_around_summit.bed clean.aydin_0h_atac_2.sam.processed.bam.bdg clean.aydin_12h_ascl1_atac_rep2.sam.processed.bam.bdg
  bedmap_chromatin aydin_neurog2_eb_treatment_48h_confidence_1000bp_around_summit.bed clean.aydin_0h_atac_2.sam.processed.bam.bdg clean.aydin_48h_neurog2_atac_rep2.sam.processed.bam.bdg
  bedmap_chromatin casey_ascl1_mesc_confidence_1000bp_around_summit.bed clean.casey_esc_uninduced_ascl1_atac.fastq.gz.sam.processed.bam.bdg clean.casey_24h_ascl1_atac.fastq.gz.sam.processed.bam.bdg
  bedmap_chromatin casey_ascl2_mesc_confidence_1000bp_around_summit.bed clean.casey_esc_uninduced_ascl2_atac.fastq.gz.sam.processed.bam.bdg clean.casey_24h_ascl2_atac.fastq.gz.sam.processed.bam.bdg
  bedmap_chromatin casey_myod1_mesc_confidence_1000bp_around_summit.bed clean.casey_esc_uninduced_myod1_atac.fastq.gz.sam.processed.bam.bdg clean.casey_24h_myod1_atac.fastq.gz.sam.processed.bam.bdg
  bedmap_chromatin lee_myod1_mef_confidence_1000bp_around_summit.bed clean.lee_mef_atac.sam.processed.bam.bdg clean.lee_myod1_48h_atac_rep1.sam.processed.bam.bdg
  bedmap_chromatin wapinski_ascl1_mef_confidence_1000bp_around_summit.bed clean.wapinski_MEF.sam.processed.bam.bdg clean.wapinski_MEF+Ascl1_48hr_B2_rep1.sam.processed.bam.bdg
  bedmap_chromatin park_ascl1_glioblastoma_treatment_confidence_1000bp_around_summit.bed clean.park_atac_control_rep1.fastq.gz.sam.processed.bam.bdg clean.park_atac_dox_rep1.fastq.sam.processed.bam.bdg
  bedmap_chromatin li_twist2_myoblast_dm_treatment_confidence_1000bp_around_summit.bed clean.li_GFP-DM_K27ac_1.fastq.gz.sam.processed.bam.bdg clean.li_Tw2-DM_K27ac_1.fastq.gz.sam.processed.bam.bdg
  bedmap_chromatin smith_neurog2_mrc5_treatment_1dpt_confidence_1000bp_around_summit.bed clean.smith_IP1_H3K27AC_MRC5_1DPT.fastq.gz.sam.processed.bam.bdg clean.smith_IP1_H3K27AC_N_1DPT.fastq.gz.sam.processed.bam.bdg
  bedmap_chromatin lin_myc_p493-6_treatment_24h_confidence_1000bp_around_summit.bed clean.lin_P493-6_T0_H3K27AC.fastq.gz.sam.processed.bam.bdg clean.lin_P493-6_T24_H3K27AC.fastq.gz.sam.processed.bam.bdg
  bedmap_chromatin zeid_mycn_shep_treatment_6h_confidence_1000bp_around_summit.bed clean.zeid_SHEP_0HR_H3K27AC.fastq.gz.sam.processed.bam.bdg clean.zeid_SHEP_6HR_H3K27AC.fastq.gz.sam.processed.bam.bdg
  bedmap_chromatin fong_neurod2_p19_confidence_1000bp_around_summit.bed clean.fong_p19_control_ach4.fastq.gz.sam.processed.bam.bdg clean.fong_p19_neurod2_induction_ach4.fastq.gz.sam.processed.bam.bdg
  bedmap_chromatin walz_myc_u2os_treatment_confidence_1000bp_around_summit.bed clean.walz_U2OS_AcH4_-Dox_ChIPseq.fastq.gz.sam.processed.bam.bdg clean.walz_U2OS_AcH4_+Dox_ChIPseq.fastq.gz.sam.processed.bam.bdg
  bedmap_chromatin barfeld_myc_lncap_treatment_confidence_1000bp_around_summit.bed clean.barfeld_H3K27ac_R1881.fastq.gz.sam.processed.bam.bdg clean.barfeld_H3K27ac_R1881_Dox.fastq.gz.sam.processed.bam.bdg
  bedmap_chromatin pereira_neurog2_astrocytes_confidence_1000bp_around_summit.bed clean.pereira_astrocytes_GSM6368722_ATACseq_GFP_rep1.sam.rmdup.bam.bdg clean.pereira_astrocytes_GSM6368725_ATACseq_Ngn2_rep1.sam.rmdup.bam.bdg
  bedmap_chromatin wang_ascl1_gi-men_confidence_1000bp_around_summit.bed clean.wang_ATAC_GIMEN_ASCL1_Neg_Rep1.sam.processed.bam.bdg clean.wang_ATAC_GIMEN_ASCL1_Pos_Rep1.sam.processed.bam.bdg
  bedmap_chromatin manandhar_myod1_fibroblasts_confidence_1000bp_around_summit.bed clean.manandhar_GSM2449222_Fb_rep1_dnase*bam.bdg clean.*manandhar_GSM2449225_Fd_rep1_dnase*bam.bdg
  bedmap_chromatin li_twist2_myoblast_gm_treatment_confidence_1000bp_around_summit.bed clean.liu_GSM3659988_GFP-GM_K27ac_1*bam.bdg clean.li_Tw2-GM_K27ac_1.fastq.gz.sam.processed.bam.bdg
  bedmap_chromatin conerly_myf5_mef_confidence_1000bp_around_summit.bed clean.mef_conerly_h4ac.sam.processed.bam.bdg clean.conerly_GSM1954035_H4Ac_Myf5*bam.bdg
  bedmap_chromatin conerly_myod1_mef_confidence_1000bp_around_summit.bed clean.mef_conerly_h4ac.sam.processed.bam.bdg clean.conerly_GSM1954034_H4Ac_MyoDs*bam.bdg
  
  


# Compute fold change
## Fold change between the maximum coverage pre induction within the peaks and the maximum coverage post induction



  cd ~/proneural
  
  for file in aydin_ascl1_eb_treatment_12h_confidence_1000bp_around_summit.bed aydin_neurog2_eb_treatment_48h_confidence_1000bp_around_summit.bed casey_ascl1_mesc_confidence_1000bp_around_summit.bed casey_ascl2_mesc_confidence_1000bp_around_summit.bed casey_myod1_mesc_confidence_1000bp_around_summit.bed lee_myod1_mef_confidence_1000bp_around_summit.bed wapinski_ascl1_mef_confidence_1000bp_around_summit.bed park_ascl1_glioblastoma_treatment_confidence_1000bp_around_summit.bed li_twist2_myoblast_dm_treatment_confidence_1000bp_around_summit.bed smith_neurog2_mrc5_treatment_1dpt_confidence_1000bp_around_summit.bed lin_myc_p493-6_treatment_24h_confidence_1000bp_around_summit.bed zeid_mycn_shep_treatment_6h_confidence_1000bp_around_summit.bed fong_neurod2_p19_confidence_1000bp_around_summit.bed walz_myc_u2os_treatment_confidence_1000bp_around_summit.bed barfeld_myc_lncap_treatment_confidence_1000bp_around_summit.bed pereira_neurog2_astrocytes_confidence_1000bp_around_summit.bed wang_ascl1_gi-men_confidence_1000bp_around_summit.bed manandhar_myod1_fibroblasts_confidence_1000bp_around_summit.bed li_twist2_myoblast_gm_treatment_confidence_1000bp_around_summit.bed conerly_myf5_mef_confidence_1000bp_around_summit.bed conerly_myod1_mef_confidence_1000bp_around_summit.bed
  do
    echo $file
    paste bedmap_pre*$file bedmap_post*$file | awk '{if ($5 == 0) print $1,$2,$3,$5,$10,"inf"; else print $1,$2,$3,$5,$10,$10/$5}' | tr " " "\t" > foldchange_induction_${file}
  done

  for file in aydin_neurog2_eb_treatment_48h_confidence_1000bp_around_summit.bed
  do
    echo $file
    paste bedmap_pre*$file bedmap_post*$file | awk '{if ($5 == 0) print $1,$2,$3,$5,$10,"inf"; else print $1,$2,$3,$5,$10,$10/$5}' | tr " " "\t" > foldchange_induction_${file}
  done




