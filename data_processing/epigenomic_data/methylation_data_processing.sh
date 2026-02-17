cd ~/proneural/meth
  
for file in *gz
do
  sbatch --partition=bigmem --wrap="gunzip $file"
done

# Change the files to be bed files with 4th column showing the percentage of CpG methylation

# bed file y la 4 columna es el porcentaje de CpGs metiladas
GSM3071998_WGBS_mm9_CpG_blk.bed GSM2800531_Microglia_Bisulfite-seq.txt GSM3526804_MCF7_WGBS.bdg GSM6506983_A1.bedGraph GSM8042633_MEF_L001_R1_CpG.bedGraph

# bed file y la 4 columna es el porcentaje de CpGs metiladas. pero start y end coords son la misma
cat GSM7300318_EP-1_deduplicated.bismark.cov | awk '{print $1,$2-1,$3,$4}' | tr ' ' '\t' > GSM7300318_EP-1_deduplicated.bismark.bed

cat GSM7841629_210408_X129_FCHYK2VDSXY_L1_CHKPE85221030252_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov | awk '{print $1,$2-1,$3,$4}' | tr ' ' '\t' > GSM7841629_210408_X129_FCHYK2VDSXY_L1_CHKPE85221030252_1_val_1_bismark_bt2_pe.deduplicated.bismark.bed


# Only end position; and in [columns] 4 and 5, the number of methylated and unmethylated CpGs.
for file in GSM5135424_GFPGR1_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz.CpG_report.txt GSM5229302_Y1_2-1063673_S1.txt.CpG_report.txt GSM7437599_22144D-01-01.CpG_report.txt
do
  awk '{if ( ($4+$5) > 0 ) print $1,$2-1,$2,$4/($4+$5)}' $file | tr ' ' '\t' > $file.bed
done

### Rename
mv GSM3071998_WGBS_mm9_CpG_blk.bed mESC_GSM3071998_WGBS_mm9_CpG_blk.bed
mv GSM5135424_GFPGR1_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz.CpG_report.txt.bed SHSY5Y_GSM5135424_GFPGR1_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz.CpG_report.txt.bed
mv GSM5229302_Y1_2-1063673_S1.txt.CpG_report.txt.bed myoblast_GSM5229302_Y1_2-1063673_S1.txt.CpG_report.txt.bed
mv GSM7437599_22144D-01-01.CpG_report.txt.bed hek_GSM7437599_22144D-01-01.CpG_report.txt.bed
mv fixed_GSM2138821_ENCFF333EXV_methylation_state_at_CpG_GRCh38.bed hESC_fixed_GSM2138821_ENCFF333EXV_methylation_state_at_CpG_GRCh38.bed
mv GSM7300318_EP-1_deduplicated.bismark.bed bjcells_GSM7300318_EP-1_deduplicated.bismark.bed
mv GSM2800531_Microglia_Bisulfite-seq.txt Microglia_GSM2800531_Bisulfite-seq.bed
mv GSM3526804_MCF7_WGBS.bdg MCF7_GSM3526804_WGBS.bed
mv GSM7841629_210408_X129_FCHYK2VDSXY_L1_CHKPE85221030252_1_val_1_bismark_bt2_pe.deduplicated.bismark.bed astrocyte_GSM7841629_210408_X129_FCHYK2VDSXY_L1_CHKPE85221030252_1_val_1_bismark_bt2_pe.deduplicated.bismark.bed
mv GSM8042633_MEF_L001_R1_CpG.bedGraph MEF_GSM8042633_L001_R1_CpG.bed



# lift over

import_shiva proneural/meth/mESC_GSM3071998_WGBS_mm9_CpG_blk.bed .
import_shiva proneural/meth/hek_GSM7437599_22144D-01-01.CpG_report.txt.bed .
import_shiva proneural/meth/MEF_GSM8042633_L001_R1_CpG.bed .
import_shiva proneural/meth/MCF7_GSM3526804_WGBS.bed .


liftOver mESC_GSM3071998_WGBS_mm9_CpG_blk.bed ~/genomes/mm9ToMm10.over.chain.gz mm10.mESC_GSM3071998_WGBS_mm9_CpG_blk.bed unmapped

for file in hek_GSM7437599_22144D-01-01.CpG_report.txt.bed MCF7_GSM3526804_WGBS.bed
do 
  liftOver $file ~/genomes/hg19ToHg38.over.chain.gz hg38.$file unmapped
done

cat MEF_GSM8042633_L001_R1_CpG.bed | grep -v bedGraph > mef_GSM8042633_L001_R1_CpG.bed
liftOver mef_GSM8042633_L001_R1_CpG.bed ~/genomes/mm39ToMm10.over.chain.gz mm10.MEF_GSM8042633_L001_R1_CpG.bed unmapped

export_shiva mm10\*bed proneural/meth
export_shiva hg38\*bed proneural/meth




module load BEDTools

# function to get CpG methylation levels in quantiles
methylation_in_quantiles ( ) {
  sort -k1,1 -k2,2n ../meth/$2 > sorted.meth.bed
  bedtools map -c 4 -o mean -a $1 -b sorted.meth.bed > $1.methylation.bed
}


cd ~/proneural/chip-seq

# microglia
methylation_in_quantiles mouse/matsuda_neurod1_microglia_confidence.bed Microglia_GSM2800531_Bisulfite-seq.bed

# astrocytes
methylation_in_quantiles mouse/pereira_neurog2_astrocytes_confidence.bed astrocyte_GSM7841629_210408_X129_FCHYK2VDSXY_L1_CHKPE85221030252_1_val_1_bismark_bt2_pe.deduplicated.bismark.bed

# bj
for file in human/manandhar_myod1_fibroblasts_confidence.bed human/soufi_myc_bj_treatment_confidence.bed
do
  methylation_in_quantiles $file bjcells_GSM7300318_EP-1_deduplicated.bismark.bed
done

# mesc
for file in human/wapinski_ascl1_mesc_confidence.bed  mouse/letourneau_sim2_mesc_b8_confidence.bed mouse/vainorius_ascl1_mesc_confidence.bed  mouse/casey_ascl1_mesc_confidence.bed   mouse/letourneau_sim2_mesc_c4_confidence.bed    mouse/vainorius_neurog2_mesc_confidence.bed mouse/casey_ascl2_mesc_confidence.bed   mouse/letourneau_sim2_mesc_eb3_confidence.bed   mouse/weber_hey1_mesc_control_confidence.bed    mouse/casey_myod1_mesc_confidence.bed   mouse/lin_mesp1_mesc_treatment_12h_confidence.bed   mouse/weber_hey1_mesc_treatment_confidence.bed  mouse/costa_atoh1_mesc_confidence.bed   mouse/lin_mesp1_mesc_treatment_24h_confidence.bed   mouse/weber_hey2_mesc_control_confidence.bed    mouse/lee_ascl1_mesc_confidence.bed mouse/mazzoni_olig2_mesc_confidence.bed mouse/weber_hey2_mesc_treatment_confidence.bed  mouse/lee_myod1_mesc_confidence.bed mouse/pataskar_neurod1_mesc_confidence.bed  mouse/letourneau_sim2_mesc_a6_confidence.bed    mouse/shang_hes1_mesc_confidence.bed 
do
  methylation_in_quantiles $file mm10.mESC_GSM3071998_WGBS_mm9_CpG_blk.bed
done

# mef
for file in mouse/fong_myod1_mef_confidence.bed human/sabo_myc_mef_treatment_confidence.bed mouse/depretis_myc_mef_treatment_10min_confidence.bed   mouse/fong_neurod2_mef_confidence.bed   mouse/GSE211864_hershbach_myod1_mef_confidence.bed  mouse/depretis_myc_mef_treatment_20min_confidence.bed   mouse/hershbach_myod1_mef_confidence.bed    mouse/conerly_myf5_mef_confidence.bed   mouse/depretis_myc_mef_treatment_2h_confidence.bed  mouse/lee_ascl1_mef_confidence.bed  mouse/conerly_myod1_mef_confidence.bed  mouse/depretis_myc_mef_treatment_30min_confidence.bed   mouse/lee_myod1_mef_confidence.bed  mouse/croci_myc_mef_confidence.bed  mouse/depretis_myc_mef_treatment_4h_confidence.bed  mouse/wapinski_ascl1_mef_confidence.bed
do
  methylation_in_quantiles $file mm10.MEF_GSM8042633_L001_R1_CpG.bed
done

# hek
for file in human/kim_neurod1_hek_treatment_pcag_confidence.bed   human/kim_twist1_hek_treatment_pcdna_confidence.bed   human/kim_twist1_hek_treatment_pcag_confidence.bed
do
  methylation_in_quantiles $file hg38.hek_GSM7437599_22144D-01-01.CpG_report.txt.bed 
done




# myblasts (gm)
methylation_in_quantiles mouse/li_twist2_myoblast_gm_treatment_confidence.bed myoblast_GSM5229302_Y1_2-1063673_S1.txt.CpG_report.txt.bed

# SH-SY5Y
for file in human/woods_ascl1_neuroblastoma_wt31_treatment_confidence.bed human/woods_ascl1_neuroblastoma_wt51_treatment_confidence.bed
do
 methylation_in_quantiles $file SHSY5Y_GSM5135424_GFPGR1_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz.CpG_report.txt.bed
done


