six_replicatesinter () {
    bedtools multiinter -i $1*bed $2*bed $3*bed $4*bed $5*bed $6*bed | cut -f1,2,3,4 | bedtools merge -i - -c 4 -o max > multiinter.bed
    cat multiinter.bed | awk '{if ($4>1) print }' | tr ' ' '\t' > twointer.bed
    bedtools intersect -wa -a $1*bed -b twointer.bed | uniq > first.bed
    bedtools intersect -v -a twointer.bed -b $1*bed | uniq > notfirst.bed
    bedtools intersect -wa -a $2*bed -b notfirst.bed | uniq > second.bed
    bedtools intersect -v -a notfirst.bed -b $2*bed | uniq > notsecond.bed
    bedtools intersect -wa -a $3*bed -b notsecond.bed | uniq > third.bed
    bedtools intersect -v -a notsecond.bed -b $3*bed | uniq > notthird.bed
    bedtools intersect -wa -a $4*bed -b notthird.bed | uniq > fourth.bed
    bedtools intersect -v -a notthird.bed -b $4*bed | uniq > notfourth.bed
    bedtools intersect -wa -a $5*bed -b notfourth.bed | uniq > fifth.bed
    cat first.bed second.bed third.bed fourth.bed fifth.bed | sort -k1,1 -k2,2n | bedtools intersect -v -a - \
    -b ~/n/genomes/$7-blacklist.v2.bed | grep -v ^chr._ | grep -v ^chr.._ | grep -v chrEBV | uniq > ${8}_confidence.bed
}


# 4 distinct replicates
four_replicatesinter () {
  bedtools multiinter -i $1*bed $2*bed $3*bed $4*bed | cut -f1,2,3,4 | bedtools merge -i - -c 4 -o max > multiinter.bed
  # hacer un bedtools merge del multiinter para juntar todos los trocitos. y luego hacer el max sobre la columna 4. para ver el trocito q tiene maximo overlap entre replicates
  cat multiinter.bed | awk '{if ($4>1) print }' | tr ' ' '\t' > twointer.bed
  bedtools intersect -wa -a $1*bed -b twointer.bed | uniq > first.bed
  bedtools intersect -v -a twointer.bed -b $1*bed | uniq > notfirst.bed
  bedtools intersect -wa -a $2*bed -b notfirst.bed | uniq > second.bed
  bedtools intersect -v -a notfirst.bed -b $2*bed | uniq > notsecond.bed
  bedtools intersect -wa -a $3*bed -b notsecond.bed | uniq > third.bed
  cat first.bed second.bed third.bed | sort -k1,1 -k2,2n | bedtools intersect -v -a - \
  -b ~/n/genomes/$5-blacklist.v2.bed | grep -v ^chr._ | grep -v ^chr.._ | grep -v chrEBV | uniq > ${6}_confidence.bed
}

three_replicatesinter () {
  bedtools intersect -wa -a $1*bed -b $3*bed > jj2
  bedtools intersect -wa -a $2*bed -b $3*bed > jj3
  bedtools intersect -v -a jj2 -b jj3 > jj4
  cat jj4 jj3 | sort -k1,1 -k2,2n | bedtools intersect -v -a - \
  -b ~/n/genomes/$4-blacklist.v2.bed | grep -v ^chr._ | grep -v ^chr.._ | grep -v chrEBV> ${5}_confidence.bed
}

two_replicatesinter () {
  bedtools intersect -wa -a $1*bed -b $2*bed | bedtools intersect -v -a - \
  -b ~/n/genomes/$3-blacklist.v2.bed | grep -v ^chr._ | grep -v ^chr.._ | grep -v chrEBV > ${4}_confidence.bed
}

one_replicate () {
  bedtools intersect -v -a $1*bed \
  -b ~/n/genomes/$2-blacklist.v2.bed | grep -v ^chr._ | grep -v ^chr.._ | grep -v chrEBV > ${3}_confidence.bed
}

# aqui lo he dejado, para correr esto
cd ~/proneural/chip-seq
one_replicate GSM2482829 hg38 bosse_mycn_kelly
one_replicate GSM2482830 hg38 bosse_mycn_nb1643
one_replicate GSM2482828 hg38 bosse_mycn_ngp
one_replicate GSM2576033 mm10 casey_ascl1_mesc
one_replicate GSM2576036 mm10 casey_ascl2_mesc
one_replicate GSM2576039 mm10 casey_myod1_mesc
one_replicate GSM1332768 mm10 chalamalasetty_msgn1_eb
one_replicate GSM1954032 mm10 conerly_myf5_mef
one_replicate GSM1954031 mm10 conerly_myod1_mef
one_replicate GSM2220066 mm10 croci_myc_liver_control
one_replicate GSM2220064 mm10 croci_myc_liver_treatment
one_replicate GSM2220154 mm10 croci_myc_mef
one_replicate GSM2595124 mm10 depretis_myc_mef_control
one_replicate GSM2595125 mm10 depretis_myc_mef_treatment_10min
one_replicate GSM2595126 mm10 depretis_myc_mef_treatment_20min
one_replicate GSM2595128 mm10 depretis_myc_mef_treatment_2h
one_replicate GSM2595127 mm10 depretis_myc_mef_treatment_30min
one_replicate GSM2595129 mm10 depretis_myc_mef_treatment_4h
one_replicate GSM3488959 hg38 emming_bhlhe40_jurkat_control_clone1
one_replicate GSM3488960 hg38 emming_bhlhe40_jurkat_control_clone2
one_replicate GSM3488957 hg38 emming_bhlhe40_jurkat_treatment_clone1
one_replicate GSM3488958 hg38 emming_bhlhe40_jurkat_treatment_clone2
one_replicate GSM857389 mm10 fong_neurod2_p19
one_replicate GSM857390 mm10 fong_myod1_p19
one_replicate GSM857391 mm10 fong_myod1_mef
one_replicate GSM857392 mm10 fong_neurod2_mef
one_replicate GSM3044605 hg38 herold_mycn_shep_control
one_replicate GSM3044606 hg38 herold_mycn_shep_treatment
one_replicate GSM1715545 hg38 jaenicke_myc_imec
one_replicate GSM3375345 hg38 jain_tcf4_sudhl
one_replicate GSM3375347 hg38 jain_tcf4_tmd8
one_replicate GSM6674604 hg38 joung_ascl1_hesc
one_replicate GSM2049318 hg38 jung_myc_u2os_control
one_replicate GSM2049319 hg38 jung_myc_u2os_treatment
one_replicate GSM7213735 hg38 kim_neurod1_hek_treatment_pcag
one_replicate GSM7213605 hg38 kim_twist1_hek_treatment_pcag
one_replicate GSM4160355 mm10 lee_ascl1_mesc
one_replicate GSM4160356 mm10 lee_myod1_mesc
one_replicate GSM6828924 hg38 lee_twist1_bmsc_control
one_replicate GSM6828926 hg38 lee_twist1_bmsc_treatment
one_replicate GSM1436059 mm10 letourneau_sim2_mesc_a6
one_replicate GSM1436061 mm10 letourneau_sim2_mesc_b8
one_replicate GSM1436063 mm10 letourneau_sim2_mesc_c4
one_replicate GSM1436065 mm10 letourneau_sim2_mesc_eb3
one_replicate GSM3073948 hg38 liang_myc_h2171
one_replicate GSM3073949 hg38 liang_myc_sw1271
one_replicate GSM5026170 mm10 lin_mesp1_mesc_treatment_12h
one_replicate GSM5026171 mm10 lin_mesp1_mesc_treatment_24h
one_replicate GSM894059 hg38 lin_myc_p493-6_control
one_replicate GSM894060 hg38 lin_myc_p493-6_treatment_1h
one_replicate GSM894058 hg38 lin_myc_p493-6_treatment_24h
one_replicate GSM3659993 mm10 li_twist2_myoblast_dm_control
one_replicate GSM1276937 mm10 liu_ascl2_cd4_control
one_replicate GSM1276938 mm10 liu_ascl2_cd4_treatment
one_replicate GSM5712744 hg38 liu_atoh8_lungcarcinoma
one_replicate GSM2050640 hg38 lorenzin_myc_u2os_control
one_replicate GSM2050639 hg38 lorenzin_myc_u2os_treatment
one_replicate GSM2800529 mm10 matsuda_neurod1_microglia
one_replicate GSM766059 mm10 mazzoni_olig2_mesc
one_replicate GSM1423725 hg38 muthalagu_myc_mcf10a
one_replicate GSM5814092 hg38 neikes_mycmax_mcf7_treatment_0001nm
one_replicate GSM5814093 hg38 neikes_mycmax_mcf7_treatment_0010nm
one_replicate GSM5814094 hg38 neikes_mycmax_mcf7_treatment_0100nm
one_replicate GSM5814095 hg38 neikes_mycmax_mcf7_treatment_1000nm
one_replicate ERR453843 mm10 raposo_ascl1_ns5
one_replicate GSM1234499 hg38 sabo_myc_lymphoblastoma_control
one_replicate GSM1234501 hg38 sabo_myc_lymphoblastoma_highmyc
one_replicate GSM1234500 hg38 sabo_myc_lymphoblastoma_lowmyc
one_replicate GSM1386342 hg38 sabo_myc_lymphoblastoma_treatment_1h
one_replicate GSM1386343 hg38 sabo_myc_lymphoblastoma_treatment_24h
one_replicate GSM1234508 hg38 sabo_myc_mef_control
one_replicate GSM1234509 hg38 sabo_myc_mef_treatment
one_replicate GSM2049387 mm10 shang_hes1_mesc
one_replicate GSM1000576 hg38 soufi_myc_bj_control
one_replicate GSM896988 hg38 soufi_myc_bj_treatment
one_replicate GSM1463456 hg38 sugita_hes1_chrondrocytes
one_replicate GSM4104624 hg38 upton_mycn_cogn415
one_replicate GSM4104628 hg38 upton_mycn_lan5
one_replicate GSM4104630 hg38 upton_mycn_nb1643
one_replicate GSM1231598 hg38 walz_myc_u2os_control
one_replicate GSM1231597 hg38 walz_myc_u2os_treatment
one_replicate GSM6616013 hg38 wang_ascl1_gi-men
one_replicate GSM1187221 mm10 wapinski_ascl1_mef
one_replicate GSM1187228 mm10 wapinski_ascl1_npc
one_replicate GSM1486068 mm10 weber_hey1_cardiomyocytes_control
one_replicate GSM1486069 mm10 weber_hey1_cardiomyocytes_treatment
one_replicate GSM1486064 mm10 weber_hey1_mesc_control
one_replicate GSM1486065 mm10 weber_hey1_mesc_treatment
one_replicate GSM1486070 mm10 weber_hey2_cardiomyocytes_control
one_replicate GSM1486071 mm10 weber_hey2_cardiomyocytes_treatment
one_replicate GSM1486066 mm10 weber_hey2_mesc_control
one_replicate GSM1486067 mm10 weber_hey2_mesc_treatment
one_replicate GSM5708358 mm10 yu_twist1_mcf7_control
one_replicate GSM5708359 mm10 yu_twist1_mcf7_treatment
one_replicate GSM2572343 hg38 zeid_mycn_shep_control
one_replicate GSM2572346 hg38 zeid_mycn_shep_treatment_2h
one_replicate GSM2572349 hg38 zeid_mycn_shep_treatment_6h
two_replicatesinter GSM3136839 GSM3136840 mm10 aydin_ascl1_eb_treatment_48h
two_replicatesinter GSM3136849 GSM3136850 mm10 aydin_neurog2_eb_treatment_48h
two_replicatesinter GSM1907203 GSM1907204 hg38 barfeld_myc_lncap_control
two_replicatesinter GSM1907205 GSM1907206 hg38 barfeld_myc_lncap_treatment
two_replicatesinter GSM4059027 GSM4059005 mm10 costa_atoh1_mesc
two_replicatesinter GSM3734509 GSM3734510 mm10 fan_twist1_eb_control
two_replicatesinter GSM1575852 GSM1575853 mm10 fong_neurod2_p19_control
two_replicatesinter GSM1575850 GSM1575851 mm10 fong_neurod2_p19_treatment
two_replicatesinter GSM2589820 GSM2589821 hg38 gambardella_tfeb_hek_treatment_18h
two_replicatesinter GSM2589822 GSM2589823 hg38 gambardella_tfeb_hek_treatment_36h
two_replicatesinter GSM2589824 GSM2589825 hg38 gambardella_tfeb_hek_treatment_90h
two_replicatesinter GSM3660010 GSM3660011 mm10 li_twist2_myoblast_dm_treatment
two_replicatesinter GSM3659985 GSM3659986 mm10 li_twist2_myoblast_gm_control
two_replicatesinter GSM3660001 GSM3660002 mm10 li_twist2_myoblast_gm_treatment
two_replicatesinter GSM4080875 GSM4080876 hg38 louphrasitthiphol_mitf_melanoma_control
two_replicatesinter GSM4080877 GSM4080878 hg38 louphrasitthiphol_mitf_melanoma_treatment
two_replicatesinter GSM1586779 GSM1586780 mm10 pataskar_neurod1_mesc
two_replicatesinter GSM8242622 GSM8242623 mm10 pereira_neurog2_astrocytes_yy-
two_replicatesinter GSM8242624 GSM8242625 mm10 pereira_neurog2_astrocytes_yy+
two_replicatesinter GSM5018587 GSM5018588 hg38 see_myc_u2os_highmyc
two_replicatesinter GSM5018585 GSM5018586 hg38 see_myc_u2os_lowmyc
two_replicatesinter GSM1970150 GSM1970151 hg38 smith_neurog2_mrc5_treatment_05dpt
two_replicatesinter GSM1970152 GSM1970153 hg38 smith_neurog2_mrc5_treatment_1dpt
two_replicatesinter GSM1970154 GSM1970155 hg38 smith_neurog2_mrc5_treatment_2dpt
two_replicatesinter GSM4012629 GSM4012630 hg38 tai-nagara_tfe3_hk-2
two_replicatesinter GSM3289768 GSM3289769 hg38 tameire_myc_colorectaladenocarcinoma
two_replicatesinter GSM1493021 GSM1493018 hg38 thomas_myc_hek_treatment
two_replicatesinter GSM2124210 GSM2124211 mm10 velasco_neurog2_eb_treatment_12h
three_replicatesinter GSM3136828 GSM3136829 GSM3136830 mm10 aydin_ascl1_eb_treatment_12h
three_replicatesinter GSM3136832 GSM3136833 GSM3136834 mm10 aydin_neurog2_eb_treatment_12h
three_replicatesinter GSM3598002 GSM3598003 GSM3598004 mm10 lee_ascl1_mef
three_replicatesinter GSM3598008 GSM3598009 GSM3598010 mm10 lee_myod1_mef
three_replicatesinter GSM2449232 GSM2449233 GSM2449233 hg38 manandhar_myod1_fibroblasts
three_replicatesinter GSM2335528 GSM2335529 GSM2335530 hg38 park_ascl1_glioblastoma_control
three_replicatesinter GSM2335531 GSM2335532 GSM2335533 hg38 park_ascl1_glioblastoma_treatment
three_replicatesinter GSM6368738 GSM6368739 GSM6368740 mm10 pereira_neurog2_astrocytes
three_replicatesinter GSM1493025 GSM1493020 GSM1493017 hg38 thomas_myc_hek_control
three_replicatesinter GSM4654885 GSM4654886 GSM4654887 hg38 woods_ascl1_neuroblastoma_wt51_control
four_replicatesinter GSM3734511 GSM3734512 GSM3734513 GSM3734514 mm10 fan_twist1_eb_treatment
four_replicatesinter GSM7213721 GSM7213733 GSM7857134 GSM7857162 hg38 kim_twist1_hek_treatment_pcdna
four_replicatesinter GSM1973344 GSM1973345 GSM1973346 GSM1973347 mm10 kress_myc_liver_control
four_replicatesinter GSM6266261 GSM6266262 GSM6266264 GSM6266265 mm10 vainorius_ascl1_mesc
four_replicatesinter GSM6266267 GSM6266268 GSM6266270 GSM6266271 mm10 vainorius_neurog2_mesc
four_replicatesinter GSM4654869 GSM4654870 GSM4654871 GSM4654872 hg38 woods_ascl1_neuroblastoma_wt31_treatment
four_replicatesinter GSM4654873 GSM4654874 GSM4654875 GSM4654876 hg38 woods_ascl1_neuroblastoma_wt51_treatment
six_replicatesinter GSM1973348 GSM1973349 GSM1973350 GSM1973351 GSM1973352 GSM1973353 mm10 kress_myc_liver_treatment




one_replicate ERR453843 mm10 ERR453843_raposo_ascl1_ns5



bedtools intersect -wa -a ~/proneural/cuts/*/*/macs2.narrow/GSM8242622_FLAGNgn2_CUTRUN_Yy1plus_astrocytes_2dpneurog2_rep1_peaks.narrowPeak -b ~/proneural/cuts/*/*/macs2.narrow/GSM8242623_FLAGNgn2_CUTRUN_Yy1plus_astrocytes_2dpneurog2_rep2_peaks.narrowPeak | bedtools intersect -v -a - -b ~/n/genomes/mm10-blacklist.v2.bed | grep -v ^chr._ | grep -v ^chr.._ | grep -v chrEBV > pereira_neurog2_astrocytes_yy+_confidence.bed
bedtools intersect -wa -a ~/proneural/cuts/*/*/macs2.narrow/GSM8242624_FLAGNgn2_CUTRUN_Yy1minus_astrocytes_2dpneurog2_rep1_peaks.narrowPeak -b ~/proneural/cuts/*/*/macs2.narrow/GSM8242625_FLAGNgn2_CUTRUN_Yy1minus_astrocytes_2dpneurog2_rep2_peaks.narrowPeak | bedtools intersect -v -a - -b ~/n/genomes/mm10-blacklist.v2.bed | grep -v ^chr._ | grep -v ^chr.._ | grep -v chrEBV > pereira_neurog2_astrocytes_yy-_confidence.bed

bedtools intersect -wa -a ~/proneural/cuts/replicate_1_myod1_peaks.narrowPeak -b ~/proneural/cuts/replicate_2_myod1_peaks.narrowPeak | bedtools intersect -v -a - -b ~/n/genomes/mm10-blacklist.v2.bed | grep -v ^chr._ | grep -v ^chr.._ | grep -v chrEBV > hershbach_myod1_mef_confidence.bed




