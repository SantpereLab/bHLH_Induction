    mapp () {
      sbatch --partition=bigmem --mem=300000 --wrap="bedtools map -o max -c 4 -a */$1 -b $2 > chromatin_${1}"
    }

    module load BEDTools
    


    cd ~/proneural/chip-seq
    
    # Microglia (mouse)
    mapp matsuda_neurod1_microglia_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.*goselin*bdg
    
    # Astrocytes
    mapp pereira_neurog2_astrocytes_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.pereira_astrocytes_GSM6368722_ATACseq_GFP.sam.rmdup.bam.bdg
    mapp pereira_neurog2_astrocytes_yy-_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.pereira_astrocytes_GSM6368748_ATACseq_yy1ko.sam.rmdup.bam.bdg
    mapp pereira_neurog2_astrocytes_yy+_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.pereira_astrocytes_GSM6368746_ATACseq_yy1wt.sam.rmdup.bam.bdg
    
    # BJ
    mapp manandhar_myod1_fibroblasts_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.*risca*atac.sam*bdg
    
    # MCF10a
    mapp jaenicke_myc_imec_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.liu_MCF10a.sam.processed.bam.bdg
    
    # P19
    mapp fong_myod1_p19_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.*salbert*bdg
    mapp fong_neurod2_p19_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.*salbert*bdg
    mapp fong_neurod2_p19_treatment_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.*salbert*bdg


    
    # EB
    mapp aydin_ascl1_eb_treatment_12h_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.aydin_0h_atac_2.sam.processed.bam.bdg
    mapp chalamalasetty_msgn1_eb_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.velasco_0h_EBs_ATAC.sam.processed.bam.bdg
    mapp fan_twist1_eb_treatment_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.velasco_0h_EBs_ATAC.sam.processed.bam.bdg
    mapp aydin_neurog2_eb_treatment_48h_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.velasco_0h_EBs_ATAC.sam.processed.bam.bdg
    mapp velasco_neurog2_eb_treatment_12h_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.velasco_0h_EBs_ATAC.sam.processed.bam.bdg


    # MRC-5
    mapp smith_neurog2_mrc5_treatment_1dpt_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.*smith_ATAC_MRC5_4DPT*bdg
    
    # Primary myoblast (differentiation medium)
    mapp li_twist2_myoblast_dm_treatment_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.li_GFP*DM_K27ac_1*.bdg


    # Glioblastoma
    mapp park_ascl1_glioblastoma_treatment_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.*park_atac_control*.bdg
    
    # mESCs
    mapp casey_ascl1_mesc_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.*casey_esc_uninduced_ascl1_atac*bdg
    mapp casey_ascl2_mesc_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.*casey_esc_uninduced_ascl2_atac*bdg
    mapp casey_myod1_mesc_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.*casey_esc_uninduced_myod1_atac*bdg
    
    for file in pataskar_neurod1_mesc_confidence_1000bp_around_summit.bed vainorius_ascl1_mesc_confidence_1000bp_around_summit.bed vainorius_neurog2_mesc_confidence_1000bp_around_summit.bed mazzoni_olig2_mesc_confidence_1000bp_around_summit.bed lee_ascl1_mesc_confidence_1000bp_around_summit.bed lee_myod1_mesc_confidence_1000bp_around_summit.bed weber_hey1_mesc_treatment_confidence_1000bp_around_summit.bed weber_hey2_mesc_treatment_confidence_1000bp_around_summit.bed
    do
      mapp $file ~/proneural/bdgs/clean.*casey_esc_uninduced_ascl1_atac*bdg
    done


    # MEF
    mapp wapinski_ascl1_mef_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.wapinski_MEF.sam.processed.bam.bdg
    
    for file in wapinski_ascl1_mef_confidence_1000bp_around_summit.bed croci_myc_mef_confidence_1000bp_around_summit.bed depretis_myc_mef_treatment_2h_confidence_1000bp_around_summit.bed fong_myod1_mef_confidence_1000bp_around_summit.bed fong_neurod2_mef_confidence_1000bp_around_summit.bed hershbach_myod1_mef_confidence_1000bp_around_summit.bed lee_ascl1_mef_confidence_1000bp_around_summit.bed lee_myod1_mef_confidence_1000bp_around_summit.bed
    do
      mapp $file ~/proneural/bdgs/clean.lee_mef_atac.sam.processed.bam.bdg
    done
    
    mapp conerly_myod1_mef_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.lee_mef_atac.sam.processed.bam.bdg
    mapp conerly_myf5_mef_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.lee_mef_atac.sam.processed.bam.bdg
    
    # U2OS osteosarcoma
    for file in lorenzin_myc_u2os_treatment_confidence_1000bp_around_summit.bed walz_myc_u2os_treatment_confidence_1000bp_around_summit.bed jung_myc_u2os_treatment_confidence_1000bp_around_summit.bed
    do
      mapp $file ~/proneural/bdgs/clean.niu_GSM7996359_ATAC_U2OS_WT_biol.sam.processed.bam.bdg
    done
    
    # P493
    mapp lin_myc_p493-6_treatment_24h_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.lin_P493-6_T0_H3K27AC.fastq.gz.sam.processed.bam.bdg
    mapp sabo_myc_lymphoblastoma_highmyc_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.lin_P493-6_T0_H3K27AC.fastq.gz.sam.processed.bam.bdg 
    mapp sabo_myc_lymphoblastoma_treatment_24h_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.lin_P493-6_T0_H3K27AC.fastq.gz.sam.processed.bam.bdg


    # Mouse pluripotent stem cells
    mapp lin_mesp1_mesc_treatment_24h_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.lin_GSM5026158_ATACseq_WT0h.sam.rmdup.bam.bdg
    
    # BJ 
    mapp soufi_myc_bj_treatment_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.*risca*atac.sam*bdg
    
    
    # NPC (mouse)
    mapp raposo_ascl1_ns5_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.*carter*atac*bdg
    mapp wapinski_ascl1_npc_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.*carter*atac*bdg


    # HEK cells
    mapp kim_neurod1_hek_treatment_pcag_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.dong_HEK293T_DMSO_ATAC-seq*bdg
    mapp kim_twist1_hek_treatment_pcag_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.dong_HEK293T_DMSO_ATAC-seq*bdg
    mapp kim_twist1_hek_treatment_pcdna_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.dong_HEK293T_DMSO_ATAC-seq*bdg
    
    # Primary myoblast (growth medium)
    mapp li_twist2_myoblast_gm_treatment_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.li_Tw2*GM_K27ac_1.fastq.gz.sam.processed.bam.bdg
    
    
    # CD4+ T cells (mouse)
    mapp liu_ascl2_cd4_treatment_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.miraldi_cd4tcells_ATAC*bdg
    
    # 501 mel
    mapp louphrasitthiphol_mitf_melanoma_treatment_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.cigrang_GSM8085812_501mel_ATAC.sam.processed.bam.bdg


    # GI-MEN
    mapp wang_ascl1_gi-men_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.wang_ATAC_GIMEN_ASCL1_Neg.sam.processed.bam.bdg
    
    # SH-SY5Y
    mapp woods_ascl1_neuroblastoma_wt31_treatment_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.zimmerman_SHSY5Y_ATAC.fastq.gz.sam.processed.bam.bdg
    mapp woods_ascl1_neuroblastoma_wt51_treatment_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.zimmerman_SHSY5Y_ATAC.fastq.gz.sam.processed.bam.bdg


    # COGN415
    mapp upton_mycn_cogn415_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.upton_NB-1643_ATAC-Seq.sam.processed.bam.bdg
    
    # LAN5
    mapp upton_mycn_lan5_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.upton_LA-N-5_ATAC.sam.processed.bam.bdg
    
    # NB1643
    mapp upton_mycn_nb1643_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.upton_NB-1643_ATAC-Seq.sam.processed.bam.bdg
    mapp bosse_mycn_nb1643_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.upton_NB-1643_ATAC-Seq.sam.processed.bam.bdg


    # NGP
    mapp bosse_mycn_ngp_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.upton_NGP_ATAC-Seq.sam.processed.bam.bdg
    
    # Kelly
    mapp bosse_mycn_kelly_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.upton_KELLY_ATAC-Seq.sam.processed.bam.bdg
    
    
    ## SUDHL
    mapp jain_tcf4_sudhl_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.grubert_GSM4130913_ATACSeq_SUDHL.sam.processed.bam.bdg
    
    # TMD8
    mapp jain_tcf4_tmd8_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.bal_GSM5525194_TMD8_H3K27Ac.sam.processed.bam.bdg

    # LNCaP
    mapp barfeld_myc_lncap_treatment_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.liu-DMSO-ATAC-seq.sam.processed.bam.bdg
    
    # YB5 (colorrectal adenocarcinoma, like ATCC CCL-221)
    mapp tameire_myc_colorectaladenocarcinoma_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.liver_zhang_ATAC*bdg
    
    # MCF10a (mammary epithelial cell line, IMEC)
    mapp muthalagu_myc_mcf10a_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.liu_MCF10a.sam.processed.bam.bdg
    
    # CORL88 (small cell lung cancer, for H2171)
    mapp liang_myc_h2171_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.patel_ATAC_CORL88.sam.processed.bam.bdg
    
    # UCLA2 cells (human embryonic ste cells, like HUES66)
    mapp joung_ascl1_hesc_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.chen_ATACseqUCLA2hESCrep1.sam.processed.bam.bdg
    
    # Cardiomyocytes (mouse)
    mapp weber_hey1_cardiomyocytes_treatment_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.hill_GSM3504371_control_CM_ATACseq.sam.processed.bam.bdg
    mapp weber_hey2_cardiomyocytes_treatment_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.hill_GSM3504371_control_CM_ATACseq.sam.processed.bam.bdg
    
    # Liver
    mapp croci_myc_liver_treatment_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.liver_cusanovich*bdg
    mapp kress_myc_liver_treatment_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.liver_cusanovich*bdg
    
    # SHEP
    mapp zeid_mycn_shep_treatment_6h_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.*zeid_SHEP_0HR_H3K27AC.*.bdg
    mapp herold_mycn_shep_treatment_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.*zeid_SHEP_0HR_H3K27AC.*.bdg


    # Jurkat
    mapp emming_bhlhe40_jurkat_treatment_clone1_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.brignall_jurkat_ATAC.sam.rmdup.bam.bdg
    mapp emming_bhlhe40_jurkat_treatment_clone2_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.brignall_jurkat_ATAC.sam.rmdup.bam.bdg

    # MCF7
    mapp neikes_mycmax_mcf7_treatment_1000nm_confidence_1000bp_around_summit.bed ~/proneural/bdgs/clean.mcf7_hafner_h3k27ac.fastq.gz.sam.processed.bam.bdg
