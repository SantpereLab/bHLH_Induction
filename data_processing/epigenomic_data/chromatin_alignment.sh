#!/usr/bin/bash

## DESCARGAR FASTQS
######################
######################
#####################
#####################
#####################
####################
####################

mkdir p/proneural; cd p/proneural
module load SRA-Toolkit

#!/usr/bin/bash

# Microglia (mouse)
sbatch --partition=bigmem --wrap="prefetch SRR5617659; fasterq-dump SRR5617659 -o microglia_goselin_atac.fastq"


# P19
sbatch --partition=bigmem --wrap="prefetch SRR3633688; fasterq-dump SRR3633688 -o p19_salbert_h3k27ac.fastq"
sbatch --partition=bigmem -n 50 --wrap="prefetch SRR396792; fasterq-dump -e 50 SRR396792 -o fong_p19_control_ach4.fastq"
sbatch --partition=bigmem -n 50 --wrap="prefetch SRR396793; fasterq-dump -e 50 SRR396793 -o fong_p19_neurod2_induction_ach4.fastq"

# EB (Embryoid Bodies)

sbatch --partition=bigmem --wrap="prefetch SRR7133003; fasterq-dump SRR7133003 -o aydin_0h_atac_2.fastq"
sbatch --partition=bigmem -n 16 --wrap="prefetch SRR7132991; fasterq-dump -e 16 SRR7132991 -o aydin_12h_ascl1_atac_rep2.fastq"
sbatch --partition=bigmem -n 16 --wrap="prefetch SRR7132999; fasterq-dump -e 16 SRR7132999 -o aydin_48h_neurog2_atac_rep2.fastq"
sbatch --partition=bigmem -n 16 --wrap="prefetch SRR3407192; fasterq-dump -e 16 SRR3407192 -o velasco_0h_EBs_ATAC_rep2.fastq"



# MRC-5
sbatch --partition=bigmem -n 50 --wrap="prefetch SRR3438293; fasterq-dump -e 50 SRR3438293 -o smith_ATAC_MRC5_4DPT.fastq"


# Primary myoblast (DM - Differentiation Medium)
sbatch --partition=normal -n 16 --wrap="prefetch SRR8695970; fasterq-dump -e 16 SRR8695970 -o ~/proneural/mouse/single/li_GFP-DM_K27ac_1.fastq"
sbatch --partition=normal -n 16 --wrap="prefetch SRR8695988; fasterq-dump -e 16 SRR8695988 -o ~/proneural/mouse/single/li_Tw2-DM_K27ac_1.fastq"


# Glioblastoma
sbatch --partition=bigmem --wrap="prefetch SRR5278783; fasterq-dump SRR5278783 -o glioblastoma_wang_h3k27ac.fastq"


## B-cell lymphoma (P493-6)
sbatch --partition=normal -n 16 --wrap='prefetch SRR444433; fasterq-dump -e 16 SRR444433 -o lin_P493-6_T0_H3K27AC.fastq'
sbatch --partition=normal -n 16 --wrap='prefetch SRR444435; fasterq-dump -e 16 SRR444435 -o lin_P493-6_T24_H3K27AC.fastq'


# SHEP
sbatch --partition=normal -n 16 --wrap='prefetch SRR5441926; fasterq-dump -e 16 SRR5441926 -o zeid_SHEP_0HR_H3K27AC.fastq'
sbatch --partition=normal -n 16 --wrap='prefetch SRR5441932; fasterq-dump -e 16 SRR5441932 -o zeid_SHEP_6HR_H3K27AC.fastq'

# mESCs (mouse Embryonic Stem Cells)
sbatch --partition=normal -n 16 --wrap="prefetch SRR5445251; fasterq-dump -e 16 SRR5445251 -o casey_esc_uninduced_ascl1_atac.fastq"
sbatch --partition=normal -n 16 --wrap="prefetch SRR5445253; fasterq-dump -e 16 SRR5445253 -o casey_esc_uninduced_ascl2_atac.fastq"
sbatch --partition=normal -n 16 --wrap="prefetch SRR5445255; fasterq-dump -e 16 SRR5445255 -o casey_esc_uninduced_myod1_atac.fastq"
sbatch --partition=bigmem -n 16 --wrap="prefetch SRR5445252; fasterq-dump -e 16 SRR5445252 -o casey_24h_ascl1_atac.fastq"
sbatch --partition=bigmem -n 16 --wrap="prefetch SRR5445254; fasterq-dump -e 16 SRR5445254 -o casey_24h_ascl2_atac.fastq"
sbatch --partition=bigmem -n 16 --wrap="prefetch SRR5445256; fasterq-dump -e 16 SRR5445256 -o casey_24h_myod1_atac.fastq"


# MEF (Mouse Embryonic Fibroblasts)
sbatch --partition=normal -n 16 --wrap="prefetch SRR5822253; fasterq-dump -e 16 SRR5822253 -o lee_mef_atac.fastq"
sbatch --partition=normal -n 16 --wrap="prefetch SRR5822268; fasterq-dump -e 16 SRR5822268 -o lee_ascl1_48h_atac_rep1.fastq"
sbatch --partition=normal -n 16 --wrap="prefetch SRR8559020; fasterq-dump -e 16 SRR8559020 -o lee_myod1_48h_atac_rep1.fastq"

sbatch --partition=normal -n 16 --wrap='prefetch SRR5822253; fasterq-dump -e 16 SRR5822253 -o wapinski_MEF.fastq'
sbatch --partition=normal -n 16 --wrap='prefetch SRR5822271; fasterq-dump -e 16 SRR5822271 -o wapinski_MEF+Ascl1_48hr_B2_rep2.fastq'


# U2OS osteosarcoma
sbatch --partition=normal -n 16 --wrap='prefetch SRR1204519; fasterq-dump -e 16 SRR1204519 -o walz_U2OS_AcH4_+Dox_ChIPseq.fastq'
sbatch --partition=normal -n 16 --wrap='prefetch SRR1204520; fasterq-dump -e 16 SRR1204520 -o walz_U2OS_AcH4_-Dox_ChIPseq.fastq'

# BJ
sbatch --partition=bigmem --wrap="prefetch SRR3571007; fasterq-dump SRR3571007 -o bj_risca_atac_run1.fastq"
sbatch --partition=bigmem --wrap="prefetch SRR3571008; fasterq-dump SRR3571008 -o bj_risca_atac_run2.fastq"


# MCF-7
sbatch --partition=bigmem --wrap="prefetch SRR5638551; fasterq-dump SRR5638551 -o mcf7_nair_atac.fastq"

# NPC (mouse)
sbatch --partition=bigmem --wrap="prefetch SRR3933615; fasterq-dump SRR3933615 -o npc_carter_atac_run1.fastq"
sbatch --partition=bigmem --wrap="prefetch SRR3933616; fasterq-dump SRR3933616 -o npc_carter_atac_run2.fastq"

# HEK cells
sbatch --partition=normal -n 16 --wrap='prefetch SRR28879545; fasterq-dump -e 16 SRR28879545 -o dong_HEK293T_DMSO_ATAC-seq_rep1.fastq'

# Primary myoblast (GM - Growth Medium)
sbatch --partition=normal -n 16 --wrap="prefetch SRR8695962; fasterq-dump -e 16 SRR8695962 -o li_GFP-GM_K27ac_1.fastq"
sbatch --partition=normal -n 16 --wrap="prefetch SRR8695979; fasterq-dump -e 16 SRR8695979 -o li_Tw2-GM_K27ac_1.fastq"


# CD4+ T cells (mouse)
sbatch --partition=bigmem -n 50 --wrap="prefetch SRR7071512; fasterq-dump SRR7071512 -e 50 -o run1_miraldi_cd4tcells_ATAC.fastq"
sbatch --partition=bigmem -n 50 --wrap="prefetch SRR7071513; fasterq-dump SRR7071513 -e 50 -o run2_miraldi_cd4tcells_ATAC.fastq"
sbatch --partition=normal --wrap="cat run1_miraldi_cd4tcells_ATAC_1.fastq run2_miraldi_cd4tcells_ATAC_1.fastq > miraldi_cd4tcells_ATAC_1.fastq"
sbatch --partition=normal --wrap="cat run1_miraldi_cd4tcells_ATAC_2.fastq run2_miraldi_cd4tcells_ATAC_2.fastq > miraldi_cd4tcells_ATAC_2.fastq"

# Jurkat T lymphocytes (human)

sbatch --partition=normal -n 50 --wrap="prefetch SRR5063990; fasterq-dump SRR5063990 -e 50 -o brignall_jurkat_ATAC.fastq"

# Human fibroblasts / MRC-5 etc (incluido arriba en MRC-5 o BJ)

# GI-MEN
sbatch --partition=normal -n 16 --wrap='prefetch SRR21803001; fasterq-dump -e 16 SRR21803001 -o wang_ATAC_GIMEN_ASCL1_Neg_Rep1.fastq'
sbatch --partition=normal -n 16 --wrap='prefetch SRR21803000; fasterq-dump -e 16 SRR21803000 -o wang_ATAC_GIMEN_ASCL1_Pos_Rep1.fastq'

# SH-SY5Y
sbatch --partition=normal -n 16 --wrap='prefetch SRR5819663; fasterq-dump -e 16 SRR5819663 -o zimmerman_SHSY5Y_ATAC_run1.fastq'
sbatch --partition=normal -n 16 --wrap='prefetch SRR5819664; fasterq-dump -e 16 SRR5819664 -o zimmerman_SHSY5Y_ATAC_run2.fastq'

# COGN415
sbatch --partition=normal -n 16 --wrap='prefetch SRR10215662; fasterq-dump -e 16 SRR10215662 -o upton_COG-N-415_ATAC-Seq_REP1.fastq'

# LAN5
sbatch --partition=normal -n 16 --wrap='prefetch SRR10215670; fasterq-dump -e 16 SRR10215670 -o upton_LA-N-5_ATAC_REP1.fastq'

# NB1643
sbatch --partition=normal -n 16 --wrap='prefetch SRR10215674; fasterq-dump -e 16 SRR10215674 -o upton_NB-1643_ATAC-Seq_REP1.fastq'

# NGP
sbatch --partition=normal -n 16 --wrap='prefetch SRR10215680; fasterq-dump -e 16 SRR10215680 -o upton_NGP_ATAC-Seq_REP1.fastq'

# Kelly
sbatch --partition=normal -n 16 --wrap='prefetch SRR10215668; fasterq-dump -e 16 SRR10215668 -o upton_KELLY_ATAC-Seq_REP1.fastq'

# LNCaP
sbatch --partition=normal -n 16 --wrap='prefetch SRR2646275; fasterq-dump -e 16 SRR2646275 -o barfeld_H3K27ac_R1881.fastq'
sbatch --partition=normal -n 16 --wrap='prefetch SRR2646276; fasterq-dump -e 16 SRR2646276 -o barfeld_H3K27ac_R1881_Dox.fastq'
sbatch --partition=normal -n 16 --wrap='prefetch SRR3622817; fasterq-dump -e 16 SRR3622817 -o liu_Rep1-DMSO-ATAC-seq.fastq'

# YB5 (colorrectal adenocarcinoma)
sbatch --partition=normal -n 16 --wrap='prefetch SRR7060329; fasterq-dump -e 16 SRR7060329 -o zhang_yb5_atac_control_2.fastq'


# HT29 (colorrectal adenocarcinoma, like ATCC CCL-221)
sbatch --partition=normal -n 16 --wrap='prefetch SRR3110399; fasterq-dump -e 16 SRR3110399 -o savic_H3K27ac_DMSO_Rep1.fastq'


# MCF10a (mammary epithelial cell line)
sbatch --partition=normal -n 16 --wrap='prefetch SRR4436440; fasterq-dump -e 16 SRR4436440 -o liu_MCF10a_rep5.fastq'

# IMEC
sbatch --partition=normal -n 16 --wrap='prefetch SRR4125262; fasterq-dump -e 16 SRR4125262 -o fagnocchi_IMEC_WT_H3K27ac.fastq'

# UCLA2 cells (human embryonic stem cells)
sbatch --partition=normal -n 16 --wrap='prefetch SRR7942726; fasterq-dump -e 16 SRR7942726 -o chen_ATACseqUCLA2hESCrep1.fastq'


# Cardiomyocytes (mouse)
sbatch --partition=normal -n 16 --wrap='prefetch SRR5318031; fasterq-dump -e 16 SRR5318031 -o ~/proneural/mouse/single/ziemann_cardiomyocyte_atac_P1_rep1.fastq'

# HK-2
sbatch --partition=normal -n 16 --wrap='prefetch SRR16303905; fasterq-dump -e 16 SRR16303905 -o patel_ATAC-seq_HK2_EV_Rep1.fastq'

# Astrocytes (mouse)
sbatch --partition=normal -n 90 --wrap="prefetch SRR20341400; fasterq-dump SRR20341400 -e 90 -o pereira_astrocytes_GSM6368722_ATACseq_GFP_rep1.fastq"
sbatch --partition=normal -n 90 --wrap="prefetch SRR20341403; fasterq-dump SRR20341403 -e 90 -o pereira_astrocytes_GSM6368725_ATACseq_Ngn2_rep1.fastq"


# Mouse pluripotent stem cells
sbatch --partition=normal -n 80 --wrap="prefetch SRR13483382; fasterq-dump SRR13483382 -e 90 -o lin_GSM5026158_ATACseq_WT0h_rep1.fastq"


# Melanoma
sbatch --partition=normal -n 90 --wrap='prefetch SRR5228489; fasterq-dump -e 90 SRR5228489 -o fontanals_501MEL_atac.fastq'

# SUDHLH2 (OCI-LY7 / B-cell)
sbatch --partition=normal -n 16 --wrap='prefetch SRR4235647; fasterq-dump -e 16 SRR4235647 -o encode_H3K27ac_ChIP-seq_from_OCI-LY7_run1.fastq'








#juntar cuando hay mas de una run
cd ~/proneural/human/single
ls *run[1234].fastq | sed -e 's/_run[1234].fastq//g' | uniq > runs_single.txt
while read line
do
  sbatch --partition=bigmem --wrap="cat $line*.fastq > ${line}.fastq "
done<runs_single.txt

cd ~/proneural

rm human/single/*run*

#gzippear todos los fastqs
for file in */*/*fastq
do
  sbatch --wrap="gzip $file"
done


##### ALIGNMENT
#####################
#####################
#####################
#####################
#####################
####################
####################


cd ~/proneural

# single-end
###################


cd ~/proneural/mouse/single
for file in *gz
do sbatch -n 50 --partition=short --wrap="bwa mem -t 50 ~/n/genomes/mm10.fa $file > $file.sam"; done


cd ~/proneural/human/single
for file in *gz
do sbatch -n 40 --partition=long --wrap="bwa mem -t 40 ~/n/genomes/hg38.fa $file > $file.sam"; done


## paired-end
#############

#mouse
ls *_[12].fastq.gz | sed -e 's/_[12].fastq.gz//g' > paired_mouse_samples.txt
while read line
do
  sbatch -n 40 -J $line.align --partition=bigmem --wrap="bwa mem -t 40 ~/p/genomes/mm10.fa ${line}_1.fastq.gz ${line}_2.fastq.gz > $line.sam "; 
done<paired_mouse_samples.txt

#human
ls *_[12].fastq.gz | sed -e 's/_[12].fastq.gz//g' > paired_human_samples.txt
while read line
do 
  sbatch -n 100 -J $line.align --partition=bigmem --mem=200000 --wrap="bwa mem -t 100 ~/p/genomes/hg38.fa ${line}_1.fastq.gz ${line}_2.fastq.gz > $line.sam "; 
done<paired_human_samples.txt




####### PROCESSING #########
############################

#single-end


## las nuevas
cd ~/proneural/human/paired
for file in *sam
do
  sbatch --partition=bigmem --wrap="samtools rmdup $file - | samtools view -h > $file.rmdup.sam;
  cat $file.rmdup.sam | grep -v -e 'XA:Z:' -e 'SA:Z:' > $file.uniqmap; cat $file.uniqmap | samtools view -bS - > $file.processed.bam"
done

cd ~/proneural/mouse/single
for file in *sam
do
  sbatch --wrap="samtools rmdup -s $file - | samtools view -h > $file.rmdup.sam;
  cat $file.rmdup.sam | grep -v -e 'XA:Z:' -e 'SA:Z:' > $file.uniqmap; cat $file.uniqmap | samtools view -bS - > $file.processed.bam"
done

cd ~/proneural/human/single
for file in *sam
do
  sbatch --wrap="samtools rmdup -s $file - | samtools view -h > $file.rmdup.sam;
  cat $file.rmdup.sam | grep -v -e 'XA:Z:' -e 'SA:Z:' > $file.uniqmap; cat $file.uniqmap | samtools view -bS - > $file.processed.bam"
done

cd ~/proneural/mouse/paired
for file in *sam
do
  sbatch --partition=bigmem --wrap="samtools view -bS $file | samtools rmdup - $file.rmdup.bam"
done








rm *uniqmap
rm *rmdup*

###### OBTAIN COVERAGE FILES
###################################################

module load Miniconda3/4.12.0
module load SAMtools

cd ~/proneural/human/paired
for file in *bam
do 
  sbatch -o $file.out --partition=short --wrap=" samtools sort -o sorted.$file $file; wait; samtools index sorted.$file; wait; \
  bamCoverage --bam sorted.$file --outFileName $file.bdg --outFileFormat bedgraph --normalizeUsing RPKM; bamCoverage --bam sorted.$file --outFileName $file.bw --normalizeUsing RPKM"
done

cd ~/proneural/human/single
for file in *bam
do 
  sbatch -o $file.out --partition=bigmem --wrap=" samtools sort -o sorted.$file $file; wait; samtools index sorted.$file; wait; \
  bamCoverage --bam sorted.$file --outFileName $file.bdg --outFileFormat bedgraph --normalizeUsing RPKM; bamCoverage --bam sorted.$file --outFileName $file.bw --normalizeUsing RPKM"
done

cd ~/proneural/mouse/single

for file in *bam
do 
  sbatch -o $file.out --partition=bigmem --wrap=" samtools sort -o sorted.$file $file; wait; samtools index sorted.$file; wait; \
  bamCoverage --bam sorted.$file --outFileName $file.bdg --outFileFormat bedgraph --normalizeUsing RPKM; bamCoverage --bam sorted.$file --outFileName $file.bw --normalizeUsing RPKM"
done


cd ~/proneural/mouse/paired
for file in *bam
do 
  sbatch -o $file.out --partition=bigmem --wrap=" samtools sort -o sorted.$file $file; wait; samtools index sorted.$file; wait; \
  bamCoverage --bam sorted.$file --outFileName $file.bdg --outFileFormat bedgraph --normalizeUsing RPKM; bamCoverage --bam sorted.$file --outFileName $file.bw --normalizeUsing RPKM"
done


cd ~/proneural/mouse/single

for file in *bdg
do
  sbatch --wrap=" cat $file | grep -v "_" > clean.$file "
done




cd ~/n/proneural
for file in *bam
do
  sbatch -o ~/proneural/$file.out --partition=bigmem --wrap=" samtools sort -o ~/proneural/sorted.$file $file; wait; samtools index ~/proneural/sorted.$file; wait; bamCoverage --bam ~/proneural/sorted.$file --outFileName ~/proneural/$file.bdg --outFileFormat bedgraph --normalizeUsing RPKM "
done

for file in *bdg
do 
  sbatch --wrap=" cat $file | sort -k1,1 -k2,2n | grep -v "_" > clean.$file "
done







