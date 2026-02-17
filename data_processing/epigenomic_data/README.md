
# bHLH Induction Analysis Pipeline: Chromatin and Methylation Processing

This repository contains the computational workflow for processing multi-omic data—including ATAC-seq, Histone ChIP-seq, and CpG Methylation—to characterize the regulatory landscape and pioneer activity of bHLH transcription factors. These scripts automate the alignment, normalization, and integration of epigenomic signals with transcription factor binding sites as described in de Martin et al. (2026).

---

## 🛠 Script Overview

### 1. Pre-induction Chromatin Processing

* 
**`chromatin_alignment.sh`**: Downloads raw FASTQ files from the SRA using `fasterq-dump` and performs read alignment to the `hg38` or `mm10` reference genomes using `bwa-mem`.


* 
**`chromatin_alignment.sh`**: Filters data by removing duplicate reads with `samtools rmdup` and retaining only uniquely mapped reads.


* 
**`chromatin_mapping.sh`**: Converts aligned BAM files to RPKM-normalized bedGraph/bigWig format using `bamCoverage`.


* 
**`chromatin_mapping.sh`**: Assigns a pre-induction accessibility value to each ChIP-seq peak by taking the maximum signal within a 1,000 bp window centered on the peak summit.



### 2. Post-induction and Dynamics

* 
**`map_post_induction_chromatin_chipseq_peaks.sh`**: Maps post-induction accessibility signals directly onto the 1,000 bp windows of TF-bound regions to allow for direct comparison between conditions.


* 
**`fold_change_chromatin_data_processing.sh`**: Calculates the  fold-change in accessibility between post- and pre-induction states.


* 
**`fold_change_chromatin_data_processing.sh`**: Stratifies peaks by genomic location into promoters ( kb of TSS) and distal enhancers ( kb from TSS) for specialized statistical modeling.



### 3. DNA Methylation

* 
**`methylation_data_processing.sh`**: Processes cell-type-matched CpG methylation data, performing coordinate liftovers (e.g., `mm9` to `mm10`) where necessary using `liftOver`.


* 
**`methylation_data_processing.sh`**: Computes the average CpG methylation percentage within each ChIP-seq peak coordinate using `bedtools map`.



---

## 🧪 Methodology Summary

The pipeline facilitates the integrative multi-omic framework detailed in the study:

* 
**Chromatin States**: Accessibility data (ATAC-seq prioritized) is used to partition peaks into quartiles based on their pre-induction "openness".


* 
**Methylation Integration**: The processing reveals how high CpG methylation levels can redirect certain factors, such as MYC and HEY, toward atypical CAT-CAC E-box variants.


* 
**Pioneer Activity**: By comparing accessibility fold-changes, the pipeline identifies the specific capacity of CAT- and CAG-preferring bHLHs to remodel closed chromatin.



---

## 🚀 Quick Start

1. 
**Alignment**: Execute `bash chromatin_alignment.sh` to process raw reads into filtered BAM files.


2. 
**Normalization**: Use `bash chromatin_mapping.sh` to generate RPKM-normalized signal tracks and map baseline accessibility to summits.


3. 
**Dynamics**: Run `bash fold_change_chromatin_data_processing.sh` to quantify chromatin remodeling across enhancers and promoters.


4. 
**Methylation**: Incorporate CpG methylation levels using `bash methylation_data_processing.sh` for downstream motif-shift analysis.



---

Would you like me to help you configure the SRA accessions for the alignment script or provide the multivariate regression code used to analyze the resulting peak tables?
