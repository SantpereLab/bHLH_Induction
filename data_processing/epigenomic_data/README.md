This repository contains the processing pipelines and analysis scripts for the multi-omic study of **basic helix-loop-helix (bHLH)** transcription factors. This framework integrates 74 induction ChIP-seq datasets across 17 different factors to decode the "motif grammar" that governs how these proteins interact with accessible and closed chromatin.

---

# bHLH Induction Analysis Framework

This project provides a unified pipeline to process ChIP-seq, ATAC-seq, and WGBS data, specifically designed to investigate the pioneer activity and DNA-binding syntax of bHLH factors such as **ASCL1, NEUROG2, and MYC**.

## 🛠 Prerequisites

The scripts assume a Linux environment with the following tools installed and available in your `$PATH`:

* **Alignment & Processing:** `SRA-Toolkit`, `BWA`, `SAMtools`.


* 
**Genomic Analysis:** `BEDTools`, `UCSC liftOver`.


* 
**Peak Calling & Coverage:** `MACS2`, `deepTools` (specifically `bamCoverage`).



---

## 📜 Script Guide

### 1. `chromatin_alignment.sh`

The "heavy lifter" for raw data. This script automates the retrieval and initial processing of sequencing reads:

* 
**Data Acquisition:** Uses `fasterq-dump` to download datasets from SRA.


* 
**Alignment:** Maps reads to **mm10** (mouse) or **hg38** (human) using `bwa mem`.


* 
**Post-processing:** Performs duplicate removal (filtering out `XA` and `SA` tags for unique mapping), sorting, and indexing.


* 
**Coverage:** Generates RPKM-normalized **BedGraph** and **BigWig** files via `bamCoverage`.




### 2. `chromatin_mapping.sh`

The primary mapping utility for spatial analysis:

* 
**Windowing:** Resizes ChIP-seq peaks to **1000bp** windows centered on the summit.


* 
**Signal Integration:** Maps maximum accessibility signals from BedGraph files onto these windows across dozens of cell types (e.g., astrocytes, mESCs, fibroblasts).



### 3. `map_post_induction_chromatin_chipseq_peaks.sh`

A specialized analysis script focused on the co-dependence of binding and accessibility:

* **Post-Induction Mapping:** Specifically maps ATAC-seq/H3K27ac signals from the "after-induction" state onto the original bHLH binding sites.
* 
**Induction Comparison:** Calculates the direct relationship between factor binding and the subsequent remodeling of the local chromatin landscape.





### 4. `fold_change_chromatin_data_processing.sh`

Dedicated to identifying dynamic chromatin regions:

* 
**Peak Calling:** Uses `macs2` for both narrow (ATAC-seq) and broad (histone PTM) peaks.


* 
**Intersection:** Merges pre- and post-induction peak sets to identify *de novo* binding sites.


* 
**Fold Change:** Calculates the  fold-change of accessibility at peak summits between conditions.



### 5. `methylation_data_processing.sh`

Handles the integration of CpG methylation data across various experimental formats:

* **Standardization:** Converts diverse file types into a unified BED format where the 4th column represents the **% of CpG methylation**.
* **Coordinate Management:** Includes logic to adjust 0-based vs 1-based start/end coordinates.
* 
**Genome LiftOver:** Updates older methylation datasets (e.g., mm9 to mm10 or hg19 to hg38) to ensure compatibility with recent accessibility maps.


* 
**Quantile Calculation:** Extracts methylation levels within specific genomic windows using `bedtools map`.



