
# bHLH Motif Analysis Toolkit

## 📜 Script Overview

### 1. `find_motifs.sh`

The primary utility for identifying E-box motifs within ChIP-seq peaks.

* 
**Core Functionality:** Utilizes R and `bedtools` to scan peak sequences for E-box patterns (`CANNTG`) and common degenerate variants.


* 
**Coordinate Management:** Automatically handles genome assembly differences (mm10 for mouse, hg38 for human) and centers analysis on peak summits.


* 
**Data Normalization:** Groups dinucleotide matches into biologically relevant E-box classes (e.g., CAT-CAT, CAG-CAG, CAC-CAC) based on the central dinucleotide.


* 
**Background Scanning:** Includes logic to generate "half-flank" peaks (background regions adjacent to actual peaks) to calculate motif enrichment indices.



### 2. `find_motifs_flanking.sh`

A specialized extension of the motif finder focused on **flanking nucleotide preferences**.

* 
**Sequence Context:** Extracts the 8-bp sequence (N-CANNTG-N) to determine how the nucleotides immediately preceding and following the core E-box influence binding specificity.


* 
**Strand Orientation:** Converts raw sequences into strand-oriented half-sites (e.g., transforming `ACATATGA` into a standardized `ACAT-TCAT` notation).


* 
**Biological Insight:** Used to validate findings such as the preference for a 'C' preceding CAT half-sites in proneural factors like NEUROD2 and NEUROG2.



### 3. `intersect_vierstra_archetypes.sh`

Integrates spatial motif data with **chromatin accessibility fold-change** data.

* 
**Vierstra Archetypes:** Maps peaks against the Vierstra transcription factor motif archetype database (v1.0) to identify secondary motifs enriched near bHLH binding sites.


* 
**Spatial Analysis:** Calculates the exact distance between the center of a motif and the peak summit to quantify **motif centrality**.


* 
**Quantile Stratification:** Splits ChIP-seq data into quartiles based on chromatin remodeling magnitude (log2 fold-change) to determine which motifs drive the most significant epigenetic changes.



---

## 🛠 Prerequisites

The scripts assume a Linux HPC environment with **Slurm** (for `sbatch` job submission) and the following modules/tools:

* 
**Alignment & Processing:** `BEDTools`.


* 
**Statistical Analysis:** `R` (with `dplyr`, `valr`, and `pracma` libraries).


* 
**Genomes:** Path access to `mm10.fa` and `hg38.fa` reference genomes.


