
# bHLH Induction Data Processing

This directory contains the core computational framework for the multi-omic study of **basic helix-loop-helix (bHLH)** transcription factors. The pipeline is designed to integrate ChIP-seq, ATAC-seq, WGBS, and transcriptomic data to decode the "motif grammar" and chromatin remodeling capacities of 17 different factors.

## 📁 Pipeline Modules

### 1. [tf_chipseq_peaks]()

**Goal:** Acquisition and refinement of high-confidence binding sites.

* Automates the download of raw peak data and maps SRX accessions to GSM IDs.
* Filters peaks to generate "high-confidence" sets (present in at least 2 replicates) and removes ENCODE blacklisted regions.

### 2. [find_motifs]()

**Goal:** Decoding sequence syntax.

* Scans for canonical and non-canonical E-box variants (`CANNTG`).
* Analyzes flanking nucleotide preferences and intersects binding sites with the **Vierstra TF archetype** database to quantify motif centrality.



### 3. [epigenomic_data]()

**Goal:** Processing of chromatin state data.

* Handles alignment of raw reads to **hg38** and **mm10**.
* Maps accessibility signals (ATAC-seq/H3K27ac) onto ChIP-seq summits.
* Standardizes CpG methylation data and performs genome liftOver (e.g., mm9 to mm10).

### 4. [RNAseq_and_microarray_processing]()

**Goal:** Transcriptomic reprogramming analysis.

* A complete pipeline from raw SRA download to **DESeq2** and **Limma** differential expression.
* Includes specialized workflows for single-replicate studies and multi-study integration.


### 5. [htselex]()

**Goal:** Biochemical validation.

* Processes high-throughput SELEX and methyl-HT-SELEX data.
* Quantifies motif enrichment across successive rounds to identify intrinsic binding preferences and methylation-dependent shifts.

### 6. [gene_orthologs_mouse_human]()

**Goal:** Cross-species integration.

* Provides utilities for mapping GEO/SRA accessions and converting gene identifiers between mouse and human orthologs for comparative analysis.

---

## 🛠 Usage & Environment

Most scripts are optimized for a **Linux HPC environment** using the **Slurm** scheduler. Key dependencies include:

* **Bioinformatics:** `BWA`, `STAR`, `BEDTools`, `SAMtools`, `MACS2`, `featureCounts`.
* **Languages:** `R` (DESeq2, GEOquery, valr) and `Bash`.
* **Genomes:** Reference fasta files for **mm10** and **hg38**.
