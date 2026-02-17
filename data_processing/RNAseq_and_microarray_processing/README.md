# Transcriptomic Pre-processing Pipeline

This pipeline standardizes the acquisition and quality control of RNA-seq data from diverse sources. It is designed to run in a Linux environment with Slurm-based HPC job submission.

## 📜 Script Guide

### 1. `01.GSM_to_SRA.r`

The entry point for data retrieval.

* **Role:** A script that converts GEO Sample Accessions (**GSM**) into Sequence Read Archive (**SRA**) identifiers.


* **Function:** It queries metadata using the `rentrez` package to resolve the underlying raw data links required for download.



### 2. `02.extract_sra_id_study_list.sh`

Metadata management and organization.

* **Role:** Parses metadata to create clean, structured lists of SRA IDs.


* **Function:** Organizes these IDs into study-specific groups to ensure that downstream processing maintains the biological context of each experiment.



### 3. `03.SRAtoFASTQ.sh`

Raw data acquisition.

* **Role:** The "heavy lifter" for data retrieval.


* **Function:** Uses the `SRA-Toolkit` (specifically `fasterq-dump`) to download raw sequencing reads and convert them into **FASTQ** format.



### 4. `04.fastp_preprocessing.sh`

Quality control and adapter trimming.

* **Role:** Prepares raw reads for alignment.


* **Function:** Uses `fastp` to perform ultra-fast quality filtering, base correction, and adapter trimming. It generates reports used to exclude samples with very low read counts (e.g., <1 million reads).



### 5. `05.multiQC.sh`

Aggregated quality reporting.

* 
**Role:** Batch visualization of QC metrics.


* **Function:** Runs `MultiQC` to combine the outputs of `fastp` and `FASTQC` across all processed samples into a single, interactive HTML report.



### 6. `06.STAR_genome_indexing.sh`

Genomic infrastructure preparation.

* **Role:** Prepares the reference framework for read mapping.


* **Function:** Uses the `STAR` aligner to generate a genomic index for **hg38** (human) or **mm10** (mouse) based on **GENCODE** annotations.



### 7. `07.STAR_alignment.sh`

Read mapping to reference genomes.

* **Role:** The primary alignment utility.


* **Function:** Uses the **STAR** universal aligner to map both single-end and paired-end reads to **hg38** or **mm10**, producing coordinate-sorted BAM files.



### 8. `08.featureCounts.sh`

Transcript quantification.

* **Role:** Quantifies genomic features.


* **Function:** Uses **featureCounts** to count mapped reads or fragments against GENCODE annotations, automatically detecting library types (SE vs. PE).



### 9. `09A.DESeq2_order.sh` & `10A.DESeq2_analysis.r`

Standard differential expression for studies with replicates.

* **Role:** Orchestrates standard differential expression analysis using **DESeq2**.


* **Function:** Employs a paired design to compare induced vs. uninduced samples and performs automated ENSEMBL-to-Symbol annotation mapping.



### 10. `09B.DESeq2_noReps_order.sh` & `10B.DESeq2_noReps_analysis.r`

Specialized workflow for single-replicate studies.

* **Role:** Analysis of experiments lacking biological replicates.


* **Function:** Implements a critical manual dispersion setting (e.g., 0.1) in **DESeq2** to enable statistical testing on single-sample datasets.



### 11. `11.microarray_full_processing.r`

Alternative platform processing.

* **Role:** A complete pipeline for legacy **Affymetrix microarray** data.


* **Function:** Uses the `affy` package for RMA normalization and the **Limma** package for differential expression analysis.



### 12. `12.merge_foldchanges_all.sh`

Final data integration.

* **Role:** Standardizes results across all platforms and species.


* **Function:** Iterates through all studies to create unified tables containing Gene Symbols, log2FoldChanges, and p-values for cross-factor comparison.
