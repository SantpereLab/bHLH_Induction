# Transcriptomic Pre-processing Pipeline

This pipeline standardizes the acquisition and quality control of RNA-seq data from diverse sources. It is designed to run in a Linux environment with Slurm-based HPC job submission.

## 📜 Script Guide

### 1. `01.GSM_to_SRA.r`

The entry point for data retrieval.

* **Role:** A script that converts GEO Sample Accessions (**GSM**) into Sequence Read Archive (**SRA**) identifiers.
* 
**Function:** It queries metadata to resolve the underlying raw data links required for download.



### 2. `02.extract_sra_id_study_list.sh`

Metadata management and organization.

* 
**Role:** Parses the output from the previous step to create clean, structured lists of SRA IDs.


* 
**Function:** Organizes these IDs into study-specific groups to ensure that downstream processing maintains the biological context of each experiment.



### 3. `03.SRAtoFASTQ.sh`

Raw data acquisition.

* **Role:** The "heavy lifter" for data retrieval.
* 
**Function:** Uses the `SRA-Toolkit` (specifically `fasterq-dump`) to download raw sequencing reads and convert them into **FASTQ** format.



### 4. `04.fastp_preprocessing.sh`

Quality control and adapter trimming.

* 
**Role:** Prepares raw reads for alignment.


* 
**Function:** Uses `fastp` to perform ultra-fast quality filtering, base correction, and adapter trimming. It generates summary reports used to exclude samples with very low read counts (e.g., <1 million reads).



### 5. `05.multiQC.sh`

Aggregated quality reporting.

* **Role:** Batch visualization of QC metrics.
* 
**Function:** Runs `MultiQC` to combine the outputs of `fastp` and `FASTQC` across all processed samples into a single, interactive HTML report.



### 6. `06.STAR_genome_indexing.sh`

Genomic infrastructure preparation.

* 
**Role:** Prepares the reference framework for read mapping.


* 
**Function:** Uses the `STAR` aligner to generate a genomic index for **hg38** (human) or **mm10** (mouse) based on **GENCODE** annotations. This is a prerequisite for the high-speed universal alignment required for multi-omic integration.



This section continues the transcriptomic pipeline, covering read alignment, gene quantification, differential expression analysis, and final data integration across RNA-seq and microarray platforms.

---

## 7. Read Mapping & Quantification

After quality control, processed reads are aligned to reference genomes to quantify transcriptional activity.

* **`07.STAR_alignment.sh`**: The primary alignment utility. It uses the **STAR** universal aligner to map both single-end and paired-end reads to **hg38** or **mm10**. It is optimized for high-performance clusters using parallel jobs and produces coordinate-sorted BAM files.


* **`08.featureCounts.sh`**: Quantifies genomic features. This script uses **featureCounts** (part of the Subread package) to count mapped reads or fragments against GENCODE annotations. It automatically detects library types (single-end vs. paired-end) to ensure accurate quantification across diverse study designs.



## 8. Differential Expression Analysis (RNA-seq)

The framework employs specialized scripts to handle varying experimental designs, particularly for studies with and without biological replicates.

* 
**`09A.DESeq2_order.sh` & `10A.DESeq2_analysis.r**`: These scripts orchestrate standard differential expression for studies with replicates using **DESeq2**. The analysis uses a paired design when possible to compare induced vs. uninduced samples within the same replicate.


* **`09B.DESeq2_noReps_order.sh` & `10B.DESeq2_noReps_analysis.r**`: A specialized workflow for single-replicate studies. Because standard DESeq2 requires replicates to estimate dispersion, script **10B** implements a critical manual dispersion setting (e.g., 0.1) to allow statistical testing on these datasets.
* **Key Features**: Both R scripts include automated annotation mapping using `org.Hs.eg.db` and `org.Mm.eg.db` to link ENSEMBL IDs to gene symbols for human and mouse samples.

## 9. Alternative Platform Processing (Microarray)

To broaden the study's scope, the pipeline incorporates older datasets that utilized different genomic technologies.

* **`11.microarray_full_processing.r`**: A complete pipeline for **Affymetrix microarray** data. It uses the `affy` package for RMA normalization and the **Limma** package for differential expression analysis. It identifies the specific GSM accessions needed and generates normalized expression sets for the final integration.



## 10. Multi-Study Data Integration

The final stage merges results from all previous steps into a unified format for cross-factor comparison.

* **`12.merge_foldchanges_all.sh`**: The "master merger" script. It iterates through all processed mouse and human studies, standardizing the output columns (Gene Symbol, log2FoldChange, p-value).
* 
**Methods**: It accounts for the differences between RNA-seq (which includes ENSEMBL IDs) and microarray outputs to generate unified species-specific tables. These tables serve as the foundation for the shared ortholog mapping and Pearson correlation clustering described in the paper.







  


