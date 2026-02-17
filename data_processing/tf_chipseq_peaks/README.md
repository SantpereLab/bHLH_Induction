# bHLH Induction ChIP-seq Processing Pipeline

This repository contains a suite of Bash scripts designed to automate the acquisition, quality control, and refinement of ChIP-seq datasets for basic helix-loop-helix (bHLH) transcription factors. The pipeline is optimized for handling induction experiments across diverse cell types and species (human/mouse).

---

## 🛠 Script Descriptions

### 1. `obtain_peaks.sh`
Handles the initial data retrieval and naming conventions:
* **Automated Download**: Retrieves `.bed` files from the ChIP-Atlas database for both `hg38` (human) and `mm10` (mouse) assemblies.
* **Accession Mapping**: Utilizes a mapping file (`srx_accessions_df.txt`) to rename raw SRX accessions to their corresponding GSM (GEO) identifiers for consistency.
* **Parallel Execution**: Uses `sbatch` to manage concurrent downloads across a compute cluster.

### 2. `high_confidence_peaks.sh`
A library of functions used to generate robust peak sets from biological replicates:
* **Replicate Intersection**: Includes functions (e.g., `two_replicatesinter` through `six_replicatesinter`) that define high-confidence peaks as those detected in at least two replicates.
* **Genomic Refinement**: Automatically excludes peaks intersecting with ENCODE blacklist regions to eliminate technical noise.
* **Chromosome Cleaning**: Removes non-canonical and non-nuclear sequences, such as `chrEBV` and `chr._`.
* **Peak Merging**: Employs `bedtools merge` to consolidate overlapping fragments and select the maximum overlap signal.

### 3. `count_peaks.sh`
Performs final quantification and dataset organization:
* **Peak Quantification**: Generates a summary table (`peak_counts.tsv`) containing GSE accessions and peak counts for every replicate and processed high-confidence set.

