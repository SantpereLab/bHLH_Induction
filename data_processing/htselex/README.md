# HT-SELEX Motif Enrichment Analysis

This section of the pipeline processes high-throughput systematic evolution of ligands by exponential enrichment (HT-SELEX) data to validate the intrinsic biochemical binding preferences of bHLH factors.

## 📜 Script Guide

### `htselex_motif_quantification.sh`

Processing and quantification of SELEX and methyl-HT-SELEX datasets.

* **Role:** Extracts raw motif frequency data from sequential SELEX cycles to determine high-affinity sequence variants and methylation-dependent shifts.


* **Function:**
* **Target Selection:** Filters the Lambert transcription factor table to identify human bHLH members and maps them to ENA accessions.


* **Data Acquisition:** Downloads raw FASTQ files from the Yin et al. (2017) study via the ENA, specifically prioritizing **methyl-HT-SELEX** over bisulfite sequencing for streamlined analysis of methylated DNA interactions.


* **Motif Scanning:** Uses `grep` to perform exact hexanucleotide matching for 16 canonical and non-canonical E-box variants (`CANNTG`) across every processed FASTQ file.


* **Read Normalization:** Calculates the total number of sequences per file to allow for the calculation of motif proportions per round.


* **Consolidation:** Aggregates motif counts across rounds 1 through 4 into a unified summary table for each factor, enabling the observation of enrichment trajectories and methylation-dependent shifts.
