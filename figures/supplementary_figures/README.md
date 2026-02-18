## 📂 Repository Structure: Supplementary Figures

| Script | Figures | Key Analysis |
| --- | --- | --- |
| `SF1.R` | **SF1a, c, d** | Canonical E-box filtering, motif archetype centrality, and peak overlap (control vs. treatment). |
| `SF2.R` | **SF2** | Comparative distribution analysis of canonical vs. non-canonical E-boxes relative to peak summits. |
| `SF3.R` | **SF3a, b** | Background sequence modeling: E-box counts and spacing distributions in flanking genomic regions. |
| `SF4.R` | **SF4a, b, c** | Quantile-based accessibility analysis, TSS distance classification, and spacing proportions for Twist factors. |
| `SF5.R` | **SF5a, b** | DNA methylation density across chromatin accessibility quantiles and methylation-windowed balloon plots. |
| `SF6.R` | **SF6a, b, c** | Fold-change distributions for bound vs. unbound peaks and motif archetype centrality across ATAC-seq peaks. |
| `SF7.R` | **SF7a, b, c** | Multivariate linear models for accessibility/expression and study-to-study Pearson correlation matrices. |
| `SF8.R` | **SF8** | GSEA-based Gene Ontology (GO) enrichment pipelines for mouse and human datasets. |


---

## 🔬 Script Overview

### **SF1.R: Motif Filtering and Peak Characterization**

* **Canonical Filtering**: Defines and extracts canonical E-boxes (sequences matching the **CA..TG** motif) and non-canonical variants from ChIP-seq data.
* **Centrality Testing**: Performs Kolmogorov-Smirnov (KS) tests against a uniform distribution to identify studies with centrally enriched motifs.
* **Archetype Intersections**: Uses `bedtools` to intersect peak summits with **Vierstra motif archetypes**, calculating the median distance of each archetype from the summit.
* **Overlap Analysis**: Quantifies the percentage of peak overlap between control (pre-induction) and treatment (post-induction) conditions.

### **SF2.R: Motif Distribution Comparison**

* **Positional Analysis**: Summarizes the frequency of canonical vs. non-canonical E-boxes in 5bp and 20bp windows around the peak summit.
* **Ridge Visualization**: Generates large-scale ridge/area plots showing the high-density "summit-hugging" behavior of canonical E-boxes compared to the flat distribution of non-canonical controls.

### **SF3.R: Background and Spacing Modeling**

* **Flanking Control**: Calculates motif counts in the 200bp upstream and downstream flanking regions of peaks to serve as a genomic background.
* **Sequence Generation**: Implements a `biasaway` pipeline to generate 100 sets of GC-content and k-mer matched background sequences for every study.
* **Spacing Heatmaps**: Analyzes the distribution of distances (1–15 bp) between adjacent E-boxes in these background sequences to establish a baseline for non-functional clustering.

### **SF4.R: Genomic Features and Accessibility Quantiles**

* **TSS Classification**: Categorizes peaks into **Promoters** (≤2kb from TSS), **Proximal Enhancers** (2-4kb), and **Distal Enhancers** (>4kb) using `ChIPseeker`.
* **Accessibility Heatmaps**: Visualizes the proportion of genomic region types across 100 accessibility quantiles, highlighting the prevalence of promoters in highly accessible sites.
* **Twist-Specific Spacing**: Detailed analysis of E-box spacing proportions specifically for Twist-related factors across accessibility states.

### **SF5.R: DNA Methylation Landscapes**

* **Methylation Integration**: Merges ChIP-seq peak coordinates with DNA methylation (`.bed`) files to calculate the median **% mCpG** per peak.
* **Density Ridges**: Generates plots showing how methylation levels shift across different accessibility quantiles.
* **Windowed Balloon Plots**: Visualizes E-box centrality and frequency within specific methylation ranges (e.g., 0-25%, 25-50%, etc.) to detect methylation-dependent binding patterns.

### **SF6.R: Functional Binding and Archetype Centrality**

* **Fold-Change Boxplots**: Compares the log-fold change of accessibility in regions bound by a TF (intersect "yes") versus similar regions that remain unbound (intersect "no").
* **ATAC-Summit Centrality**: Analyzes the centrality of various motif archetypes relative to **ATAC-seq summits** (rather than ChIP-seq summits) in the top 4th quantile of accessibility.

### **SF7.R: Multivariate Modeling and Cross-Study Correlation**

* **Refined LM Fits**: Performs multivariate linear regression for accessibility fold-change vs. motif enrichment, specifically controlling for initial chromatin accessibility.
* **Ortholog Correlation**: Uses `biomaRt` to map human and mouse genes to common orthologs, allowing for the creation of a Pearson correlation matrix that compares gene expression programs across all 37+ studies.

### **SF8.R: Gene Ontology Enrichment**

* **GSEA Pipeline**: Runs `gseGO` (Biological Process) on human and mouse differential expression data using `clusterProfiler`.
* **Cluster Visualization**: Generates the master dendrogram and balloon plots for all GO terms, summarizing the biological pathways (e.g., "Neural development," "Skeletal muscle development") activated by each factor.

---

## 🛠 Usage and System Requirements

* **External Tools**: Several scripts (`SF1.R`, `SF3.R`) require **BEDTools** and **Miniconda** (for `biasaway`) to be available in the system path.
* **RDS Files**: Ensure that `study_mapping.RDS` and `orden.RDS` are present in your working directory to maintain consistent naming and factor ordering.
* **Output**: Plots are exported as high-resolution PDFs to the `~/Dropbox/induction/` directory structure.

**Would you like me to generate a dependency map showing which scripts must be run in sequence to produce the final `combined_data.RDS` used by these figures?**
