
## 📂 Repository Structure

The repository consists of five primary scripts, each corresponding to a main figure in the manuscript:

| Script | Figure | Key Analysis |
| --- | --- | --- |
| `MF1.R` | **Figure 1** | E-box motif centrality, flanking nucleotide preferences, and motif clustering/spacing. |
| `MF2.R` | **Figure 2** | Chromatin accessibility heatmaps and multivariate modeling of accessibility vs. motif enrichment. |
| `MF3.R` | **Figure 3** | SELEX vs. Methyl-SELEX affinity comparisons and *in vivo* DNA methylation levels across binding sites. |
| `MF4.R` | **Figure 4** | Nucleosome positioning analysis (TF summit distributions relative to MNase-seq dyads). |
| `MF5.R` | **Figure 5** | Functional remodeling: Linear models of accessibility/expression fold-change and grouped Gene Ontology (GO) enrichment. |

---

## 🛠 Prerequisites

### R Packages

To run these scripts, you will need the following libraries installed:

* **Data Manipulation:** `dplyr`, `tidyr`, `reshape2`, `broom`, `valr`
* **Visualization:** `ggplot2`, `pheatmap`, `viridis`, `ggridges`, `patchwork`, `gridExtra`, `tidyheatmaps`
* **Stats/Modeling:** `scales`, `glmnet`, `zoo`, `signal`

### Data Dependencies

The scripts rely on several pre-processed `.RDS` files and genomic `.bed` files. Ensure your working directory contains the following (or update the paths in the scripts):

* `combined_canonical_eboxes.RDS`: Master motif database.
* `study_mapping.RDS`: Metadata for renaming and ordering datasets.
* `filtered_studies.RDS`: List of studies passing QC.
* `combined_data.RDS`: Integrated motif and accessibility data.

---

## 🔬 Script Details

### **MF1.R: The E-box Landscape**

This script characterizes the "anatomy" of TF peaks. It calculates the **centrality** (distance to peak summit) and **motif rate** (frequency) for different E-box variants (e.g., CAT-CAT vs CAT-CAG). It also generates the "balloon plots" that visualize how motif preferences shift across different proneural factors.

### **MF2.R: Chromatin Context**

Focuses on the relationship between DNA sequence and **pre-existing accessibility**. It includes multivariate linear models to determine which E-box variants are most predictive of low-accessibility binding sites.

### **MF3.R: DNA Methylation**

Integrates *in vitro* SELEX data with *in vivo* bisulfite sequencing. It compares how methylation-sensitive factors (like HEY1/2) behave in Methyl-HT-SELEX vs. how they bind to methylated CpG sites in different cell types.

### **MF4.R: Nucleosome Positioning**

Uses chemical  data to map **dyad positions**. It calculates the distribution of ChIP-seq summits around these dyads to determine if TFs prefer "translational" positions (on the nucleosome) or "linker" DNA.

### **MF5.R: Functional Remodeling & GO**

The "business end" of the pipeline. It models how motif enrichment predicts **log2 fold changes** in accessibility and gene expression during cellular induction. It concludes with a grouped GO analysis to show the biological processes regulated by different TF clusters.

---

## 🚀 Usage

1. **Set Working Directory:** Most scripts assume a structure relative to `~/proneural/chip-seq`. You may need to run `setwd()` accordingly.
2. **Load Order:** The scripts often load `orden.RDS` to ensure consistent factor ordering (e.g., NEUROD1 -> NEUROG2 -> ASCL1) across all plots.
3. **Execution:** Scripts are designed to be run top-to-bottom. They automatically save high-resolution PDFs of the figures.

---

> **Note:** If you are looking to replicate the multivariate linear models in `MF5.R`, ensure you have the fold-change accessibility data (`foldchange_accessibility_with_tss_distance.RDS`) available, as these models require the comparison between pre- and post-induction states.

**Would you like me to help you refine the path-handling logic in these scripts to make them more portable for other users?**
