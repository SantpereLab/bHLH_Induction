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

## 🔬 Script Details

### **SF1.R & SF2.R: Motif Refinement**

* **SF1.R** defines the canonical E-box as **CA..TG** and splits the motif database into canonical and non-canonical sets. It also uses `bedtools` to intersect ChIP-seq summits with Vierstra motif archetypes.
* **SF2.R** generates high-density facet plots comparing the "summit-hugging" distribution of canonical variants against the flat distribution of non-canonical control motifs.

### **SF3.R: Background Modeling**

* To prove enrichment, this script calculates the frequency and spacing of E-boxes in the 200bp **flanking regions** (upstream/downstream) of actual peaks.
* It utilizes `biasaway` to generate 100 iterations of k-mer and GC-content matched background sequences for every study.

### **SF4.R & SF5.R: Chromatin & Methylation States**

* **SF4.R** classifies peaks into **Promoters**, **Proximal Enhancers**, and **Distal Enhancers** using `ChIPseeker` and genomic annotations (hg38/mm10).
* **SF5.R** processes methylation `.bed` files to calculate the median **% mCpG** per peak, subsequently visualizing how methylation density shifts across accessibility quantiles.

### **SF6.R & SF7.R: Comparative Dynamics**

* **SF6.R** evaluates "functional binding" by comparing the accessibility fold-change of regions bound by TFs vs. those that remain unbound after induction.
* **SF7.R** computes a global correlation matrix using **orthologous genes** between human and mouse studies to identify conserved regulatory programs.

### **SF8.R: Biological Processes (GSEA)**

* This script runs `gseGO` (Biological Process) for every study using `clusterProfiler`.
* It generates the master **GO Balloon Plot**, where the size represents  and the color represents the **Normalized Enrichment Score (NES)**.

---

## 🛠 Usage Notes

1. **System Dependencies:** `SF1.R` and `SF3.R` utilize `system()` calls to external tools like **BEDTools**, **biasaway**, and **SLURM** (`sbatch`).
2. **Mapping:** Most scripts require `study_mapping.RDS` and `orden.RDS` to synchronize naming conventions across figures.
3. **Genomes:** Pathing in `SF4.R` points to Gencode v47 (hg38) and vM24 (mm10).

