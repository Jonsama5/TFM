# Microbiome Analysis for Master's Thesis

This repository contains the scripts used for the analysis of microbiome data related to my Master's thesis, which focuses on the metagenomic analysis of *Anopheles gambiae* microbiota from the 1000 Genomes Project.

---

## Pipeline Overview

* **Data source:**  
  Sequences were downloaded from the *Anopheles gambiae* dataset provided by phase 3 of the 1000 Genomes Project (MalariaGEN).

* **Quality control:**  
  Initial quality assessment was performed using **FastQC (v0.12.1)**. When necessary, reads were trimmed with Trimmomatic.

* **Host read removal:**  
  Reads belonging to *Anopheles gambiae* and *Homo sapiens* were removed using the aligners **BWA MEM (v0.7.17)** and **BWA SAMSE (v0.7.17)**, combined with **SAMtools (v1.12)** to extract unmapped reads.

* **Microbiome profiling:**  
  The filtered reads, presumed to represent microbiome content, were analyzed using the **KMA aligner (v1.3.15)** integrated into the **CCMetagen pipeline** for taxonomic classification against the NCBI NT database. Output tables with abundance and taxonomy for bacteria, viruses, protists, and fungi were generated and manually curated to include *Plasmodium falciparum*.

* **File formats:**  
  Results were exported as CSV files and imported into R for downstream analysis.

---

## Data Analysis Workflow in R

The data analysis pipeline was implemented in R (version 4.4.2) using the following key packages, with versions noted:

- **phyloseq (v1.5)** — handling microbiome data objects and analyses  
- **microbiome (v1.28)** — microbiome-specific metrics and visualizations  
- **vegan (v2.6.1)** — ecological and diversity statistics  
- **DESeq2 (v1.46)** — differential abundance analysis  
- **mvabund (v4.2.1)** — multivariate abundance modeling  
- **gt (v0.11.1)** — elegant tabular summaries  
- **ggplot2 (v3.5.2)** — graphics and visualization  
- **ggpubr (v0.4.0)** — publication-ready plots  
- **FSA (v0.9.6)** — post hoc tests including Dunn’s test  
- **readxl (v1.4.5)** — reading Excel files (for metadata or data imports)  
- **dplyr (v1.1.4)**, **tibble (v3.2.1)**, **tidyr (v1.3.1)** — data wrangling and reshaping  
- **scales (v1.4)** and **ggrepel (v0.9.6)** — enhanced plot labeling and scaling

---

### Steps Performed:

1. **Data Import and Wrangling:**  
   CSV files from CCMetagen were imported using `readxl` and transformed into abundance, taxonomy, and metadata tables suitable for `phyloseq` objects via `dplyr`, `tibble`, and `tidyr`.

2. **Exploratory Analysis:**  
   - Summaries of read counts and taxa were generated using `gt` and visualized with `ggplot2`.  
   - Diversity indices (Shannon, Simpson for alpha diversity; Bray-Curtis, Jaccard for beta diversity) were calculated using `phyloseq` and `vegan`.  
   - Multivariate analyses including PCA and NMDS helped visualize community structure.

3. **Statistical Testing:**  
   - Differential abundance analyses were conducted with `DESeq2`.  
   - Multivariate generalized linear models were fit with `mvabund` using negative binomial distributions to assess the impact of metadata variables.  
   - Correlations between taxa abundances and *Plasmodium* parasite load were analyzed using Spearman’s rank correlation.  
   - Non-parametric tests (Kruskal-Wallis and Dunn’s post hoc) were used to assess diversity differences, adjusted for multiple testing using FDR.

4. **Visualization and Reporting:**  
   - Final plots were generated using `ggplot2` and extensions (`ggpubr`, `ggrepel`, `scales`).  
   - Tables summarizing statistical results and metadata were prepared with `gt`.

---

## Notes

- Protist and fungal reads were not detected in the dataset except for *Plasmodium*, which was manually included in analyses.  
- The pipeline and analysis scripts are organized to facilitate reproducibility.  (TO BE ADDED)
- The entire workflow, including data preprocessing and statistical scripts, can be found in this repository. (TO BE ADDED)

---

## Repository structure

- `scripts/` — R scripts used for data import, cleaning, analysis, and visualization  
- `data/` — CSV files with taxonomic and abundance tables  
- `results/` — output figures and summary tables  
- `README.md` — this documentation  

---

## Contact & Further Information

For detailed methods and additional context, please refer to the thesis document and associated materials.

## License
This repository is licensed under the [MIT License](LICENSE) for its original code.

This project makes use of the [CCMetagen pipeline](https://github.com/vrmarcelino/CCMetagen), which is licensed under the [GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.html). CCMetagen is used as an external tool in this pipeline to process sequencing data, and is not modified or included directly in this repository.
