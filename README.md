# Microbiome Analysis for Master's Thesis
This repository contains the scripts used for the analysis of microbiome data related to my Master's thesis, which focuses on the metagenomic analysis of Anopheles gambiae microbiota from the 1000 Genomes Project.

-- Pipeline Overview --
*Data source: Sequences were downloaded from the Anopheles gambiae dataset provided by the phase 3 of the 1000 Genomes Project (MalariaGEN).
Quality control: Initial quality assessment was performed using FastQC (v0.12.1). When necessary, reads were trimmed with Trimmomatic.
Host read removal: Reads belonging to Anopheles gambiae and Homo sapiens were removed using the aligners BWA MEM (v0.7.17) and BWA SAMSE (v0.7.17), combined with SAMtools (v1.12) to extract unmapped reads.
Microbiome profiling: The filtered reads, presumed to represent microbiome content, were analyzed using the KMA aligner (v1.3.15) integrated into the CCMetagen pipeline for taxonomic classification against the NCBI NT database. Output tables with abundance and taxonomy for bacteria, viruses, protists, and fungi were generated and manually curated to include Plasmodium falciparum.
File formats: Results were exported as CSV files and imported into R for downstream analysis.
-- Data analysis --
After an initial review of the data, reads were not detected for the protist and fungal kingdoms, leaving only bacterial and viral datasets for further study. The following steps of the thesis are documented in this repository:
1. Data preprocessing: Required libraries are installed, and datasets are converted into a suitable format for analysis.
2. Exploratory analysis: Techniques such as Principal Component Analysis (PCA) and Non-metric Multidimensional Scaling (NMDS) are used to visualize the data and detect patterns or other phenomena.
3. Characterization: Metrics (e.g., average read count, sparsity, singletons) and compositional plots at various taxonomic ranks are generated to visualize the structure of the microbiome.
4. Diversity: Common diversity metrics are calculated and visualized.
5. ANCOM and network analysis: Differential abundance of specific taxa is evaluated using ANCOM, and network plots are generated to explore the dynamics of microbial communities.
6. Hypothesis testing: Statistical tests and visualizations are used to examine associations between taxa and the presence of key microorganisms.

## License
This repository is licensed under the [MIT License](LICENSE) for its original code.

This project makes use of the [CCMetagen pipeline](https://github.com/vrmarcelino/CCMetagen), which is licensed under the [GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.html). CCMetagen is used as an external tool in this pipeline to process sequencing data, and is not modified or included directly in this repository.
