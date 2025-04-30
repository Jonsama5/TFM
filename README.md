# TFM
This repository contains the scripts used for the analysis of microbiome data related to my Master's thesis.

-- Pipeline --
Sequences were downloaded from the 1000 Genomes Project Anopheles gambiae repository and subjected to quality control using FastQC, followed by Trimmomatic when necessary. Subsequently, a series of alignments were performed to remove Anopheles and Homo sapiens reads from the samples using BWA MEM and BWA SAMSE. The remaining sequences — corresponding to microbiome content — were analyzed using the CCMetagen pipeline. Sequences were aligned to the NCBI NT database, and abundance and taxonomic classification data for bacteria, viruses, protists, and fungi were extracted and converted into comma-separated values (CSV) for further analysis in R.

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
