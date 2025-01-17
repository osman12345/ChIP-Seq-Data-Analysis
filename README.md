# ChIP-Seq Analysis Script

## Overview

This repository contains an R script for analyzing ChIP-Seq data, specifically for histone modifications such as H3K27Ac, H3K27me3, and H3K4me1. The script leverages popular bioinformatics libraries like `ChIPseeker`, `clusterProfiler`, and `ReactomePA` to perform data visualization, functional enrichment analysis, and peak annotation.

---

## Features

### 1. **Peak Analysis**
- Visualizes ChIP-Seq peaks across chromosomes.
- Generates heatmaps and average profile plots.

### 2. **Peak Annotation**
- Annotates peaks relative to genomic regions such as TSS (transcription start sites).
- Produces pie charts, bar plots, and upset plots for genomic annotation.

### 3. **Functional Enrichment**
- Enriches pathways using the Reactome database.
- Visualizes pathway analysis results with dot plots and enrichment maps.

### 4. **Data Comparison**
- Compares ChIP-Seq datasets.
- Profiles multiple datasets with heatmaps, bar plots, and TSS binding visualizations.

### 5. **Statistical Testing**
- Conducts statistical testing for ChIP-Seq overlap enrichment.

---

## Dependencies

The following R packages are required to run the script:

- [`conflicted`](https://cran.r-project.org/package=conflicted)
- [`RColorBrewer`](https://cran.r-project.org/package=RColorBrewer)
- [`ChIPseeker`](https://bioconductor.org/packages/ChIPseeker)
- [`clusterProfiler`](https://bioconductor.org/packages/clusterProfiler)
- [`TxDb.Mmusculus.UCSC.mm10.knownGene`](https://bioconductor.org/packages/TxDb.Mmusculus.UCSC.mm10.knownGene)
- [`ReactomePA`](https://bioconductor.org/packages/ReactomePA)
- [`devtools`](https://cran.r-project.org/package=devtools)

To install these packages, run:
```R
install.packages(c("conflicted", "RColorBrewer", "devtools"))
BiocManager::install(c("ChIPseeker", "clusterProfiler", "TxDb.Mmusculus.UCSC.mm10.knownGene", "ReactomePA"))
