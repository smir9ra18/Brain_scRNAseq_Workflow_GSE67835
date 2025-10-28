# Brain_scRNAseq_Workflow_GSE67835
Single-cell RNA-seq analysis workflow for human brain (GSE67835).

A comprehensive single-cell RNA sequencing analysis pipeline for brain tissue using the GSE67835 dataset from NCBI GEO.

## 🧬 Overview

This repository contains a complete, reproducible workflow for analyzing single-cell RNA-seq data from mouse brain tissue (cortex and hippocampus). The pipeline performs quality control, normalization, dimensionality reduction, clustering, and cell type annotation to identify major brain cell populations.

## 🎯 Objectives

- **Download and preprocess** scRNA-seq data from NCBI GEO (GSE67835)
- **Perform quality control** to filter low-quality cells and genes
- **Normalize and identify** highly variable genes for downstream analysis
- **Reduce dimensionality** using PCA and UMAP for visualization
- **Cluster cells** using the Leiden algorithm to identify distinct populations
- **Identify marker genes** that distinguish different cell clusters
- **Annotate cell types** based on known brain cell markers
- **Generate publication-quality** visualizations and save results

## 📊 Dataset

**GEO Accession**: [GSE67835](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67835)

**Title**: Single-cell RNA-seq of mouse brain cells

**Organism**: *Mus musculus* (Mouse)

**Tissue**: Adult cortex and hippocampus

**Platform**: Single-cell RNA sequencing

## 🔧 Tools and Technologies

### Core Libraries
- **Scanpy** - Single-cell analysis in Python
- **GEOparse** - Downloading data from NCBI GEO
- **Pandas** & **NumPy** - Data manipulation
- **Matplotlib** & **Seaborn** - Data visualization

### Key Algorithms
- **Leiden Algorithm** - Community detection for cell clustering
- **UMAP** - Dimensionality reduction for visualization
- **PCA** - Linear dimensionality reduction
- **Wilcoxon Rank-Sum Test** - Differential expression analysis

## 🚀 Quick Start

### Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/brain-scrna-seq-workflow.git
cd brain-scrna-seq-workflow

# Install required packages
pip install scanpy GEOparse leidenalg python-igraph matplotlib seaborn pandas numpy scipy
```

### Run the Analysis

```bash
# Run the Jupyter notebook
jupyter notebook brain_scrna_seq_analysis.ipynb

# Or run as Python script
python brain_scrna_seq_analysis.py
```

## 📁 Repository Structure

```
brain-scrna-seq-workflow/
├── README.md
├── brain_scrna_seq_analysis.ipynb    # Main analysis notebook
├── data/                              # Downloaded raw data (auto-created)
└── results/
    ├── figures/                       # Generated plots
    │   ├── 01_qc_metrics.png
    │   ├── 02_highly_variable_genes.png
    │   ├── 03_pca_variance.png
    │   ├── 04_umap_clusters.png
    │   ├── 05_marker_heatmap.png
    │   ├── 06_marker_dotplot.png
    │   ├── 07_umap_cell_types.png
    │   ├── 08_marker_expression_umap.png
    │   ├── 09_marker_violin.png
    │   └── 10_complete_analysis_overview.png
    ├── tables/                        # Data tables and annotations
    │   ├── marker_genes_top10.csv
    │   ├── cell_type_annotations.csv
    │   ├── analysis_summary.csv
    │   └── cell_metadata.csv
    └── processed_adata.h5ad           # Processed AnnData object
```

## 🔬 Analysis Workflow

### Step 1: Data Acquisition
- Download GSE67835 dataset from NCBI GEO using GEOparse
- Load expression matrix and metadata

### Step 2: Quality Control
- Calculate QC metrics (total counts, genes per cell, mitochondrial %)
- Filter low-quality cells and genes
- Visualize QC distributions

### Step 3: Normalization & Preprocessing
- Normalize counts to 10,000 per cell
- Log-transform expression values
- Identify highly variable genes (HVGs)
- Scale data for downstream analysis

### Step 4: Dimensionality Reduction
- Perform PCA on HVGs
- Compute neighborhood graph
- Generate UMAP embedding for visualization

### Step 5: Cell Clustering
- Apply Leiden algorithm for community detection
- Identify distinct cell populations
- Visualize clusters on UMAP

### Step 6: Marker Gene Identification
- Compute differential expression between clusters
- Identify top marker genes per cluster
- Create heatmaps and dot plots

### Step 7: Cell Type Annotation
- Score cells for known brain cell type markers:
  - **Neurons** (Snap25, Slc17a7)
  - **GABAergic neurons** (Gad1)
  - **Astrocytes** (Aqp4, Gfap, Aldoc)
  - **Oligodendrocytes** (Mog, Mbp, Plp1)
  - **Microglia** (Cx3cr1)
  - **Endothelial cells** (Cldn5, Pecam1)
- Assign cell types based on highest marker scores

### Step 8: Visualization & Results
- Generate comprehensive plots
- Export annotated data and marker genes
- Save processed AnnData object

## 📈 Key Results

The pipeline generates:

- **10 publication-quality figures** showing QC metrics, clustering, cell types, and marker expression
- **Cell type annotations** for all analyzed cells
- **Top marker genes** for each cluster
- **Summary statistics** of the analysis
- **Processed data** ready for further downstream analysis

## 🧠 Identified Cell Types

The analysis identifies major brain cell populations:

1. **Excitatory Neurons** - Express Snap25, Slc17a7
2. **GABAergic Neurons** - Express Gad1
3. **Astrocytes** - Express Aqp4, Gfap, Aldoc
4. **Oligodendrocytes** - Express Mog, Mbp, Plp1
5. **Microglia** - Express Cx3cr1
6. **Endothelial Cells** - Express Cldn5, Pecam1

## 💻 System Requirements

- **Python**: 3.7 or higher
- **Memory**: 8GB RAM minimum (16GB recommended)
- **Storage**: ~2GB for data and results
- **OS**: Windows, macOS, or Linux

## 📝 Notes

- The notebook is optimized for laptop compatibility with subsampled data
- Full dataset analysis may require more computational resources
- Processing time: ~5-15 minutes depending on system specifications
- Results are reproducible with the same random seed

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## 📄 License

This project is open source and available under the MIT License.

## 📚 References

- **Scanpy**: Wolf, F. A. et al. (2018). SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology*, 19:15.
- **Leiden Algorithm**: Traag, V. A. et al. (2019). From Louvain to Leiden: guaranteeing well-connected communities. *Scientific Reports*, 9:5233.
- **UMAP**: McInnes, L. et al. (2018). UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction. *ArXiv*.

## 📧 Contact

For questions or feedback, please open an issue in this repository.

---

**Happy analyzing! 🔬✨**
