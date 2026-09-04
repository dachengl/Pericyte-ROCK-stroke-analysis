# Pericyte-ROCK-stroke-analysis

This repository contains the analysis code for the manuscript:
**"Single-Cell Transcriptomic Analysis Reveals Pericyte Subtype-Specific ROCK Signaling Suppression and Phenotypic Transition in Ischemic Stroke"**

## Dependencies

All analyses were performed using R (v4.2.0) with the following packages:

- Seurat (v4.3.0)
- CellChat (v1.6.1)
- Monocle (v2.26.0)
- SCENIC (v1.3.1)
- clusterProfiler (v4.4.4)
- GSVA (v1.44.5)

## File description

| File | Description |
|------|-------------|
| `r.01_Singlecell.R` | Single-cell RNA-seq data processing, quality control, clustering, and cell type annotation using Seurat |
| `r.02_Expression.R` | Differential expression analysis and ROCK pathway activity scoring |
| `r.03_GSEA.R` | Gene Ontology (GO) and KEGG enrichment analysis, and Gene Set Enrichment Analysis (GSEA) |
| `r.04_Immune.R` | Immune infiltration analysis (ssGSEA) and correlation analysis |

## Data availability

The raw sequencing data are publicly available in the NCBI Gene Expression Omnibus (GEO) under accession numbers:
- **GSE174574** (single-cell RNA-seq)
- **GSE104036** (bulk RNA-seq)

## Citation

If you find this code useful, please cite our manuscript (pending publication).
