# Analysis Code for HGSOC CTRC-CELA3A-CELA3B Spatial Multi-Omics Study

This repository contains the R scripts used for the spatial multi-omics analysis described in the manuscript: "A Pancreatic-Like CTRC–CELA3A–CELA3B Program at the Tumor-Stroma Interface of High-Grade Serous Ovarian Cancer: A Spatial Multi-Omics Study".

## Contents
The R script includes the complete pipeline for:
- WGCNA co-expression network analysis
- Single-cell RNA-seq processing (Seurat)
- Spatial transcriptomics analysis (Visium)
- CODEX spatial proteomics processing
- decoupleR transcription factor activity inference
- Mediation analysis
- Virtual knockout (GENIE3/bulkKnk) and sRACIPE simulations
- Molecular docking analysis
- Exploratory clinical boundary analyses

## Data Availability
Due to size limits, raw datasets (e.g., TCGA-OV, GSE154600, GSE288483) are not included in this repository. Please download them from their respective public databases (UCSC Xena, GEO) as described in the Methods section of the manuscript. 

## Requirements
R version 4.5.3 or higher. All required packages (e.g., `survival`, `Seurat`, `WGCNA`, `decoupleR`, `mediation`, `GENIE3`, `bulkKnk`) are loaded at the beginning of the script.
