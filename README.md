# Analysis Code for HGSOC CTRC-CELA3A-CELA3B Module

This repository contains the R script used for the multi-omics exploratory analysis in the manuscript.

## Requirements
* **R version**: 4.x (Tested on R 4.5.3).
* **Required R packages**: All necessary packages are loaded at the beginning of the script (e.g., `survival`, `Seurat`, `WGCNA`, `TwoSampleMR`, `oncoPredict`, `bulkKnk`, `GENIE3`, etc.). Please ensure these are installed before running.

## Data Preparation
The raw data are publicly available. Due to size limits, they are not provided in this repository. Please download the datasets from their respective public databases (as described in the manuscript) and place them in the **same working directory** as the `HGSOC_Analysis_Code.R` file.

Please ensure the downloaded files are renamed exactly as follows (as referenced in the code):
* `TCGA_OV_CBS.rds` (Preprocessed TCGA-OV expression matrix from Mendeley Data)
* `survival_OV_survival.txt` (TCGA-OV clinical survival data)
* `GSE49997_expr_hgsoc.rds` and `GSE49997_survival.csv` (Validation cohort)
* `CTRC_triple_score_pc1.rds` and `ssgsea_scores.rds` (Intermediate data generated in previous steps or provided in the `data/` folder)
* `predicted_IC50_GDSC2_all_drugs_optimized.rds` (Output of oncoPredict)

## Usage Instructions
1. Open the `HGSOC_Analysis_Code.R` script in RStudio.
2. Set your working directory to the folder containing the script and the downloaded data (`Session -> Set Working Directory -> To Source File Location`).
3. **Important**: Before running the Mendelian Randomization section (Part 2), you **MUST** configure your OpenGWAS JWT token in your local environment (e.g., `Sys.setenv(OPENGWAS_JWT = "your_token_here")`). The token is not provided in this repository for security reasons.
4. Run the script sequentially from top to bottom. Some steps (like GENIE3 network inference or bulkKnk virtual knockout) may take several hours depending on your computer's performance.
