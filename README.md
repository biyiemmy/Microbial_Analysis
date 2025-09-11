# RNA-Seq Analysis of Staphylococcus aureus under Sodium Propionate Treatment

## Project Overview
This repository contains the bioinformatics pipeline and analysis for a transcriptomics study comparing gene expression in *Staphylococcus aureus* under sodium propionate treatment vs. control conditions.

## Project Structure
Microbial_EnvRNASeq/
├── scripts/ # Analysis scripts (Python/R)
├── notes/ # Any project notes
├── sra_list.txt # SRA accessions used
├── summary.txt # Analysis summary
└── README.md


## Analysis Workflow
1. **Data Acquisition:** SRA data downloaded using `fasterq-dump`
2. **Quality Control:** `FastQC` & `fastp`
3. **Alignment:** `Bowtie2` to *S. aureus* reference genome
4. **Quantification:** `featureCounts`
5. **Differential Expression:** `limma-voom` in R
6. **Functional Annotation:** `Biopython` & NCBI Entrez queries

## Key Results
- Top DEGs: tRNA genes (tRNA-Ser, tRNA-Leu, tRNA-His) were significantly upregulated
- Full results: `results/NCBI_Final_Annotations.csv`

## How to Run
1. Install dependencies: `conda env create -f environment.yml`
2. Run RNA-Seq pipeline: `bash scripts/run_pipeline.sh`
3. Generate figures: `Rscript scripts/differential_expression.R`