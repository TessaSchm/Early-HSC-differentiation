# Early-HSC-differentiation
The code in this repository is relevant to the combined analysis of the data associated with the following studies:

[00_CSV files](folder100_CSV_files/) 
- This folder contains the csv source files used in this study. The accompanying excel sheet provides comprehensive information and descriptions regarding the utilized BM samples.

[01_Data processing and transformation](01_Data_processing_and_transformation) 
- This code contains: 
01. Generating seurat objects from csv files. 
02. Automated annotation of each cartridge using Simon Haas WTA reference.
03. Combination of all cartridges into one seurat object. 
04. Batch effect correction for mrna analyses corrected by CCA method for combined data and by age group.  

[02_Data analysis and visualization](02_Data_analysis_and_visualization) 
- This code contains:
05. Analysis and visualization of CD34+ HSPCs. 
06. Re-clustering of immature cell subset for in-depth analysis. 
07. Annotation and visualization of immature clusters.
08. Manual gating to confirm cluster annotation.
09. Volcanoplots comparing protein (abseq) expression of HSC1vsHSC2 subset.
10. Pseudobulk analysis to compare different age groups. 
11. Heatmaps for visualizing differently expressed gene subsets.
12. Violinplots comparing different age groups within the HSC-1 cluster. 
13. Slingshot and tradeseq analysis to gain insights into the cell fate decisions. 
