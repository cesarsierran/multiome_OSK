# Multiome OSK
Code used for the manuscript "Cognitive rejuvenation through partial reprogramming of engram cells"

- ```0_multiome_objs.R```: create a Seurat objects from Cell Ranger ARC outputs of different samples.
- ```1_merge.R```: merge all samples into a single Seurat object
- ```2_QC_labelling.Rmd```: RNA and ATAC QC filtering. Dimensionality reduction and clustering. Labelling of engram+ cells. Save RNA and ATAC objects.
- ```3_ATAC_1_process.Rmd```: ATAC data processing, dimensionality reduction and clustering.
- ```3_ATAC_2_peak.R```: ATAC peak calling for each annotated cell type and group using MACS2
- ```3_ATAC_3_1_DA.R```: Differential Accessibility Analysis.
- ```3_ATAC_3_2_DA.Rmd```: Differential Accessibility Analysis plotting and downstream analyses
- ```3_ATAC_4_min_footprint.R```: ATAC footprinting.
- ```3_ATAC_5_min_motifs.Rmd```: ATAC motif enrichment analysis and plotting.
- ```3_ATAC_6_LinkPeaks.R```: Signac Peak linking to genes in the dataset
- ```3_RNA_1_composition.Rmd```: cell types visualization and compositional analysis
- ```3_RNA_2_DE.Rmd```: Differential expression analysis
- ```3_RNA_3_GSEA.Rmd```: Gene Set Enrichment Analysis
- ```3_RNA_3_identity.Rmd```: Identity score comparison between groups
- ```3_RNA_4_pseudotime.Rmd```: Monocle pseudotime analysis based on engram and AD DEGs
- ```3_RNA_5_CHEA3.Rmd```: TF enrichment analysis of downregulated genes
- ```4_integration.Rmd```: Integration analyses, including RNA-ATAC correlations
- ```scCODA.py```: compositional analysis using ths scCODA package

The transcriptional landscapes from Figure S7 can be visualized [here]([https://pages.github.com/](https://cesarsierra.shinyapps.io/multiome_shiny/)).
