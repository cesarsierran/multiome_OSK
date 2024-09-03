# multiome_OSK
Code used for the manuscript "Cognitive rejuvenation through partial reprogramming of engram cells" 

0_multiome_objs: create a Seurat objects from Cell Ranger ARC outputs of different samples.
1_merge: merge all samples into a single Seurat object
2_QC_labelling: RNA and ATAC QC filtering. Dimensionality reduction and clustering. Labelling of engram+ cells. Save RNA and ATAC objects.
3_ATAC_1_process: ATAC data processing, dimensionality reduction and clustering.
3_ATAC_2_peak_geneact: ATAC peak calling for each annotated cell type and group using MACS2
3_ATAC_3_1_DA: Differential Accessibility Analysis.
3_ATAC_3_2_DA: Differential Accessibility Analysis plotting and downstream analyses
3_ATAC_4_min_footprint: ATAC footprinting.
3_ATAC_5_min_motifs: ATAC motif enrichment analysis and plotting.
3_ATAC_6_LinkPeaks: Signac Peak linking to genes in the dataset
3_RNA_1_composition: cell types visualization and compositional analysis
3_RNA_2_DE: Differential expression analysis
3_RNA_3_GSEA: GSEA analysis
3_RNA_3_identity: Identity score comparison between groups
4_integration: Integration analyses, including RNA-ATAC correlations
scCODA: compositional analysis using ths scCODA package
