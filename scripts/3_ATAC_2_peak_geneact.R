library(Seurat)
library(dplyr)
library(Signac)
library(irlba)
library(argparse)
library(ggplot2)
library(IRanges)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(ensembldb)

#Call peaks per cell type and condition
#Load object
load(file="/scratch/csierra/multiome_osk/objects/3_ATAC_1_obj_all.rds")

## Re-call peaks for each annotated cell type and condition using MACS2
peaks <- CallPeaks(obj_atac, assay = 'ATAC',
                    macs2.path = '/home/csierra/miniconda3/envs/sc/bin/macs2',
                    group.by = c('celltype',"Group"),
                    outdir = 'MACS2_output',
                    fragment.tempdir = 'MACS2_output',
                    cleanup = FALSE)

 save(peaks,file="/scratch/csierra/multiome_osk/objects/peaks_celltype_group.rds")
 # load(file="/scratch/csierra/multiome_osk/objects/peaks.rds")

# Remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = 'coarse')
peaks <- subsetByOverlaps(x = peaks,
                           ranges = blacklist_mm10,
                           invert = TRUE)

# Quantify counts in each peak
 DefaultAssay(obj_atac) <- 'ATAC'
 frags <- Fragments(obj_atac)

 macs_count <- FeatureMatrix(
   fragments = frags,
   features = peaks,
   cells = colnames(obj_atac)
 )
save(macs_count,file="/scratch/csierra/multiome_osk/objects/macs_count_celltypes_group.rds")

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
annotations <- renameSeqlevels(annotations, mapSeqlevels(seqlevels(annotations), "UCSC"))
genome(annotations) <- "GRCm38"

obj_atac[['peaks_celltypes_group']] <- CreateChromatinAssay(
   counts = macs_count,
   sep = c(":", "-"),
   fragments = frags,
   annotation = annotations
 )

DefaultAssay(obj_atac) <- 'peaks_celltypes_group'
obj_atac <- RunTFIDF(obj_atac, method = 3)
obj_atac <- FindTopFeatures(obj_atac, min.cutoff = 'q75')


#save object
save(obj_atac,file="/scratch/csierra/multiome_osk/objects/3_ATAC_2_peaks_geneact.rds")
