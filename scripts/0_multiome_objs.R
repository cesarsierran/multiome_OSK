#SBATCH -J ${file}
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cesar.sierra@epfl.ch
#SBATCH --export=ALL
#SBATCH --nodes=1 --ntasks-per-node=48
#SBATCH --signal=2
#SBATCH --no-requeue
#SBATCH --mem=320GB

#Reference paper: https://www.science.org/doi/10.1126/sciadv.adg3754

setwd("/scratch/csierra/multiome_osk")

library(Seurat)
library(dplyr)
library(Signac)
library(irlba)
library(argparse)
library(hdf5r)
library(ggplot2)
library(IRanges)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(ensembldb)

# create parser object
parser <- ArgumentParser()
parser$add_argument("--input",type="character", default="input")
args <- parser$parse_args()
file <- args$input

# load the RNA data
counts <- Read10X_h5(paste0("samples/",file,"/outs/filtered_feature_bc_matrix.h5"))
metadata<-read.csv(paste0("samples/",file,"/outs/per_barcode_metrics.csv"), header = TRUE, row.names = 1, stringsAsFactors = FALSE)


#ATAC data
fragpath <- paste0("samples/",file,"/outs/atac_fragments.tsv.gz")

# create Seurat objects containing the RNA adata
seu <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA",
  project = file,
  meta.data = metadata
)

#ATAC
# Only use peaks in standard chromosomes
atac_counts <- counts$Peaks
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]

# Get gene annotations for GRCm38 - This part is done locally because it was not working in cluster
# Still not working: maybe contact helpdesk
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
annotations <- renameSeqlevels(annotations, mapSeqlevels(seqlevels(annotations), "UCSC"))
genome(annotations) <- "GRCm38"

# create ATAC assay and add it to the object
seu[["ATAC"]] <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  fragments = fragpath,
  min.cells = 10,
  annotation = annotations
)

##Quality control

# ATAC - calculate parameters
DefaultAssay(seu) <- "ATAC"
seu <- NucleosomeSignal(seu)
seu <- TSSEnrichment(seu)
seu$pct_reads_in_peaks <- seu$atac_peak_region_fragments / seu$atac_fragments * 100
seu$blacklist_fraction <- FractionCountsInRegion(seu, assay = 'ATAC', regions = blacklist_hg38_unified)


## Peak calling using MACS2

peaks <- CallPeaks(seu, 
                   assay = 'ATAC',
                   macs2.path = '/home/csierra/miniconda3/envs/sc/bin/macs2')

peaks <- keepStandardChromosomes(peaks, pruning.mode = 'coarse')
peaks <- subsetByOverlaps(x = peaks,
                          ranges = blacklist_mm10, 
                          invert = TRUE)

# quantify counts in each peak
macs_count <- FeatureMatrix(fragments = Fragments(seu),
                            features = peaks,
                            cells = colnames(seu))

seu[['ATAC']] <- CreateChromatinAssay(
  counts = macs_count,
  sep = c(":", "-"),
  fragments = fragpath,
  min.cells = 10,
  annotation = annotations
)

## Save the sample-specific Seurat object 
save(seu, file = paste0(file, '_processed.rda'))

