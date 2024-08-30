library(tidyselect)
library(fastmap)
library(Seurat)
library(dplyr)
library(Signac)
library(irlba)
library(argparse)
library(stringr)
library(cowplot)
library(ChIPseeker)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(patchwork)
library(future)

#Enable parallelization
plan("multicore", workers = 10)
plan()
options(future.globals.maxSize = 320 * 1024 ^ 3)

#Load object with peaks called per cell type and condition
load(file="/scratch/csierra/multiome_osk/objects/3_ATAC_4_min_2.rds")
DefaultAssay(obj_atac_min)<-"peaks_celltypes_group"

###Motif enrichment
# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
obj_atac_min <- AddMotifs(
  object = obj_atac_min,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)

#Footprinting
obj_atac_min <- Footprint(
     object = obj_atac_min,
     motif.name = c("KLF4", "SOX2", "POU5F1","Pou5f1::Sox2"),
     genome = BSgenome.Mmusculus.UCSC.mm10,
     in.peaks=TRUE
)
  
#save object with footprinting done
save(obj_atac_min,file="/scratch/csierra/multiome_osk/objects/3_ATAC_5_min_footprint_3.rds")