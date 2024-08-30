#Reference paper: https://www.science.org/doi/10.1126/sciadv.adg3754

setwd("/scratch/csierra/multiome_osk")

library(Seurat)
library(dplyr)
library(Signac)
library(irlba)
library(argparse)
library(future)

#Enable parallelization
plan("multicore", workers = 10)
plan()
options(future.globals.maxSize = 320 * 1024 ^ 3)

#Load object with peaks called per cell type and condition
load(file="/scratch/csierra/multiome_osk/objects/3_ATAC_2_peaks_geneact.rds")
DefaultAssay(obj_atac)<-"peaks_celltypes_group"

#Define functions
################################################################################
## Differential expression by cluster
################################################################################

DAbyCluster <- function(clusters,seuData, Group1, Group2) {
  ## Define comparisons
  obj_atac.ident <- subset(x = seuData, idents = clusters)
  Idents(object = obj_atac.ident) <- obj_atac.ident$Group
  ##DA analysis  
  FindMarkers(obj_atac.ident,
              ident.1 = Group1,
              ident.2 = Group2,
              logfc.threshold = 0,
              min.pct=0.05,
              test.use="LR",
              latent.vars=c("Batch", "nFeature_ATAC"))
}

#1. AD_GFP v CON_GFP
##Cell types
################################################################################
## Run DE for each cell type
################################################################################
## Define clusters to be compared
# Idents(obj_atac) <- obj_atac$celltype
# clusters <- levels(Idents(obj_atac))
#  
# wil <- lapply(clusters,DAbyCluster,obj_atac, Group1='AD_GFP', Group2='CON_GFP')
# names(wil)  <- clusters
#  
# save(wil,file="/scratch/csierra/multiome_osk/Analyses/ATAC/1_AD_GFPvCON_GFP/AD_GFPvCON_GFP_celltype.rds")
#  
# ##Engram cells
# ## Filter infected cells
# idx <-which(obj_atac$infected == 'Pos')
# obj_atac <- obj_atac[,idx]
# # # ################################################################################
# # # # ## Run DE for each cell type
# # # # ################################################################################
# ## Define clusters to be compared
# Idents(obj_atac) <- obj_atac$celltype
# clusters <- levels(Idents(obj_atac))
# 
# wil <- lapply(clusters,DAbyCluster,obj_atac, "AD_GFP", "CON_GFP")
# names(wil)  <- clusters
# names(wil) <- gsub("/","-",names(wil))
# 
# save(wil,file="/scratch/csierra/multiome_osk/Analyses/ATAC/1_AD_GFPvCON_GFP/AD_GFPvCON_GFP_engrams.rds")
# 

# #2. AD_OSK v AD_GFP 
# ##Cell types
# ################################################################################
# ## Run DE for each cell type
# ################################################################################
## Define clusters to be compared
# Idents(obj_atac) <- obj_atac$celltype
# clusters <- levels(Idents(obj_atac))
#  
# wil <- lapply(clusters,DAbyCluster,obj_atac, "AD_OSK", "AD_GFP")
# names(wil)  <- clusters
#   
# save(wil,file="/scratch/csierra/multiome_osk/Analyses/ATAC/2_AD_OSKvAD_GFP/AD_OSKvAD_GFP_celltype.rds")
#   
#  ##Engram cells
#  ## Filter infected cells
# idx <-which(obj_atac$infected == 'Pos')
# obj_atac_pos <- obj_atac[,idx]
# # # # ################################################################################
# # # # ## Run DE for each cell type
# # # # ################################################################################
# # ## Define clusters to be compared
# Idents(obj_atac_pos) <- obj_atac_pos$celltype
# clusters <- levels(Idents(obj_atac_pos))
# 
# wil <- lapply(clusters,DAbyCluster,obj_atac_pos, "AD_OSK", "AD_GFP")
# names(wil)  <- clusters
# names(wil) <- gsub("/","-",names(wil))
# #
# save(wil,file="/scratch/csierra/multiome_osk/Analyses/ATAC/2_AD_OSKvAD_GFP/AD_OSKvAD_GFP_engrams.rds")
# #
# #3. AD_OSK v CON_GFP
# ##Cell types
# ################################################################################
# ## Run DE for each cell type
# ################################################################################
# ## Define clusters to be compared
Idents(obj_atac) <- obj_atac$celltype
clusters <- levels(Idents(obj_atac))

wil <- lapply(clusters,DAbyCluster,obj_atac, "AD_OSK", "CON_GFP")
names(wil)  <- clusters

save(wil,file="/scratch/csierra/multiome_osk/Analyses/ATAC/3_AD_OSKvCON_GFP/AD_OSKvCON_GFP_da_celltypes.rds")

##Engram cells
## Filter infected cells
idx <-which(obj_atac$infected == 'Pos')
obj_atac_pos <- obj_atac[,idx]
# # ################################################################################
# # ## Run DE for each cell type
# # ################################################################################
# # ## Define clusters to be compared
Idents(obj_atac_pos) <- obj_atac_pos$celltype
clusters <- levels(Idents(obj_atac_pos))

wil <- lapply(clusters,DAbyCluster,obj_atac_pos, "AD_OSK", "CON_GFP")
names(wil)  <- clusters
names(wil) <- gsub("/","-",names(wil))

save(wil,file="/scratch/csierra/multiome_osk/Analyses/ATAC/3_AD_OSKvCON_GFP/AD_OSKvCON_GFP_da_engrams.rds")
