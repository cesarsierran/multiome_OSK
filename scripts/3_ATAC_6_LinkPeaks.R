#Reference paper: https://www.science.org/doi/10.1126/sciadv.adg3754

setwd("/scratch/csierra/multiome_osk")

library(tidyselect, lib.loc = "/home/csierra/R/x86_64-pc-linux-gnu-library/4.1")
library(fastmap)
library(Seurat)
library(dplyr)
library(Signac)
library(irlba)
library(argparse)
library(stringr)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(ChIPseeker)
library(BSgenome.Mmusculus.UCSC.mm10)
library(future)
library(readxl)

#Enable parallelization
plan("multicore", workers = 10)
plan()
options(future.globals.maxSize = 320 * 1024 ^ 3)

#Load object with peaks called per condition
load(file="/scratch/csierra/multiome_osk/objects/3_ATAC_2_peaks_geneact.rds")
DefaultAssay(obj_atac) <- "RNA"
obj_atac <- JoinLayers(obj_atac)  

#Get list of genes from universe
uni<-read_excel("./Analyses/RNA/2_AD_OSKvAD_GFP/AD_OSKvAD_GFP_universe_exc.xlsx")
uni_genes<-uni$`rownames(wil.de$exc)`

#Linking peaks
DefaultAssay(obj_atac) <- "peaks_celltypes_group"

# # first compute the GC content for each peak
obj_atac <- RegionStats(obj_atac, genome = BSgenome.Mmusculus.UCSC.mm10)
 
# # link peaks to genes
obj_atac <- LinkPeaks(
   object = obj_atac,
   peak.assay = "peaks_celltypes_group",
   expression.assay = "RNA",
   genes.use = uni_genes,
   distance = 1.5e+05
 )
# 
# #Save object with linked peaks
 save(obj_atac,file="/scratch/csierra/multiome_osk/objects/3_ATAC_X_LinkedPeaks_ALL.rds")
# #load("/scratch/csierra/multiome_osk/objects/3_ATAC_X_LinkedPeaks.rds")
# 
# 
# 
#Get list of all links
# Create an empty dataframe
df <- tibble(index = integer(), char_string = character())

for (i in uni_genes) {
   
   peaks<-GetLinkedPeaks(obj_atac, i, assay = NULL, min.abs.score = 0.1)  
   
   for (peak in peaks) {
   # Create a new row as a tibble
     new_row <- tibble(x = i, char_string = peak)
     
     # Append the new row to the dataframe
     df <- bind_rows(df, new_row)
   }
 }
 
linked<-df[,2:3]
colnames(linked)<-c("query_region", "gene")
 
write.csv(linked,file="/scratch/csierra/multiome_osk/LinkedPeaks_ALL.csv")