#SBATCH -J ${file}
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cesar.sierra@epfl.ch
#SBATCH --export=ALL
#SBATCH --nodes=1 --ntasks-per-node=48
#SBATCH --signal=2
#SBATCH --no-requeue
#SBATCH --mem=320GB

setwd("/scratch/csierra/multiome_osk")

library(Seurat)
library(dplyr)
library(Signac)
library(irlba)
library(hdf5r)

######## Combine data for all the samples (n=12) to a single Seurat object ########
load("objects/CON_GFP1_processed.rda")
load("objects/CON_GFP2_processed.rda")
load("objects/CON_GFP3_processed.rda")
load("objects/AD_GFP1_processed.rda")
load("objects/AD_GFP2_processed.rda")
load("objects/AD_GFP3_processed.rda")
load("objects/AD_OSK1_processed.rda")
load("objects/AD_OSK2_processed.rda")
load("objects/AD_OSK3_processed.rda")
load("objects/AD_OSK4_processed.rda")


obj_all <- merge(AD_OSK4, y = c(CON_GFP1, CON_GFP2, CON_GFP3, AD_GFP1, AD_GFP2, AD_GFP3, AD_OSK1, AD_OSK2, AD_OSK3), 
                 add.cell.id = c("AD_OSK4","CON_GFP1","CON_GFP2","CON_GFP3", "AD_GFP1", "AD_GFP2","AD_GFP3", "AD_OSK1", "AD_OSK2", "AD_OSK3"), 
                 project = "brain_all")

save(obj_all, file = 'objects/1_obj_all.rda')
