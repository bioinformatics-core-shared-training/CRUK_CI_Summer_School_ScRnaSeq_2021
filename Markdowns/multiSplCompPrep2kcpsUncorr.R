# short script used to:
# - read in the 'uncorrected' SCE object made with all cells in the data sets
# - downsample to 2000 cells per sample, choosing cells randomly 
# - TODO: maybe make new SCE from counts to start from scratch and avoid dragging incorrect QC metrics.

setwd("/home/ubuntu/Course_Materials/scRNAseq/Markdowns")

library(scater)
library(scran)

# load the uncorrected data
uncorrected <- readRDS("../Robjects/caron_sce_nz_postDeconv_allCells_dsi_PBMMC_ETV6-RUNX1_uncorr.Rds")
uncorrected$SampleName <- uncorrected$Sample.Name2
uncorrected$SampleGroup <- uncorrected$source_name

# downsample to 2000 cells per sample

#--- feature-annotation ---#
rownames(uncorrected) <- uniquifyFeatureNames(
  rowData(uncorrected)$ensembl_gene_id,
  rowData(uncorrected)$Symbol)

#--- quality-control ---#
# done already

# downsample
sampleNameToGet <- colData(uncorrected)[,c("SampleName", "SampleGroup")] %>%
  unique() %>%
  data.frame %>%
  pull(SampleName)
bcToKeep <- lapply(sampleNameToGet, function(x){ 
  colData(uncorrected) %>%
    data.frame %>%
    dplyr::filter(SampleName == x) %>%
    sample_n(., 1200) %>%
    pull(Barcode)
}) %>% unlist
colnames(uncorrected) <- uncorrected$Barcode
uncorrected <- uncorrected[,bcToKeep]

# only keep genes detected in at 20 cells (about 0.15% of 14,000 cells) 
newFeatQc <- perFeatureQCMetrics(uncorrected)
featToKeep <- newFeatQc$detected > 0.15
uncorrected <- uncorrected[featToKeep,]

#--- normalization ---#
set.seed(100) # clusters with PCA from irlba with approximation
library(BiocParallel)
# prepare 'cluster' for multicore processing
# use 7 workers
bpp <- MulticoreParam(7)
clust <- quickCluster(uncorrected, BPPARAM=bpp) # slow with all cells.
table(clust)
uncorrected <- computeSumFactors(uncorrected,
                                 cluster = clust,
                                 min.mean = 0.1,
                                 BPPARAM = bpp)
uncorrected <- logNormCounts(uncorrected) # adds logcounts
# check list of assays stored:
print(assayNames(uncorrected))

#--- variance-modelling ---#
dec.uncorr <- modelGeneVar(uncorrected, block=uncorrected$SampleName)
chosen.hvgs <- dec.uncorr$bio > 0

# write SCE object to file:
#saveRDS(uncorrected, file = "../Robjects/caron_sce_nz_postDeconv_2kcps_dsi_PBMMC_ETV6-RUNX1_uncorr.Rds")
saveRDS(uncorrected, file = "../Robjects/caron_sce_nz_postDeconv_1p2kcps_dsi_PBMMC_ETV6-RUNX1_uncorr.Rds")
