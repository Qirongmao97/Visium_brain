# Mouse brain Visium data preparation
#
# Purpose
#   Load 12 Visium Seurat objects from two processing batches, assign genotype
#   labels, and create a merged Seurat object with joined assay layers.
#   Use batch 1 annotations to map fine brain-region labels to broader regions
#   in each individual sample.
#
# Inputs
#   The 12 RDS files specified in sample_manifest. Place them in data/visium_rds
#   relative to the repository root, or set VISIUM_RDS_DIR to their directory.
#   Batch 1 objects must contain both annotation and annotation_hlevel metadata.
#
# Outputs
#   In-memory objects only. Source this script from the repository root to
#   access merged, mapping_table, seurat_list_mapped, and the named samples.
#   No files are written.
#
# Example
#   source("visium_brain_preprocess.R")
#
# Requirements
#   Seurat v5 for JoinLayers. The original package imports are retained for
#   compatibility with subsequent analysis code.

# Load packages ---------------------------------------------------------------

library(SingleCellExperiment)
library(Seurat)
library(patchwork)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(limma)
library(muscat)
library(purrr)
library(scater)
library(EnhancedVolcano)
library(ComplexHeatmap)
library(openxlsx)
library(readxl)
library(tidyr)
library(ggVennDiagram)

# Load samples and assign genotype labels --------------------------------------

data_dir <- Sys.getenv("VISIUM_RDS_DIR", unset = file.path("data", "visium_rds"))

# Row order preserves the original merge order.
sample_manifest <- data.frame(
  sample_id = c(
    "wildtype_1", "wildtype_2", "wildtype_3", "wildtype_4",
    "mdxbl6_1", "mdxbl6_2", "mdxbl6_3",
    "mdx4cv_1", "mdx4cv_2", "mdx4cv_3",
    "mdx52_1", "mdx52_2"
  ),
  genotype = c(rep("wildtype", 4), rep("mdxbl6", 3), rep("mdx4cv", 3), rep("mdx52", 2)),
  batch = c(1L, 2L, 2L, 2L, 1L, 2L, 2L, 1L, 2L, 2L, 1L, 2L),
  stringsAsFactors = FALSE
)

sampleNames <- sample_manifest$sample_id
sample_manifest$filename <- ifelse(
  sample_manifest$batch == 1L,
  paste0(sampleNames, "_longshort.rds"),
  paste0(sampleNames, ".rds")
)

rds_paths <- file.path(data_dir, sample_manifest$filename)
missing_files <- sample_manifest$filename[!file.exists(rds_paths)]
if (length(missing_files) > 0L) {
  stop(
    "Missing input RDS files. Set VISIUM_RDS_DIR or check data/visium_rds.\n",
    paste(missing_files, collapse = "\n"),
    call. = FALSE
  )
}

seuratObjList <- purrr::map2(rds_paths, sample_manifest$genotype, function(path, genotype) {
  obj <- readRDS(path)
  obj$Group <- genotype
  obj
})

# Merge Seurat objects and join assay layers -----------------------------------

# As in the original workflow, the merge precedes annotation harmonization.
merged <- Reduce(function(x, y) merge(x, y), seuratObjList)
merged <- JoinLayers(merged)

# Harmonize broad brain-region labels in individual samples ------------------

reference_indices <- which(sample_manifest$batch == 1L)
mapping_table <- seuratObjList[reference_indices] %>%
  purrr::map(~ dplyr::select(.x@meta.data, annotation, annotation_hlevel)) %>%
  dplyr::bind_rows() %>%
  dplyr::distinct() %>%
  dplyr::arrange(annotation_hlevel, annotation)

mapping_vec <- setNames(mapping_table$annotation_hlevel, mapping_table$annotation)

add_hlevel <- function(obj, mapping_vec, fine_col = "annotation", high_col = "annotation_hlevel") {
  md <- obj@meta.data
  if (!fine_col %in% names(md)) {
    message("Skipping: no ", fine_col, " column")
    return(obj)
  }

  # Look up labels by name, including when the fine labels are factors.
  # Fine labels absent from the reference mapping receive NA.
  md[[high_col]] <- unname(mapping_vec[as.character(md[[fine_col]])])
  obj@meta.data <- md
  obj
}

seurat_list_mapped <- purrr::map(seuratObjList, ~ add_hlevel(.x, mapping_vec))
names(seurat_list_mapped) <- sampleNames

# Keep individual object names available to downstream analysis scripts.
for (nm in sampleNames) {
  assign(nm, seurat_list_mapped[[nm]])
}
