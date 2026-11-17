# R analysis scripts

These scripts prepare and filter annotated mouse brain Visium sections, summarize precomputed regional differential expression results, and generate manuscript and supporting figures. See the [repository overview](../README.md) for the publication and dataset.

Run shell commands from the repository root so paths such as `data/` resolve correctly. The preprocessing, QC, and correlation steps use the **same R session**. The MAST-based figure scripts can be run as separate `Rscript` commands. The filenames below match the scripts in this folder.

## R environment

Use Seurat v5 for the Seurat workflows. The scripts list or import their required packages. To install the dependencies for all scripts in this folder:

~~~r
install.packages(c(
  "Seurat", "BiocManager", "dplyr", "tidyr", "purrr", "tibble", "stringr",
  "data.table", "ggplot2", "ggrepel", "patchwork", "gridExtra",
  "RColorBrewer", "scales", "openxlsx", "readxl", "ggVennDiagram", "UpSetR",
  "svglite"
))

BiocManager::install(c(
  "SingleCellExperiment", "limma", "muscat", "scater", "EnhancedVolcano",
  "ComplexHeatmap", "AnnotationDbi", "org.Mm.eg.db", "clusterProfiler",
  "enrichplot", "decoupleR"
))
~~~

Some plots use Cairo PDF. Check `capabilities("cairo")` if PDF export fails. For a reproducible rerun, keep a record of `sessionInfo()`, the mouse annotation database version, and the CollecTRI network used.

## Seurat preprocessing and spot QC

Start with 12 annotated Seurat RDS objects in `data/visium_rds/` or set `VISIUM_RDS_DIR`. The four reference samples use the filename pattern `<sample>_1_longshort.rds`, where `<sample>` is `wildtype`, `mdxbl6`, `mdx4cv`, or `mdx52`. The remaining eight use `<sample>_<number>.rds`. See `sample_manifest` in [visium_brain_preprocess.R](visium_brain_preprocess.R) for the exact names.

~~~r
source("R/visium_brain_preprocess.R")
source("R/visium_brain_spot_qc.R")
~~~

| Script | What it does | Output |
| --- | --- | --- |
| [visium_brain_preprocess.R](visium_brain_preprocess%281%29.R) | Loads 12 samples and uses reference labels to map fine annotations to broad brain regions in the individual objects. | Named sample objects and `seurat_list_mapped` in memory. The separate `merged` object is created before broad labels are harmonized. |
| [visium_brain_spot_qc.R](visium_brain_spot_qc%281%29.R) | Applies hippocampal and non-hippocampal QC thresholds, then filters features and log-normalizes the Spatial assay. The feature filter does not prune ONT features. | Named `*_filtered` objects, `seuratObjList_filtered`, and `results/visium_qc/summary_table.xlsx`. |

If you rerun Python label transfer, use `C_scANVI` from its per-sample CSVs when preparing the `annotation` metadata for the query Seurat objects. Keep the original reference annotations. Match by sample and original barcode, not by row order. The scripts here do not perform the H5AD-to-Seurat conversion or import the annotation CSVs automatically.

## Spatial and module-score analyses

| Script | Required input | Output |
| --- | --- | --- |
| [visium_myelin_mito_correlation.R](visium_myelin_mito_correlation%281%29.R) | Ten filtered WT, mdxbl6, and mdx4cv objects plus the original gene vectors `complex_I` through `complex_V`. Source it after QC once those vectors are defined. | `visium_merged`, `cor_summary`, and `p_cor_heatmap` in memory. Correlations are descriptive and calculated across spots within each genotype and region. |
| [ETC_spatial_pathway_scores.R](ETC_spatial_pathway_scores.R) | A merged, normalized `visium_merged` Seurat object with spatial images and a joined Spatial assay. Source the script with the object in memory, or set `VISIUM_MERGED_RDS` to use a saved object. | PDF maps for respiratory-chain complexes I–V and GO:0006007 in `results/visium_spatial_modules/`. The display range is fixed at 0–1.2 and does not truncate the stored values. |

## Regional MAST results and figures

Unless otherwise noted, the scripts below read `data/ALL_MAST_results.RData`. This file must contain `result_list` with named regional and contrast-specific `$deg` and/or `$diag` tables. The exact requirements are in each script header. Set `MAST_RESULTS_FILE` if the file is elsewhere. These scripts consume existing model results rather than refitting the models.

| Script | Analysis and outputs | Default location |
| --- | --- | --- |
| [multi_region_mast_analysis.R](multi_region_mast_analysis%281%29.R) | Regional volcano plots, heatmaps and UpSet plots for shared genes, GO enrichment, and GSEA ridgeplots. | `results/multi_region_mast/` |
| [ETC_complex_heatmap.R](ETC_complex_heatmap.R) | Mean and gene-level coefficients for respiratory-chain complexes I–V and GO:0006007. Exports PDFs and source tables. | `results/mast_pathway_panels/` |
| [mast_targeted_gene_heatmaps.R](mast_targeted_gene_heatmaps%281%29.R) | Targeted heatmaps for nitric oxide synthases, astroglial channels, and dystrophin-associated genes. Exports SVGs, source CSVs, and a table-matching diagnostic. | `results/multi_region_mast/targeted_gene_heatmaps/` |
| [mast_tf_activity_analysis.R](mast_tf_activity_analysis%281%29.R) | Uses decoupleR with mouse CollecTRI for transcription-factor activity inference. Exports regional plots and tables and, where results are available, a clustered status heatmap and GO enrichment. | `results/multi_region_mast/tf_activity/` |
| [regional_gsea_enrichment_curves.R](regional_gsea_enrichment_curves%281%29.R) | Overlays brain region-specific running GSEA enrichment curves for GO:0006007 in each genotype contrast. | `results/multi_region_mast/gsea_curves/` |
| [cross_contrast_gsea_dotplot.R](cross_contrast_gsea_dotplot%281%29.R) | Reads precomputed `GSEA_Master_MultiDB_<contrast>.csv` inputs and plots selected GO Biological Process terms across contrasts and regions. The other scripts listed here do not create those CSVs. | `results/multi_region_mast/` |

For example, to regenerate the targeted heatmaps with a stored MAST results file:

~~~bash
MAST_RESULTS_FILE=data/ALL_MAST_results.RData \
  Rscript "R/mast_targeted_gene_heatmaps.R"
~~~

## Optional path settings

| Environment variable | Purpose |
| --- | --- |
| `VISIUM_RDS_DIR`, `VISIUM_QC_OUTPUT_DIR` | Seurat RDS input directory and QC workbook output directory. |
| `VISIUM_MERGED_RDS`, `VISIUM_SPATIAL_OUTPUT_DIR` | Saved merged Seurat input and spatial module plot destination. |
| `MAST_RESULTS_FILE`, `MAST_OUTPUT_DIR` | Precomputed MAST input and main multi-region output directory. |
| `MAST_PATHWAY_OUTPUT_DIR`, `MAST_HEATMAP_OUTPUT_DIR` | Output directories for pathway panels and targeted heatmaps. |
| `MAST_HEATMAP_REGIONS_FILE` | Optional text file with one region per line to reproduce the reference sample's region selection. |
| `TF_OUTPUT_DIR`, `TF_NETWORK_CSV` | TF analysis destination and optional saved mouse network. |
| `GSEA_CURVE_OUTPUT_DIR` | Running enrichment curve destination. |
| `GSEA_INPUT_DIR`, `GSEA_OUTPUT_DIR` | Cross-contrast GSEA CSV input and dot plot output directories. |
