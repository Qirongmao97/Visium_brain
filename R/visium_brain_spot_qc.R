# Region-specific spot QC and feature filtering for mouse brain Visium
#
# Purpose
#   Apply hippocampus-specific spot QC to 12 annotated Visium samples, then
#   remove ribosomal, hemoglobin, and mitochondrial features only from the
#   Spatial assay and log-normalize that assay. ONT features are not pruned.
#
# Inputs
#   In-memory Seurat objects created by visium_brain_preprocess.R (or equivalent):
#   wildtype_1–4, mdxbl6_1–3, mdx4cv_1–3, and mdx52_1–2.
#   Each object needs Spatial counts, nCount_Spatial, nFeature_Spatial, and
#   annotation_hlevel. Spots not labeled as Hippocampus Formation, including
#   unannotated spots, use the non-hippocampal thresholds.
#
# Outputs
#   summary_table.xlsx reports spot retention and overlapping failure counts
#   for each criterion. It is written to results/visium_qc by default. Set
#   VISIUM_QC_OUTPUT_DIR to use another location.
#   Filtered, normalized Seurat objects are returned in memory as
#   seuratObjList_filtered and <sample>_filtered. No RDS files are written.
#
# Usage
#   With the annotated objects in memory, run from the repository root:
#   source("visium_brain_spot_qc.R")
#
# Requirements
#   Seurat v5 and openxlsx.

suppressPackageStartupMessages(library(Seurat))

# Sample IDs and thresholds ----------------------------------------------------

sampleNames <- c(
  paste0("wildtype_", seq_len(4)),
  paste0("mdxbl6_", seq_len(3)),
  paste0("mdx4cv_", seq_len(3)),
  paste0("mdx52_", seq_len(2))
)
script_env <- environment()
missing_samples <- sampleNames[
  !vapply(sampleNames, exists, logical(1), envir = script_env, inherits = TRUE)
]
if (length(missing_samples) > 0L) {
  stop("Missing annotated sample objects: ", paste(missing_samples, collapse = ", "),
       call. = FALSE)
}
seuratObjList <- mget(sampleNames, envir = script_env, inherits = TRUE)

max_hbb <- 5
max_counts <- 35000
max_features <- 6500

max_mt_hipp <- 40
min_counts_hipp <- 200
min_features_hipp <- 100

max_mt_outside <- 25
min_counts_outside <- 500
min_features_outside <- 300

OUT_DIR <- Sys.getenv(
  "VISIUM_QC_OUTPUT_DIR",
  unset = file.path("results", "visium_qc")
)

# Apply spot-level QC ----------------------------------------------------------

filter_spatial_object <- function(obj, obj_name = NULL, region_col = "annotation_hlevel",
                                  verbose = TRUE) {
  if (!inherits(obj, "Seurat") || !"Spatial" %in% Assays(obj)) {
    stop("A Seurat object with a Spatial assay is required: ", obj_name, call. = FALSE)
  }
  if (is.null(obj_name)) obj_name <- Project(obj)

  md <- obj[[]]
  required_cols <- c(region_col, "nCount_Spatial", "nFeature_Spatial")
  missing_cols <- setdiff(required_cols, names(md))
  if (length(missing_cols) > 0L) {
    stop(obj_name, " is missing metadata: ", paste(missing_cols, collapse = ", "),
         call. = FALSE)
  }

  obj <- PercentageFeatureSet(obj, pattern = "^mt-", col.name = "percent.mt", assay = "Spatial")
  obj <- PercentageFeatureSet(obj, pattern = "^Hb[ab]", col.name = "percent.hbb", assay = "Spatial")
  md <- obj[[]]
  pct_mt <- md$percent.mt
  pct_hbb <- md$percent.hbb
  region <- as.character(md[[region_col]])
  is_hipp <- !is.na(region) & grepl("Hippocampus Formation", region, ignore.case = TRUE)

  flags <- data.frame(
    fail_hbb = pct_hbb >= max_hbb,
    fail_maxcount = md$nCount_Spatial >= max_counts,
    fail_maxfeat = md$nFeature_Spatial >= max_features,
    fail_mt = pct_mt >= ifelse(is_hipp, max_mt_hipp, max_mt_outside),
    fail_mincount = md$nCount_Spatial <= ifelse(is_hipp, min_counts_hipp, min_counts_outside),
    fail_minfeat = md$nFeature_Spatial <= ifelse(is_hipp, min_features_hipp, min_features_outside),
    row.names = rownames(md)
  )
  # Match the original rule: missing measurements do not trigger a failure.
  flags[] <- lapply(flags, function(flag) replace(flag, is.na(flag), FALSE))
  flags$fail_any <- Reduce(`|`, flags)

  kept_spots <- sum(!flags$fail_any)
  if (kept_spots == 0L) {
    stop("No spots remain after QC for ", obj_name, ".", call. = FALSE)
  }
  total_spots <- ncol(obj)
  summary_row <- data.frame(
    Sample = obj_name,
    Total_Spots = total_spots,
    Kept_Spots = kept_spots,
    Kept_Percent = round(100 * kept_spots / total_spots, 1),
    Filtered_Spots = total_spots - kept_spots,
    HBB_Fails = sum(flags$fail_hbb),
    MaxCount_Fails = sum(flags$fail_maxcount),
    MinCount_Fails = sum(flags$fail_mincount),
    MaxFeature_Fails = sum(flags$fail_maxfeat),
    MinFeature_Fails = sum(flags$fail_minfeat),
    Hipp_MT_Fails = sum(flags$fail_mt & is_hipp),
    Outside_MT_Fails = sum(flags$fail_mt & !is_hipp)
  )

  qc_metadata <- data.frame(is_hipp = is_hipp, flags, row.names = rownames(md))
  obj <- AddMetaData(obj, metadata = qc_metadata)
  obj_filtered <- subset(obj, cells = rownames(md)[!flags$fail_any])

  if (verbose) {
    message(obj_name, ": retained ", kept_spots, "/", total_spots, " spots (",
            summary_row$Kept_Percent, "%).")
  }
  list(filtered = obj_filtered, summary = summary_row)
}

results <- Map(function(obj, nm) filter_spatial_object(obj, obj_name = nm),
               seuratObjList, sampleNames)
summary_table <- do.call(rbind, lapply(results, function(result) result$summary))
rownames(summary_table) <- NULL

if (!dir.exists(OUT_DIR)) {
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
}
if (!dir.exists(OUT_DIR)) {
  stop("Cannot create the output directory. Check VISIUM_QC_OUTPUT_DIR.", call. = FALSE)
}
summary_path <- file.path(OUT_DIR, "summary_table.xlsx")
openxlsx::write.xlsx(summary_table, summary_path, rowNames = FALSE)
message("Wrote QC summary to ", summary_path, ".")

# Remove excluded features from Spatial only and normalize ---------------------

seuratObjList_filtered <- lapply(results, function(result) result$filtered)
names(seuratObjList_filtered) <- sampleNames

for (nm in sampleNames) {
  obj <- seuratObjList_filtered[[nm]]
  spatial_assay <- obj[["Spatial"]]
  spatial_features <- rownames(spatial_assay)
  exclude_genes <- grep(
    "^Rps|^Rpl|^Hba|^Hbb|^mt-", spatial_features, ignore.case = TRUE, value = TRUE
  )
  spatial_keep <- setdiff(spatial_features, exclude_genes)
  if (length(spatial_keep) == 0L) {
    stop("No Spatial features remain after filtering ", nm, ".", call. = FALSE)
  }

  # Subset the assay rather than using a union of Spatial and ONT feature IDs.
  # This also removes a matching Spatial feature when its ID exists in ONT.
  if (length(exclude_genes) > 0L) {
    obj[["Spatial"]] <- subset(spatial_assay, features = spatial_keep)
  }
  obj <- NormalizeData(obj, assay = "Spatial", normalization.method = "LogNormalize")
  seuratObjList_filtered[[nm]] <- obj
  assign(paste0(nm, "_filtered"), obj, envir = script_env)
  message(nm, ": removed ", length(exclude_genes), " Spatial features.")
}
