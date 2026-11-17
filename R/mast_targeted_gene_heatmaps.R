# Targeted mouse brain MAST heatmaps
#
# Purpose
#   Plot precomputed MAST coefficients for nitric oxide synthase genes,
#   astroglial water and potassium channels, and dystrophin-associated genes
#   across brain regions and genotype contrasts. No MAST models are fitted.
#
# Inputs
#   data/ALL_MAST_results.RData, or a file specified by MAST_RESULTS_FILE.
#   The file must contain result_list with a named deg list. Table keys must
#   identify a region and an ordered contrast, for example
#   Hippocampus_Formation_wildtype_vs_mdxbl6. Each table must provide gene
#   symbols, adjusted p-values, and coefficients.
#   The script reads region labels from the MAST table keys by default.
#   To use the original reference-sample region set, source it with
#   wildtype_2_filtered in memory or set MAST_HEATMAP_REGIONS_FILE to a text
#   file containing one region label per line.
#
# Outputs
#   Three SVG heatmaps, their source CSVs, and a CSV documenting unmatched
#   and ambiguous MAST keys. By default, outputs are written to
#   results/multi_region_mast/targeted_gene_heatmaps. Set MAST_HEATMAP_OUTPUT_DIR
#   to override this location. Missing coefficients appear as gray tiles.
#   The fill scale is capped at +/-0.25, but CSV coefficients are not capped.
#   An asterisk marks FDR < 0.05 and an absolute coefficient > 0.1.
#
# Usage from the repository root
#   Rscript mast_targeted_gene_heatmaps.R
#
# Requirements
#   ggplot2, dplyr, tidyr, scales, and svglite.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# Paths and analysis settings -------------------------------------------------

MAST_RESULTS_FILE <- Sys.getenv(
  "MAST_RESULTS_FILE",
  unset = file.path("data", "ALL_MAST_results.RData")
)
OUTPUT_DIR <- Sys.getenv(
  "MAST_HEATMAP_OUTPUT_DIR",
  unset = file.path("results", "multi_region_mast", "targeted_gene_heatmaps")
)
REGIONS_FILE <- Sys.getenv("MAST_HEATMAP_REGIONS_FILE", unset = "")

if (!file.exists(MAST_RESULTS_FILE)) {
  stop("MAST results file not found. Check MAST_RESULTS_FILE.", call. = FALSE)
}
mast_input <- new.env(parent = emptyenv())
load(MAST_RESULTS_FILE, envir = mast_input)
if (!exists("result_list", envir = mast_input, inherits = FALSE)) {
  stop("The MAST results file must contain an object named result_list.", call. = FALSE)
}
MAST_RESULT_OBJECT <- mast_input[["result_list"]]
if (!is.list(MAST_RESULT_OBJECT) || !is.list(MAST_RESULT_OBJECT$deg) ||
    length(MAST_RESULT_OBJECT$deg) == 0L) {
  stop("result_list$deg must be a non-empty, named list of MAST tables.", call. = FALSE)
}
DEG_LIST <- MAST_RESULT_OBJECT$deg
if (is.null(names(DEG_LIST)) || anyNA(names(DEG_LIST)) || any(!nzchar(names(DEG_LIST)))) {
  stop("MAST table keys must be non-empty.", call. = FALSE)
}

EXCLUDED_REGIONS <- c("Hindbrain", "Meninge", "Ventricle")
DESIRED_REGION_ORDER <- c(
  "Isocortex", "Hippocampus Formation", "Thalamus", "Amygdala",
  "Central subpallium", "Fiber tract", "Unknown"
)
DESIRED_ORDER <- c("WT vs mdxbl6", "WT vs mdx4cv", "mdxbl6 vs mdx4cv")
FDR_CUTOFF <- 0.05
COEF_CUTOFF <- 0.1
COLOR_LIMIT <- 0.25

# Preserve the original panel identifiers for the output filenames.
GENE_SETS <- list(
  Q1_NOS = c("Nos1", "Nos2", "Nos3"),
  Q5_AQP4_Kir41 = c("Aqp4", "Kcnj10"),
  Q6_DGC = c("Dmd", "Dag1", "Dtna", "Dtnb", "Snta1", "Sntb1", "Sntb2", "Utrn")
)

FIG_TITLES <- c(
  Q1_NOS = "Q1 - Nitric oxide synthase genes",
  Q5_AQP4_Kir41 = "Q5 - Astroglial water and potassium channels",
  Q6_DGC = "Q6 - Dystrophin-associated complex genes"
)

# Match ordered contrasts and full region labels -----------------------------

normalize_key <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("wild[_. -]*type", "wt", x)
  x <- gsub("mdx[_. -]*bl6", "mdxbl6", x)
  x <- gsub("mdx[_. -]*4cv", "mdx4cv", x)
  trimws(gsub("[^a-z0-9]+", " ", x))
}

build_key_index <- function(deg_names, comparisons) {
  key_norm <- normalize_key(deg_names)
  dplyr::bind_rows(lapply(comparisons, function(comparison) {
    # Support a contrast at either edge, including older keys without "vs".
    contrast <- gsub(" vs ", " (vs )?", normalize_key(comparison), fixed = TRUE)
    suffix <- paste0(" ", contrast, "$")
    prefix <- paste0("^", contrast, " ")
    suffix_hit <- grepl(suffix, key_norm)
    prefix_hit <- !suffix_hit & grepl(prefix, key_norm)
    hits <- which(suffix_hit | prefix_hit)
    region <- as.character(ifelse(
      suffix_hit[hits],
      sub(suffix, "", key_norm[hits]),
      sub(prefix, "", key_norm[hits])
    ))
    data.frame(
      index = hits, region = region, comparison = rep(comparison, length(hits)),
      key = deg_names[hits], stringsAsFactors = FALSE
    )
  }))
}

key_index <- build_key_index(names(DEG_LIST), DESIRED_ORDER)

if (nzchar(REGIONS_FILE)) {
  if (!file.exists(REGIONS_FILE)) {
    stop("Region labels file not found. Check MAST_HEATMAP_REGIONS_FILE.", call. = FALSE)
  }
  available_regions <- readLines(REGIONS_FILE, warn = FALSE)
} else if (exists("wildtype_2_filtered", inherits = TRUE)) {
  # Keep the original reference-sample region selection when available.
  available_regions <- wildtype_2_filtered$annotation_hlevel
  if (is.null(available_regions)) {
    stop("wildtype_2_filtered must contain annotation_hlevel metadata.", call. = FALSE)
  }
} else {
  region_keys <- unique(key_index$region)
  known_regions <- c(DESIRED_REGION_ORDER, EXCLUDED_REGIONS)
  idx <- match(region_keys, normalize_key(known_regions))
  available_regions <- tools::toTitleCase(region_keys)
  available_regions[!is.na(idx)] <- known_regions[idx[!is.na(idx)]]
}

available_regions <- unique(trimws(as.character(available_regions)))
available_regions <- available_regions[!is.na(available_regions) & nzchar(available_regions)]
target_candidates <- available_regions[!normalize_key(available_regions) %in% normalize_key(EXCLUDED_REGIONS)]
TARGET_REGIONS <- c(
  intersect(DESIRED_REGION_ORDER, target_candidates),
  setdiff(target_candidates, DESIRED_REGION_ORDER)
)
if (length(TARGET_REGIONS) == 0L) {
  stop("No brain regions remain. Check the MAST keys or provide MAST_HEATMAP_REGIONS_FILE.",
       call. = FALSE)
}
message("Regions: ", paste(TARGET_REGIONS, collapse = ", "))

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
}
if (!dir.exists(OUTPUT_DIR)) {
  stop("Cannot create the output directory. Check MAST_HEATMAP_OUTPUT_DIR.", call. = FALSE)
}

find_deg_key <- function(region, comparison, key_index) {
  hits <- key_index[key_index$region == normalize_key(region) &
                      key_index$comparison == comparison, , drop = FALSE]
  if (nrow(hits) == 1L) {
    return(list(index = hits$index, status = "matched", names = hits$key))
  }
  if (nrow(hits) == 0L) {
    return(list(index = NA_integer_, status = "no_match", names = ""))
  }
  list(index = NA_integer_, status = "ambiguous", names = paste(hits$key, collapse = " | "))
}

find_column <- function(df, candidates) {
  nm <- names(df)
  idx <- match(tolower(candidates), tolower(nm))
  idx <- idx[!is.na(idx)]
  if (length(idx)) nm[idx[1]] else NA_character_
}

extract_target_stats <- function(df, genes, key_name) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  gene_col <- find_column(df, c("primerid", "gene", "symbol", "rowname", "id"))
  fdr_col  <- find_column(df, c("fdr", "p_val_adj", "padj", "p.adj", "adjp", "p.adjust"))
  coef_col <- find_column(df, c("coef", "logfc", "avg_log2fc", "log2fc", "estimate"))
  if (is.na(gene_col)) {
    rn <- rownames(df)
    if (!is.null(rn) && !identical(rn, as.character(seq_len(nrow(df))))) {
      df$.gene_from_rownames <- rn
      gene_col <- ".gene_from_rownames"
    }
  }
  if (anyNA(c(gene_col, fdr_col, coef_col))) {
    stop("Missing required gene / FDR / coefficient column(s) in: ", key_name,
         ". Actual columns: ", paste(names(df), collapse = ", "))
  }

  comp_col <- find_column(df, "component")
  if (!is.na(comp_col)) {
    component <- as.character(df[[comp_col]])
    # Retain the component preference used in the original analysis.
    if (any(component == "H", na.rm = TRUE)) {
      df <- df[component %in% "H", , drop = FALSE]
    } else if (any(component == "logFC", na.rm = TRUE)) {
      df <- df[component %in% "logFC", , drop = FALSE]
    } else {
      stop("No H or logFC component found in: ", key_name,
           ". Verify the MAST output before choosing C/D coefficients.")
    }
  }

  out <- data.frame(
    Gene = as.character(df[[gene_col]]),
    Coefficient = suppressWarnings(as.numeric(as.character(df[[coef_col]]))),
    FDR = suppressWarnings(as.numeric(as.character(df[[fdr_col]]))),
    stringsAsFactors = FALSE
  )
  out <- out[out$Gene %in% genes, , drop = FALSE]
  if (anyDuplicated(out$Gene)) {
    stop("Duplicate genes remain in ", key_name,
         ". Check the MAST component and gene columns; do not collapse arbitrarily.")
  }
  out
}

# Extract coefficients and write matching diagnostics ------------------------
ALL_TARGET_GENES <- unique(unlist(GENE_SETS, use.names = FALSE))
extracted <- list()
audit <- list()
for (comparison in DESIRED_ORDER) {
  for (region in TARGET_REGIONS) {
    match_info <- find_deg_key(region, comparison, key_index)
    audit[[length(audit) + 1L]] <- data.frame(
      Region = region,
      Comparison = comparison,
      Status = match_info$status,
      DEG_table = match_info$names,
      stringsAsFactors = FALSE
    )
    if (match_info$status != "matched") next
    idx <- match_info$index
    piece <- extract_target_stats(DEG_LIST[[idx]], ALL_TARGET_GENES, names(DEG_LIST)[idx])
    if (!nrow(piece)) next
    piece$Region <- region
    piece$Comparison <- comparison
    extracted[[length(extracted) + 1L]] <- piece
  }
}

audit_df <- dplyr::bind_rows(audit)
write.csv(audit_df, file.path(OUTPUT_DIR, "00_DEG_table_matching_diagnostics.csv"), row.names = FALSE)
if (any(audit_df$Status != "matched")) {
  warning("Some DEG tables were unmatched or ambiguous. Review 00_DEG_table_matching_diagnostics.csv; missing tiles will be gray.")
}

stats_df <- dplyr::bind_rows(extracted)
if (!nrow(stats_df)) {
  stop("No target genes recovered from MAST DEG tables. Inspect the names of the DEG tables and the matching diagnostics CSV.")
}

# Export gene-set heatmaps and source tables --------------------------------
plot_set <- function(set_name, genes) {
  grid <- tidyr::expand_grid(
    Gene = genes,
    Region = TARGET_REGIONS,
    Comparison = DESIRED_ORDER
  )
  # Preserve unavailable/not-tested combinations as NA (gray tiles).
  plot_df <- dplyr::left_join(grid, stats_df, by = c("Gene", "Region", "Comparison"))
  plot_df$Significance <- ifelse(
    !is.na(plot_df$FDR) & !is.na(plot_df$Coefficient) &
      plot_df$FDR < FDR_CUTOFF & abs(plot_df$Coefficient) > COEF_CUTOFF,
    "*", ""
  )
  plot_df$Comparison <- factor(plot_df$Comparison, levels = DESIRED_ORDER)
  plot_df$Region <- factor(plot_df$Region, levels = TARGET_REGIONS)
  plot_df$Gene <- factor(plot_df$Gene, levels = rev(genes))

  # The actual MAST coefficient remains unchanged in the CSV.
  # Only values outside the color scale are visually capped.
  p <- ggplot(plot_df, aes(x = Region, y = Gene, fill = Coefficient)) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_text(aes(label = Significance), size = 2.8, color = "black") +
    scale_fill_gradient2(
      low = "#313695", mid = "white", high = "#a50026", midpoint = 0,
      limits = c(-COLOR_LIMIT, COLOR_LIMIT), oob = scales::squish,
      breaks = c(-0.25, -0.125, 0, 0.125, 0.25),
      labels = c("-0.25", "-0.125", "0", "0.125", "0.25"),
      na.value = "grey85", name = "MAST coefficient",
      guide = guide_colorbar(title.position = "top", title.hjust = 0.5,
                             barwidth = grid::unit(5, "cm"),
                             barheight = grid::unit(0.3, "cm"))
    ) +
    facet_grid(~ Comparison, drop = FALSE) +
    labs(
      title = unname(FIG_TITLES[set_name]),
      x = NULL, y = NULL,
      caption = "Color capped at +/-0.25; gray = unavailable/not tested. * FDR < 0.05 and |coefficient| > 0.1."
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                                 face = "bold", size = 8),
      axis.text.y = element_text(face = "italic", size = 9),
      strip.background = element_rect(fill = "#f0f0f0", linewidth = 0.35),
      strip.text = element_text(face = "bold", size = 9),
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9),
      plot.caption = element_text(size = 7.5, hjust = 0),
      plot.margin = margin(8, 8, 8, 8)
    )

  svg_file <- file.path(OUTPUT_DIR, paste0(set_name, "_MAST_heatmap.svg"))
  csv_file <- file.path(OUTPUT_DIR, paste0(set_name, "_MAST_results.csv"))
  figure_height <- max(3.8, 2.6 + 0.32 * length(genes))
  ggsave(svg_file, p, width = 12.5, height = figure_height,
         units = "in", device = svglite::svglite, bg = "white")
  write.csv(plot_df, csv_file, row.names = FALSE)
  message("Saved: ", svg_file)
  invisible(p)
}

for (nm in names(GENE_SETS)) plot_set(nm, GENE_SETS[[nm]])
message("All three SVG heatmaps saved to: ", OUTPUT_DIR)
