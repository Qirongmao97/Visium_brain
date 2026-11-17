# Regional GSEA running-score curves for a selected Gene Ontology term
#
# Purpose
#   Use precomputed mouse brain MAST log fold changes to generate a running
#   GSEA enrichment curve for each brain region and genotype contrast. The
#   default target is GO:0006007 (glucose catabolic process). This script
#   does not fit MAST models or select genes by differential expression FDR.
#
# Inputs
#   data/ALL_MAST_results.RData, or a file specified by MAST_RESULTS_FILE.
#   The file must contain result_list with a named diag list, deg list, or
#   both. Keys use <region>_<group1>_vs_<group2>. The parser uses matched
#   hurdle and logFC rows from diag when possible, otherwise primerid and
#   coef from deg. WT is accepted as an alias for wildtype.
#   GO gene membership is read from org.Mm.eg.db using GOALL. Keep the
#   package version fixed if exact GO membership must be reproduced.
#
# Ranking
#   Rank genes by descending MAST coefficient, keeping the coefficient
#   with largest absolute magnitude for duplicate cleaned gene symbols.
#   GSEA uses a p-value cutoff of 1 so curves can be plotted even when
#   enrichment is not significant.
#
# Outputs
#   One 3 x 2 inch GSEA_Curve_<GO_ID>_<contrast>.pdf per available contrast,
#   with region-specific curves overlaid. Outputs go to
#   results/multi_region_mast/gsea_curves by default. Set GSEA_CURVE_OUTPUT_DIR
#   to override the output location.
#
# Usage from the repository root
#   Rscript regional_gsea_enrichment_curves.R
#
# Requirements
#   dplyr, data.table, ggplot2, clusterProfiler, AnnotationDbi,
#   org.Mm.eg.db, and enrichplot.

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(ggplot2)
})

# Settings --------------------------------------------------------------------

MAST_RESULTS_FILE <- Sys.getenv(
  "MAST_RESULTS_FILE",
  unset = file.path("data", "ALL_MAST_results.RData")
)
OUT_DIR <- Sys.getenv(
  "GSEA_CURVE_OUTPUT_DIR",
  unset = file.path("results", "multi_region_mast", "gsea_curves")
)

TARGET_CONTRASTS <- c("wildtype vs mdxbl6", "wildtype vs mdx4cv", "mdxbl6 vs mdx4cv")
REGIONS_TO_INCLUDE <- c("Isocortex", "Hippocampus", "Thalamus", "Amygdala", "Central", "Fiber")
TARGET_GO_ID <- "GO:0006007"
TARGET_GO_NAME <- "glucose catabolic process"
PLOT_WIDTH_IN <- 3
PLOT_HEIGHT_IN <- 2

region_color_palette <- c(
  Isocortex = "#E495A5", Hippocampus = "#86B875", Thalamus = "#7DB0DD",
  Amygdala = "#FFB347", Central = "#B39EB5", Fiber = "#CFCFCF"
)

# Load precomputed MAST results -----------------------------------------------

if (!file.exists(MAST_RESULTS_FILE)) {
  stop("MAST results file not found. Check MAST_RESULTS_FILE.", call. = FALSE)
}
mast_input <- new.env(parent = emptyenv())
load(MAST_RESULTS_FILE, envir = mast_input)
if (!exists("result_list", envir = mast_input, inherits = FALSE)) {
  stop("The MAST results file must contain an object named result_list.", call. = FALSE)
}
res_model <- mast_input[["result_list"]]
if (!is.list(res_model) || !any(c("diag", "deg") %in% names(res_model))) {
  stop("result_list must contain diag and/or deg results.", call. = FALSE)
}

# Parse region and contrast labels without truncating multiword regions.
region_label <- function(region) {
  region <- trimws(gsub("[_[:space:]./]+", " ", region))
  aliases <- c(
    "hippocampus formation" = "Hippocampus",
    "central subpallium" = "Central",
    "fiber tract" = "Fiber",
    "fiber tracts" = "Fiber"
  )
  j <- match(tolower(region), names(aliases))
  if (!is.na(j)) return(unname(aliases[[j]]))
  j <- match(tolower(region), tolower(REGIONS_TO_INCLUDE))
  if (is.na(j)) region else REGIONS_TO_INCLUDE[[j]]
}

group_label <- function(group) {
  group <- tolower(group)
  if (group == "wt") "wildtype" else group
}

get_diagnostic_coefficients <- function(diagnostic, group2) {
  required <- c("component", "primerid", "contrast", "Pr(>Chisq)", "coef", "ci.hi", "ci.lo")
  if (!is.data.frame(diagnostic) || !all(required %in% names(diagnostic))) {
    return(NULL)
  }
  dt <- as.data.table(diagnostic)
  contrast_labels <- unique(as.character(dt$contrast[!is.na(dt$contrast)]))
  expected <- paste0("Group", group2)
  if (expected %in% contrast_labels) {
    selected <- expected
  } else {
    candidates <- contrast_labels[endsWith(contrast_labels, group2)]
    if (length(candidates) == 0L) return(NULL)
    selected <- candidates[[1]]
  }

  hurdle <- dt[contrast == selected & component == "H", .(primerid)]
  effects <- dt[contrast == selected & component == "logFC", .(primerid, coef)]
  if (nrow(hurdle) == 0L || nrow(effects) == 0L) return(NULL)
  merge(hurdle, effects, by = "primerid", all = FALSE)
}

parse_mast_results <- function(res_model) {
  keys <- unique(c(names(res_model$diag), names(res_model$deg)))
  key_pattern <- "^(.+)_([^_[:space:]]+)[_[:space:]]+vs[_[:space:]]+([^_[:space:]]+)$"
  records <- regmatches(keys, regexec(key_pattern, keys, ignore.case = TRUE))

  out <- lapply(seq_along(keys), function(i) {
    parts <- records[[i]]
    if (length(parts) != 4L) return(NULL)

    k <- keys[[i]]
    region <- region_label(parts[2])
    comparison <- paste(group_label(parts[3]), "vs", group_label(parts[4]))
    if (!region %in% REGIONS_TO_INCLUDE || !comparison %in% TARGET_CONTRASTS) {
      return(NULL)
    }

    fc <- get_diagnostic_coefficients(res_model$diag[[k]], group2 = parts[4])
    if (is.null(fc) || nrow(fc) == 0L) {
      dg <- res_model$deg[[k]]
      if (!is.data.frame(dg) || !all(c("primerid", "coef") %in% names(dg))) return(NULL)
      fc <- as.data.table(dg)[, .(primerid, coef)]
    }
    if (nrow(fc) == 0L) return(NULL)
    fc[, .(key = k, Region = region, Contrast = comparison, primerid, coef)]
  })
  rbindlist(out)
}

message("Parsing MAST coefficients.")
mast_table <- parse_mast_results(res_model)
if (nrow(mast_table) == 0L) {
  stop("No MAST coefficients matched the configured regions and contrasts. Check the table keys.",
       call. = FALSE)
}
duplicate_groups <- mast_table[, .(n_keys = uniqueN(key)), by = .(Region, Contrast)][n_keys > 1L]
if (nrow(duplicate_groups) > 0L) {
  stop(
    "Multiple MAST keys map to the same region and contrast. Check region and genotype aliases.",
    call. = FALSE
  )
}

# Obtain the selected GO gene set ----------------------------------------------

go_genes <- suppressMessages(AnnotationDbi::select(
  org.Mm.eg.db::org.Mm.eg.db,
  keys = TARGET_GO_ID,
  columns = "SYMBOL",
  keytype = "GOALL"
))
go_symbols <- unique(go_genes$SYMBOL[!is.na(go_genes$SYMBOL)])
if (length(go_symbols) == 0L) {
  stop("No mouse gene symbols were found for ", TARGET_GO_ID, ".", call. = FALSE)
}
custom_t2g <- data.frame(term = TARGET_GO_NAME, gene = go_symbols)

# Run GSEA and collect curves -------------------------------------------------

all_curves <- list()

for (cc in TARGET_CONTRASTS) {
  for (reg in REGIONS_TO_INCLUDE) {
    region_data <- mast_table %>%
      filter(Contrast == cc, Region == reg, is.finite(coef), !is.na(primerid)) %>%
      mutate(clean_gene = sub("\\..*$", "", trimws(as.character(primerid)))) %>%
      filter(nzchar(clean_gene)) %>%
      arrange(desc(abs(coef))) %>%
      distinct(clean_gene, .keep_all = TRUE)

    if (nrow(region_data) < 10L) next

    gene_list <- sort(setNames(region_data$coef, region_data$clean_gene), decreasing = TRUE)
    if (length(intersect(names(gene_list), go_symbols)) < 5L) next

    set.seed(42)
    gsea_res <- clusterProfiler::GSEA(
      geneList = gene_list,
      TERM2GENE = custom_t2g,
      pvalueCutoff = 1.0,
      minGSSize = 5,
      maxGSSize = 1000,
      verbose = FALSE
    )

    if (is.null(gsea_res) || nrow(as.data.frame(gsea_res)) == 0L) next

    # Use the running-score data behind enrichplot's public plotting function.
    curve <- enrichplot::gseaplot(gsea_res, geneSetID = TARGET_GO_NAME, by = "runningScore")
    if (!all(c("x", "runningScore") %in% names(curve$data))) {
      stop("enrichplot did not return the expected running-score columns.", call. = FALSE)
    }
    all_curves[[paste(cc, reg)]] <- data.frame(
      x = curve$data$x, runningScore = curve$data$runningScore, Contrast = cc, Region = reg
    )
  }
}

if (length(all_curves) == 0L) {
  stop("No GSEA curves could be calculated. Check gene-set overlap and the input table keys.",
       call. = FALSE)
}
df_plot <- bind_rows(all_curves)
df_plot$Region <- factor(df_plot$Region, levels = REGIONS_TO_INCLUDE)

# Export one PDF per contrast -------------------------------------------------

if (!dir.exists(OUT_DIR)) {
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
}
if (!dir.exists(OUT_DIR)) {
  stop("Cannot create the output directory. Check GSEA_CURVE_OUTPUT_DIR.", call. = FALSE)
}

for (cc in unique(df_plot$Contrast)) {
  df_subset <- df_plot %>% filter(Contrast == cc)
  p_gsea <- ggplot(df_subset, aes(x = x, y = runningScore, color = Region)) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
    geom_line(linewidth = 0.6) +
    scale_color_manual(values = region_color_palette) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 7),
      plot.subtitle = element_text(hjust = 0.5, size = 6, color = "grey40"),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 5),
      axis.title = element_text(face = "bold", size = 6),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 6),
      legend.text = element_text(size = 5),
      legend.key.size = grid::unit(0.3, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      plot.margin = margin(2, 4, 2, 2)
    ) +
    labs(
      title = paste("Enrichment Score:", TARGET_GO_NAME),
      subtitle = paste("Contrast:", cc),
      x = "Rank in Ordered Dataset",
      y = "Enrichment score (ES)",
      color = "Region"
    )

  safe_cc <- gsub(" ", "_", cc, fixed = TRUE)
  safe_term <- gsub(":", "_", TARGET_GO_ID, fixed = TRUE)
  plot_path <- file.path(OUT_DIR, paste0("GSEA_Curve_", safe_term, "_", safe_cc, ".pdf"))
  ggsave(plot_path, plot = p_gsea, width = PLOT_WIDTH_IN, height = PLOT_HEIGHT_IN, units = "in")
  message("Saved ", plot_path, ".")
}
