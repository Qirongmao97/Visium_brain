# Regional transcription factor activity inference from mouse brain MAST results
#
# Purpose
#   Infer transcription factor (TF) activities separately for each region and
#   genotype contrast using decoupleR's univariate linear model (ULM) and the
#   mouse CollecTRI network. Compare WT vs mdxbl6 with WT vs mdx4cv in each
#   region, cluster TFs by their regional classification patterns, and test
#   the resulting TF groups for GO Biological Process enrichment.
#
# Inputs
#   data/ALL_MAST_results.RData, or a file specified by MAST_RESULTS_FILE.
#   The RData must contain a result_list with a named deg list. Use keys of
#   the form <region>_<group1>_vs_<group2>. Spaces around "vs" are supported.
#   WT and wildtype are accepted as wild-type group names.
#   Tables require primerid, coef, and FDR (or Pr(>Chisq) for BH correction).
#   By default the mouse CollecTRI network is retrieved by decoupleR. To reuse
#   a fixed network, set TF_NETWORK_CSV to a CSV with source, target,
#   and mor (or weight) columns.
#
# Analysis
#   The ULM input is sign(coef) * -log10(FDR), not a measured t statistic.
#   Zero FDR is represented by -log10(FDR) = 300. At least 50 input rows
#   per region/contrast and five overlapping targets per TF are required.
#   TF categories use nominal ULM p < 0.05. "Exclusive" means that only one
#   contrast met the threshold. It is not a direct test of activity differences
#   between the two mutant genotypes.
#
# Outputs
#   TF activity barplots, regional comparison scatterplots and status CSVs,
#   a clustered TF status heatmap, and a GO dotplot and results CSV when
#   enrichment is available. The default output is results/multi_region_mast/
#   tf_activity. Set TF_OUTPUT_DIR to choose a different directory.
#
# Usage from the repository root
#   Rscript mast_tf_activity_analysis.R
#
# Requirements
#   decoupleR, dplyr, ggplot2, data.table, tibble, RColorBrewer, ggrepel,
#   tidyr, stringr, clusterProfiler, org.Mm.eg.db, enrichplot, and Cairo PDF.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(data.table)
})

# Paths and settings ----------------------------------------------------------

MAST_RESULTS_FILE <- Sys.getenv(
  "MAST_RESULTS_FILE",
  unset = file.path("data", "ALL_MAST_results.RData")
)
OUT_DIR <- Sys.getenv(
  "TF_OUTPUT_DIR",
  unset = file.path("results", "multi_region_mast", "tf_activity")
)
TF_NETWORK_CSV <- Sys.getenv("TF_NETWORK_CSV", unset = "")

SCATTER_WIDTH <- 8
SCATTER_HEIGHT <- 8
HEATMAP_WIDTH <- 4.0
HEATMAP_HEIGHT <- 5.5
LEGEND_WRAP <- 18
N_TF_GROUPS <- 4L

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
    length(MAST_RESULT_OBJECT$deg) == 0L || is.null(names(MAST_RESULT_OBJECT$deg))) {
  stop("result_list must contain a non-empty, named deg list.", call. = FALSE)
}

if (!dir.exists(OUT_DIR)) {
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
}
if (!dir.exists(OUT_DIR)) {
  stop("Cannot create the output directory. Check TF_OUTPUT_DIR.", call. = FALSE)
}

# Parse MAST results ----------------------------------------------------------

build_mast_table <- function(deg_list) {
  keys <- names(deg_list)
  if (anyDuplicated(keys)) {
    stop("The MAST deg list contains duplicated contrast keys.", call. = FALSE)
  }
  # Anchor the contrast suffix so a region may contain more than one word.
  key_pattern <- "^(.+)_([^_[:space:]]+)[_[:space:]]+vs[_[:space:]]+([^_[:space:]]+)$"
  key_parts <- regmatches(keys, regexec(key_pattern, keys, ignore.case = TRUE))
  if (any(lengths(key_parts) != 4L)) {
    malformed <- keys[lengths(key_parts) != 4L]
    stop(
      "Expected <region>_<group1>_vs_<group2> MAST keys. Examples of nonmatching keys: ",
      paste(head(malformed, 3L), collapse = ", "),
      call. = FALSE
    )
  }
  out <- lapply(seq_along(keys), function(i) {
    parts <- key_parts[[i]]
    dg <- as.data.table(deg_list[[i]])
    if (nrow(dg) == 0L) return(NULL)
    if (!all(c("primerid", "coef") %in% names(dg))) {
      stop("Missing primerid or coef in MAST table: ", keys[[i]], call. = FALSE)
    }
    if (!"FDR" %in% names(dg)) {
      if (!"Pr(>Chisq)" %in% names(dg)) {
        stop("Missing FDR and Pr(>Chisq) in MAST table: ", keys[[i]], call. = FALSE)
      }
      dg[, FDR := p.adjust(get("Pr(>Chisq)"), method = "BH")]
    }
    dg[, annotation_hlevel := parts[2]]
    dg[, universal_contrast := paste(parts[3], "vs", parts[4], sep = "_")]
    dg[, .(annotation_hlevel, universal_contrast, primerid, coef, FDR)]
  })
  rbindlist(out)
}

message("Parsing MAST results.")
dt_filtered <- build_mast_table(MAST_RESULT_OBJECT$deg)
if (nrow(dt_filtered) == 0L) {
  stop("No non-empty MAST result tables were found.", call. = FALSE)
}

dt_filtered <- dt_filtered[!annotation_hlevel %in% c("Hindbrain", "Meninge", "Ventricle")]
dt_filtered[, neg_log_fdr := -log10(FDR)]
dt_filtered[is.infinite(neg_log_fdr), neg_log_fdr := 300]
dt_filtered <- dt_filtered[is.finite(coef) & is.finite(neg_log_fdr)]

# Obtain the mouse TF network -------------------------------------------------

if (nzchar(TF_NETWORK_CSV)) {
  if (!file.exists(TF_NETWORK_CSV)) {
    stop("TF network CSV not found. Check TF_NETWORK_CSV.", call. = FALSE)
  }
  raw_net <- read.csv(TF_NETWORK_CSV, stringsAsFactors = FALSE)
} else {
  message("Retrieving the mouse CollecTRI network.")
  raw_net <- decoupleR::get_collectri(organism = "mouse", split_complexes = FALSE)
}
if (!all(c("source", "target") %in% names(raw_net)) ||
    !any(c("weight", "mor") %in% names(raw_net))) {
  stop("TF network requires source, target, and weight or mor.", call. = FALSE)
}
regulation <- if ("mor" %in% names(raw_net)) raw_net$mor else rep(NA_real_, nrow(raw_net))
weight <- if ("weight" %in% names(raw_net)) raw_net$weight else rep(NA_real_, nrow(raw_net))
net <- dplyr::distinct(data.frame(
  source = raw_net$source,
  target = raw_net$target,
  mor = dplyr::coalesce(weight, regulation)
))
if (nrow(net) == 0L || !is.numeric(net$mor) || any(!is.finite(net$mor))) {
  stop("TF network must have finite numeric regulation weights.", call. = FALSE)
}

# Infer TF activities per region and contrast ---------------------------------
dt_filtered[, combined_condition := paste(annotation_hlevel, universal_contrast, sep = " | ")]
all_conditions <- unique(dt_filtered$combined_condition)

# We will store all results here for the overlap analysis in Step 4
master_tf_list <- list()

cat(">>> Running TF Activity Inference on", length(all_conditions), "conditions...\n")

for (target_condition in all_conditions) {

  deg_sub <- dt_filtered[combined_condition == target_condition]
  if(nrow(deg_sub) < 50) next

  # Calculate the signed log-FDR score for each gene.
  deg <- deg_sub %>%
    dplyr::mutate(signed_log_fdr = sign(coef) * neg_log_fdr) %>%
    dplyr::filter(!is.na(primerid), nzchar(primerid), is.finite(signed_log_fdr)) %>%
    dplyr::select(primerid, signed_log_fdr) %>%
    dplyr::distinct(primerid, .keep_all = TRUE) %>%
    tibble::column_to_rownames('primerid')

  if (nrow(deg) == 0L) next

  # Run ULM.
  contrast_acts <- decoupleR::run_ulm(
    mat = as.matrix(deg), net = net, .source = 'source', .target = 'target', .mor = 'mor', minsize = 5
  )

  if(nrow(contrast_acts) == 0) next

  # Save to master list for Step 4
  region_name <- unique(deg_sub$annotation_hlevel)[1]
  contrast_name <- unique(deg_sub$universal_contrast)[1]

  df_save <- contrast_acts %>% dplyr::mutate(Region = region_name, Contrast = contrast_name)
  master_tf_list[[target_condition]] <- df_save

  # Print Individual Barplot for top 15 TFs
  top_tfs <- contrast_acts %>% dplyr::arrange(desc(abs(score))) %>% head(15)
  colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])

  p_bar <- ggplot(top_tfs, aes(x = stats::reorder(source, score), y = score, fill = score)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.5) +
    scale_fill_gradient2(low = colors[1], mid = "whitesmoke", high = colors[2], midpoint = 0) +
    coord_flip() + theme_minimal() +
    labs(title = "Top TF Activity", subtitle = target_condition, x = "TF", y = "Activity Score (ULM)")

  safe_name <- gsub(" |/|\\|", "_", target_condition)
  ggsave(file.path(OUT_DIR, paste0("TF_Barplot_", safe_name, ".pdf")), plot = p_bar, width = 6, height = 5, device = grDevices::cairo_pdf)
}

# Combine all results into one big dataframe
master_tf_df <- dplyr::bind_rows(master_tf_list)


if (nrow(master_tf_df) == 0L) {
  stop("No TF activity results were produced. Check gene IDs and network overlap.", call. = FALSE)
}

# Compare TF activity between the two mutant contrasts -----------------------

cross_reg_list <- list()
all_regions <- unique(master_tf_df$Region)

for (reg in all_regions) {

  df_reg <- master_tf_df %>% dplyr::filter(Region == reg)
  df_bl6 <- df_reg %>% dplyr::filter(grepl("^(WT|wildtype)_vs_mdxbl6$", Contrast, ignore.case = TRUE))
  df_4cv <- df_reg %>% dplyr::filter(grepl("^(WT|wildtype)_vs_mdx4cv$", Contrast, ignore.case = TRUE))

  if(nrow(df_bl6) == 0 || nrow(df_4cv) == 0) next

  if (anyDuplicated(df_bl6$source) || anyDuplicated(df_4cv$source)) {
    stop("Duplicate TF activities for a contrast in region ", reg, ".", call. = FALSE)
  }
  comp_df <- dplyr::full_join(
    df_bl6 %>% dplyr::select(source, score, p_value),
    df_4cv %>% dplyr::select(source, score, p_value),
    by = "source", suffix = c("_bl6", "_4cv")
  )

  # Clean up missing data
  comp_df <- comp_df %>%
    dplyr::mutate(
      p_value_bl6 = ifelse(is.na(p_value_bl6), 1, p_value_bl6),
      p_value_4cv = ifelse(is.na(p_value_4cv), 1, p_value_4cv),
      score_bl6   = ifelse(is.na(score_bl6), 0, score_bl6),
      score_4cv   = ifelse(is.na(score_4cv), 0, score_4cv)
    )

  # Categorize every TF based on its significance
  comp_df <- comp_df %>%
    dplyr::mutate(
      TF_Status = dplyr::case_when(
        p_value_bl6 < 0.05 & p_value_4cv < 0.05 & sign(score_bl6) == sign(score_4cv) & score_bl6 > 0 ~ "Shared: Activated",
        p_value_bl6 < 0.05 & p_value_4cv < 0.05 & sign(score_bl6) == sign(score_4cv) & score_bl6 < 0 ~ "Shared: Repressed",
        p_value_bl6 < 0.05 & p_value_4cv < 0.05 & sign(score_bl6) != sign(score_4cv) ~ "Discordant (Opposite)",
        p_value_bl6 < 0.05 & p_value_4cv >= 0.05 ~ "Exclusive to WT vs mdxbl6",
        p_value_bl6 >= 0.05 & p_value_4cv < 0.05 ~ "Exclusive to WT vs mdx4cv",
        TRUE ~ "Not Significant"
      )
    )

  # Retain all classifications for the cross-region clustering.
  cross_reg_list[[reg]] <- dplyr::mutate(comp_df, Region = reg)

  # Count the groups for the plot subtitle
  counts <- table(comp_df$TF_Status)

  safe_count <- function(cat_name) {
    val <- as.numeric(counts[cat_name])
    if(is.na(val)) return(0) else return(val)
  }

  n_shared_up <- safe_count("Shared: Activated")
  n_shared_dn <- safe_count("Shared: Repressed")
  n_ex_bl6    <- safe_count("Exclusive to WT vs mdxbl6")
  n_ex_4cv    <- safe_count("Exclusive to WT vs mdx4cv")

  subtitle_text <- sprintf("Shared: %d Up, %d Down | Exclusive: %d bl6, %d 4cv",
                           n_shared_up, n_shared_dn, n_ex_bl6, n_ex_4cv)

  # Custom Colors
  status_colors <- c(
    "Shared: Activated" = "#B31B21",
    "Shared: Repressed" = "#1465AC",
    "Exclusive to WT vs mdxbl6" = "#FF7F00",
    "Exclusive to WT vs mdx4cv" = "#984EA3",
    "Discordant (Opposite)" = "#4DAF4A",
    "Not Significant" = "grey85"
  )

  # Generate Scatter Plot
  p_scatter <- ggplot(comp_df, aes(x = score_bl6, y = score_4cv, color = TF_Status)) +
    geom_point(alpha = 0.8, size = 1) +
    scale_color_manual(values = status_colors) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +

    ggrepel::geom_text_repel(
      data = dplyr::filter(comp_df, TF_Status != "Not Significant"),
      aes(label = source),
      size = 6,
      fontface = "bold",
      family = "sans",
      show.legend = FALSE,
      max.overlaps = 8,
      segment.color = "grey50"
    ) +
    theme_bw(base_size = 12, base_family = "sans") +
    labs(
      title = paste("TF Activity Comparison:", reg),
      subtitle = subtitle_text,
      x = "TF Activity Score (WT vs mdxbl6)",
      y = "TF Activity Score (WT vs mdx4cv)",
      color = "TF Classification"
    ) +
    theme(
      legend.position = "bottom",
      legend.direction = "vertical",
      plot.title = element_text(face = "bold", size = 14)
    ) +
    guides(color = guide_legend(ncol = 2))

  if (interactive()) print(p_scatter)

  # Save plots and CSVs
  safe_reg <- gsub(" |/", "_", reg)
  out_scatter <- file.path(OUT_DIR, paste0("TF_Comprehensive_Scatter_", safe_reg, ".pdf"))

  ggsave(
    filename = out_scatter,
    plot = p_scatter,
    width = SCATTER_WIDTH,
    height = SCATTER_HEIGHT,
    device = "pdf",
    useDingbats = FALSE
  )

  out_csv <- file.path(OUT_DIR, paste0("TF_Status_Table_", safe_reg, ".csv"))
  write.csv(comp_df %>% dplyr::filter(TF_Status != "Not Significant") %>% dplyr::arrange(TF_Status, desc(abs(score_bl6))),
            file = out_csv, row.names = FALSE)
}

cross_reg_df <- dplyr::bind_rows(cross_reg_list)

# Cluster classification patterns and run GO enrichment ----------------------

sig_tf_names <- if (nrow(cross_reg_df) > 0L) {
  unique(cross_reg_df$source[cross_reg_df$TF_Status != "Not Significant"])
} else {
  character()
}

if (length(sig_tf_names) == 0L) {
  message("No significant TF classifications for clustering or GO enrichment.")
} else {
  # Create binary categorical matrix for clustering
  categorical_matrix <- cross_reg_df %>%
    dplyr::filter(source %in% sig_tf_names) %>%
    dplyr::select(source, Region, TF_Status) %>%
    dplyr::mutate(value = 1) %>%
    tidyr::pivot_wider(names_from = c(Region, TF_Status), values_from = value, values_fill = 0)

  score_matrix <- as.matrix(categorical_matrix %>% dplyr::select(-source))
  rownames(score_matrix) <- categorical_matrix$source

  # Cluster based on visual presence/absence of categories
  set.seed(42)
  if (nrow(score_matrix) > 1L) {
    dist_matrix <- dist(score_matrix, method = "binary")
    hc <- hclust(dist_matrix, method = "ward.D2")
    tf_factors <- cutree(hc, k = min(N_TF_GROUPS, nrow(score_matrix)))
    tf_order <- hc$labels[hc$order]
  } else {
    tf_order <- rownames(score_matrix)
    tf_factors <- setNames(1L, tf_order)
  }

  factor_df <- data.frame(
    source = names(tf_factors),
    TF_Factor = paste("Factor", tf_factors)
  )

  plot_df <- cross_reg_df %>%
    dplyr::filter(source %in% sig_tf_names) %>%
    dplyr::left_join(factor_df, by = "source")

  plot_df$source <- factor(plot_df$source, levels = tf_order)

  # Generate & Save the Factor-Grouped Heatmap
  status_colors <- c(
    "Shared: Activated" = "#B31B21",
    "Shared: Repressed" = "#1465AC",
    "Exclusive to WT vs mdxbl6" = "#FF7F00",
    "Exclusive to WT vs mdx4cv" = "#984EA3",
    "Discordant (Opposite)" = "#4DAF4A",
    "Not Significant" = "grey93"
  )

  # Calculate a safe Y-axis font size to prevent overlapping TF names
  safe_y_text_size <- max(4, (HEATMAP_HEIGHT / length(sig_tf_names)) * 60)

  p_factor_heat <- ggplot(plot_df, aes(x = Region, y = source, fill = TF_Status)) +
    geom_tile(color = "white", linewidth = 0.5) +

    # Automatically wrap the extremely long legend labels
    scale_fill_manual(
      values = status_colors,
      labels = function(x) stringr::str_wrap(x, width = LEGEND_WRAP)
    ) +
    facet_grid(TF_Factor ~ ., scales = "free_y", space = "free_y") +
    theme_minimal(base_size = 8) +

    # Strip main titles and axis labels to preserve drawing space
    labs(
      title = NULL,
      x = NULL,
      y = NULL,
      fill = "TF Classification"
    ) +
    theme(
      # X-axis: Brain Regions (Angled and shrunken)
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black", size = 6),

      # Y-axis: Apply the dynamic text size calculation
      axis.text.y = element_text(size = safe_y_text_size, face = "bold", color = "black"),

      panel.grid = element_blank(),

      # Legend formatting for narrow plots
      legend.position = "bottom",
      legend.title = element_text(size = 6, face = "bold"),
      legend.text = element_text(size = 5),
      legend.key.size = grid::unit(0.3, "cm"),
      legend.margin = margin(t = 2, r = 0, b = 0, l = 0),

      # Facet strip formatting
      strip.text.y = element_text(angle = 0, face = "bold", size = 6),
      strip.background = element_rect(fill = "grey85", color = "white"),

      # Squeeze the margins tight
      plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "pt")
    ) +
    # Force the legend into two neat columns so it doesn't stretch off the page
    guides(fill = guide_legend(ncol = 2, title.position = "top", title.hjust = 0.5))

  out_heat <- file.path(OUT_DIR, "TF_Clustered_Factors_Heatmap.pdf")
  ggsave(filename = out_heat, plot = p_factor_heat, width = HEATMAP_WIDTH, height = HEATMAP_HEIGHT, limitsize = FALSE, device = grDevices::cairo_pdf)

  cat(">>> Saved customizable heatmap to:", out_heat, "\n")


  cat("\n>>> Running Gene Ontology (GO) Enrichment on TF Factors...\n")

  # Run comparative GO Enrichment (Biological Process)
  go_compare <- tryCatch({
    clusterProfiler::compareCluster(
      geneClusters = source ~ TF_Factor,
      data = factor_df,
      fun = "enrichGO",
      OrgDb = org.Mm.eg.db::org.Mm.eg.db,
      keyType = "SYMBOL",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2
    )
  }, error = function(e) {
    message("  -> GO enrichment failed or found no terms: ", e$message)
    return(NULL)
  })

  if(!is.null(go_compare) && nrow(as.data.frame(go_compare)) > 0) {

    p_go <- enrichplot::dotplot(go_compare, showCategory = 5) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::labs(
        title = "Biological Functions of TF Co-Regulation Factors",
        x = "TF Factor Group"
      ) +
      ggplot2::theme(
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black"),
        axis.text.y = element_text(size = 10, face = "bold", color = "black")
      )

    out_go_plot <- file.path(OUT_DIR, "TF_Factor_GO_Enrichment_Dotplot.pdf")
    ggsave(filename = out_go_plot, plot = p_go, width = 10, height = 8, device = grDevices::cairo_pdf)

    write.csv(as.data.frame(go_compare), file = file.path(OUT_DIR, "TF_Factor_GO_Enrichment_Results.csv"), row.names = FALSE)
    cat("  -> Saved GO enrichment dotplot to:", out_go_plot, "\n")
  }

}

message("Saved available TF analysis results to ", OUT_DIR, ".")
