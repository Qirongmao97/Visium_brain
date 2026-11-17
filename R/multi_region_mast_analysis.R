# Multi-region analysis of mouse brain MAST results
#
# Purpose
#   Summarize precomputed pairwise MAST differential expression results across
#   brain regions. Plot regional effects, identify genes shared across regions,
#   and compare GO Biological Process enrichment across genotype contrasts.
#   This script does not fit MAST models.
#
# Input
#   data/ALL_MAST_results.RData, or a file specified by MAST_RESULTS_FILE.
#   The file must contain a list named result_list with contrast tables in
#   $diag and/or $deg. Contrast keys follow <region>_<group1>_vs_<group2>
#   (region names contain no underscores). primerid contains mouse gene symbols.
#
# Outputs
#   Regional volcano PDFs, contrast-level heatmap PDFs, an UpSet SVG, and GO
#   enrichment dotplot and GSEA ridgeplot PDFs when results pass the filters.
#   Default output directory: results/multi_region_mast (override with
#   MAST_OUTPUT_DIR).
#
# Notes
#   Heatmaps and GO over-representation analysis use genes significant in
#   at least three regions with consistent direction. The UpSet analysis uses
#   the same region count without requiring directional consistency. GSEA ranks
#   all available genes by their mean coefficient across retained regions.
#   The heatmaps and ridgeplots use Cairo PDF. Set MAST_PLOT_FONT to override
#   the default sans font, for example to Arial if it is installed.
#
# Usage from the repository root
#   Rscript multi_region_mast_analysis.R

# Setup -----------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(UpSetR)
  library(svglite)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
})

MAST_RESULTS_FILE <- Sys.getenv(
  "MAST_RESULTS_FILE",
  unset = file.path("data", "ALL_MAST_results.RData")
)
OUT_DIR <- Sys.getenv(
  "MAST_OUTPUT_DIR",
  unset = file.path("results", "multi_region_mast")
)
PLOT_FONT <- Sys.getenv("MAST_PLOT_FONT", unset = "sans")

if (!file.exists(MAST_RESULTS_FILE)) {
  stop(
    "MAST results file not found. Place ALL_MAST_results.RData in data/ ",
    "or set MAST_RESULTS_FILE.",
    call. = FALSE
  )
}

# Load the expected object without importing unrelated objects into the workspace.
mast_input <- new.env(parent = emptyenv())
load(MAST_RESULTS_FILE, envir = mast_input)
if (!exists("result_list", envir = mast_input, inherits = FALSE)) {
  stop("The MAST results file must contain an object named result_list.", call. = FALSE)
}
MAST_RESULT_OBJECT <- mast_input[["result_list"]]
if (!is.list(MAST_RESULT_OBJECT) ||
    !any(c("diag", "deg") %in% names(MAST_RESULT_OBJECT))) {
  stop("result_list must be a list with diag and/or deg elements.", call. = FALSE)
}

if (!dir.exists(OUT_DIR)) {
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
}
if (!dir.exists(OUT_DIR)) {
  stop("Cannot create the output directory. Check MAST_OUTPUT_DIR.", call. = FALSE)
}

# --- FILTERS ---
# Regions to Remove
REGIONS_TO_REMOVE <- c("Hindbrain", "Ventricle", "Meninge")

# Groups to Remove
GROUPS_TO_REMOVE  <- c("mdx52")

# Statistical Thresholds
FDR_CUTOFF <- 0.05
LFC_CUTOFF <- 0.1
MIN_REGIONS_SHARED <- 3  # For Core Genes: Must be significant in at least this many regions
# -------------------------------------------------------------------------


#########################################################################################
#                                 Step 1: Function Definitions                          #
#########################################################################################

# --- Function 1: Parse MAST results into a list ---
build_fcHurdle_for_model_autodetect <- function(res_model, use_metric = c("FDR","p"), cutoff = 0.05) {
  use_metric <- match.arg(use_metric)
  keys <- unique(c(names(res_model$diag), names(res_model$deg)))
  out_fc <- vector("list", length(keys)); names(out_fc) <- keys

  for (k in keys) {
    region <- sub("_.*$", "", k)
    tail   <- sub(paste0("^", region, "_"), "", k)
    parts  <- strsplit(tail, "_vs_", fixed = TRUE)[[1]]
    if (length(parts) != 2) next
    g1 <- parts[1]; g2 <- parts[2]

    fc <- NULL
    # Path A: Try $diag
    dd <- res_model$diag[[k]]
    if (is.data.frame(dd) && all(c("component","primerid") %in% colnames(dd))) {
      dt <- as.data.table(dd)
      chosen_contrast <- NULL
      if ("contrast" %in% names(dt)) {
        contrs <- unique(dt$contrast)
        cand1  <- paste0("Group", g2)
        if (cand1 %in% contrs) chosen_contrast <- cand1 else {
          hits <- contrs[grepl(paste0(g2, "$"), contrs)]
          if (length(hits) >= 1) chosen_contrast <- hits[1]
        }
      }
      needed_cols <- c("Pr(>Chisq)","coef","ci.hi","ci.lo")
      if (!is.null(chosen_contrast) && all(needed_cols %in% names(dt))) {
        h  <- dt[contrast == chosen_contrast & component == "H", .(primerid, `Pr(>Chisq)`)]
        lf <- dt[contrast == chosen_contrast & component == "logFC", .(primerid, coef, ci.hi, ci.lo)]
        if (nrow(h) && nrow(lf)) {
          fc <- merge(h, lf, by = "primerid", all = FALSE)
          fc[, FDR := p.adjust(`Pr(>Chisq)`, method = "BH")]
          fc[, `:=`(annotation = region, group1 = g1, group2 = g2)]
        }
      }
    }
    # Path B: Try $deg
    if (is.null(fc) || !nrow(fc)) {
      dg <- res_model$deg[[k]]
      if (is.data.frame(dg) && "primerid" %in% names(dg)) {
        dt <- as.data.table(dg)
        if (!"FDR" %in% names(dt) && "Pr(>Chisq)" %in% names(dt)) dt[, FDR := p.adjust(`Pr(>Chisq)`, "BH")]
        for (nm in c("ci.hi","ci.lo")) if (!(nm %in% names(dt))) dt[, (nm) := NA_real_]
        dt[, `:=`(annotation = region, group1 = g1, group2 = g2)]
        fc <- dt
      }
    }
    out_fc[[k]] <- fc
  }
  return(out_fc)
}

#########################################################################################
#                                 Step 2: Data Cleaning & Filtering                     #
#########################################################################################

cat(">>> Parsing MAST results...\n")
raw_list <- build_fcHurdle_for_model_autodetect(res_model = MAST_RESULT_OBJECT, use_metric = "FDR")

cat(">>> Merging and cleaning...\n")
dt_all <- rbindlist(raw_list, idcol = "key")

# Clean group names and create contrast column
clean_group_func <- function(x) sapply(strsplit(as.character(x), "_"), tail, 1)
dt_all[, clean_g1 := clean_group_func(group1)]
dt_all[, clean_g2 := clean_group_func(group2)]

# -----------------------------------------------------------------------
# FILTERING (Regions AND Groups)
# -----------------------------------------------------------------------
if(length(REGIONS_TO_REMOVE) > 0) {
  cat(paste0(">>> Removing regions: ", paste(REGIONS_TO_REMOVE, collapse=", "), "\n"))
  dt_all <- dt_all[!annotation %in% REGIONS_TO_REMOVE]
}

if(length(GROUPS_TO_REMOVE) > 0) {
  cat(paste0(">>> Removing groups: ", paste(GROUPS_TO_REMOVE, collapse=", "), "\n"))
  dt_all <- dt_all[!clean_g1 %in% GROUPS_TO_REMOVE & !clean_g2 %in% GROUPS_TO_REMOVE]
}

dt_all[, universal_contrast := paste0(clean_g1, " vs ", clean_g2)]

# Handle numeric issues
min_nonzero <- min(dt_all[FDR > 0]$FDR, na.rm = TRUE)
dt_all[FDR == 0, FDR := min_nonzero]
dt_all[, neg_log_fdr := -log10(FDR)]
max_val <- max(dt_all[is.finite(neg_log_fdr)]$neg_log_fdr)
dt_all[!is.finite(neg_log_fdr), neg_log_fdr := max_val]

# Define significance
dt_all[, diff_type := "NS"]
dt_all[FDR < FDR_CUTOFF & coef > LFC_CUTOFF, diff_type := "UP"]
dt_all[FDR < FDR_CUTOFF & coef < -LFC_CUTOFF, diff_type := "DOWN"]

dt_filtered <- dt_all
cat("    Retained regions:", paste(unique(dt_filtered$annotation), collapse=", "), "\n")

if(nrow(dt_filtered) == 0) stop("Error: No data left after filtering.")

#########################################################################################
#                                 Step 3: Volcano Plots                                 #
#########################################################################################

# *** DEFINE YOUR ORDER HERE ***
CONTRAST_ORDER <- c("wildtype vs mdxbl6", "wildtype vs mdx4cv", "mdxbl6 vs mdx4cv")

cat("\n>>> Generating Volcano Plots...\n")

existing_contrasts <- unique(dt_filtered$universal_contrast)
valid_levels <- intersect(CONTRAST_ORDER, existing_contrasts)
final_levels <- c(valid_levels, setdiff(existing_contrasts, valid_levels))

dt_filtered$universal_contrast <- factor(dt_filtered$universal_contrast, levels = final_levels)

unique_regions <- unique(dt_filtered$annotation)

for (region in unique_regions) {
  plot_data <- dt_filtered[annotation == region]
  top_genes <- plot_data[diff_type != "NS"][order(FDR), .SD[1:10], by = universal_contrast]

  p <- ggplot(plot_data, aes(x = coef, y = neg_log_fdr, color = diff_type)) +
    geom_point(alpha = 0.6, size = 1) +
    facet_wrap(~ universal_contrast, scales = "fixed", nrow = 1) +
    scale_color_manual(values = c("UP"="#B31B21", "DOWN"="#1465AC", "NS"="grey80")) +
    geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype="dashed", color="grey30") +
    geom_hline(yintercept = -log10(FDR_CUTOFF), linetype="dashed", color="grey30") +
    geom_text_repel(data = top_genes, aes(label = primerid), size = 3, color = "black", max.overlaps = 15) +
    theme_bw() +
    labs(title = region, x = "Log2FC", y = "-log10(FDR)") +
    theme(legend.position = "none")

  if (interactive()) print(p)
  ggsave(file.path(OUT_DIR, paste0("Volcano_", region, ".pdf")), p, width = 12, height = 5)
}

#########################################################################################
#                                 Step 4: Heatmap Loop Analysis                         #
#########################################################################################


cat("\n>>> Starting Heatmap Analysis for all contrasts...\n")
all_contrasts <- unique(dt_filtered$universal_contrast)

for (current_contrast in all_contrasts) {

  cat(paste0("\n=== Analyzing Contrast: ", current_contrast, " ===\n"))

  # 1. Identify Core Genes & Calculate Stats
  sig_dt <- dt_filtered[universal_contrast == current_contrast & FDR < FDR_CUTOFF & abs(coef) > LFC_CUTOFF]

  gene_stats <- sig_dt[, .(
    n_regions = uniqueN(annotation),
    mean_lfc  = mean(coef),
    is_consistent_up = all(coef > 0),
    is_consistent_down = all(coef < 0)
  ), by = primerid]

  # 2. Filter & Sort UP/DOWN Candidates
  cands_up   <- gene_stats[n_regions >= MIN_REGIONS_SHARED & is_consistent_up == TRUE]
  cands_down <- gene_stats[n_regions >= MIN_REGIONS_SHARED & is_consistent_down == TRUE]

  if(nrow(cands_up) > 0) cands_up <- cands_up[order(-n_regions, -abs(mean_lfc))]
  if(nrow(cands_down) > 0) cands_down <- cands_down[order(-n_regions, -abs(mean_lfc))]

  cat(paste0("  Core UP: ", nrow(cands_up), " | Core DOWN: ", nrow(cands_down), "\n"))

  # 3. Select Top Genes for Heatmap
  top_up   <- head(cands_up$primerid, 20)
  top_down <- head(cands_down$primerid, 20)
  top_genes <- c(top_up, top_down)

  # 4. Generate Heatmap
  if(length(top_genes) > 0) {
    # Pull data for top genes (including non-significant rows so the grid is full)
    plot_data_hm <- dt_filtered[universal_contrast == current_contrast & primerid %in% top_genes]

    # Mark cells that meet the FDR and absolute effect-size thresholds.
    plot_data_hm[, sig_label := ifelse(FDR < FDR_CUTOFF & abs(coef) >= LFC_CUTOFF, "*", "")]

    # Order Facets
    plot_data_hm[, category := ifelse(primerid %in% cands_up$primerid, "UP", "DOWN")]
    plot_data_hm[, category := factor(category, levels = c("UP", "DOWN"))]

    # Order Y-axis genes
    final_levels <- unique(c(top_down, rev(top_up)))
    plot_data_hm$primerid <- factor(plot_data_hm$primerid, levels = final_levels)

    # Order X-axis (Annotations)
    region_order <- c("Isocortex", "Hippocampus", "Thalamus","Amygdala","Central","Fiber","Unknown")
    plot_data_hm[, annotation := factor(annotation, levels = intersect(region_order, unique(plot_data_hm$annotation)))]

    # Build Plot
    p_hm <- ggplot(plot_data_hm, aes(x = annotation, y = primerid, fill = coef)) +
      geom_tile(color = "white") +
      geom_text(aes(label = sig_label), color = "black", size = 5, vjust = 0.8, family = PLOT_FONT)  +
      scale_fill_gradient2(
        low = "#1465AC", mid = "white", high = "#B31B21", midpoint = 0,
        limits = c(-1, 1),
        oob = scales::squish
      ) +
      facet_grid(category ~ ., scales = "free_y", space = "free_y") +
      labs(
        x = NULL, y = NULL,
        fill = "Coefficient"
      ) +
      theme_minimal(base_family = PLOT_FONT) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, face="bold", family = PLOT_FONT),
        strip.background = element_rect(fill="grey95", color=NA),
        strip.text = element_text(size = 20, face = "bold", family = PLOT_FONT)
      )

    if (interactive()) print(p_hm)

    # Export a vector PDF with the selected font.
    ggsave(
      filename = file.path(OUT_DIR, paste0("Heatmap_", gsub(" ", "_", current_contrast), ".pdf")),
      plot = p_hm,
      width = 8,
      height = 10,
      device = cairo_pdf,
      family = PLOT_FONT
    )
  }
}
#########################################################################################
#                                 Step 5: UpSet Plot Analysis                           #
#########################################################################################

cat("\n>>> Generating UpSet Plot for Core DEGs across Contrasts...\n")

# 1. First, ensure upset_list is generated based on your Core Genes logic
upset_list <- list()
all_contrasts <- unique(dt_filtered$universal_contrast)

for (current_contrast in all_contrasts) {
  sig_dt <- dt_filtered[universal_contrast == current_contrast & FDR < FDR_CUTOFF & abs(coef) > LFC_CUTOFF]

  # Calculate how many regions each gene is significant in
  gene_stats <- sig_dt[, .(n_regions = uniqueN(annotation)), by = primerid]

  # Filter for genes meeting the MIN_REGIONS_SHARED threshold
  core_genes <- gene_stats[n_regions >= MIN_REGIONS_SHARED]$primerid

  if(length(core_genes) > 0) {
    upset_list[[current_contrast]] <- core_genes
  }
}

# 2. Now run the plotting logic
if(length(upset_list) > 0 && max(lengths(upset_list)) > 0) {

  # -------------------------------------------------------------------------
  # EXPORT USING SVGLITE (Clean vector paths, no clipping masks)
  # -------------------------------------------------------------------------
  svg_filename <- file.path(OUT_DIR, "UpSet_Plot_Core_DEGs_Scaled.svg")

  svglite(svg_filename,
          width = 5.2,
          height = 5.4,
          fix_text_size = FALSE) # Helps Illustrator read the fonts correctly

  # ASSIGN the plot to an object
  my_upset_plot <- UpSetR::upset(
    fromList(upset_list),
    nsets = length(upset_list),
    order.by = "freq",
    keep.order = TRUE,
    mainbar.y.max = 220,
    main.bar.color = "#2E3440",
    sets.bar.color = "#1465AC",
    matrix.color = "#B31B21",
    set_size.show = TRUE,
    text.scale = 1.5,
    point.size = 3.5,
    line.size = 1.2,
    mb.ratio = c(0.60, 0.40)
  )

  # EXPLICITLY PRINT IT
  print(my_upset_plot)

  dev.off() # Close the SVG device
  # -------------------------------------------------------------------------

  cat(paste0("    UpSet plot saved successfully as: ", svg_filename, "\n"))

} else {
  cat("    Not enough overlapping data across contrasts to generate an UpSet plot.\n")
}

#########################################################################################
#                                 Step 6: Core DEG GO Enrichment                        #
#########################################################################################

cat("\n>>> Extracting Multi-Region (Core) Genes per Contrast...\n")

all_contrasts <- unique(dt_filtered$universal_contrast)
core_genes_list <- list()

for (current_contrast in all_contrasts) {

  # 1. Identify Significant Genes for the specific contrast
  sig_dt <- dt_filtered[universal_contrast == current_contrast & FDR < FDR_CUTOFF & abs(coef) > LFC_CUTOFF]

  # 2. Calculate Stats (Number of regions and directional consistency)
  gene_stats <- sig_dt[, .(
    n_regions = uniqueN(annotation),
    is_consistent_up = all(coef > 0),
    is_consistent_down = all(coef < 0)
  ), by = primerid]

  # 3. Filter for Core UP/DOWN Candidates based on MIN_REGIONS_SHARED
  cands_up   <- gene_stats[n_regions >= MIN_REGIONS_SHARED & is_consistent_up == TRUE]
  cands_down <- gene_stats[n_regions >= MIN_REGIONS_SHARED & is_consistent_down == TRUE]

  # 4. Combine into a single vector of core genes for this contrast
  core_genes <- unique(c(cands_up$primerid, cands_down$primerid))

  if (length(core_genes) > 0) {
    core_genes_list[[current_contrast]] <- core_genes
    cat(paste0("  -> ", current_contrast, ": ", length(core_genes), " core genes extracted.\n"))
  } else {
    cat(paste0("  -> ", current_contrast, ": 0 core genes found. Skipping.\n"))
  }
}

# -------------------------------------------------------------------------
# RUN GO ENRICHMENT & PLOT
# -------------------------------------------------------------------------

if(length(core_genes_list) > 0) {
  cat("\n>>> Running GO Enrichment via compareCluster...\n")

  # Note: Adjust 'universe' parameter here if you have a defined background gene list
  go_res <- tryCatch({
    compareCluster(
      geneCluster = core_genes_list,
      fun = "enrichGO",
      OrgDb = org.Mm.eg.db,
      keyType = "SYMBOL",  # Change to ENSEMBL or ENTREZID if your primerids are not gene symbols
      ont = "BP",          # Biological Process
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2
    )
  }, error = function(e) {
    cat("  -> Error running compareCluster: ", e$message, "\n")
    return(NULL)
  })

  if(!is.null(go_res) && nrow(as.data.frame(go_res)) > 0) {

    cat("\n>>> Simplifying terms and generating Dotplot...\n")

    # Remove highly redundant GO terms
    go_res_simplified <- simplify(go_res, cutoff = 0.7, by = "p.adjust", select_fun = min)

    # Generate the Dotplot
    p_dot <- enrichplot::dotplot(go_res_simplified,
                     showCategory = 10,   # Top 10 terms per contrast
                     font.size = 10,
                     title = paste0("GO Enrichment: Core DEGs (\u2265 ", MIN_REGIONS_SHARED, " shared regions)")) +
      theme_bw(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black"),
        axis.text.y = element_text(color = "black", lineheight = 0.8),
        plot.title = element_text(face = "bold", hjust = 0.5)
      ) +
      # Wrap long GO terms so they don't stretch the plot
      scale_y_discrete(labels = scales::label_wrap(width = 45))

    if (interactive()) print(p_dot)

    # Save the plot
    out_file <- file.path(OUT_DIR, "Dotplot_Core_DEGs_GO.pdf")
    ggsave(out_file, plot = p_dot, width = 8, height = 7, device = "pdf", useDingbats = FALSE)
    cat(paste0("  -> Saved Dotplot successfully to: ", out_file, "\n"))

  } else {
    cat("  -> No significant GO terms found across the extracted core gene sets.\n")
  }
} else {
  cat("\n>>> Enrichment skipped. No core genes met the threshold across any contrast.\n")
}


#########################################################################################
#                                 Step 7: GSEA Ridgeplots per Contrast                  #
#########################################################################################

cat("\n>>> Running Global GSEA and generating Ridgeplots per Contrast...\n")

all_contrasts <- unique(dt_filtered$universal_contrast)

for (current_contrast in all_contrasts) {

  cat(paste0("\n=== Processing GSEA for: ", current_contrast, " ===\n"))

  # 1. Extract ALL genes for this contrast (not just significant ones!)
  contrast_dt <- dt_filtered[universal_contrast == current_contrast]

  # 2. Average the Log2FC (coef) across all regions to get a single global ranking per gene
  ranked_dt <- contrast_dt[, .(mean_lfc = mean(coef, na.rm = TRUE)), by = primerid]

  # Remove any NAs just to be safe
  ranked_dt <- ranked_dt[!is.na(mean_lfc)]

  # 3. Create the named numeric vector required by clusterProfiler for GSEA
  gene_list <- ranked_dt$mean_lfc
  names(gene_list) <- ranked_dt$primerid

  # 4. Sort the list in decreasing order (CRITICAL FOR GSEA)
  gene_list <- sort(gene_list, decreasing = TRUE)

  # 5. Run GSEA
  gsea_res <- tryCatch({
    gseGO(
      geneList      = gene_list,
      OrgDb         = org.Mm.eg.db,
      keyType       = "SYMBOL",
      ont           = "BP",
      minGSSize     = 10,
      maxGSSize     = 500,
      pvalueCutoff  = 0.05,
      pAdjustMethod = "BH",
      verbose       = FALSE
    )
  }, error = function(e) {
    cat("  -> Error running gseGO: ", e$message, "\n")
    return(NULL)
  })

  # 6. Plotting
  if (!is.null(gsea_res) && nrow(as.data.frame(gsea_res)) > 0) {

    cat("  -> Significant GSEA terms found. Simplifying and generating ridgeplot...\n")

    # Simplify redundant GO terms for a cleaner plot
    gsea_res_simp <- simplify(gsea_res, cutoff = 0.5, by = "p.adjust", select_fun = min)

    # Check again if terms remain after simplification
    if (nrow(as.data.frame(gsea_res_simp)) > 0) {

      # Generate the Ridgeplot
      p_ridge <- enrichplot::ridgeplot(gsea_res_simp, showCategory = 20, fill = "p.adjust") +
        labs(
          title = paste0("GSEA Ridgeplot: ", current_contrast),
          x = "Global Enrichment Distribution (Log2FC)",
          y = NULL
        ) +
        theme_bw(base_size = 12) +
        theme(
          axis.text.y = element_text(color = "black", size = 11, lineheight = 0.8, face = "bold"),
          axis.text.x = element_text(color = "black", size = 10),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
          panel.grid.minor = element_blank()
        ) +
        # Wrap long GO terms automatically
        scale_y_discrete(labels = scales::label_wrap(width = 45))

      # Save the plot (Using cairo_pdf to prevent unicode warning issues!)
      safe_name <- gsub(" ", "_", current_contrast)
      out_file <- file.path(OUT_DIR, paste0("Ridgeplot_GSEA_", safe_name, ".pdf"))

      ggsave(out_file, plot = p_ridge, width = 10, height = 8, device = cairo_pdf)
      cat(paste0("  -> Saved Ridgeplot successfully to: ", out_file, "\n"))

    } else {
      cat("  -> No terms remained after simplification.\n")
    }
  } else {
    cat("  -> No significant GSEA terms found for this contrast.\n")
  }
}

cat("\n>>> GSEA Analysis Complete!\n")
