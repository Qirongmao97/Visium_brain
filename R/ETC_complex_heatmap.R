# Regional expression panels for respiratory-chain and glucose pathway genes
#
# Purpose
#   Extract gene-level coefficients and adjusted p-values from precomputed
#   mouse brain MAST contrasts for respiratory-chain complexes I–V and genes
#   annotated to GO:0006007. This script does not refit MAST models.
#
# Inputs
#   data/ALL_MAST_results.RData, or a file specified by MAST_RESULTS_FILE.
#   The file must contain result_list with a named $deg list. Use keys of
#   the form <region>_<group1>_vs_<group2>. Spaces, underscores, periods, and
#   slashes are treated as equivalent separators, and WT can replace wildtype.
#   Each table must contain a gene symbol, coefficient, and adjusted p-value
#   column. Supported column names are listed in get_gene_stats().
#   Complex gene lists are embedded. GO:0006007 membership is retrieved from
#   the installed org.Mm.eg.db annotation package via GOALL. Use the same
#   package version to reproduce the exact GO membership.
#
# Outputs
#   A CSV and heatmap of mean gene coefficients without a significance
#   filter, a boxplot and dot plot matrix colored by the FDR and effect-size
#   thresholds, and an Excel table of extracted gene-level results.
#   All outputs go to results/mast_pathway_panels unless overridden with
#   MAST_PATHWAY_OUTPUT_DIR.
#
# Usage from the repository root
#   Rscript mast_pathway_expression_panels.R

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(openxlsx)
})

# Settings --------------------------------------------------------------------

MAST_RESULTS_FILE <- Sys.getenv(
  "MAST_RESULTS_FILE",
  unset = file.path("data", "ALL_MAST_results.RData")
)
OUT_DIR <- Sys.getenv(
  "MAST_PATHWAY_OUTPUT_DIR",
  unset = file.path("results", "mast_pathway_panels")
)

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

# Use the study's panel regions directly, without loading a Seurat sample.
DESIRED_REGION_ORDER <- c(
  "Isocortex", "Hippocampus Formation", "Thalamus", "Amygdala",
  "Central subpallium", "Fiber tract"
)
TARGET_REGIONS <- DESIRED_REGION_ORDER
DESIRED_ORDER <- c("WT vs mdxbl6", "WT vs mdx4cv", "mdxbl6 vs mdx4cv")

FDR_CUTOFF <- 0.05
LFC_CUTOFF <- 0.1
USER_HEATMAP_WIDTH <- 2.62
USER_HEATMAP_HEIGHT <- 3.81
USER_MATRIX_WIDTH <- 16
USER_MATRIX_HEIGHT <- 8

# Match a region and ordered contrast to its MAST table -----------------------

normalize_key <- function(x) {
  tolower(gsub("[[:space:]_./]+", "_", trimws(x)))
}

group_aliases <- function(group) {
  if (group == "WT") c("WT", "wildtype") else group
}

find_deg_table <- function(region, comparison, deg_list) {
  groups <- strsplit(comparison, " vs ", fixed = TRUE)[[1]]
  contrast_keys <- outer(
    group_aliases(groups[1]), group_aliases(groups[2]),
    function(a, b) paste(a, "vs", b, sep = "_")
  )
  keys <- normalize_key(paste(region, as.vector(contrast_keys), sep = "_"))
  hits <- which(normalize_key(names(deg_list)) %in% keys)
  if (length(hits) > 1L) {
    stop(
      "Multiple MAST tables match ", region, " / ", comparison, ": ",
      paste(names(deg_list)[hits], collapse = ", "),
      call. = FALSE
    )
  }
  hits
}

get_gene_stats <- function(reg, comp, df, target_genes) {
  cols_lower <- tolower(names(df))
  gene_col <- names(df)[cols_lower %in% c("primerid", "gene", "rowname", "symbol", "id")][1]
  fdr_col <- names(df)[cols_lower %in% c("fdr", "adjp", "p_val_adj", "p.adj", "padj", "p.adjust")][1]
  lfc_col <- names(df)[cols_lower %in% c("coef", "logfc", "avg_log2fc")][1]

  if (anyNA(c(gene_col, fdr_col, lfc_col))) {
    stop(
      "Missing gene, FDR, or coefficient column for ", reg, " / ", comp, ".",
      call. = FALSE
    )
  }

  if ("component" %in% cols_lower) {
    comp_col <- names(df)[cols_lower == "component"][1]
    df <- df[df[[comp_col]] %in% c("H", "logFC"), , drop = FALSE]
  }

  df_sub <- df[df[[gene_col]] %in% target_genes, , drop = FALSE]
  if (nrow(df_sub) == 0L) return(NULL)

  as_stat <- function(x) {
    if (is.factor(x)) x <- as.character(x)
    as.numeric(x)
  }

  data.frame(
    Gene = df_sub[[gene_col]],
    Region = reg,
    Comparison = comp,
    LFC = as_stat(df_sub[[lfc_col]]),
    FDR = as_stat(df_sub[[fdr_col]])
  )
}

# Define gene sets ------------------------------------------------------------
my_gene_lists <- list(
  "Complex I" = c("Acad9", "Aifm1", "Dmac1", "Dmac2", "Ecsit", "Foxred1", "Lyrm2", "mt-Nd1", "mt-Nd2", "mt-Nd3", "mt-Nd4", "mt-Nd4l", "mt-Nd5", "mt-Nd6", "Ndufa1", "Ndufa10", "Ndufa11", "Ndufa12", "Ndufa13", "Ndufa2", "Ndufa3", "Ndufa5", "Ndufa6", "Ndufa7", "Ndufa8", "Ndufa9", "Ndufab1", "Ndufaf1", "Ndufaf2", "Ndufaf3", "Ndufaf4", "Ndufaf5", "Ndufaf6", "Ndufaf7", "Ndufaf8", "Ndufb10", "Ndufb11", "Ndufb2", "Ndufb3", "Ndufb4", "Ndufb5", "Ndufb6", "Ndufb7", "Ndufb8", "Ndufb9", "Ndufc1", "Ndufc2", "Ndufs1", "Ndufs2", "Ndufs3", "Ndufs4", "Ndufs5", "Ndufs6", "Ndufs7", "Ndufs8", "Ndufv1", "Ndufv2", "Ndufv3", "Nubpl", "Timmdc1", "Tmem126a", "Tmem126b", "Tmem186", "Tmem70"),
  "Complex II" = c("Sdha", "Sdhaf1", "Sdhaf2", "Sdhaf3", "Sdhaf4", "Sdhb", "Sdhc", "Sdhd"),
  "Complex III" = c("Bcs1l", "Cyc1", "Lyrm7", "mt-Cytb", "Ttc19", "Uqcc1", "Uqcc2", "Uqcc3", "Uqcr10", "Uqcr11", "Uqcrb", "Uqcrc1", "Uqcrc2", "Uqcrfs1", "Uqcrh", "Uqcrq"),
  "Complex IV" = c("Cep89", "Cmc1", "Cmc2", "Coa3", "Coa4", "Coa5", "Coa6", "Coa7", "Coa8", "Cox10", "Cox11", "Cox14", "Cox15", "Cox16", "Cox17", "Cox18", "Cox19", "Cox20", "Cox4i1", "Cox4i2", "Cox5a", "Cox5b", "Cox6a1", "Cox6a2", "Cox6b1", "Cox6b2", "Cox6c", "Cox7a1", "Cox7a2", "Cox7a2l", "Cox7b", "Cox7c", "Cox8a", "Cox8b", "Cox8c", "Higd1a", "mt-Co1", "mt-Co2", "mt-Co3", "Ndufa4", "Pet100", "Pet117", "Pnkd", "Sco1", "Sco2", "Smim20", "Surf1", "Taco1", "Timm21", "Tmem177"),
  "Complex V" = c("Atp5a1", "Atp5b", "Atp5c1", "Atp5d", "Atp5e", "Atp5g1", "Atp5g2", "Atp5g3", "Atp5h", "Atp5j", "Atp5j2", "Atp5k", "Atp5l", "Atp5md", "Atp5mpl", "Atp5o", "Atp5pb", "Atpaf1", "Atpaf2", "Atpif1", "Atpsckmt", "Dmac2l", "Fmc1", "mt-Atp6", "mt-Atp8", "Tmem70"))

# Dynamically fetch genes for GO:0006007
message("\n>>> Fetching genes for GO:0006007 from org.Mm.eg.db...")
go_query <- AnnotationDbi::select(
  org.Mm.eg.db::org.Mm.eg.db,
  keys = "GO:0006007",
  columns = "SYMBOL",
  keytype = "GOALL"
)
go_symbols <- unique(go_query$SYMBOL[!is.na(go_query$SYMBOL)])
my_gene_lists[["GO:0006007 (Glucose Metab)"]] <- go_symbols
message(">>> Successfully mapped ", length(go_symbols), " genes to GO:0006007.")

# Extract gene-level results ---------------------------------------------------

raw_data_list <- list()
deg_list <- MAST_RESULT_OBJECT$deg

for (comp in DESIRED_ORDER) {
  for (reg in TARGET_REGIONS) {
    idx <- find_deg_table(reg, comp, deg_list)
    if (length(idx) == 0L) {
      message("No MAST table for ", reg, " / ", comp, ".")
      next
    }
    df <- as.data.frame(deg_list[[idx]])
    for (complex_name in names(my_gene_lists)) {
      stats <- get_gene_stats(reg, comp, df, my_gene_lists[[complex_name]])
      if (!is.null(stats)) {
        stats$Complex <- complex_name
        raw_data_list[[paste(complex_name, reg, comp)]] <- stats
      }
    }
  }
}

if (length(raw_data_list) == 0L) {
  stop("No targeted genes found. Check the MAST table keys, regions, and gene symbols.",
       call. = FALSE)
}
output_order <- unlist(lapply(names(my_gene_lists), function(complex_name) {
  unlist(lapply(DESIRED_ORDER, function(comp) {
    paste(complex_name, TARGET_REGIONS, comp)
  }), use.names = FALSE)
}), use.names = FALSE)
raw_df <- bind_rows(raw_data_list[intersect(output_order, names(raw_data_list))])

if (!dir.exists(OUT_DIR)) {
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
}
if (!dir.exists(OUT_DIR)) {
  stop("Cannot create the output directory. Check MAST_PATHWAY_OUTPUT_DIR.", call. = FALSE)
}

# Set universal factor levels to enforce plotting orders
raw_df$Comparison <- factor(raw_df$Comparison, levels = DESIRED_ORDER)
raw_df$Region <- factor(raw_df$Region, levels = DESIRED_REGION_ORDER)
raw_df$Complex <- factor(raw_df$Complex, levels = names(my_gene_lists))

# Calculate unfiltered averages and plot the heatmap --------------------------
heatmap_df <- raw_df %>%
  group_by(Complex, Comparison, Region) %>%
  summarise(
    Avg_LFC = mean(LFC, na.rm = TRUE),
    Total_Genes_Found = n(),
    .groups = "drop"
  )

write.csv(heatmap_df, file.path(OUT_DIR, "Pathway_Global_Averages.csv"), row.names = FALSE)

# Reverse complex order for Heatmap Y-Axis
heatmap_df$Complex <- factor(heatmap_df$Complex, levels = rev(names(my_gene_lists)))

# Generate Heatmap stacked vertically to fit 2.62in x 3.81in
p_heat <- ggplot(heatmap_df, aes(x = Region, y = Complex, fill = Avg_LFC)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#313695", mid = "white", high = "#a50026",
    midpoint = 0, name = "Avg LFC" # Shortened legend name for small width
  ) +
  # Stack contrasts vertically for the narrow panel.
  facet_grid(Comparison ~ ., scales = "free_y", space = "free") +
  theme_bw(base_size = 7) + # Scaled down base font size drastically for small canvas
  theme(
    plot.title = element_text(face = "bold", size = 8, hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 6),
    axis.text.y = element_text(face = "bold", size = 6),
    axis.title = element_blank(), # Removed axis titles to save space on a 2.6in plot
    strip.background = element_rect(fill = "#f0f0f0", color = "black", linewidth = 0.3),
    strip.text = element_text(face = "bold", size = 6, margin = margin(2, 2, 2, 2)),
    panel.spacing = grid::unit(0, "lines"),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.key.width = grid::unit(0.4, "cm"), # Smaller legend key
    legend.key.height = grid::unit(0.2, "cm"),
    legend.title = element_text(size = 6, vjust = 0.8),
    legend.text = element_text(size = 5),
    legend.margin = margin(t = 2),
    plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "pt")
  ) +
  labs(title = "Average LFC per Pathway")

out_heat <- file.path(OUT_DIR, "Panel_Pathway_Heatmap_Vertical_Small.pdf")
ggsave(out_heat, plot = p_heat, width = USER_HEATMAP_WIDTH, height = USER_HEATMAP_HEIGHT, device = "pdf")
message(">>> Generated Vertical Heatmap: ", out_heat)


# Plot gene-level effects, colored by significance -----------------------------

raw_df$Significance <- ifelse(
  raw_df$FDR < FDR_CUTOFF & abs(raw_df$LFC) > LFC_CUTOFF,
  "Significant",
  "Not Significant"
)
raw_df$Significance <- factor(raw_df$Significance, levels = c("Significant", "Not Significant"))

p_box_dot <- ggplot(raw_df, aes(x = Region, y = LFC)) +
  geom_boxplot(fill = "white", color = "black", outlier.shape = NA, alpha = 0.5) +
  geom_jitter(aes(color = Significance), width = 0.2, size = 1.0, alpha = 0.7) +
  scale_color_manual(values = c("Significant" = "#a50026", "Not Significant" = "grey75")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray", linewidth = 0.5) +
  facet_grid(Comparison ~ Complex, scales = "free_y") +
  theme_bw(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "#f0f0f0", color = "black", linewidth = 0.5),
    strip.text = element_text(face = "bold", size = 9),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt")
  ) +
  labs(
    title = "Pathway Gene Expression Profiles (Colored by Significance)",
    y = "Log2 Fold Change",
    x = "Brain Region"
  )

out_matrix <- file.path(OUT_DIR, "Panel_Pathway_Boxplot_Dotplot_Matrix_Colored_Significance.pdf")
ggsave(out_matrix, plot = p_box_dot, width = USER_MATRIX_WIDTH, height = USER_MATRIX_HEIGHT, device = "pdf")
message(">>> Generated Boxplot Matrix: ", out_matrix)

# Export supporting gene-level results ----------------------------------------

out_xlsx <- file.path(OUT_DIR, "Supplementary_Table_Raw_Gene_Stats.xlsx")
openxlsx::write.xlsx(raw_df, out_xlsx, rowNames = FALSE)
message("Exported raw gene statistics to: ", out_xlsx)
