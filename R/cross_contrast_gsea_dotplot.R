# Cross-contrast GO Biological Process GSEA dot plot
#
# Purpose
#   Compare precomputed GSEA results across genotype contrasts and brain
#   regions. This script does not run GSEA or fit differential expression models.
#
# Inputs
#   GSEA_Master_MultiDB_<contrast>.csv files, with spaces in the contrast
#   replaced by underscores. Each file must contain Database, Cluster,
#   Description, NES, and p.adjust columns.
#   By default, CSVs are read from results/multi_region_mast relative to the
#   repository root. Set GSEA_INPUT_DIR to override this location.
#
# Term selection
#   Keep GO_BP entries and exclude Unknown regions. Select the five terms per
#   contrast with the smallest adjusted p-values, using the largest absolute
#   NES to break ties. Plot the union of these terms across all available
#   contrasts and the six listed brain regions. Color indicates the signed
#   normalized enrichment score and size indicates -log10(adjusted p-value).
#
# Output
#   Cross_Contrast_GSEA_Comparison.pdf in the input directory by default.
#   Set GSEA_OUTPUT_DIR to write the PDF elsewhere.
#
# Usage from the repository root
#   Rscript cross_contrast_gsea_dotplot.R

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

# Settings --------------------------------------------------------------------

INPUT_DIR <- Sys.getenv(
  "GSEA_INPUT_DIR",
  unset = file.path("results", "multi_region_mast")
)
OUT_DIR <- Sys.getenv("GSEA_OUTPUT_DIR", unset = INPUT_DIR)

TARGET_CONTRASTS <- c("wildtype vs mdxbl6", "wildtype vs mdx4cv", "mdxbl6 vs mdx4cv")
CONTRAST_LABELS <- c("WT vs mdxbl6", "WT vs mdx4cv", "mdxbl6 vs mdx4cv")
REGION_ORDER <- c(
  "Isocortex", "Hippocampus", "Thalamus", "Amygdala", "Central subpallium", "Fiber"
)
TOP_N_TERMS <- 5
WRAP_WIDTH <- 35

# Read and combine available contrasts ---------------------------------------

message("Loading and combining GSEA results.")
required_columns <- c("Database", "Cluster", "Description", "NES", "p.adjust")
all_gsea_data <- list()

for (contrast in TARGET_CONTRASTS) {
  csv_name <- paste0(
    "GSEA_Master_MultiDB_", gsub(" ", "_", contrast, fixed = TRUE), ".csv"
  )
  csv_path <- file.path(INPUT_DIR, csv_name)
  if (!file.exists(csv_path)) {
    message("Skipping missing GSEA input: ", csv_name)
    next
  }

  contrast_data <- read.csv(csv_path, stringsAsFactors = FALSE)
  missing_columns <- setdiff(required_columns, names(contrast_data))
  if (length(missing_columns) > 0L) {
    stop(
      "Missing columns in ", csv_name, ": ", paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  contrast_data$universal_contrast <- rep(contrast, nrow(contrast_data))
  all_gsea_data[[contrast]] <- contrast_data
}

if (length(all_gsea_data) == 0L) {
  stop("No GSEA CSVs found. Check GSEA_INPUT_DIR.", call. = FALSE)
}

df_combined <- bind_rows(all_gsea_data)

# Filter to GO Biological Process and select top terms ------------------------

df_go <- df_combined %>%
  filter(Database == "GO_BP", tolower(Cluster) != "unknown")

top_terms_per_contrast <- df_go %>%
  group_by(universal_contrast, Description) %>%
  summarise(
    max_abs_NES = max(abs(NES), na.rm = TRUE),
    min_p = min(p.adjust, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(universal_contrast) %>%
  arrange(min_p, desc(max_abs_NES)) %>%
  slice_head(n = TOP_N_TERMS) %>%
  pull(Description) %>%
  unique()

plot_data <- df_go %>%
  filter(Description %in% top_terms_per_contrast)

# Order terms by contrast sharing, then by adjusted p-value -------------------

term_sharing <- plot_data %>%
  group_by(Description) %>%
  summarise(
    n_contrasts_shared = n_distinct(universal_contrast),
    overall_min_p = min(p.adjust),
    .groups = "drop"
  ) %>%
  arrange(desc(n_contrasts_shared), overall_min_p)

plot_data <- plot_data %>%
  mutate(
    Description = str_wrap(Description, width = WRAP_WIDTH),
    Description = factor(
      Description,
      levels = rev(str_wrap(term_sharing$Description, width = WRAP_WIDTH))
    ),
    universal_contrast = factor(
      universal_contrast, levels = TARGET_CONTRASTS, labels = CONTRAST_LABELS
    ),
    Cluster = factor(Cluster, levels = REGION_ORDER)
  ) %>%
  filter(!is.na(Cluster))

if (nrow(plot_data) == 0L) {
  stop("No selected GO_BP terms remain in the configured brain regions.", call. = FALSE)
}

# Plot and save ---------------------------------------------------------------

message("Generating the GSEA comparison plot.")
p_cross_contrast <- ggplot(
  plot_data,
  aes(x = universal_contrast, y = Description, size = -log10(p.adjust), color = NES)
) +
  geom_point() +
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0, name = "NES\n(Direction)"
  ) +
  facet_grid(. ~ Cluster) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 8),
    axis.text.y = element_text(size = 7, lineheight = 0.8),
    strip.background = element_rect(fill = "lightgrey", color = "black"),
    strip.text = element_text(face = "bold", size = 9),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.key.size = grid::unit(0.4, "cm")
  ) +
  labs(x = "Contrast Comparison", y = "GO Biological Process", size = "-log10(p.adj)")

if (interactive()) print(p_cross_contrast)

if (!dir.exists(OUT_DIR)) {
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
}
if (!dir.exists(OUT_DIR)) {
  stop("Cannot create the output directory. Check GSEA_OUTPUT_DIR.", call. = FALSE)
}

plot_path <- file.path(OUT_DIR, "Cross_Contrast_GSEA_Comparison.pdf")
ggsave(plot_path, plot = p_cross_contrast, width = 8.6, height = 4)
message("Saved GSEA plot to: ", plot_path)
