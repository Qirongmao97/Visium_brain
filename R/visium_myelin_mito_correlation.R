# Myelin and mitochondrial respiratory-chain module scores in mouse brain Visium
#
# Purpose
#   Merge filtered Visium samples, score a five-gene myelin module and
#   respiratory-chain complexes I–V, then correlate the scores across spots
#   within each genotype and brain region.
#
# Inputs
#   Run after filtering and annotating the individual Seurat objects:
#   wildtype_1_filtered through wildtype_4_filtered, mdxbl6_1_filtered through
#   mdxbl6_3_filtered, and mdx4cv_1_filtered through mdx4cv_3_filtered.
#   Each object must have a Spatial assay and annotation_hlevel metadata.
#   Define complex_I through complex_V as character vectors of mouse gene
#   symbols matching the Spatial assay. Supply the original study gene lists.
#
# Outputs
#   visium_merged contains the module scores. cor_summary contains the
#   spot-level Spearman correlations, and p_cor_heatmap is the faceted plot.
#   These are descriptive, pooled spot-level correlations, not replicate-level
#   estimates. No files are written.
#
# Usage
#   With the required inputs in memory, run:
#   source("visium_myelin_mito_correlation.R")
#
# Requirements
#   Seurat v5 (for JoinLayers), dplyr, tidyr, and ggplot2.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# Check the required objects and gene sets ------------------------------------

sample_ids <- c(
  paste0("wildtype_", seq_len(4)),
  paste0("mdxbl6_", seq_len(3)),
  paste0("mdx4cv_", seq_len(3))
)
sample_objects <- paste0(sample_ids, "_filtered")
complex_ids <- paste0("complex_", c("I", "II", "III", "IV", "V"))
required_inputs <- c(sample_objects, complex_ids)
script_env <- environment()
missing_inputs <- required_inputs[
  !vapply(required_inputs, exists, logical(1), envir = script_env, inherits = TRUE)
]
if (length(missing_inputs) > 0L) {
  stop("Missing in-memory inputs: ", paste(missing_inputs, collapse = ", "), call. = FALSE)
}

filtered_objects <- mget(sample_objects, envir = script_env, inherits = TRUE)
complex_gene_sets <- mget(complex_ids, envir = script_env, inherits = TRUE)
for (nm in names(filtered_objects)) {
  obj <- filtered_objects[[nm]]
  if (!inherits(obj, "Seurat") ||
      !"Spatial" %in% Assays(obj) ||
      !"annotation_hlevel" %in% names(obj[[]])) {
    stop(nm, " must be a Seurat object with a Spatial assay and annotation_hlevel metadata.",
         call. = FALSE)
  }
}
for (nm in complex_ids) {
  genes <- complex_gene_sets[[nm]]
  if (!is.character(genes) || length(genes) == 0L || anyNA(genes) || any(!nzchar(genes))) {
    stop(nm, " must be a non-empty character vector of gene symbols.", call. = FALSE)
  }
}

# Label and merge samples -----------------------------------------------------

add_meta <- function(obj, group_name, sample_id) {
  obj$Group <- group_name
  obj$Sample <- sample_id
  obj
}

group_labels <- rep(c("wildtype", "mdxbl6", "mdx4cv"), times = c(4, 3, 3))
all_objects <- Map(add_meta, filtered_objects, group_labels, sample_ids)
names(all_objects) <- sample_ids

# Prefix spot barcodes with the sample IDs after adding metadata.
visium_merged <- merge(
  x = all_objects[[1]],
  y = all_objects[-1],
  add.cell.ids = sample_ids
)

# Normalize and score the Spatial assay ---------------------------------------

visium_merged <- JoinLayers(visium_merged, assay = "Spatial")
visium_merged <- NormalizeData(
  visium_merged, assay = "Spatial", normalization.method = "LogNormalize"
)

myelin_genes <- c("Mbp", "Plp1", "Mobp", "Cnp", "Mag")
score_prefixes <- paste0("Complex_", c("I", "II", "III", "IV", "V"), "_Score")
module_gene_sets <- c(list(Myelin_Score = myelin_genes), setNames(complex_gene_sets, score_prefixes))

for (score_prefix in names(module_gene_sets)) {
  visium_merged <- AddModuleScore(
    visium_merged,
    features = list(module_gene_sets[[score_prefix]]),
    assay = "Spatial",
    name = score_prefix
  )
}

# Reshape scores and calculate Spearman correlations --------------------------

REGIONS_TO_REMOVE <- c("Hindbrain", "Ventricle", "Meninge")
data_df <- visium_merged[[]] %>%
  dplyr::select(
    Group, Sample, annotation_hlevel,
    Myelin_Score1,
    Complex_I_Score1, Complex_II_Score1, Complex_III_Score1,
    Complex_IV_Score1, Complex_V_Score1
  ) %>%
  filter(!annotation_hlevel %in% REGIONS_TO_REMOVE)

data_long <- data_df %>%
  pivot_longer(
    cols = starts_with("Complex"),
    names_to = "Complex",
    values_to = "Mito_Score"
  ) %>%
  mutate(
    Complex = sub("_Score1$", "", Complex),
    Complex = gsub("_", " ", Complex, fixed = TRUE)
  )

spearman_complete <- function(x, y) {
  keep <- complete.cases(x, y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < 2L || length(unique(x)) < 2L || length(unique(y)) < 2L) {
    return(NA_real_)
  }
  cor(x, y, method = "spearman")
}

cor_summary <- data_long %>%
  group_by(Group, annotation_hlevel, Complex) %>%
  filter(n() > 10L) %>%
  summarise(
    correlation = spearman_complete(Myelin_Score1, Mito_Score),
    .groups = "drop"
  )

# Plot one heatmap panel per respiratory-chain complex ------------------------

p_cor_heatmap <- ggplot(cor_summary, aes(x = Group, y = annotation_hlevel, fill = correlation)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient2(
    low = "#2166ac", mid = "white", high = "#b2182b",
    midpoint = 0, limits = c(-1, 1), space = "Lab",
    name = "Spearman\nCorrelation"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", size = 12)
  ) +
  facet_wrap(~ Complex, nrow = 1) +
  labs(
    title = "Myelin and respiratory-chain module score correlations",
    subtitle = "Spot-level Spearman correlations within genotypes and brain regions",
    x = "Mouse Model",
    y = "Spatial Annotation"
  )

if (interactive() && nrow(cor_summary) > 0L) print(p_cor_heatmap)
