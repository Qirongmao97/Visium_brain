# Spatial maps of respiratory-chain and glucose pathway module scores
#
# Purpose
#   Score respiratory-chain complexes I–V and genes annotated to GO:0006007
#   in a merged mouse brain Visium dataset, then map each score across
#   the tissue sections.
#
# Inputs
#   A merged Seurat object with spatial images and a joined, log-normalized
#   Spatial assay. Either provide visium_merged in memory or set
#   VISIUM_MERGED_RDS to the path of a saved Seurat object. The script does
#   not merge or normalize samples.
#   Complex gene sets are embedded. GO:0006007 membership is read from the
#   installed org.Mm.eg.db via GOALL. Use the same annotation package version
#   to reproduce the exact GO membership.
#
# Outputs
#   Six Spatial_Distribution_<module>_ModuleScore.pdf files by default in
#   results/visium_spatial_modules. Set VISIUM_SPATIAL_OUTPUT_DIR to override.
#   When sourced, visium_merged also retains the score columns in memory.
#   The color scale is fixed at 0–1.2 for all panels. Negative scores and scores
#   above 1.2 are clipped for display only. The stored scores are unchanged.
#
# Usage from the repository root
#   With visium_merged in memory: source("visium_spatial_pathway_scores.R")
#   With VISIUM_MERGED_RDS set: Rscript visium_spatial_pathway_scores.R

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

# Input and figure settings ----------------------------------------------------

INPUT_RDS <- Sys.getenv("VISIUM_MERGED_RDS", unset = "")
OUT_DIR <- Sys.getenv(
  "VISIUM_SPATIAL_OUTPUT_DIR",
  unset = file.path("results", "visium_spatial_modules")
)

SPATIAL_WIDTH <- 20
SPATIAL_HEIGHT <- 8
N_COLUMNS <- 5
SCORE_LIMITS <- c(0, 1.2)
SCORE_BREAKS <- c(0, 0.4, 0.8, 1.2)

if (nzchar(INPUT_RDS)) {
  if (!file.exists(INPUT_RDS)) {
    stop("Merged Visium RDS file not found. Check VISIUM_MERGED_RDS.", call. = FALSE)
  }
  visium_merged <- readRDS(INPUT_RDS)
} else if (!exists("visium_merged", inherits = TRUE)) {
  stop("Provide visium_merged in memory or set VISIUM_MERGED_RDS.", call. = FALSE)
}

if (!inherits(visium_merged, "Seurat") || !"Spatial" %in% Assays(visium_merged)) {
  stop("visium_merged must be a Seurat object with a Spatial assay.", call. = FALSE)
}
if (!"data" %in% Layers(visium_merged[["Spatial"]])) {
  stop("The Spatial assay must contain a joined, normalized data layer.", call. = FALSE)
}
if (length(Images(visium_merged)) == 0L) {
  stop("visium_merged has no spatial images to plot.", call. = FALSE)
}

# Define gene sets ------------------------------------------------------------

my_gene_lists <- list(
  "Complex I" = c("Acad9", "Aifm1", "Dmac1", "Dmac2", "Ecsit", "Foxred1", "Lyrm2", "mt-Nd1", "mt-Nd2", "mt-Nd3", "mt-Nd4", "mt-Nd4l", "mt-Nd5", "mt-Nd6", "Ndufa1", "Ndufa10", "Ndufa11", "Ndufa12", "Ndufa13", "Ndufa2", "Ndufa3", "Ndufa5", "Ndufa6", "Ndufa7", "Ndufa8", "Ndufa9", "Ndufab1", "Ndufaf1", "Ndufaf2", "Ndufaf3", "Ndufaf4", "Ndufaf5", "Ndufaf6", "Ndufaf7", "Ndufaf8", "Ndufb10", "Ndufb11", "Ndufb2", "Ndufb3", "Ndufb4", "Ndufb5", "Ndufb6", "Ndufb7", "Ndufb8", "Ndufb9", "Ndufc1", "Ndufc2", "Ndufs1", "Ndufs2", "Ndufs3", "Ndufs4", "Ndufs5", "Ndufs6", "Ndufs7", "Ndufs8", "Ndufv1", "Ndufv2", "Ndufv3", "Nubpl", "Timmdc1", "Tmem126a", "Tmem126b", "Tmem186", "Tmem70"),
  "Complex II" = c("Sdha", "Sdhaf1", "Sdhaf2", "Sdhaf3", "Sdhaf4", "Sdhb", "Sdhc", "Sdhd"),
  "Complex III" = c("Bcs1l", "Cyc1", "Lyrm7", "mt-Cytb", "Ttc19", "Uqcc1", "Uqcc2", "Uqcc3", "Uqcr10", "Uqcr11", "Uqcrb", "Uqcrc1", "Uqcrc2", "Uqcrfs1", "Uqcrh", "Uqcrq"),
  "Complex IV" = c("Cep89", "Cmc1", "Cmc2", "Coa3", "Coa4", "Coa5", "Coa6", "Coa7", "Coa8", "Cox10", "Cox11", "Cox14", "Cox15", "Cox16", "Cox17", "Cox18", "Cox19", "Cox20", "Cox4i1", "Cox4i2", "Cox5a", "Cox5b", "Cox6a1", "Cox6a2", "Cox6b1", "Cox6b2", "Cox6c", "Cox7a1", "Cox7a2", "Cox7a2l", "Cox7b", "Cox7c", "Cox8a", "Cox8b", "Cox8c", "Higd1a", "mt-Co1", "mt-Co2", "mt-Co3", "Ndufa4", "Pet100", "Pet117", "Pnkd", "Sco1", "Sco2", "Smim20", "Surf1", "Taco1", "Timm21", "Tmem177"),
  "Complex V" = c("Atp5a1", "Atp5b", "Atp5c1", "Atp5d", "Atp5e", "Atp5g1", "Atp5g2", "Atp5g3", "Atp5h", "Atp5j", "Atp5j2", "Atp5k", "Atp5l", "Atp5md", "Atp5mpl", "Atp5o", "Atp5pb", "Atpaf1", "Atpaf2", "Atpif1", "Atpsckmt", "Dmac2l", "Fmc1", "mt-Atp6", "mt-Atp8", "Tmem70")
)

message("Fetching genes annotated to GO:0006007.")
go_query <- AnnotationDbi::select(
  org.Mm.eg.db::org.Mm.eg.db,
  keys = "GO:0006007",
  columns = "SYMBOL",
  keytype = "GOALL"
)
go_symbols <- unique(go_query$SYMBOL[!is.na(go_query$SYMBOL)])
if (length(go_symbols) == 0L) {
  stop("No mouse gene symbols were found for GO:0006007.", call. = FALSE)
}
my_gene_lists[["Glucose_Metabolism"]] <- go_symbols

# Calculate module scores -----------------------------------------------------

message("Calculating pathway module scores.")
score_columns <- gsub("[^A-Za-z0-9_]", "", gsub(" ", "_", names(my_gene_lists)))

for (i in seq_along(my_gene_lists)) {
  clean_name <- score_columns[[i]]
  visium_merged <- AddModuleScore(
    object = visium_merged,
    features = list(my_gene_lists[[i]]),
    assay = "Spatial",
    name = clean_name,
    ctrl = 100
  )

  # AddModuleScore appends 1 for a single gene set. Keep the unnumbered column.
  default_name <- paste0(clean_name, "1")
  visium_merged[[clean_name]] <- visium_merged[[]][[default_name]]
  visium_merged[[default_name]] <- NULL
  message("Added module score: ", clean_name)
}

# Visualize module scores across tissue sections -------------------------------

if (!dir.exists(OUT_DIR)) {
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
}
if (!dir.exists(OUT_DIR)) {
  stop("Cannot create the output directory. Check VISIUM_SPATIAL_OUTPUT_DIR.", call. = FALSE)
}

score_palette <- rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral"))

for (score_col in score_columns) {
  message("Plotting ", score_col, ".")
  p_spatial <- SpatialFeaturePlot(
    visium_merged,
    features = score_col,
    ncol = N_COLUMNS,
    pt.size.factor = 1.3,
    image.alpha = 0,
    alpha = c(0.05, 1),
    min.cutoff = SCORE_LIMITS[1],
    max.cutoff = SCORE_LIMITS[2],
    stroke = 0
  ) +
    plot_annotation(
      title = paste("Spatial Distribution of", gsub("_", " ", score_col), "Module Score"),
      theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
    ) &
    scale_fill_gradientn(
      colors = score_palette,
      breaks = SCORE_BREAKS,
      limits = SCORE_LIMITS,
      na.value = "white"
    )

  out_file <- file.path(OUT_DIR, paste0("Spatial_Distribution_", score_col, "_ModuleScore.pdf"))
  ggsave(
    out_file,
    plot = p_spatial,
    width = SPATIAL_WIDTH,
    height = SPATIAL_HEIGHT,
    limitsize = FALSE
  )
}

message("Saved spatial plots to ", OUT_DIR, ".")
