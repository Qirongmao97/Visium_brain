# --- 1. SETTINGS ---
library(dplyr)
library(ggplot2)
library(openxlsx)
# --- 1. SETTINGS ---
load("/exports/archive/hg-exon-skip/Qirong/DEG_pairwise_MAST/With_cngeneson/ALL/ALL_MAST_results.RData") 
# Assuming your data is loaded into DEG_with_cngenson
DEG_with_cngenson <- result_list
MAST_RESULT_OBJECT <- DEG_with_cngenson 

# List ALL the regions present in your dataset here
ALL_AVAILABLE_REGIONS <- unique(wildtype_2_filtered$annotation_hlevel)

# List the exact names of any regions you want to EXCLUDE
EXCLUDED_REGIONS <- c("Hindbrain","Meninge","Ventricle","Unknown") 

# Add your custom X-axis order here
DESIRED_REGION_ORDER <- c("Isocortex", "Hippocampus Formation", "Thalamus", "Amygdala", "Central subpallium", "Fiber tract")

TARGET_REGIONS <- setdiff(ALL_AVAILABLE_REGIONS, EXCLUDED_REGIONS)
message("Regions to be analyzed: ", paste(TARGET_REGIONS, collapse = ", "))

# Standard Settings 
DESIRED_ORDER <- c("WT vs mdxbl6", "WT vs mdx4cv", "mdxbl6 vs mdx4cv")

# Added cutoffs back to determine Significance in the dotplot
FDR_CUTOFF <- 0.05
LFC_CUTOFF <- 0.1

# CRITICAL FIX: EXACT HEATMAP DIMENSIONS
USER_HEATMAP_WIDTH <- 2.62
USER_HEATMAP_HEIGHT <- 3.81

USER_MATRIX_WIDTH <- 16 
USER_MATRIX_HEIGHT <- 8

# --- 2. REQUIRED HELPER FUNCTIONS ---
get_gene_stats <- function(reg, comp, mast_data, target_genes) {
  deg_list <- mast_data$deg
  p1 <- ifelse(strsplit(comp, " vs ")[[1]][1] == "WT", "wildtype|WT", strsplit(comp, " vs ")[[1]][1])
  p2 <- ifelse(strsplit(comp, " vs ")[[1]][2] == "WT", "wildtype|WT", strsplit(comp, " vs ")[[1]][2])
  
  safe_reg <- gsub("[ /]", "[_\\\\.]", reg)
  reg_pattern <- paste0("(^|[_\\\\.])", safe_reg, "([_\\\\.]|$)")
  
  idx <- which(grepl(reg_pattern, names(deg_list)) & grepl(p1, names(deg_list)) & grepl(p2, names(deg_list)))
  if(length(idx) == 0) return(NULL)
  
  df <- as.data.frame(deg_list[[idx[1]]])
  cols_lower <- tolower(colnames(df))
  gene_col <- colnames(df)[cols_lower %in% c("primerid", "gene", "rowname", "symbol", "id")][1]
  fdr_col  <- colnames(df)[cols_lower %in% c("fdr", "adjp", "p_val_adj", "p.adj", "padj")][1]
  lfc_col  <- colnames(df)[cols_lower %in% c("coef", "logfc", "avg_log2fc")][1]
  
  if ("component" %in% cols_lower) {
    comp_col <- colnames(df)[cols_lower == "component"][1]
    df <- df[df[[comp_col]] %in% c("H", "logFC"), , drop = FALSE]
  }
  
  df_sub <- df[df[[gene_col]] %in% target_genes, ]
  if (nrow(df_sub) == 0) return(NULL)
  
  res <- data.frame(
    Gene = df_sub[[gene_col]],
    Region = reg,
    Comparison = comp,
    LFC = as.numeric(df_sub[[lfc_col]]),
    FDR = as.numeric(df_sub[[fdr_col]])
  )
  return(res)
}

# --- 3. GENE LIST SETUP (ETC + GO TERM) ---
my_gene_lists <- list(
  "Complex I" = c("Acad9", "Aifm1", "Dmac1", "Dmac2", "Ecsit", "Foxred1", "Lyrm2", "mt-Nd1", "mt-Nd2", "mt-Nd3", "mt-Nd4", "mt-Nd4l", "mt-Nd5", "mt-Nd6", "Ndufa1", "Ndufa10", "Ndufa11", "Ndufa12", "Ndufa13", "Ndufa2", "Ndufa3", "Ndufa5", "Ndufa6", "Ndufa7", "Ndufa8", "Ndufa9", "Ndufab1", "Ndufaf1", "Ndufaf2", "Ndufaf3", "Ndufaf4", "Ndufaf5", "Ndufaf6", "Ndufaf7", "Ndufaf8", "Ndufb10", "Ndufb11", "Ndufb2", "Ndufb3", "Ndufb4", "Ndufb5", "Ndufb6", "Ndufb7", "Ndufb8", "Ndufb9", "Ndufc1", "Ndufc2", "Ndufs1", "Ndufs2", "Ndufs3", "Ndufs4", "Ndufs5", "Ndufs6", "Ndufs7", "Ndufs8", "Ndufv1", "Ndufv2", "Ndufv3", "Nubpl", "Timmdc1", "Tmem126a", "Tmem126b", "Tmem186", "Tmem70"),
  "Complex II" = c("Sdha", "Sdhaf1", "Sdhaf2", "Sdhaf3", "Sdhaf4", "Sdhb", "Sdhc", "Sdhd"),
  "Complex III" = c("Bcs1l", "Cyc1", "Lyrm7", "mt-Cytb", "Ttc19", "Uqcc1", "Uqcc2", "Uqcc3", "Uqcr10", "Uqcr11", "Uqcrb", "Uqcrc1", "Uqcrc2", "Uqcrfs1", "Uqcrh", "Uqcrq"),
  "Complex IV" = c("Cep89", "Cmc1", "Cmc2", "Coa3", "Coa4", "Coa5", "Coa6", "Coa7", "Coa8", "Cox10", "Cox11", "Cox14", "Cox15", "Cox16", "Cox17", "Cox18", "Cox19", "Cox20", "Cox4i1", "Cox4i2", "Cox5a", "Cox5b", "Cox6a1", "Cox6a2", "Cox6b1", "Cox6b2", "Cox6c", "Cox7a1", "Cox7a2", "Cox7a2l", "Cox7b", "Cox7c", "Cox8a", "Cox8b", "Cox8c", "Higd1a", "mt-Co1", "mt-Co2", "mt-Co3", "Ndufa4", "Pet100", "Pet117", "Pnkd", "Sco1", "Sco2", "Smim20", "Surf1", "Taco1", "Timm21", "Tmem177"),
  "Complex V" = c("Atp5a1", "Atp5b", "Atp5c1", "Atp5d", "Atp5e", "Atp5g1", "Atp5g2", "Atp5g3", "Atp5h", "Atp5j", "Atp5j2", "Atp5k", "Atp5l", "Atp5md", "Atp5mpl", "Atp5o", "Atp5pb", "Atpaf1", "Atpaf2", "Atpif1", "Atpsckmt", "Dmac2l", "Fmc1", "mt-Atp6", "mt-Atp8", "Tmem70"))

# Dynamically fetch genes for GO:0006007
message("\n>>> Fetching genes for GO:0006007 from org.Mm.eg.db...")
go_query <- AnnotationDbi::select(org.Mm.eg.db, keys = "GO:0006007", columns = "SYMBOL", keytype = "GOALL")
go_symbols <- unique(go_query$SYMBOL[!is.na(go_query$SYMBOL)])
my_gene_lists[["GO:0006007 (Glucose Metab)"]] <- go_symbols
message(">>> Successfully mapped ", length(go_symbols), " genes to GO:0006007.")

# --- 4. EXTRACT ALL RAW DATA ---
message("\n-----------------------------------------")
message("Extracting all available gene stats...")

raw_data_list <- list()

for (complex_name in names(my_gene_lists)) {
  current_genes <- my_gene_lists[[complex_name]]
  for (comp in DESIRED_ORDER) {
    for (reg in TARGET_REGIONS) {
      stats <- get_gene_stats(reg, comp, MAST_RESULT_OBJECT, current_genes)
      if (!is.null(stats)) {
        stats$Complex <- complex_name
        raw_data_list[[paste(complex_name, reg, comp)]] <- stats
      }
    }
  }
}

raw_df <- bind_rows(raw_data_list)

if (nrow(raw_df) > 0) {
  
  # Set universal factor levels to enforce plotting orders
  raw_df$Comparison <- factor(raw_df$Comparison, levels = DESIRED_ORDER)
  raw_df$Region <- factor(raw_df$Region, levels = DESIRED_REGION_ORDER)
  raw_df$Complex <- factor(raw_df$Complex, levels = names(my_gene_lists))
  
  # --- 5A. CALCULATE UNFILTERED AVERAGES FOR HEATMAP ---
  heatmap_df <- raw_df %>%
    group_by(Complex, Comparison, Region) %>%
    summarise(
      Avg_LFC = mean(LFC, na.rm = TRUE),
      Total_Genes_Found = n(),
      .groups = "drop"
    )
  
  write.csv(heatmap_df, "Pathway_Global_Averages.csv", row.names = FALSE)
  
  # Reverse complex order for Heatmap Y-Axis
  heatmap_df$Complex <- factor(heatmap_df$Complex, levels = rev(names(my_gene_lists)))
  
  # Generate Heatmap stacked vertically to fit 2.62in x 3.81in
  p_heat <- ggplot(heatmap_df, aes(x = Region, y = Complex, fill = Avg_LFC)) +
    geom_tile() + 
    scale_fill_gradient2(
      low = "#313695", mid = "white", high = "#a50026", 
      midpoint = 0, name = "Avg LFC" # Shortened legend name for small width
    ) +
    # CRITICAL FIX: Changed to vertical stacking (Comparison ~ .) to fit taller dimensions
    facet_grid(Comparison ~ ., scales = "free_y", space = "free") +
    theme_bw(base_size = 7) + # Scaled down base font size drastically for small canvas
    theme(
      plot.title = element_text(face = "bold", size = 8, hjust = 0.5), 
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 6), 
      axis.text.y = element_text(face = "bold", size = 6),
      axis.title = element_blank(), # Removed axis titles to save space on a 2.6in plot
      strip.background = element_rect(fill = "#f0f0f0", color = "black", linewidth = 0.3),
      strip.text = element_text(face = "bold", size = 6, margin = margin(2, 2, 2, 2)), 
      panel.spacing = unit(0, "lines"), 
      panel.grid = element_blank(),
      legend.position = "bottom", 
      legend.key.width = unit(0.4, "cm"), # Smaller legend key
      legend.key.height = unit(0.2, "cm"),
      legend.title = element_text(size = 6, vjust = 0.8),
      legend.text = element_text(size = 5),
      legend.margin = margin(t = 2), 
      plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "pt") 
    ) +
    labs(title = "Average LFC per Pathway")
  
  out_heat <- "Panel_Pathway_Heatmap_Vertical_Small.pdf"
  ggsave(out_heat, plot = p_heat, width = USER_HEATMAP_WIDTH, height = USER_HEATMAP_HEIGHT, device = "pdf")
  message(">>> Generated Vertical Heatmap: ", out_heat)
  
  
  # --- 5B. GENERATE BOXPLOT + DOTPLOT MATRIX (ALL GENES COLORED BY SIGNIFICANCE) ---
  
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
  
  out_matrix <- "Panel_Pathway_Boxplot_Dotplot_Matrix_Colored_Significance.pdf"
  ggsave(out_matrix, plot = p_box_dot, width = USER_MATRIX_WIDTH, height = USER_MATRIX_HEIGHT, device = "pdf")
  message(">>> Generated Boxplot Matrix: ", out_matrix)
  
} else {
  message("No targeted genes found in the MAST results.")
}

# --- 6. EXPORT RAW DATA TO EXCEL (SUPPLEMENTARY) ---


message("\n-----------------------------------------")
message("Exporting raw_df to Excel for supplementary material...")

# Define the output file name
out_xlsx <- "Supplementary_Table_Raw_Gene_Stats.xlsx"

# Write the dataframe to an Excel file
write.xlsx(raw_df, out_xlsx, rowNames = FALSE)

message(">>> Successfully exported raw data to: ", out_xlsx)
