# =============================================================================
# SpxA1 and SpxA2 Function as a Stoichiometry-Dependent Regulatory Rheostat
# Governing Virulence Gene Expression in Group A Streptococcus
#
# Flores Streptococcal Laboratory
# Division of Pediatric Infectious Diseases
# Vanderbilt University Medical Center
#
# Authors:  Anthony R. Flores, Matthew A. Sanson, Luis A. Vega [et al.]
# Contact:  flores-lab.org
# Journal:  mBio (in submission)
# GEO:      PRJNA1472884 [accession pending]
# GitHub:   https://github.com/flores-lab/SpxA1-SpxA2-MultiOmic-GAS
#
# SAMPLE NAMING CONVENTIONS USED THROUGHOUT:
#   WT          = wild-type MGAS10870 (emm3)
#   spxA1KO     = isogenic ΔspxA1 deletion mutant
#   spxA2KO     = isogenic ΔspxA2 deletion mutant
#   spxA2OE     = LiaSQ146A constitutively active LiaS mutant
#                 (drives high SpxA2; used as SpxA2 overexpression surrogate)
#   UI          = uninduced (vehicle control)
#   BAC         = sub-MIC bacitracin
#   hNP1        = sub-MIC human neutrophil peptide-1
#
# WORKING DIRECTORY:
#   Set BASE_DIR to the root of the repository clone before running.
#   All input and output paths are relative to BASE_DIR.
# =============================================================================

# =============================================================================
# MODULE 1: NanoString nCounter Transcriptional Profiling
# =============================================================================
# Description:
#   - Import nSolver-normalized counts (119 probes: 116 endogenous + 3 HK)
#   - PCA quality control; exclusion of spxA1KO hNP1 replicate 1 (PCA outlier)
#   - limma empirical Bayes DE across 12 contrasts
#     (6 within-strain: BAC/UI, hNP1/UI per strain)
#     (6 cross-strain: spxA1KO/WT and spxA2KO/WT per condition)
#   - Hierarchical clustering -> k=4 gene modules
#   - Heatmap and module profile plot generation
#   - Export: supplemental tables, Prism-ready data, R workspace
#
# Input files (place in data/ subdirectory):
#   Combo_counts-norm_05-2026.csv     — nSolver-normalized counts (macrophage
#                                       columns present but excluded in script)
#
# Output files:
#   ST1_DE_tables_noMo.xlsx           — all 12 contrast DE tables
#   NanoString_normalized_counts.txt  — GEO submission file (27 samples)
#   NanoString_gene_modules.txt       — GEO submission file (k=4 modules)
#   S3_Module_heatmap_k4_noMo.pdf     — Figure
#   S4_Module_profiles_noMo.pdf       — Figure
#   NanoString_workspace_noMo.RData   — workspace for Module 2
# =============================================================================

library(limma)
library(readxl)
library(openxlsx)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(tidyr)

# =============================================================================
# SECTION 0: PATHS — set BASE_DIR to repository root before running
# =============================================================================
BASE_DIR  <- "."                                   # change if needed
DATA_DIR  <- file.path(BASE_DIR, "data")
OUT_DIR   <- file.path(BASE_DIR, "NanoString_RNAseq", "output")
dir.create(OUT_DIR, showWarnings=FALSE, recursive=TRUE)


# ── PATHS — EDIT THESE ────────────────────────────────────────────────────────
data_dir <- "/path/to/your/data"   # directory containing input files
outdir   <- file.path(data_dir, "new_noMo")
dir.create(outdir, showWarnings = FALSE)

# ── COLOUR PALETTE ────────────────────────────────────────────────────────────
strain_colors <- c(
  WT   = "#2166AC",
  A1KO = "#D6604D",
  A2KO = "#4DAC26"
)
module_colors <- c(
  M1 = "#E41A1C",
  M2 = "#377EB8",
  M3 = "#4DAC26",
  M4 = "#984EA3"
)

# =============================================================================
# MODULE 1.1 — IMPORT AND PARSE NORMALIZED COUNTS
# =============================================================================

counts_raw <- read.csv(
  file.path(data_dir, file.path(DATA_DIR, "Combo_counts-norm_05-2026.csv")),
  row.names = 1,
  check.names = FALSE
)

cat("Raw count matrix:", nrow(counts_raw), "probes x",
    ncol(counts_raw), "samples\n")

# Build sample metadata from column names
# Expected format: Strain_Condition_Rep (e.g. WT_UI_A, A2KO_BAC_B)
metadata <- data.frame(
  SampleID  = colnames(counts_raw),
  stringsAsFactors = FALSE
)
metadata$Strain    <- sub("_(UI|BAC|hNP|Mo)_[ABC]$", "", metadata$SampleID)
metadata$Condition <- sub(".*_(UI|BAC|hNP|Mo)_[ABC]$", "\\1", metadata$SampleID)
metadata$Replicate <- sub(".*_([ABC])$", "\\1", metadata$SampleID)

cat("Samples by strain:\n")
print(table(metadata$Strain))
cat("Samples by condition:\n")
print(table(metadata$Condition))

# ── PROBE ANNOTATION ──────────────────────────────────────────────────────────
# Annotation stored in count file or separately; adapt as needed.
# Expected columns: ProbeID, Gene, Label, Annotation, M1_specific
# M1_specific = TRUE for probes targeting emm1-specific genes absent in emm3

# Minimum required: ProbeID and Label (display name)
# Build from rownames if no annotation file available:
annot <- data.frame(
  ProbeID    = rownames(counts_raw),
  Gene       = rownames(counts_raw),
  Label      = rownames(counts_raw),
  Annotation = NA_character_,
  M1_specific = FALSE,
  stringsAsFactors = FALSE
)
# NOTE: Replace with your actual annotation table.
# Known M1-specific probes (absent in MGAS10870 emm3 background):
m1_specific_probes <- c(
  "M5005_Spy_0106c",  # rofA
  "M5005_Spy_0107",   # cpa
  "M5005_Spy_0109",   # prtF
  "M5005_Spy_0114",   # srtB
  "M5005_Spy_0356c",  # speJ
  "M5005_Spy_0561",   # epf
  "M5005_Spy_0805",   # srtK
  "M5005_Spy_1169",   # spd3
  "M5005_Spy_1415c",  # sdaD2
  "M5005_Spy_1702",   # smeZ
  "M5005_Spy_1718c"   # sic1
)
annot$M1_specific[annot$ProbeID %in% m1_specific_probes] <- TRUE

# Probes targeting spxA1 and spxA2 themselves — excluded from DE analysis
ko_probes <- c("M5005_Spy_0959c",   # spxA1
               "M5005_Spy_1798c")   # spxA2

# =============================================================================
# MODULE 1.2 — EXCLUDE MACROPHAGE CO-CULTURE SAMPLES
# =============================================================================

metadata_no_mo <- metadata[metadata$Condition != "Mo", ]
counts_analysis_nomo <- counts_raw[, metadata_no_mo$SampleID]

cat("\nAfter Mo exclusion:", ncol(counts_analysis_nomo), "samples\n")
cat("Conditions retained:",
    paste(unique(metadata_no_mo$Condition), collapse = ", "), "\n")

# =============================================================================
# MODULE 1.3 — QUALITY CONTROL
# =============================================================================

# Log2 transform for QC (nSolver counts already normalized)
counts_filt      <- counts_analysis_nomo
counts_filt_log2 <- log2(counts_filt + 1)

# ── Identify zero-variance probes ────────────────────────────────────────────
gene_vars  <- apply(counts_filt_log2, 1, var)
zero_var   <- names(gene_vars[gene_vars == 0])
cat("\nZero-variance probes (excluded from PCA only):\n")
print(zero_var)
# Note: spd3 (M5005_Spy_1169) shows zero variance in UI/BAC/hNP dataset —
# expression is Mo-specific; probe retained in DE analysis matrix.

# ── PCA ──────────────────────────────────────────────────────────────────────
pca_data   <- t(counts_filt_log2[!rownames(counts_filt_log2) %in% zero_var, ])
pca_result <- prcomp(pca_data, center = TRUE, scale. = FALSE)
pca_var    <- summary(pca_result)$importance
pc1_pct    <- round(pca_var[2, 1] * 100, 1)
pc2_pct    <- round(pca_var[2, 2] * 100, 1)

pca_df           <- as.data.frame(pca_result$x[, 1:3])
pca_df$SampleID  <- rownames(pca_df)
pca_df           <- merge(pca_df, metadata_no_mo, by = "SampleID")
pca_df$Strain    <- factor(pca_df$Strain,
                            levels = c("WT", "A1KO", "A2KO"))
pca_df$Condition <- factor(pca_df$Condition,
                            levels = c("UI", "BAC", "hNP"))
pca_df$Shape     <- c(UI = 16, BAC = 17, hNP = 15)[pca_df$Condition]

p_pca <- ggplot(pca_df,
    aes(x = PC1, y = PC2,
        color = Strain,
        shape = Condition,
        label = SampleID)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text_repel(size = 2.8, max.overlaps = 20,
                   show.legend = FALSE) +
  scale_color_manual(values = strain_colors) +
  scale_shape_manual(values = c(UI = 16, BAC = 17, hNP = 15)) +
  labs(
    title = "NanoString PCA: no-Mo dataset (UI/BAC/hNP)",
    x     = paste0("PC1 (", pc1_pct, "%)"),
    y     = paste0("PC2 (", pc2_pct, "%)")
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(OUT_DIR, "QC_PCA_noMo.pdf"),
       plot = p_pca, width = 8, height = 6,
       useDingbats = FALSE)

# ── OUTLIER EXCLUSION ─────────────────────────────────────────────────────────
# ΔspxA1 hNP-1 replicate A identified as global expression outlier by PCA.
# No technical QC flags; excluded on biological grounds (global separation).
outlier_sample <- "spxA1KO_hNP_A"   # adjust to match your SampleID convention

if (outlier_sample %in% metadata_no_mo$SampleID) {
  metadata_filt <- metadata_no_mo[
    metadata_no_mo$SampleID != outlier_sample, ]
  cat("\nOutlier excluded:", outlier_sample, "\n")
} else {
  metadata_filt <- metadata_no_mo
  warning("Outlier sample not found — check SampleID convention.")
}

counts_filt      <- counts_analysis_nomo[, metadata_filt$SampleID]
counts_filt_log2 <- log2(counts_filt + 1)
cat("Final sample count:", nrow(metadata_filt), "\n")

# Replicate correlation heatmap
cor_mat <- cor(counts_filt_log2, method = "pearson")
pdf(file.path(OUT_DIR, "QC_replicate_correlations_noMo.pdf"),
    width = 10, height = 9)
pheatmap(cor_mat,
  color             = colorRampPalette(
    rev(brewer.pal(9, "RdBu")))(100),
  breaks            = seq(0.85, 1.0, length.out = 101),
  display_numbers   = TRUE,
  number_format     = "%.3f",
  fontsize_number   = 7,
  clustering_method = "complete",
  main              = "Pearson correlations — normalized log2 counts (no Mo)")
dev.off()

# =============================================================================
# MODULE 1.4 — LIMMA EMPIRICAL BAYES DIFFERENTIAL EXPRESSION
# =============================================================================

metadata_filt$Group <- paste(
  metadata_filt$Strain, metadata_filt$Condition, sep = "_")
metadata_filt$Group <- factor(metadata_filt$Group)

design <- model.matrix(~ 0 + Group, data = metadata_filt)
colnames(design) <- gsub("Group", "", colnames(design))

cat("\nDesign matrix columns:\n")
print(colnames(design))

fit <- lmFit(counts_filt_log2, design)

# ── CONTRAST MATRIX ───────────────────────────────────────────────────────────
# 12 valid contrasts (no Mo):
# Within-strain: BAC vs UI and hNP vs UI for each strain
# Cross-strain:  each KO vs WT under UI, BAC, hNP

within_strain_contrasts <- c(
  "WT_BAC_vs_UI"   = "WT_BAC - WT_UI",
  "WT_hNP_vs_UI"   = "WT_hNP - WT_UI",
  "spxspxA1KO_BAC_vs_UI" = "spxA1KO_BAC - A1KO_UI",
  "spxspxA1KO_hNP1_vs_UI" = "spxA1KO_hNP - A1KO_UI",
  "spxspxA2KO_BAC_vs_UI" = "spxA2KO_BAC - A2KO_UI",
  "spxspxA2KO_hNP1_vs_UI" = "spxA2KO_hNP - A2KO_UI"
)

cross_strain_contrasts <- c(
  "spxspxA1KO_vs_WT_UI"  = "spxA1KO_UI  - WT_UI",
  "spxspxA2KO_vs_WT_UI"  = "spxA2KO_UI  - WT_UI",
  "spxspxA1KO_vs_WT_BAC" = "spxA1KO_BAC - WT_BAC",
  "spxspxA2KO_vs_WT_BAC" = "spxA2KO_BAC - WT_BAC",
  "spxspxA1KO_vs_WT_hNP1" = "spxA1KO_hNP - WT_hNP",
  "spxspxA2KO_vs_WT_hNP1" = "spxA2KO_hNP - WT_hNP"
)

all_contrasts   <- c(within_strain_contrasts, cross_strain_contrasts)
contrast_matrix <- makeContrasts(
  contrasts = all_contrasts,
  levels    = design
)
colnames(contrast_matrix) <- names(all_contrasts)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Column labels for module analysis (12 contrasts in order)
col_labels <- c(
  "WT BAC/UI", "WT hNP/UI",
  "A1KO BAC/UI", "A1KO hNP/UI",
  "A2KO BAC/UI", "A2KO hNP/UI",
  "A1KO/WT UI", "A2KO/WT UI",
  "A1KO/WT BAC", "A2KO/WT BAC",
  "A1KO/WT hNP", "A2KO/WT hNP"
)

# ── HELPER: extract DE table for one contrast ─────────────────────────────────
get_full_table <- function(fit2_obj, contrast_name, annot_df,
                           exclude_probes = NULL) {
  tt <- topTable(fit2_obj,
                 coef      = contrast_name,
                 number    = Inf,
                 sort.by   = "P")
  tt$ProbeID   <- rownames(tt)
  tt           <- merge(tt,
                        annot_df[, c("ProbeID", "Label",
                                     "Gene", "Annotation")],
                        by = "ProbeID", all.x = TRUE)
  tt$Sig_p05   <- tt$P.Value < 0.05 &
                  abs(tt$logFC) > log2(1.5)
  tt$Sig_FDR05 <- tt$adj.P.Val < 0.05 &
                  abs(tt$logFC) > log2(1.5)
  tt$Direction <- ifelse(tt$logFC > 0, "Up", "Down")
  tt$Contrast  <- contrast_name
  if (!is.null(exclude_probes)) {
    tt <- tt[!tt$ProbeID %in% exclude_probes, ]
  }
  tt[order(tt$P.Value), ]
}

# =============================================================================
# MODULE 1.5 — SAVE ST1: DIFFERENTIAL EXPRESSION TABLES
# =============================================================================

wb_st1    <- createWorkbook()
header_style <- createStyle(
  fontName  = "Arial",
  fontSize  = 10,
  textDecoration = "bold",
  fgFill    = "#2E4057",
  fontColour = "#FFFFFF",
  halign    = "center"
)

for (cname in names(all_contrasts)) {
  tt <- get_full_table(fit2, cname, annot, ko_probes)
  addWorksheet(wb_st1, cname)
  writeData(wb_st1, cname,
    tt[, c("ProbeID", "Label", "Gene", "Annotation",
           "logFC", "AveExpr", "t", "P.Value",
           "adj.P.Val", "Sig_p05", "Sig_FDR05",
           "Direction", "Contrast")],
    headerStyle = header_style
  )
  freezePane(wb_st1, cname, firstRow = TRUE)
  setColWidths(wb_st1, cname, cols = 1:13, widths = "auto")
}

saveWorkbook(wb_st1,
  file.path(OUT_DIR, "ST1_DE_tables_noMo.xlsx"),
  overwrite = TRUE)
message("ST1 saved: ST1_DE_tables_noMo.xlsx")

# ── SUMMARISE DE COUNTS ───────────────────────────────────────────────────────
de_counts <- sapply(names(all_contrasts), function(cname) {
  tt <- get_full_table(fit2, cname, annot, ko_probes)
  sum(tt$Sig_p05)
})
cat("\nDE gene counts (|log2FC|>log2(1.5), p<0.05):\n")
print(de_counts)

# =============================================================================
# MODULE 1.6 — OVERLAP ANALYSIS (BAC+hNP)
# =============================================================================

get_sig_labels <- function(contrast_name) {
  tt <- get_full_table(fit2, contrast_name, annot, ko_probes)
  tt$Label[tt$Sig_p05]
}

sig_sets <- list(
  WT_BAC   = get_sig_labels("WT_BAC_vs_UI"),
  WT_hNP   = get_sig_labels("WT_hNP_vs_UI"),
  A1KO_BAC = get_sig_labels("spxspxA1KO_BAC_vs_UI"),
  A1KO_hNP = get_sig_labels("spxspxA1KO_hNP1_vs_UI"),
  A2KO_BAC = get_sig_labels("spxspxA2KO_BAC_vs_UI"),
  A2KO_hNP = get_sig_labels("spxspxA2KO_hNP1_vs_UI")
)

overlap_df <- data.frame(
  Strain   = c("WT", "A1KO", "A2KO"),
  BAC_only = c(
    length(setdiff(sig_sets$WT_BAC,   sig_sets$WT_hNP)),
    length(setdiff(sig_sets$A1KO_BAC, sig_sets$A1KO_hNP)),
    length(setdiff(sig_sets$A2KO_BAC, sig_sets$A2KO_hNP))
  ),
  Both     = c(
    length(intersect(sig_sets$WT_BAC,   sig_sets$WT_hNP)),
    length(intersect(sig_sets$A1KO_BAC, sig_sets$A1KO_hNP)),
    length(intersect(sig_sets$A2KO_BAC, sig_sets$A2KO_hNP))
  ),
  hNP_only = c(
    length(setdiff(sig_sets$WT_hNP,   sig_sets$WT_BAC)),
    length(setdiff(sig_sets$A1KO_hNP, sig_sets$A1KO_BAC)),
    length(setdiff(sig_sets$A2KO_hNP, sig_sets$A2KO_BAC))
  )
)
cat("\nBAC / hNP overlap by strain:\n")
print(overlap_df)

# =============================================================================
# MODULE 1.7 — GENE MODULE ANALYSIS (k=4 HIERARCHICAL CLUSTERING)
# =============================================================================

# Build logFC matrix: one row per probe, one column per contrast
logfc_matrix <- do.call(cbind, lapply(names(all_contrasts), function(cname) {
  tt <- get_full_table(fit2, cname, annot, ko_probes)
  setNames(tt$logFC, tt$Label)
}))
colnames(logfc_matrix) <- col_labels

# Filter: keep genes with |logFC| > log2(1.5) in >= 2 contrasts
n_sig          <- rowSums(abs(logfc_matrix) > log2(1.5))
logfc_filtered <- logfc_matrix[n_sig >= 2, ]

cat("\nGenes passing module filter:", nrow(logfc_filtered), "\n")

# Fix dagger characters for PDF rendering
rownames(logfc_filtered) <- gsub("\u2020", "*",
                                  rownames(logfc_filtered))

# Hierarchical clustering — Euclidean distance, complete linkage
hc_genes <- hclust(
  dist(logfc_filtered, method = "euclidean"),
  method = "complete"
)

# k=4 selected based on biological coherence;
# k=5 evaluated and rejected (magnitude gradient split only within M2)
set.seed(42)
gene_clusters_k4 <- cutree(hc_genes, k = 4)

cat("\nModule sizes (k=4):\n")
print(table(paste0("M", gene_clusters_k4)))

# Build cluster summary data frame
cluster_df_k4 <- data.frame(
  Label   = names(gene_clusters_k4),
  Cluster = paste0("M", gene_clusters_k4),
  logfc_filtered,
  check.names = FALSE
)

# Module mean profiles
cluster_means_k4 <- cluster_df_k4 %>%
  group_by(Cluster) %>%
  summarise(across(all_of(col_labels), mean), .groups = "drop")

cat("\nModule mean logFC profiles:\n")
print(as.data.frame(cluster_means_k4))

# Module lookup: Label → module
module_lookup <- setNames(
  paste0("M", gene_clusters_k4),
  names(gene_clusters_k4)
)

# ── ANNOTATION OBJECTS FOR HEATMAP ────────────────────────────────────────────
anno_row_k4 <- data.frame(
  Cluster   = paste0("M", gene_clusters_k4),
  row.names = names(gene_clusters_k4)
)
rownames(anno_row_k4) <- gsub("\u2020", "*", rownames(anno_row_k4))

anno_col_module <- data.frame(
  Condition = gsub(".*/UI$|.*/BAC$|.*/hNP$",
                   function(x) sub(".*/", "", x),
                   col_labels),
  Type      = ifelse(grepl("/WT", col_labels),
                      "Cross-strain", "Within-strain"),
  row.names = col_labels
)

anno_colors_module <- list(
  Cluster   = module_colors,
  Type      = c(`Cross-strain`  = "#5C3A1E",
                `Within-strain` = "#2C5F8A")
)

anno_colors_k4 <- anno_colors_module

# ── HEATMAP: k=4 ─────────────────────────────────────────────────────────────
pdf(file.path(OUT_DIR, "S3_Module_heatmap_k4_noMo.pdf"),
    width = 11, height = 13, useDingbats = FALSE)
pheatmap(logfc_filtered,
  annotation_col    = anno_col_module,
  annotation_row    = anno_row_k4,
  annotation_colors = anno_colors_k4,
  color             = colorRampPalette(
    rev(brewer.pal(9, "RdBu")))(100),
  breaks            = seq(-3, 3, length.out = 101),
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method        = "complete",
  show_rownames = TRUE, show_colnames = TRUE,
  fontsize = 8, fontsize_row = 6, fontsize_col = 8,
  main = "Gene modules: logFC profiles (k=4, no Mo)")
dev.off()

# ── MODULE PROFILE PLOT ───────────────────────────────────────────────────────
cluster_means_long <- cluster_means_k4 %>%
  pivot_longer(cols = all_of(col_labels),
               names_to  = "Contrast",
               values_to = "mean_logFC")
cluster_means_long$Contrast   <- factor(
  cluster_means_long$Contrast, levels = col_labels)
cluster_means_long$Comparison <- ifelse(
  grepl("/WT", cluster_means_long$Contrast),
  "Cross-strain", "Within-strain")

p_profiles <- ggplot(cluster_means_long,
    aes(x = Contrast, y = mean_logFC,
        color = Cluster, group = Cluster)) +
  geom_hline(yintercept = 0, color = "grey70",
             linetype = "dashed", linewidth = 0.4) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_color_manual(values = module_colors) +
  facet_wrap(~ Comparison, scales = "free_x", nrow = 1) +
  labs(title = "Gene module mean logFC profiles (k=4, no Mo)",
       x = "", y = "Mean log\u2082 Fold Change",
       color = "Module") +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position    = "bottom",
    strip.background   = element_rect(fill = "#F0F0F0")
  )

ggsave(file.path(OUT_DIR, "S4_Module_profiles_noMo.pdf"),
       plot = p_profiles, width = 10, height = 5,
       useDingbats = FALSE)

# =============================================================================
# MODULE 1.8 — WITHIN-STRAIN AND CROSS-STRAIN COMPARISON TABLES
# =============================================================================

# Comparisons for scatter plots and figure exports
build_combined_table <- function(tt_A, tt_B, suffix_A, suffix_B,
                                  by_col = "Label") {
  shared_cols <- c(by_col, "ProbeID", "logFC", "P.Value", "Sig_p05")
  A <- tt_A[, shared_cols]
  B <- tt_B[, shared_cols]
  colnames(A)[3:5] <- paste0(c("logFC_", "P.Value_", "Sig_p05_"),
                               suffix_A)
  colnames(B)[3:5] <- paste0(c("logFC_", "P.Value_", "Sig_p05_"),
                               suffix_B)
  merged <- merge(A, B, by = c(by_col, "ProbeID"))
  merged$Pattern <- case_when(
    merged[[paste0("Sig_p05_", suffix_A)]] &
      merged[[paste0("Sig_p05_", suffix_B)]]  ~ paste0("Both"),
    merged[[paste0("Sig_p05_", suffix_A)]] &
      !merged[[paste0("Sig_p05_", suffix_B)]] ~ suffix_A,
    !merged[[paste0("Sig_p05_", suffix_A)]] &
      merged[[paste0("Sig_p05_", suffix_B)]]  ~ suffix_B,
    TRUE ~ "NS"
  )
  merged
}

# BAC: WT vs A2KO (within-strain logFC for scatter Figure A panels B/C)
bac_wt   <- get_full_table(fit2, "WT_BAC_vs_UI",   annot, ko_probes)
bac_a1ko <- get_full_table(fit2, "spxspxA1KO_BAC_vs_UI", annot, ko_probes)
bac_a2ko <- get_full_table(fit2, "spxspxA2KO_BAC_vs_UI", annot, ko_probes)
bac_combined <- build_combined_table(bac_wt, bac_a2ko, "WT", "A2KO")

hnp_wt   <- get_full_table(fit2, "WT_hNP_vs_UI",   annot, ko_probes)
hnp_a1ko <- get_full_table(fit2, "spxspxA1KO_hNP1_vs_UI", annot, ko_probes)
hnp_a2ko <- get_full_table(fit2, "spxspxA2KO_hNP1_vs_UI", annot, ko_probes)
hnp_combined <- build_combined_table(hnp_wt, hnp_a2ko, "WT", "A2KO")

# Cross-strain UI comparisons
ui_a1ko <- get_full_table(fit2, "spxspxA1KO_vs_WT_UI", annot, ko_probes)
ui_a2ko <- get_full_table(fit2, "spxspxA2KO_vs_WT_UI", annot, ko_probes)
ui_combined <- build_combined_table(ui_a1ko, ui_a2ko, "A1KO", "A2KO")

# Significance sets for overlap analysis
bac_sig_any <- union(
  bac_wt$Label[bac_wt$Sig_p05],
  union(bac_a1ko$Label[bac_a1ko$Sig_p05],
        bac_a2ko$Label[bac_a2ko$Sig_p05])
)
hnp_sig_any <- union(
  hnp_wt$Label[hnp_wt$Sig_p05],
  union(hnp_a1ko$Label[hnp_a1ko$Sig_p05],
        hnp_a2ko$Label[hnp_a2ko$Sig_p05])
)

# =============================================================================
# MODULE 1.9 — PRISM DATA EXPORT
# =============================================================================

wb_prism <- createWorkbook()

# Figure A Panel A: DE counts
de_counts_out <- data.frame(
  Strain = c("WT", "A1KO", "A2KO"),
  BAC    = c(sum(bac_wt$Sig_p05),
              sum(bac_a1ko$Sig_p05),
              sum(bac_a2ko$Sig_p05)),
  hNP    = c(sum(hnp_wt$Sig_p05),
              sum(hnp_a1ko$Sig_p05),
              sum(hnp_a2ko$Sig_p05))
)
addWorksheet(wb_prism, "FigA_A_DE_counts")
writeData(wb_prism, "FigA_A_DE_counts",
          de_counts_out, headerStyle = header_style)

# Figure A Panels B/C: WT vs A2KO scatter (BAC and hNP)
bac_scatter <- merge(
  bac_combined[, c("Label", "ProbeID",
                    "logFC_WT", "logFC_A2KO",
                    "Sig_p05_WT", "Sig_p05_A2KO", "Pattern")],
  data.frame(Label = names(module_lookup),
             Module = module_lookup),
  by = "Label", all.x = TRUE
)
bac_scatter$Module[is.na(bac_scatter$Module)] <- "NS"
bac_scatter <- bac_scatter[!bac_scatter$ProbeID %in% ko_probes, ]
bac_scatter <- bac_scatter[order(bac_scatter$Module, bac_scatter$Label), ]

addWorksheet(wb_prism, "FigA_B_BAC_WT_vs_A2KO")
writeData(wb_prism, "FigA_B_BAC_WT_vs_A2KO",
          bac_scatter, headerStyle = header_style)
freezePane(wb_prism, "FigA_B_BAC_WT_vs_A2KO", firstRow = TRUE)

hnp_scatter <- merge(
  hnp_combined[, c("Label", "ProbeID",
                    "logFC_WT", "logFC_A2KO",
                    "Sig_p05_WT", "Sig_p05_A2KO", "Pattern")],
  data.frame(Label = names(module_lookup),
             Module = module_lookup),
  by = "Label", all.x = TRUE
)
hnp_scatter$Module[is.na(hnp_scatter$Module)] <- "NS"
hnp_scatter <- hnp_scatter[!hnp_scatter$ProbeID %in% ko_probes, ]
hnp_scatter <- hnp_scatter[order(hnp_scatter$Module, hnp_scatter$Label), ]

addWorksheet(wb_prism, "FigA_C_hNP_WT_vs_A2KO")
writeData(wb_prism, "FigA_C_hNP_WT_vs_A2KO",
          hnp_scatter, headerStyle = header_style)
freezePane(wb_prism, "FigA_C_hNP_WT_vs_A2KO", firstRow = TRUE)

# Figure A Panel D: overlap
addWorksheet(wb_prism, "FigA_D_overlap")
writeData(wb_prism, "FigA_D_overlap",
          overlap_df, headerStyle = header_style)

# Figure B Panels A/B: A1KO vs A2KO scatter (BAC and hNP)
for (cond in c("BAC", "hNP")) {
  tt_a1 <- get_full_table(
    fit2, paste0("spxA1KO_", cond, "_vs_UI"), annot, ko_probes)
  tt_a2 <- get_full_table(
    fit2, paste0("spxA2KO_", cond, "_vs_UI"), annot, ko_probes)
  ko_scatter <- build_combined_table(tt_a1, tt_a2, "A1KO", "A2KO")
  ko_scatter <- merge(
    ko_scatter,
    data.frame(Label = names(module_lookup), Module = module_lookup),
    by = "Label", all.x = TRUE)
  ko_scatter$Module[is.na(ko_scatter$Module)] <- "NS"
  ko_scatter <- ko_scatter[!ko_scatter$ProbeID %in% ko_probes, ]
  ko_scatter <- ko_scatter[order(ko_scatter$Pattern, ko_scatter$Label), ]
  sheet_name <- paste0("FigB_", ifelse(cond == "BAC", "A", "B"),
                        "_", cond, "_A1KO_vs_A2KO")
  addWorksheet(wb_prism, sheet_name)
  writeData(wb_prism, sheet_name, ko_scatter, headerStyle = header_style)
  freezePane(wb_prism, sheet_name, firstRow = TRUE)
}

# Figure B Panel C: CovR bridge (sagA excluded)
# CovR Class 1 and Class 2 targets — A2KO/WT BAC and hNP comparisons
covr_class1 <- c("sclA", "ska", "codY", "sclB", "prtS", "hasA")
covr_class2 <- c("mga", "emm1.0")

get_lfc <- function(gene, contrast) {
  if (gene %in% rownames(logfc_filtered))
    round(logfc_filtered[gene, contrast], 3)
  else NA
}

figB_C <- data.frame(
  Gene       = c(covr_class1, covr_class2),
  CovR_Class = c(rep("Class_1", length(covr_class1)),
                  rep("Class_2", length(covr_class2))),
  logFC_BAC  = sapply(c(covr_class1, covr_class2),
                       get_lfc, "A2KO/WT BAC"),
  logFC_hNP  = sapply(c(covr_class1, covr_class2),
                       get_lfc, "A2KO/WT hNP"),
  stringsAsFactors = FALSE
)
addWorksheet(wb_prism, "FigB_C_CovR_bridge")
writeData(wb_prism, "FigB_C_CovR_bridge",
          figB_C, headerStyle = header_style)

cat("\nCovR bridge values (sagA excluded):\n")
print(figB_C)

saveWorkbook(wb_prism,
  file.path(OUT_DIR, "Prism_data_final.xlsx"),
  overwrite = TRUE)
message("Prism data saved: Prism_data_final.xlsx")

# =============================================================================
# MODULE 1.10 — SAVE WORKSPACE
# =============================================================================

save(
  counts_raw, counts_analysis_nomo,
  counts_filt, counts_filt_log2,
  metadata, metadata_no_mo, metadata_filt,
  design, contrast_matrix, fit, fit2,
  annot, module_lookup, header_style,
  within_strain_contrasts, cross_strain_contrasts,
  all_contrasts, col_labels, ko_probes,
  logfc_matrix, logfc_filtered,
  hc_genes, gene_clusters_k4,
  cluster_df_k4, cluster_means_k4,
  anno_col_module, anno_colors_module,
  anno_row_k4, module_colors,
  bac_wt, bac_a1ko, bac_a2ko, bac_combined,
  hnp_wt, hnp_a1ko, hnp_a2ko, hnp_combined,
  ui_combined, ui_a1ko, ui_a2ko,
  bac_sig_any, hnp_sig_any,
  sig_sets, overlap_df, de_counts,
  strain_colors, outdir, data_dir,
  file = file.path(OUT_DIR, "NanoString_workspace_noMo.RData")
)
message("Workspace saved: NanoString_workspace_noMo.RData")
message("\nModule 1 complete.")

# =============================================================================
# SESSION INFO
# =============================================================================
cat("\n=== Session Info ===\n")
writeLines(capture.output(sessionInfo()),
           file.path(OUT_DIR, "NanoString_sessionInfo.txt"))
sessionInfo()
