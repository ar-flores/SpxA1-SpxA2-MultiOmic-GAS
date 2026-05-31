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
# Proteomics Analysis — DIA LFQ Differential Abundance
# =============================================================================
# Description:
#   Differential abundance analysis of DIA proteomics data comparing
#   isogenic mutant strains of S. pyogenes MGAS10870 (emm3) to WT.
#
#   Main comparisons (manuscript figures):
#     spxA1KO vs WT, spxA2KO vs WT
#   Supplemental comparisons:
#     clpX_KO vs WT, pA2_A1 vs WT, pA1_A2 vs WT
#
#   Generates: volcano plots, functional category bar charts, FC heatmap,
#              proteome comparison scatter, supplemental tables
#
# Input files (place in data/ subdirectory):
#   Strain_lfq_table.xlsx         — DIA-NN LFQ output (sheet: "DIA analysis")
#   cog_annotation_by_protname.csv — functional annotation (UniProt entry names)
#
# Output files (Proteomics/output/):
#   ST_Prot1-5 supplemental tables (.xlsx)
#   Volcano plots, heatmap, scatter figures (.pdf/.png)
#   proteomics_sessionInfo.txt
# =============================================================================

library(tidyverse)
library(readxl)
library(limma)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(writexl)
library(patchwork)


# =============================================================================
# 1. LOAD DATA
# =============================================================================

lfq_raw <- read_excel(file.path(DATA_DIR, "Strain_lfq_table.xlsx"), sheet = "DIA analysis")

# Separate annotation columns from intensity matrix
meta_cols <- c("PG.Genes", "PG.ProteinDescriptions", "PG.ProteinNames")
annot     <- lfq_raw[, meta_cols]

# Extract intensity matrix, coerce to numeric, replace 0/NaN with NA
mat_raw <- lfq_raw[, !names(lfq_raw) %in% meta_cols]
mat_raw <- as.data.frame(lapply(mat_raw,
                                function(x) as.numeric(as.character(x))))
mat_raw[mat_raw == 0]              <- NA
mat_raw[is.nan(as.matrix(mat_raw))] <- NA

# Log2-transform
mat_log2           <- log2(as.matrix(mat_raw))
rownames(mat_log2) <- annot$PG.ProteinNames   # PG.ProteinNames is the unique key

# Sample grouping — 6 groups x 3 replicates
# Verify column order matches the header comment before running
groups <- factor(
  rep(c("WT", "spxA2", "clpX", "spxA1", "pA2_A1", "pA1_A2"), each = 3),
  levels = c("WT", "spxA1", "spxA2", "clpX", "pA2_A1", "pA1_A2")
)

cat(sprintf("Proteins in input: %d\n", nrow(mat_log2)))
cat(sprintf("Samples: %d  |  Groups: %d x %d replicates\n",
            ncol(mat_log2), nlevels(groups), 3))


# =============================================================================
# 2. NORMALIZATION
# =============================================================================
# LFQ intensities are exported from DIA-NN following internal global median
# normalization across runs. No additional normalization is applied here.
# Adequacy confirmed by within-group replicate correlations (see Step 3).

gene_names_all     <- rownames(mat_log2)
mat_norm           <- mat_log2
rownames(mat_norm) <- gene_names_all


# =============================================================================
# 3. MISSING VALUE FILTER AND QUALITY CONTROL
# =============================================================================

# Retain proteins with >= 2/3 valid values in at least one group
keep <- apply(mat_norm, 1, function(x) {
  counts <- tapply(!is.na(x), groups, sum)
  any(counts >= 2)
})

mat_filt             <- mat_norm[keep, ]
annot_filt           <- annot[keep, ]
rownames(mat_filt)   <- annot_filt$PG.ProteinNames  # re-attach after subsetting

cat(sprintf("Proteins passing 2/3 filter: %d / %d\n", sum(keep), nrow(mat_norm)))
cat(sprintf("Rowname check (first 3): %s\n",
            paste(head(rownames(mat_filt), 3), collapse = ", ")))

# Guard: rownames must be gene/protein identifiers, not integers
stopifnot(
  "Rownames are integers — check data loading and column order" =
    !all(grepl("^[0-9]+$", head(rownames(mat_filt), 10)))
)

# Within-group replicate correlations (quality check)
cat("\nWithin-group replicate correlations (Pearson, log2 LFQ):\n")
for (grp in levels(groups)) {
  grp_cols    <- which(groups == grp)
  sub         <- mat_filt[, grp_cols]
  complete_rows <- rowSums(!is.na(sub)) == ncol(sub)
  if (sum(complete_rows) > 10) {
    r_mat   <- cor(sub[complete_rows, ], use = "complete.obs")
    off_diag <- r_mat[lower.tri(r_mat)]
    cat(sprintf("  %-8s: mean r = %.4f\n", grp, mean(off_diag)))
  }
}


# =============================================================================
# 4. SIGNIFICANCE THRESHOLDS, FUNCTIONAL ANNOTATION, AND LIMMA MODEL
# =============================================================================

# Significance thresholds — consistent with RNA-seq and NanoString analyses
FC_THRESH   <- log2(1.5)  # 0.585; equivalent to 1.5-fold change
PADJ_THRESH <- 0.05

# Load functional annotation keyed on PG.ProteinNames (unique UniProt entry name)
# Built from UniProt proteome UP000000564 (S. pyogenes MGAS315, accessed May 2026)
# with manual overrides for GAS-specific proteins
cog_annot  <- read_csv(file.path(DATA_DIR, "cog_annotation_by_protname.csv"), show_col_types = FALSE)

# Merge annotation into filtered protein table
annot_filt <- annot_filt %>%
  left_join(cog_annot %>%
              select(PG.ProteinNames, COG_letter, COG_label, SpyM3_locus_tag,
                     `Protein names`, `Function [CC]`, Keywords,
                     `Gene Ontology (biological process)`, Reviewed) %>%
              distinct(PG.ProteinNames, .keep_all = TRUE),
            by = "PG.ProteinNames") %>%
  mutate(
    COG_letter     = replace_na(COG_letter, "S"),
    COG_label      = replace_na(COG_label,  "Unknown function")
  )

cat(sprintf("\nCOG annotation merged: %d proteins\n", nrow(annot_filt)))
cat("COG category distribution:\n")
print(sort(table(annot_filt$COG_label), decreasing = TRUE))

# Limma linear model
design <- model.matrix(~ 0 + groups)
colnames(design) <- levels(groups)

fit <- lmFit(mat_filt, design)

contrasts_mat <- makeContrasts(
  spxA1_vs_WT  = spxA1  - WT,
  spxA2_vs_WT  = spxA2  - WT,
  clpX_vs_WT   = clpX   - WT,
  pA2A1_vs_WT  = pA2_A1 - WT,
  pA1A2_vs_WT  = pA1_A2 - WT,
  levels = design
)

fit2 <- contrasts.fit(fit, contrasts_mat)
fit2 <- eBayes(fit2)

comparison_names <- colnames(contrasts_mat)

# Extract results for all comparisons
# topTable returns PG.ProteinNames as rownames; capture via rownames_to_column
results_list <- lapply(comparison_names, function(coef) {
  topTable(fit2, coef = coef, number = Inf, sort.by = "P") %>%
    rownames_to_column("PG.ProteinNames") %>%
    { if ("ID" %in% names(.)) select(., -ID) else . } %>%
    left_join(
      annot_filt %>%
        select(PG.ProteinNames, PG.Genes, PG.ProteinDescriptions),
      by = "PG.ProteinNames"
    ) %>%
    mutate(
      comparison = coef,
      direction  = case_when(
        adj.P.Val < PADJ_THRESH & logFC >= FC_THRESH  ~ "Increased",
        adj.P.Val < PADJ_THRESH & logFC <= -FC_THRESH ~ "Decreased",
        TRUE ~ "NS"
      )
    )
})
names(results_list) <- comparison_names

# Attach COG annotation to results immediately after extraction
results_list <- lapply(results_list, function(res) {
  left_join(
    res,
    annot_filt %>%
      select(PG.ProteinNames, COG_letter, COG_label, SpyM3_locus_tag) %>%
      distinct(PG.ProteinNames, .keep_all = TRUE),
    by = "PG.ProteinNames"
  )
})

cat(sprintf("\nCOG_label present in results: %s\n",
            "COG_label" %in% names(results_list[[1]])))
cat(sprintf("PG.Genes present in results:  %s\n",
            "PG.Genes" %in% names(results_list[[1]])))


# =============================================================================
# 5. SUMMARY STATISTICS
# =============================================================================

summary_stats <- map_dfr(results_list, function(res) {
  sig <- filter(res, adj.P.Val < PADJ_THRESH, abs(logFC) >= FC_THRESH)
  tibble(
    Comparison  = unique(res$comparison),
    N_detected  = nrow(res),
    N_sig       = nrow(sig),
    N_increased = sum(sig$logFC > 0),
    N_decreased = sum(sig$logFC < 0)
  )
})

cat("\n=== SUMMARY STATISTICS ===\n")
print(summary_stats)
write_csv(summary_stats, "proteomics_summary_statistics.csv")


# =============================================================================
# 6. VOLCANO PLOTS
# =============================================================================

# Proteins highlighted from manuscript Figure 4
# SpyM3 locus tags used where standard gene name not in dataset
highlight_spxA1 <- c(
  "sodA",          # Superoxide dismutase
  "cysM",          # Cysteine synthase
  "SpyM3_0212",    # SufB — Fe-S cluster assembly
  "nox.1",         # NADH oxidase
  "dpr",           # Peroxide resistance protein
  "SpyM3_0428",    # Glutathione peroxidase (gpoA)
  "SpyM3_1197",    # CRP/FNR-type regulator (fnr)
  "rex",           # Redox-sensing transcriptional repressor
  "SpyM3_0317",    # Manganese transport regulator (mntR)
  "ahpC",          # Alkyl hydroperoxide reductase
  "nox.2"          # Thioredoxin reductase
)

highlight_spxA2 <- c(
  "speB",          # Cysteine protease
  "SpyM3_0583",    # IgG protease (ideS)
  "ska",           # Streptokinase
  "grab",          # Alpha-2-macroglobulin-binding protein
  "pepO",          # Endopeptidase O
  "sagC",          # Streptolysin S biosynthesis
  "arcA",          # Arginine deiminase
  "prtS",          # Intracellular protease
  "slo",           # Streptolysin O
  "hasA"           # Hyaluronate synthase
)

highlight_clpX <- c(
  "spx (SpxA1)", "SpyM3_1799 (SpxA2)",
  "clpP", "clpX", "clpE", "clpC", "dnaK", "groEL"
)

highlight_map <- list(
  spxA1_vs_WT = highlight_spxA1,
  spxA2_vs_WT = highlight_spxA2,
  clpX_vs_WT  = highlight_clpX,
  pA2A1_vs_WT = unique(c(highlight_spxA1, highlight_spxA2)),
  pA1A2_vs_WT = unique(c(highlight_spxA1, highlight_spxA2))
)

volcano_titles <- list(
  spxA1_vs_WT = "\u0394spxA1 vs. WT",
  spxA2_vs_WT = "\u0394spxA2 vs. WT",
  clpX_vs_WT  = "\u0394clpX vs. WT",
  pA2A1_vs_WT = "pA2_A1 vs. WT (SpxA1 under spxA2 promoter)",
  pA1A2_vs_WT = "pA1_A2 vs. WT (SpxA2 under spxA1 promoter)"
)

make_volcano <- function(res_df, highlights, title,
                         color_up = "#C0392B", color_dn = "#2980B9") {
  res_df <- res_df %>%
    mutate(
      sig   = factor(direction, levels = c("Increased", "Decreased", "NS")),
      label = ifelse(PG.Genes %in% highlights & direction != "NS",
                     PG.Genes, NA_character_)
    )
  n_up   <- sum(res_df$direction == "Increased", na.rm = TRUE)
  n_down <- sum(res_df$direction == "Decreased", na.rm = TRUE)

  ggplot(res_df,
         aes(x = logFC, y = -log10(adj.P.Val), color = sig)) +
    geom_point(data = filter(res_df, sig == "NS"),
               alpha = 0.3, size = 0.9, color = "grey70") +
    geom_point(data = filter(res_df, sig != "NS"),
               alpha = 0.75, size = 1.6) +
    geom_text_repel(aes(label = label), size = 2.8, show.legend = FALSE,
                    max.overlaps = 25, segment.size = 0.3,
                    box.padding = 0.35, point.padding = 0.2,
                    min.segment.length = 0.1) +
    geom_hline(yintercept = -log10(PADJ_THRESH),
               linetype = "dashed", color = "grey40", linewidth = 0.4) +
    geom_vline(xintercept = c(-FC_THRESH, FC_THRESH),
               linetype = "dashed", color = "grey40", linewidth = 0.4) +
    scale_color_manual(
      values = c(Increased = color_up, Decreased = color_dn, NS = "grey70"),
      labels = c(sprintf("Increased (n=%d)", n_up),
                 sprintf("Decreased (n=%d)", n_down),
                 "Not significant"),
      drop = FALSE
    ) +
    labs(
      title = title,
      x     = expression(log[2]~"Fold Change"),
      y     = expression(-log[10]~"Adjusted p-value"),
      color = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position  = "bottom",
      legend.text      = element_text(size = 9),
      panel.grid.minor = element_blank(),
      plot.title       = element_text(size = 11, face = "bold")
    )
}

# Generate all volcano plots
volcano_plots <- lapply(comparison_names, function(comp) {
  make_volcano(results_list[[comp]], highlight_map[[comp]], volcano_titles[[comp]])
})
names(volcano_plots) <- comparison_names

# Main figure — side by side
main_volcano <- volcano_plots$spxA1_vs_WT + volcano_plots$spxA2_vs_WT +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave(file.path(OUT_DIR, "Figure4_volcanos_main.pdf", main_volcano, width = 11, height = 5)

# Supplemental — three panels
supp_volcano <- (volcano_plots$clpX_vs_WT |
                   volcano_plots$pA2A1_vs_WT |
                   volcano_plots$pA1A2_vs_WT) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave("FigureS_volcanos_supplemental.pdf", supp_volcano, width = 16, height = 5)

cat("Volcano plots saved.\n")


# =============================================================================
# 7. FUNCTIONAL CATEGORY BAR CHARTS
# =============================================================================

make_category_chart <- function(comparisons, result_list, title) {
  comp_labels <- c(
    spxA1_vs_WT = "\u0394spxA1 vs WT",
    spxA2_vs_WT = "\u0394spxA2 vs WT",
    clpX_vs_WT  = "\u0394clpX vs WT",
    pA2A1_vs_WT = "pA2_A1 vs WT",
    pA1A2_vs_WT = "pA1_A2 vs WT"
  )

  plot_data <- map_dfr(comparisons, function(comp) {
    sig <- result_list[[comp]] %>%
      filter(adj.P.Val < PADJ_THRESH, abs(logFC) >= FC_THRESH)
    if (nrow(sig) == 0) return(NULL)
    sig %>%
      count(COG_label, direction) %>%
      mutate(
        comparison = comp_labels[comp],
        n_dir      = ifelse(direction == "Decreased", -n, n)
      )
  })

  if (is.null(plot_data) || nrow(plot_data) == 0) {
    warning("No significant proteins found for category chart.")
    return(NULL)
  }

  ggplot(plot_data, aes(x = COG_label, y = n_dir, fill = direction)) +
    geom_col(position = "identity", alpha = 0.85, width = 0.7) +
    facet_wrap(~ comparison, nrow = 1) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
    scale_fill_manual(
      values = c(Increased = "#C0392B", Decreased = "#2980B9"),
      labels = c("Increased abundance", "Decreased abundance")
    ) +
    scale_y_continuous(labels = abs) +
    coord_flip() +
    labs(
      title   = title,
      x       = NULL,
      y       = "Number of proteins",
      fill    = NULL,
      caption = "Bars right of 0 = increased; bars left of 0 = decreased vs. WT"
    ) +
    theme_bw(base_size = 10) +
    theme(
      legend.position  = "bottom",
      strip.text       = element_text(size = 9, face = "bold"),
      panel.grid.minor = element_blank(),
      axis.text.y      = element_text(size = 9)
    )
}

main_cat <- make_category_chart(
  c("spxA1_vs_WT", "spxA2_vs_WT"),
  results_list,
  "Functional categories — main comparisons"
)
if (!is.null(main_cat))
  ggsave("Figure4_functional_categories_main.pdf", main_cat, width = 10, height = 5)

all_cat <- make_category_chart(
  comparison_names,
  results_list,
  "Functional categories — all comparisons"
)
if (!is.null(all_cat))
  ggsave("FigureS_functional_categories_all.pdf"), all_cat, width = 18, height = 5)

cat("Functional category charts saved.\n")


# =============================================================================
# 8. LOG2FC HEATMAP OF HIGHLIGHTED PROTEINS
# =============================================================================

all_highlights <- unique(c(highlight_spxA1, highlight_spxA2))

# Identify which highlights are detected in the results
all_detected      <- unique(unlist(lapply(results_list, function(r) r$PG.Genes)))
detected_highlights <- intersect(all_highlights, all_detected)

cat(sprintf("\nHeatmap: %d / %d highlighted proteins detected\n",
            length(detected_highlights), length(all_highlights)))
missing_h <- setdiff(all_highlights, detected_highlights)
if (length(missing_h) > 0)
  cat("  Not detected:", paste(missing_h, collapse = ", "), "\n")

if (length(detected_highlights) < 2) {
  warning("Fewer than 2 highlighted proteins detected. Check gene name list.")
} else {

  # Build FC matrix: rows = detected highlights, cols = all 5 comparisons
  fc_matrix <- data.frame(PG.Genes = detected_highlights,
                          stringsAsFactors = FALSE)
  for (comp in comparison_names) {
    tmp <- results_list[[comp]] %>%
      filter(PG.Genes %in% detected_highlights) %>%
      select(PG.Genes, logFC) %>%
      rename(!!comp := logFC)
    fc_matrix <- left_join(fc_matrix, tmp, by = "PG.Genes")
  }

  fc_mat            <- as.matrix(fc_matrix[, -1])
  rownames(fc_mat)  <- fc_matrix$PG.Genes
  colnames(fc_mat)  <- c("\u0394spxA1", "\u0394spxA2", "\u0394clpX",
                          "pA2_A1",     "pA1_A2")
  fc_mat[is.na(fc_mat)] <- 0  # proteins not detected in a comparison -> 0

  # Significance asterisks
  sig_marks <- matrix("", nrow = nrow(fc_mat), ncol = ncol(fc_mat),
                      dimnames = dimnames(fc_mat))
  for (i in seq_along(comparison_names)) {
    comp     <- comparison_names[i]
    sig_hits <- results_list[[comp]] %>%
      filter(PG.Genes %in% rownames(fc_mat),
             adj.P.Val < PADJ_THRESH,
             abs(logFC) >= FC_THRESH) %>%
      pull(PG.Genes)
    sig_marks[rownames(fc_mat) %in% sig_hits, i] <- "*"
  }

  # Row annotation — functional category
  row_annot <- annot_filt %>%
    filter(PG.Genes %in% rownames(fc_mat)) %>%
    select(PG.Genes, COG_label) %>%
    distinct(PG.Genes, .keep_all = TRUE) %>%
    column_to_rownames("PG.Genes")
  row_annot <- row_annot[
    rownames(fc_mat)[rownames(fc_mat) %in% rownames(row_annot)], ,
    drop = FALSE
  ]
  colnames(row_annot) <- "Functional category"

  n_cats    <- length(unique(row_annot[["Functional category"]]))
  cat_colors <- setNames(
    colorRampPalette(brewer.pal(8, "Set2"))(n_cats),
    unique(row_annot[["Functional category"]])
  )

  # Substitute readable names for SpyM3 locus tags on y-axis
  locus_to_name <- c(
    "SpyM3_0583" = "ideS",
    "SpyM3_0428" = "gpoA",
    "SpyM3_0212" = "sufB",
    "SpyM3_0317" = "mntR",
    "SpyM3_1197" = "fnr"
  )
  display_names <- rownames(fc_mat)
  for (tag in names(locus_to_name)) {
    display_names[display_names == tag] <- locus_to_name[tag]
  }

  pdf("Figure4_heatmap_highlighted_proteins.pdf", width = 7, height = 8)
  pheatmap(
    fc_mat,
    color             = colorRampPalette(c("#2980B9", "white", "#C0392B"))(101),
    breaks            = seq(-3, 3, length.out = 102),
    display_numbers   = sig_marks,
    number_color      = "black",
    fontsize_number   = 9,
    annotation_row    = row_annot,
    annotation_colors = list(`Functional category` = cat_colors),
    cluster_cols      = FALSE,
    cluster_rows      = TRUE,
    border_color      = "grey80",
    labels_row        = display_names,
    main              = sprintf(
      "log\u2082 Fold Change vs. WT  (* padj < %.2f, |log\u2082FC| \u2265 log\u2082(1.5))",
      PADJ_THRESH),
    fontsize          = 9,
    cellwidth         = 45,
    cellheight        = 18
  )
  dev.off()
  cat("Heatmap saved.\n")
}


# =============================================================================
# 9. PROTEOME COMPARISON SCATTER PLOT
# =============================================================================
# Compares the full ∆spxA1 and ∆spxA2 fold-change distributions across all
# detected proteins to illustrate global co-regulation and paralog-specific
# divergence.

full_compare <- results_list$spxA1_vs_WT %>%
  select(PG.ProteinNames, PG.Genes, SpyM3_locus_tag, logFC, adj.P.Val,
         COG_label) %>%
  inner_join(
    results_list$spxA2_vs_WT %>%
      select(PG.ProteinNames, logFC, adj.P.Val) %>%
      rename(logFC_A2 = logFC, padj_A2 = adj.P.Val),
    by = "PG.ProteinNames"
  ) %>%
  rename(logFC_A1 = logFC, padj_A1 = adj.P.Val) %>%
  mutate(
    sig_cat = case_when(
      (padj_A1 < PADJ_THRESH & abs(logFC_A1) >= FC_THRESH) &
        (padj_A2 < PADJ_THRESH & abs(logFC_A2) >= FC_THRESH) ~ "Both",
      (padj_A1 < PADJ_THRESH & abs(logFC_A1) >= FC_THRESH) &
       !(padj_A2 < PADJ_THRESH & abs(logFC_A2) >= FC_THRESH) ~ "\u0394spxA1 only",
      !(padj_A1 < PADJ_THRESH & abs(logFC_A1) >= FC_THRESH) &
        (padj_A2 < PADJ_THRESH & abs(logFC_A2) >= FC_THRESH) ~ "\u0394spxA2 only",
      TRUE ~ "Neither"
    ),
    sig_cat = factor(sig_cat,
                     levels = c("Both", "\u0394spxA1 only",
                                "\u0394spxA2 only", "Neither"))
  )

# Save for external use
write_csv(full_compare, "proteome_comparison_full.csv")

# Correlation statistics
clean_fc <- full_compare %>% drop_na(logFC_A1, logFC_A2)
r_p  <- cor(clean_fc$logFC_A1, clean_fc$logFC_A2, method = "pearson")
r_s  <- cor(clean_fc$logFC_A1, clean_fc$logFC_A2, method = "spearman")
cat(sprintf("\nProteome comparison correlation (n=%d):\n", nrow(clean_fc)))
cat(sprintf("  Pearson r  = %.3f\n  Spearman r = %.3f\n", r_p, r_s))

# Proteins to label on scatter plot
label_genes <- c("speB", "sodA", "arcA", "slo", "ska", "hasA", "pepO",
                 "prtS", "sagC", "grab", "dpr", "nox.1", "cysM",
                 "SpyM3_0583", "SpyM3_1799 (SpxA2)", "copY")

full_compare <- full_compare %>%
  mutate(
    label = ifelse(PG.Genes %in% label_genes & sig_cat != "Neither",
                   PG.Genes, NA_character_),
    # Use readable names for locus tags
    label = case_when(
      label == "SpyM3_0583"         ~ "ideS",
      label == "SpyM3_1799 (SpxA2)" ~ "spxA2",
      TRUE                          ~ label
    )
  )

ann_text <- sprintf("Pearson r = %.3f\nSpearman r = %.3f\nn = %d",
                    r_p, r_s, nrow(clean_fc))

p_compare <- ggplot(
  full_compare %>% drop_na(logFC_A1, logFC_A2),
  aes(x = logFC_A1, y = logFC_A2, color = sig_cat)
) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_hline(yintercept = c(-FC_THRESH,  FC_THRESH),
             linetype = "dashed", color = "grey75", linewidth = 0.3) +
  geom_vline(xintercept = c(-FC_THRESH, FC_THRESH),
             linetype = "dashed", color = "grey75", linewidth = 0.3) +
  geom_point(data = . %>% filter(sig_cat == "Neither"),
             alpha = 0.25, size = 0.9) +
  geom_point(data = . %>% filter(sig_cat != "Neither"),
             alpha = 0.80, size = 1.8) +
  geom_smooth(method = "lm", se = TRUE, color = "grey30",
              linewidth = 0.6, linetype = "dashed", alpha = 0.1) +
  geom_text_repel(aes(label = label), size = 2.7, show.legend = FALSE,
                  max.overlaps = 20, segment.size = 0.3,
                  box.padding = 0.4, point.padding = 0.2) +
  scale_color_manual(
    values = c(
      "Both"          = "#8B2FC9",
      "\u0394spxA1 only" = "#C0392B",
      "\u0394spxA2 only" = "#2980B9",
      "Neither"       = "grey75"
    ),
    labels = c(
      sprintf("Both significant (n=%d)",
              sum(full_compare$sig_cat == "Both", na.rm = TRUE)),
      sprintf("\u0394spxA1 only (n=%d)",
              sum(full_compare$sig_cat == "\u0394spxA1 only", na.rm = TRUE)),
      sprintf("\u0394spxA2 only (n=%d)",
              sum(full_compare$sig_cat == "\u0394spxA2 only", na.rm = TRUE)),
      sprintf("Neither (n=%d)",
              sum(full_compare$sig_cat == "Neither", na.rm = TRUE))
    )
  ) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.5,
           label = ann_text, size = 2.8, color = "grey30") +
  labs(
    x     = expression("\u0394"*spxA1~~log[2]~"Fold Change vs. WT"),
    y     = expression("\u0394"*spxA2~~log[2]~"Fold Change vs. WT"),
    color = NULL,
    title = "\u0394spxA1 vs. \u0394spxA2 proteome comparison"
  ) +
  coord_cartesian(xlim = c(-9, 4), ylim = c(-5, 5)) +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "bottom",
    legend.text      = element_text(size = 9),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(size = 11, face = "bold")
  )

ggsave(file.path(OUT_DIR, "Figure4_proteome_comparison_scatter.pdf"),
       p_compare, width = 6, height = 6)
cat("Proteome comparison scatter saved.\n")


# =============================================================================
# 10. NARRATIVE CONSISTENCY CHECK — SUPPLEMENTAL STRAINS
# =============================================================================
# Assesses whether promoter swap and clpX-KO proteomes parallel the primary
# mutant narratives, supporting functional non-redundancy of SpxA1 and SpxA2.

cat("\n=== NARRATIVE CONSISTENCY CHECK ===\n")

sig_genes <- lapply(results_list, function(res) {
  res %>%
    filter(adj.P.Val < PADJ_THRESH, abs(logFC) >= FC_THRESH) %>%
    select(PG.Genes, logFC, direction, COG_label)
})

# pA2_A1 (SpxA1 expressed from spxA2 promoter) vs. ∆spxA1
cat("\npA2_A1 vs. \u0394spxA1 significant protein overlap:\n")
overlap_1 <- intersect(sig_genes$pA2A1_vs_WT$PG.Genes,
                       sig_genes$spxA1_vs_WT$PG.Genes)
cat(sprintf("  Shared: %d  (\u0394spxA1: %d sig; pA2_A1: %d sig)\n",
            length(overlap_1),
            nrow(sig_genes$spxA1_vs_WT),
            nrow(sig_genes$pA2A1_vs_WT)))
if (length(overlap_1) > 0) {
  conc1 <- sig_genes$pA2A1_vs_WT %>%
    filter(PG.Genes %in% overlap_1) %>%
    rename(dir_pA2A1 = direction) %>%
    left_join(sig_genes$spxA1_vs_WT %>%
                select(PG.Genes, direction) %>%
                rename(dir_spxA1 = direction),
              by = "PG.Genes") %>%
    mutate(concordant = dir_pA2A1 == dir_spxA1)
  cat(sprintf("  Direction concordance: %d/%d (%.0f%%)\n",
              sum(conc1$concordant, na.rm = TRUE),
              nrow(conc1),
              100 * mean(conc1$concordant, na.rm = TRUE)))
}

# pA1_A2 (SpxA2 expressed from spxA1 promoter) vs. ∆spxA2
cat("\npA1_A2 vs. \u0394spxA2 significant protein overlap:\n")
overlap_2 <- intersect(sig_genes$pA1A2_vs_WT$PG.Genes,
                       sig_genes$spxA2_vs_WT$PG.Genes)
cat(sprintf("  Shared: %d  (\u0394spxA2: %d sig; pA1_A2: %d sig)\n",
            length(overlap_2),
            nrow(sig_genes$spxA2_vs_WT),
            nrow(sig_genes$pA1A2_vs_WT)))
if (length(overlap_2) > 0) {
  conc2 <- sig_genes$pA1A2_vs_WT %>%
    filter(PG.Genes %in% overlap_2) %>%
    rename(dir_pA1A2 = direction) %>%
    left_join(sig_genes$spxA2_vs_WT %>%
                select(PG.Genes, direction) %>%
                rename(dir_spxA2 = direction),
              by = "PG.Genes") %>%
    mutate(concordant = dir_pA1A2 == dir_spxA2)
  cat(sprintf("  Direction concordance: %d/%d (%.0f%%)\n",
              sum(conc2$concordant, na.rm = TRUE),
              nrow(conc2),
              100 * mean(conc2$concordant, na.rm = TRUE)))
}

# ∆clpX: SpxA1/SpxA2 protein abundance
cat("\n\u0394clpX — SpxA1 and SpxA2 protein abundance:\n")
clpx_spx <- results_list$clpX_vs_WT %>%
  filter(PG.Genes %in% c("spx (SpxA1)", "SpyM3_1799 (SpxA2)", "clpP", "clpX"))
if (nrow(clpx_spx) > 0)
  print(clpx_spx %>% select(PG.Genes, logFC, adj.P.Val, direction))


# =============================================================================
# 11. EXPORT SUPPLEMENTAL TABLES
# =============================================================================

export_results <- lapply(comparison_names, function(comp) {
  results_list[[comp]] %>%
    left_join(cog_annot %>% select(PG.ProteinNames, SpyM3_locus_tag),
              by = "PG.ProteinNames") %>%
    select(PG.ProteinNames, PG.Genes, SpyM3_locus_tag, PG.ProteinDescriptions,
           logFC, AveExpr, t, P.Value, adj.P.Val, direction, COG_label) %>%
    rename(
      UniProt_ID          = PG.ProteinNames,
      Gene                = PG.Genes,
      SpyM3_Locus_Tag     = SpyM3_locus_tag,
      Description         = PG.ProteinDescriptions,
      log2FC              = logFC,
      MeanExpression      = AveExpr,
      t_statistic         = t,
      P_value             = P.Value,
      Padj_BH             = adj.P.Val,
      Direction           = direction,
      Functional_Category = COG_label
    ) %>%
    arrange(Padj_BH)
})

names(export_results) <- c(
  "spxA1_vs_WT", "spxA2_vs_WT", "clpX_vs_WT", "pA2A1_vs_WT", "pA1A2_vs_WT"
)

write_xlsx(export_results,
           "Supplemental_Proteomics_Results_AllComparisons.xlsx")

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Output files:\n")
cat("  Figure4_volcanos_main.pdf\n")
cat("  Figure4_functional_categories_main.pdf\n")
cat("  Figure4_heatmap_highlighted_proteins.pdf\n")
cat("  Figure4_proteome_comparison_scatter.pdf\n")
cat("  FigureS_volcanos_supplemental.pdf\n")
cat("  FigureS_functional_categories_all.pdf\n")
cat("  Supplemental_Proteomics_Results_AllComparisons.xlsx\n")
cat("  proteomics_summary_statistics.csv\n")
cat("  proteome_comparison_full.csv\n")

# Session information for reproducibility
cat("\nSession info:\n")
print(sessionInfo())

