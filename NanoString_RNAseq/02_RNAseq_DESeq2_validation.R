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
# MODULE 2: RNA-seq DESeq2 Analysis and NanoString Cross-Platform Validation
# =============================================================================
# Description:
#   - Loads NanoString workspace from Module 1
#   - Maps NanoString probes to MGAS10870 C9Q_ locus tags via crosswalk
#   - Imports RNA-seq raw counts from CLC output (9 samples)
#   - DESeq2: spxA1KO vs WT and spxA2KO vs WT (UI condition)
#   - Cross-platform concordance: NanoString logFC vs DESeq2 logFC
#   - Expression correlation: WT UI NanoString vs RNA-seq log2(TPM+1)
#   - QC figures: PCA, sample distances, dispersion, MA plots
#
# Prerequisite: Run 01_NanoString_analysis.R first
#
# Input files (place in data/ subdirectory):
#   NanoString_workspace_noMo.RData    — from Module 1
#   NS_crosswalk_final.csv             — from 02_NS_genbank_crosswalk.py
#   Spx_RNAseq.xlsx                    — CLC RNA-seq output (9 sheets)
#   MGAS10870_locus-tag_key.xlsx       — C9Q_ to SpyM3 locus tag mapping
#
# GEO submission files generated:
#   DESeq2_RNAseq_results.xlsx         — differential expression results
#   ST_Validation_NanoString_RNAseq.xlsx — cross-platform concordance
# =============================================================================

library(DESeq2)
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
BASE_DIR  <- "."
DATA_DIR  <- file.path(BASE_DIR, "data")
OUT_DIR   <- file.path(BASE_DIR, "NanoString_RNAseq", "output")
dir.create(OUT_DIR, showWarnings=FALSE, recursive=TRUE)


# ── LOAD NANOSTRING WORKSPACE ─────────────────────────────────────────────────
# Set paths before running — must match Module 1
data_dir <- "/path/to/your/data"
outdir   <- file.path(data_dir, "new_noMo")

load(file.path(OUT_DIR, "NanoString_workspace_noMo.RData"))
message("NanoString workspace loaded.")

# =============================================================================
# MODULE 2.1 — LOAD CROSSWALK (from NS_genbank_crosswalk.py)
# =============================================================================

crosswalk <- read.csv(
  file.path(data_dir, file.path(DATA_DIR, "NS_crosswalk_final.csv")),
  stringsAsFactors = FALSE
)

cat("Crosswalk loaded:", nrow(crosswalk), "probes\n")
cat("Match type summary:\n")
print(table(crosswalk$Match_type))

# Flag probes with clean, unambiguous C9Q_ assignment
crosswalk$Has_C9Q <- (
  !is.na(crosswalk$C9Q_tag) &
  crosswalk$C9Q_tag != "" &
  !grepl("\\|", crosswalk$C9Q_tag) &
  crosswalk$Match_type != "absent_M10870_verified" &
  crosswalk$Match_type != "not_in_MGAS5005"
)

cat("\nProbes with clean C9Q_ assignment:", sum(crosswalk$Has_C9Q), "\n")
cat("Probes excluded from validation:", sum(!crosswalk$Has_C9Q), "\n")

# =============================================================================
# MODULE 2.2 — IMPORT LOCUS TAG KEY (C9Q_ → SpyM3 → gene name)
# =============================================================================

locus_key <- read_excel(
  file.path(data_dir, file.path(DATA_DIR, "MGAS10870_locus-tag_key.xlsx")))
colnames(locus_key) <- c("Name", "C9Q_tag", "SpyM3_tag",
                           "Start", "Stop", "Annotation_RNAseq")
locus_key <- locus_key[
  !is.na(locus_key$SpyM3_tag) &
  locus_key$SpyM3_tag != "#N/A", ]

cat("Locus key rows retained:", nrow(locus_key), "\n")

# Build lookup maps
c9q_to_name <- setNames(
  as.character(locus_key$Name),
  as.character(locus_key$C9Q_tag)
)
spym3_map <- setNames(
  as.character(locus_key$SpyM3_tag),
  as.character(locus_key$C9Q_tag)
)

# =============================================================================
# MODULE 2.3 — IMPORT RNA-seq COUNTS (CLC output)
# =============================================================================

rnaseq_file <- file.path(data_dir, "Spx_RNAseq_counts.xlsx")
sheets      <- excel_sheets(rnaseq_file)
cat("RNA-seq sheets:", paste(sheets, collapse = ", "), "\n")

# Each sheet: Name, Region, TPM, RPKM, Gene_length, Unique_reads, Total_reads
import_sheet <- function(sheet_name) {
  df <- read_excel(rnaseq_file, sheet = sheet_name)
  colnames(df) <- c("Name", "Region", "TPM", "RPKM",
                     "Gene_length", "Unique_reads", "Total_reads")
  df[, c("Name", "TPM", "Total_reads")]
}

counts_list <- lapply(sheets, import_sheet)
names(counts_list) <- sheets

# Wide TPM matrix (for expression correlation with NanoString)
tpm_wide <- counts_list[[1]][, c("Name", "TPM")]
colnames(tpm_wide)[2] <- sheets[1]
for (i in 2:length(sheets)) {
  tmp <- counts_list[[i]][, c("Name", "TPM")]
  colnames(tmp)[2] <- sheets[i]
  tpm_wide <- merge(tpm_wide, tmp, by = "Name", all = TRUE)
}

# Raw count matrix for DESeq2
count_matrix <- counts_list[[1]][, c("Name", "Total_reads")]
colnames(count_matrix)[2] <- sheets[1]
for (i in 2:length(sheets)) {
  tmp <- counts_list[[i]][, c("Name", "Total_reads")]
  colnames(tmp)[2] <- sheets[i]
  count_matrix <- merge(count_matrix, tmp, by = "Name", all = TRUE)
}
rownames(count_matrix) <- count_matrix$Name
count_matrix <- as.matrix(count_matrix[, -1])
mode(count_matrix) <- "integer"
count_matrix[is.na(count_matrix)] <- 0L

cat("\nRaw count matrix:", nrow(count_matrix), "x",
    ncol(count_matrix), "\n")
cat("Library sizes (millions):\n")
print(round(colSums(count_matrix) / 1e6, 2))

# =============================================================================
# MODULE 2.4 — BUILD MATCHED NANOSTRING ↔ RNA-SEQ MATRIX
# =============================================================================

# Merge NanoString probes → C9Q_ tag → RNA-seq TPM
crosswalk_clean <- crosswalk[crosswalk$Has_C9Q, ]

# Map C9Q_ → RNA-seq row name (gene name or C9Q_ tag)
tpm_mapped <- merge(tpm_wide,
                     locus_key[, c("Name", "C9Q_tag", "Annotation_RNAseq")],
                     by = "Name", all.x = TRUE)

ns_rnaseq_map <- merge(
  crosswalk_clean[, c("NanoString_ProbeID", "M5005_gene",
                       "C9Q_tag", "C9Q_gene_name", "SpyM3_tag",
                       "Match_type", "Kmer_identity")],
  tpm_mapped[, c("C9Q_tag", sheets, "Annotation_RNAseq")],
  by = "C9Q_tag", all.x = TRUE
)

# Add NanoString gene Label
label_map <- setNames(annot$Label, annot$ProbeID)
ns_rnaseq_map$Label <- label_map[ns_rnaseq_map$NanoString_ProbeID]
ns_rnaseq_map$Label <- ifelse(
  is.na(ns_rnaseq_map$Label),
  ns_rnaseq_map$NanoString_ProbeID,
  ns_rnaseq_map$Label)

cat("\nMatched genes (NanoString + RNA-seq):", nrow(ns_rnaseq_map), "\n")

# Mean TPM per strain (UI only)
wt_sheets   <- grep("^WT",   sheets, value = TRUE)
a1ko_sheets <- grep("^A1KO", sheets, value = TRUE)
a2ko_sheets <- grep("^A2KO", sheets, value = TRUE)

for (col in c(sheets)) {
  ns_rnaseq_map[[col]] <- as.numeric(ns_rnaseq_map[[col]])
}

ns_rnaseq_map$TPM_WT_mean   <- rowMeans(ns_rnaseq_map[, wt_sheets],
                                          na.rm = TRUE)
ns_rnaseq_map$TPM_A1KO_mean <- rowMeans(ns_rnaseq_map[, a1ko_sheets],
                                          na.rm = TRUE)
ns_rnaseq_map$TPM_A2KO_mean <- rowMeans(ns_rnaseq_map[, a2ko_sheets],
                                          na.rm = TRUE)
ns_rnaseq_map$log2TPM_WT    <- log2(ns_rnaseq_map$TPM_WT_mean   + 1)
ns_rnaseq_map$log2TPM_A1KO  <- log2(ns_rnaseq_map$TPM_A1KO_mean + 1)
ns_rnaseq_map$log2TPM_A2KO  <- log2(ns_rnaseq_map$TPM_A2KO_mean + 1)

# NanoString UI mean expression per strain
wt_ui_samples   <- metadata_filt$SampleID[
  metadata_filt$Strain == "WT"   & metadata_filt$Condition == "UI"]
a1ko_ui_samples <- metadata_filt$SampleID[
  metadata_filt$Strain == "spxA1KO" & metadata_filt$Condition == "UI"]
a2ko_ui_samples <- metadata_filt$SampleID[
  metadata_filt$Strain == "spxA2KO" & metadata_filt$Condition == "UI"]

ns_ui <- data.frame(
  NanoString_ProbeID = rownames(counts_filt_log2),
  NS_WT_mean   = rowMeans(counts_filt_log2[, wt_ui_samples]),
  NS_A1KO_mean = rowMeans(counts_filt_log2[, a1ko_ui_samples]),
  NS_A2KO_mean = rowMeans(counts_filt_log2[, a2ko_ui_samples]),
  stringsAsFactors = FALSE
)

ns_rnaseq_map <- merge(ns_rnaseq_map, ns_ui,
                         by = "NanoString_ProbeID", all.x = TRUE)

# RNA-seq logFC (log2TPM difference)
ns_rnaseq_map$RNAseq_logFC_A1KO_vs_WT <-
  as.numeric(ns_rnaseq_map$log2TPM_A1KO) -
  as.numeric(ns_rnaseq_map$log2TPM_WT)
ns_rnaseq_map$RNAseq_logFC_A2KO_vs_WT <-
  as.numeric(ns_rnaseq_map$log2TPM_A2KO) -
  as.numeric(ns_rnaseq_map$log2TPM_WT)

cat("Final matched matrix:", nrow(ns_rnaseq_map), "genes\n")
cat("With both NS and RNA-seq data:",
    sum(!is.na(ns_rnaseq_map$NS_WT_mean) &
        !is.na(ns_rnaseq_map$log2TPM_WT)), "\n")

# =============================================================================
# MODULE 2.5 — EXPRESSION CORRELATION: WT UI (NanoString vs RNA-seq TPM)
# =============================================================================

cor_df <- ns_rnaseq_map[
  !is.na(ns_rnaseq_map$NS_WT_mean) &
  !is.na(ns_rnaseq_map$log2TPM_WT) &
  ns_rnaseq_map$log2TPM_WT > 0, ]

r_wt  <- cor(cor_df$NS_WT_mean, cor_df$log2TPM_WT, method = "pearson")
r2_wt <- r_wt^2
cat("\nWT UI Pearson r:", round(r_wt, 4), "\n")
cat("WT UI R²:", round(r2_wt, 4), "\n")
cat("n genes:", nrow(cor_df), "\n")

p_cor_wt <- ggplot(cor_df,
    aes(x = log2TPM_WT, y = NS_WT_mean, label = Label)) +
  geom_point(size = 2.5, alpha = 0.75, color = "#2166AC") +
  geom_smooth(method = "lm", se = TRUE,
              color = "#D6604D", linewidth = 0.8,
              fill = "#D6604D", alpha = 0.15) +
  geom_text_repel(
    data = subset(cor_df,
      NS_WT_mean   > quantile(cor_df$NS_WT_mean,   0.9) |
      log2TPM_WT   > quantile(cor_df$log2TPM_WT,   0.9)),
    size = 2.5, max.overlaps = 20, box.padding = 0.3,
    fontface = "italic", min.segment.length = 0,
    segment.color = "grey50", segment.size = 0.3
  ) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
           label = paste0("r = ", round(r_wt, 3),
                          "\nR\u00b2 = ", round(r2_wt, 3),
                          "\nn = ", nrow(cor_df)),
           size = 4, color = "#D6604D") +
  labs(
    title    = "NanoString vs RNA-seq: WT UI expression",
    subtitle = "Independent experiments, matched growth conditions",
    x        = "RNA-seq log\u2082(TPM + 1)",
    y        = "NanoString log\u2082(normalized counts + 1)"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(outdir, "Validation_NanoString_vs_RNAseq_WT_UI.pdf"),
       plot = p_cor_wt, width = 7, height = 6, useDingbats = FALSE)

# =============================================================================
# MODULE 2.6 — DESEQ2 DIFFERENTIAL EXPRESSION
# =============================================================================

# Sample metadata for DESeq2
col_data <- data.frame(
  Sample    = colnames(count_matrix),
  Strain    = case_when(
    grepl("^WT",   colnames(count_matrix)) ~ "WT",
    grepl("^A1KO", colnames(count_matrix)) ~ "spxA1KO",
    grepl("^A2KO", colnames(count_matrix)) ~ "spxA2KO"
  ),
  Replicate = sub(".*-([ABC])$", "\\1", colnames(count_matrix)),
  stringsAsFactors = FALSE
)
col_data$Strain   <- factor(col_data$Strain,
                              levels = c("WT", "spxA1KO", "spxA2KO"))
rownames(col_data) <- col_data$Sample

cat("\nDESeq2 sample metadata:\n")
print(col_data[, c("Sample", "Strain", "Replicate")])

stopifnot(all(colnames(count_matrix) == rownames(col_data)))

# Build DESeq2 object and filter low-count genes
dds  <- DESeqDataSetFromMatrix(countData = count_matrix,
                                 colData   = col_data,
                                 design    = ~ Strain)
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]
cat("\nGenes after pre-filter (>=10 total counts):", nrow(dds), "\n")

# Run DESeq2
set.seed(42)
dds <- DESeq(dds)

cat("\nSize factors:\n")
print(round(sizeFactors(dds), 3))

# Extract results
res_a1ko    <- results(dds, contrast = c("Strain", "spxA1KO", "WT"),
                        alpha = 0.05, pAdjustMethod = "BH")
res_a2ko    <- results(dds, contrast = c("Strain", "spxA2KO", "WT"),
                        alpha = 0.05, pAdjustMethod = "BH")
res_a1ko_df <- as.data.frame(res_a1ko); res_a1ko_df$Name <- rownames(res_a1ko_df)
res_a2ko_df <- as.data.frame(res_a2ko); res_a2ko_df$Name <- rownames(res_a2ko_df)

# Significance flag: padj<0.05 AND |log2FC|>log2(1.5)
res_a1ko_df$Sig <- !is.na(res_a1ko_df$padj) &
                    res_a1ko_df$padj < 0.05 &
                    abs(res_a1ko_df$log2FoldChange) > log2(1.5)
res_a2ko_df$Sig <- !is.na(res_a2ko_df$padj) &
                    res_a2ko_df$padj < 0.05 &
                    abs(res_a2ko_df$log2FoldChange) > log2(1.5)

cat("\nDESeq2 DE genes (padj<0.05, |log2FC|>log2(1.5)):\n")
cat("A1KO vs WT:", sum(res_a1ko_df$Sig, na.rm = TRUE), "\n")
cat("A2KO vs WT:", sum(res_a2ko_df$Sig, na.rm = TRUE), "\n")

# Map DESeq2 results back to NanoString probes via C9Q_ tag
name_to_c9q <- setNames(names(c9q_to_name), c9q_to_name)

add_c9q <- function(df) {
  df$C9Q_tag <- ifelse(
    df$Name %in% names(name_to_c9q),
    name_to_c9q[df$Name],
    df$Name
  )
  df
}
res_a1ko_df <- add_c9q(res_a1ko_df)
res_a2ko_df <- add_c9q(res_a2ko_df)

map_deseq <- function(res_df, suffix) {
  merged <- merge(
    crosswalk[crosswalk$Has_C9Q,
               c("NanoString_ProbeID", "C9Q_tag", "M5005_gene")],
    res_df[, c("C9Q_tag", "log2FoldChange", "pvalue", "padj", "Sig")],
    by = "C9Q_tag", all.x = TRUE
  )
  colnames(merged)[4:7] <- paste0(
    c("DESeq2_log2FC_", "DESeq2_pvalue_",
      "DESeq2_padj_", "DESeq2_Sig_"), suffix)
  merged
}
deseq_a1ko_mapped <- map_deseq(res_a1ko_df, "spxA1KO")
deseq_a2ko_mapped <- map_deseq(res_a2ko_df, "spxA2KO")

# =============================================================================
# MODULE 2.7 — CROSS-PLATFORM CONCORDANCE
# =============================================================================

# NanoString cross-strain logFC from limma
ns_a1ko_lfc <- get_full_table(fit2, "spxA1KO_vs_WT_UI", annot, ko_probes)[,
  c("ProbeID", "logFC", "Sig_p05")]
colnames(ns_a1ko_lfc)[2:3] <- c("NS_log2FC_spxA1KO_vs_WT", "NS_Sig_spxA1KO")

ns_a2ko_lfc <- get_full_table(fit2, "spxA2KO_vs_WT_UI", annot, ko_probes)[,
  c("ProbeID", "logFC", "Sig_p05")]
colnames(ns_a2ko_lfc)[2:3] <- c("NS_log2FC_spxA2KO_vs_WT", "NS_Sig_spxA2KO")

# Build concordance table
concordance_df <- merge(
  ns_rnaseq_map[, c("NanoString_ProbeID", "Label", "C9Q_tag",
                     "RNAseq_logFC_A1KO_vs_WT",
                     "RNAseq_logFC_A2KO_vs_WT")],
  ns_a1ko_lfc, by.x = "NanoString_ProbeID", by.y = "ProbeID",
  all.x = TRUE)
concordance_df <- merge(concordance_df, ns_a2ko_lfc,
  by.x = "NanoString_ProbeID", by.y = "ProbeID", all.x = TRUE)

# Add DESeq2 results
concordance_df <- merge(concordance_df,
  deseq_a1ko_mapped[, c("NanoString_ProbeID", "DESeq2_log2FC_spxA1KO",
                          "DESeq2_padj_A1KO", "DESeq2_Sig_spxA1KO")],
  by = "NanoString_ProbeID", all.x = TRUE)
concordance_df <- merge(concordance_df,
  deseq_a2ko_mapped[, c("NanoString_ProbeID", "DESeq2_log2FC_spxA2KO",
                          "DESeq2_padj_A2KO", "DESeq2_Sig_spxA2KO")],
  by = "NanoString_ProbeID", all.x = TRUE)

# Directional agreement
concordance_df$direction_concordance_spxA1KO_DESeq2 <-
  sign(concordance_df$NS_log2FC_spxA1KO_vs_WT) ==
  sign(concordance_df$DESeq2_log2FC_spxA1KO)
concordance_df$direction_concordance_spxA2KO_DESeq2 <-
  sign(concordance_df$NS_log2FC_spxA2KO_vs_WT) ==
  sign(concordance_df$DESeq2_log2FC_spxA2KO)

cat("\n=== CROSS-PLATFORM CONCORDANCE ===\n")
cat("\nDirectional agreement — all mapped genes:\n")
cat("A1KO:", sum(concordance_df$direction_concordance_spxA1KO_DESeq2, na.rm = TRUE), "/",
    sum(!is.na(concordance_df$direction_concordance_spxA1KO_DESeq2)),
    sprintf("(%.1f%%)\n",
      100 * mean(concordance_df$direction_concordance_spxA1KO_DESeq2, na.rm = TRUE)))
cat("A2KO:", sum(concordance_df$direction_concordance_spxA2KO_DESeq2, na.rm = TRUE), "/",
    sum(!is.na(concordance_df$direction_concordance_spxA2KO_DESeq2)),
    sprintf("(%.1f%%)\n",
      100 * mean(concordance_df$direction_concordance_spxA2KO_DESeq2, na.rm = TRUE)))

ns_de_a1ko <- concordance_df[
  !is.na(concordance_df$NS_Sig_spxA1KO) & concordance_df$NS_Sig_spxA1KO, ]
ns_de_a2ko <- concordance_df[
  !is.na(concordance_df$NS_Sig_spxA2KO) & concordance_df$NS_Sig_spxA2KO, ]

cat("\nAmong NanoString-DE genes:\n")
cat("A1KO (n =", nrow(ns_de_a1ko), "):",
    sum(ns_de_a1ko$direction_concordance_spxA1KO_DESeq2, na.rm = TRUE), "/",
    sum(!is.na(ns_de_a1ko$direction_concordance_spxA1KO_DESeq2)),
    sprintf("(%.1f%%)\n",
      100 * mean(ns_de_a1ko$direction_concordance_spxA1KO_DESeq2, na.rm = TRUE)))
cat("A2KO (n =", nrow(ns_de_a2ko), "):",
    sum(ns_de_a2ko$direction_concordance_spxA2KO_DESeq2, na.rm = TRUE), "/",
    sum(!is.na(ns_de_a2ko$direction_concordance_spxA2KO_DESeq2)),
    sprintf("(%.1f%%)\n",
      100 * mean(ns_de_a2ko$direction_concordance_spxA2KO_DESeq2, na.rm = TRUE)))

both_sig_a1ko <- ns_de_a1ko[
  !is.na(ns_de_a1ko$DESeq2_Sig_spxA1KO) &
  ns_de_a1ko$DESeq2_Sig_spxA1KO &
  ns_de_a1ko$direction_concordance_spxA1KO_DESeq2, ]
both_sig_a2ko <- ns_de_a2ko[
  !is.na(ns_de_a2ko$DESeq2_Sig_spxA2KO) &
  ns_de_a2ko$DESeq2_Sig_spxA2KO &
  ns_de_a2ko$direction_concordance_spxA2KO_DESeq2, ]

cat("\nNanoString DE genes also significant by DESeq2:\n")
cat("A1KO:", nrow(both_sig_a1ko), "of", nrow(ns_de_a1ko),
    sprintf("(%.1f%%)\n",
      100 * nrow(both_sig_a1ko) / nrow(ns_de_a1ko)))
cat("A2KO:", nrow(both_sig_a2ko), "of", nrow(ns_de_a2ko),
    sprintf("(%.1f%%)\n",
      100 * nrow(both_sig_a2ko) / nrow(ns_de_a2ko)))

r_deseq_a1ko <- cor(
  concordance_df$NS_log2FC_spxA1KO_vs_WT,
  concordance_df$DESeq2_log2FC_spxA1KO,
  use = "complete.obs", method = "pearson")
r_deseq_a2ko <- cor(
  concordance_df$NS_log2FC_spxA2KO_vs_WT,
  concordance_df$DESeq2_log2FC_spxA2KO,
  use = "complete.obs", method = "pearson")

cat("\nlogFC Pearson r (NanoString vs DESeq2):\n")
cat("A1KO: r =", round(r_deseq_a1ko, 3), "\n")
cat("A2KO: r =", round(r_deseq_a2ko, 3), "\n")

# ── CONCORDANCE SCATTER PLOTS ─────────────────────────────────────────────────
make_concordance_scatter <- function(df, x_col, y_col, sig_col,
                                      r_val, x_lab, y_lab, title) {
  df_plot <- df[!is.na(df[[x_col]]) & !is.na(df[[y_col]]), ]
  ggplot(df_plot,
      aes_string(x = x_col, y = y_col,
                  color = sig_col, label = "Label")) +
    geom_hline(yintercept = 0, color = "grey80", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey80", linewidth = 0.3) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", color = "grey55", linewidth = 0.4) +
    geom_point(size = 2.2, alpha = 0.85) +
    geom_text_repel(
      data        = df_plot[df_plot[[sig_col]] == TRUE, ],
      size        = 2.8, max.overlaps = 25, box.padding = 0.4,
      fontface    = "italic", min.segment.length = 0,
      segment.color = "grey50", segment.size = 0.3,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = c("TRUE" = "#D6604D", "FALSE" = "grey70"),
      labels = c("TRUE" = "DE in NanoString", "FALSE" = "Not DE"),
      name = ""
    ) +
    annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
             label = paste0("r = ", round(r_val, 3),
                            "\nn = ", nrow(df_plot)),
             size = 4, color = "#D6604D") +
    labs(title = title,
         subtitle = "UI condition — cross-platform concordance",
         x = x_lab, y = y_lab) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          legend.position  = "bottom")
}

p_conc_a2ko <- make_concordance_scatter(
  concordance_df, "DESeq2_log2FC_spxA2KO", "NS_log2FC_spxA2KO_vs_WT",
  "NS_Sig_spxA2KO", r_deseq_a2ko,
  "DESeq2 log\u2082FC (\u0394spxA2 vs WT, UI)",
  "NanoString log\u2082FC (\u0394spxA2 vs WT, UI)",
  expression(paste(Delta, "spxA2 vs WT: NanoString vs DESeq2 logFC"))
)
ggsave(file.path(outdir, "Validation_DESeq2_concordance_A2KO.pdf"),
       plot = p_conc_a2ko, width = 7, height = 6, useDingbats = FALSE)

p_conc_a1ko <- make_concordance_scatter(
  concordance_df, "DESeq2_log2FC_spxA1KO", "NS_log2FC_spxA1KO_vs_WT",
  "NS_Sig_spxA1KO", r_deseq_a1ko,
  "DESeq2 log\u2082FC (\u0394spxA1 vs WT, UI)",
  "NanoString log\u2082FC (\u0394spxA1 vs WT, UI)",
  expression(paste(Delta, "spxA1 vs WT: NanoString vs DESeq2 logFC"))
)
ggsave(file.path(outdir, "Validation_DESeq2_concordance_A1KO.pdf"),
       plot = p_conc_a1ko, width = 7, height = 6, useDingbats = FALSE)

# =============================================================================
# MODULE 2.8 — RNA-seq QC FIGURES
# =============================================================================

# Library sizes
lib_sizes <- data.frame(
  Sample     = colnames(count_matrix),
  TotalReads = colSums(count_matrix),
  stringsAsFactors = FALSE
)
lib_sizes <- merge(lib_sizes, col_data, by = "Sample")
lib_sizes$Strain <- factor(lib_sizes$Strain, levels = c("WT","spxA1KO","spxA2KO"))

p_libsize <- ggplot(lib_sizes,
    aes(x = Sample, y = TotalReads / 1e6, fill = Strain)) +
  geom_bar(stat = "identity", color = "white", linewidth = 0.3) +
  scale_fill_manual(values = strain_colors) +
  geom_text(aes(label = round(TotalReads / 1e6, 1)),
            vjust = -0.4, size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "RNA-seq library sizes", x = "",
       y = "Total reads (millions)") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor    = element_blank())
ggsave(file.path(outdir, "RNAseq_QC_library_sizes.pdf"),
       plot = p_libsize, width = 8, height = 4.5, useDingbats = FALSE)

# VST PCA
vst_data <- vst(dds, blind = TRUE)
vst_mat  <- assay(vst_data)
pca_rna  <- prcomp(t(vst_mat), center = TRUE, scale. = FALSE)
pca_var  <- summary(pca_rna)$importance
pc1_rna  <- round(pca_var[2, 1] * 100, 1)
pc2_rna  <- round(pca_var[2, 2] * 100, 1)
pc3_rna  <- round(pca_var[2, 3] * 100, 1)

pca_rna_df           <- as.data.frame(pca_rna$x[, 1:3])
pca_rna_df$Sample    <- rownames(pca_rna_df)
pca_rna_df$Strain    <- col_data[pca_rna_df$Sample, "Strain"]
pca_rna_df$Replicate <- col_data[pca_rna_df$Sample, "Replicate"]

cat("\nRNA-seq PCA variance explained:\n")
cat("PC1:", pc1_rna, "% | PC2:", pc2_rna, "% | PC3:", pc3_rna, "%\n")

p_pca_rna <- ggplot(pca_rna_df,
    aes(x = PC1, y = PC2, color = Strain, label = Sample)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text_repel(size = 3, max.overlaps = 20, show.legend = FALSE) +
  scale_color_manual(values = strain_colors) +
  labs(title = "RNA-seq PCA: PC1 vs PC2 (VST counts)",
       x = paste0("PC1 (", pc1_rna, "%)"),
       y = paste0("PC2 (", pc2_rna, "%)")) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())
ggsave(file.path(outdir, "RNAseq_QC_PCA_PC1vPC2.pdf"),
       plot = p_pca_rna, width = 7, height = 5.5, useDingbats = FALSE)

# Sample-sample distance heatmap
samp_dists    <- dist(t(vst_mat))
samp_dist_mat <- as.matrix(samp_dists)
anno_rna      <- data.frame(
  Strain    = col_data$Strain,
  Replicate = col_data$Replicate,
  row.names = col_data$Sample
)
pdf(file.path(outdir, "RNAseq_QC_sample_distance_heatmap.pdf"),
    width = 7, height = 6)
pheatmap(samp_dist_mat,
  annotation_col    = anno_rna,
  annotation_row    = anno_rna,
  annotation_colors = list(
    Strain    = strain_colors,
    Replicate = c(A = "#E8E8E8", B = "#A8A8A8", C = "#585858")
  ),
  color             = colorRampPalette(rev(brewer.pal(9, "Blues")))(100),
  clustering_distance_rows = samp_dists,
  clustering_distance_cols = samp_dists,
  clustering_method = "complete",
  fontsize = 9,
  main = "Sample-sample Euclidean distances (VST)")
dev.off()

# Dispersion plot
pdf(file.path(outdir, "RNAseq_QC_dispersion_plot.pdf"),
    width = 7, height = 5)
plotDispEsts(dds, main = "DESeq2 dispersion estimates")
dev.off()

# MA plots
pdf(file.path(outdir, "RNAseq_QC_MA_plots.pdf"),
    width = 10, height = 5)
par(mfrow = c(1, 2))
plotMA(res_a1ko, alpha = 0.05, main = "MA plot: A1KO vs WT",
       ylim = c(-5, 5))
abline(h = c(-log2(1.5), log2(1.5)), col = "grey50", lty = 2)
plotMA(res_a2ko, alpha = 0.05, main = "MA plot: A2KO vs WT",
       ylim = c(-5, 5))
abline(h = c(-log2(1.5), log2(1.5)), col = "grey50", lty = 2)
dev.off()

# Normalized count distributions
norm_counts <- counts(dds, normalized = TRUE)
norm_log    <- log2(norm_counts + 1)
norm_long   <- as.data.frame(norm_log) %>%
  tibble::rownames_to_column("Gene") %>%
  pivot_longer(cols = -Gene, names_to = "Sample",
               values_to = "log2_norm_count") %>%
  merge(col_data[, c("Sample", "Strain")], by = "Sample")
norm_long$Strain <- factor(norm_long$Strain,
                            levels = c("WT", "spxA1KO", "spxA2KO"))

p_norm <- ggplot(norm_long,
    aes(x = Sample, y = log2_norm_count, fill = Strain)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3, linewidth = 0.4) +
  scale_fill_manual(values = strain_colors) +
  labs(title = "RNA-seq normalized count distributions", x = "",
       y = "log\u2082(normalized counts + 1)") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        panel.grid.minor = element_blank())
ggsave(file.path(outdir, "RNAseq_QC_norm_distributions.pdf"),
       plot = p_norm, width = 8, height = 5, useDingbats = FALSE)

message("RNA-seq QC figures saved.")

# =============================================================================
# MODULE 2.9 — SAVE OUTPUT FILES
# =============================================================================

# ST_Validation_NanoString_RNAseq.xlsx
wb_val <- createWorkbook()

addWorksheet(wb_val, "NanoString_RNAseq_map")
writeData(wb_val, "NanoString_RNAseq_map",
          ns_rnaseq_map[order(ns_rnaseq_map$Label), ],
          headerStyle = header_style)
freezePane(wb_val, "NanoString_RNAseq_map", firstRow = TRUE)
setColWidths(wb_val, "NanoString_RNAseq_map",
             cols = 1:ncol(ns_rnaseq_map), widths = "auto")

addWorksheet(wb_val, "DE_concordance")
writeData(wb_val, "DE_concordance",
          concordance_df[order(concordance_df$Label), ],
          headerStyle = header_style)
freezePane(wb_val, "DE_concordance", firstRow = TRUE)
setColWidths(wb_val, "DE_concordance",
             cols = 1:ncol(concordance_df), widths = "auto")

addWorksheet(wb_val, "Crosswalk")
writeData(wb_val, "Crosswalk",
          crosswalk[order(crosswalk$NanoString_ProbeID), ],
          headerStyle = header_style)
freezePane(wb_val, "Crosswalk", firstRow = TRUE)
setColWidths(wb_val, "Crosswalk",
             cols = 1:ncol(crosswalk), widths = "auto")

# Validation metrics summary
val_metrics <- data.frame(
  Metric = c(
    "Probes mapped to MGAS10870",
    "Probes in expression correlation (n)",
    "WT UI expression Pearson r (NanoString vs RNA-seq TPM)",
    "WT UI R\u00b2",
    "A1KO vs WT logFC Pearson r (NanoString vs DESeq2)",
    "A2KO vs WT logFC Pearson r (NanoString vs DESeq2)",
    "A1KO directional agreement — all mapped genes",
    "A2KO directional agreement — all mapped genes",
    "A1KO directional agreement — NanoString DE genes only",
    "A2KO directional agreement — NanoString DE genes only",
    "A1KO: NanoString DE genes also significant by DESeq2",
    "A2KO: NanoString DE genes also significant by DESeq2"
  ),
  Value = c(
    paste0(sum(crosswalk$Has_C9Q), " of ", nrow(crosswalk)),
    nrow(cor_df),
    round(r_wt, 3),
    round(r2_wt, 3),
    round(r_deseq_a1ko, 3),
    round(r_deseq_a2ko, 3),
    sprintf("%.1f%% (%d/%d)",
      100 * mean(concordance_df$direction_concordance_spxA1KO_DESeq2, na.rm = TRUE),
      sum(concordance_df$direction_concordance_spxA1KO_DESeq2, na.rm = TRUE),
      sum(!is.na(concordance_df$direction_concordance_spxA1KO_DESeq2))),
    sprintf("%.1f%% (%d/%d)",
      100 * mean(concordance_df$direction_concordance_spxA2KO_DESeq2, na.rm = TRUE),
      sum(concordance_df$direction_concordance_spxA2KO_DESeq2, na.rm = TRUE),
      sum(!is.na(concordance_df$direction_concordance_spxA2KO_DESeq2))),
    sprintf("%.1f%% (%d/%d)",
      100 * mean(ns_de_a1ko$direction_concordance_spxA1KO_DESeq2, na.rm = TRUE),
      sum(ns_de_a1ko$direction_concordance_spxA1KO_DESeq2, na.rm = TRUE),
      sum(!is.na(ns_de_a1ko$direction_concordance_spxA1KO_DESeq2))),
    sprintf("%.1f%% (%d/%d)",
      100 * mean(ns_de_a2ko$direction_concordance_spxA2KO_DESeq2, na.rm = TRUE),
      sum(ns_de_a2ko$direction_concordance_spxA2KO_DESeq2, na.rm = TRUE),
      sum(!is.na(ns_de_a2ko$direction_concordance_spxA2KO_DESeq2))),
    sprintf("%.1f%% (%d/%d)",
      100 * nrow(both_sig_a1ko) / nrow(ns_de_a1ko),
      nrow(both_sig_a1ko), nrow(ns_de_a1ko)),
    sprintf("%.1f%% (%d/%d)",
      100 * nrow(both_sig_a2ko) / nrow(ns_de_a2ko),
      nrow(both_sig_a2ko), nrow(ns_de_a2ko))
  ),
  stringsAsFactors = FALSE
)
addWorksheet(wb_val, "Validation_metrics")
writeData(wb_val, "Validation_metrics",
          val_metrics, headerStyle = header_style)
setColWidths(wb_val, "Validation_metrics",
             cols = 1:2, widths = c(60, 25))

saveWorkbook(wb_val,
  file.path(outdir, "ST_Validation_NanoString_RNAseq.xlsx"),
  overwrite = TRUE)

# DESeq2_RNAseq_results.xlsx
wb_deseq <- createWorkbook()

addWorksheet(wb_deseq, "DESeq2_spxA1KO_vs_WT")
writeData(wb_deseq, "DESeq2_spxA1KO_vs_WT",
          res_a1ko_df[order(res_a1ko_df$padj), ],
          headerStyle = header_style)
freezePane(wb_deseq, "DESeq2_spxA1KO_vs_WT", firstRow = TRUE)
setColWidths(wb_deseq, "DESeq2_spxA1KO_vs_WT",
             cols = 1:ncol(res_a1ko_df), widths = "auto")

addWorksheet(wb_deseq, "DESeq2_spxA2KO_vs_WT")
writeData(wb_deseq, "DESeq2_spxA2KO_vs_WT",
          res_a2ko_df[order(res_a2ko_df$padj), ],
          headerStyle = header_style)
freezePane(wb_deseq, "DESeq2_spxA2KO_vs_WT", firstRow = TRUE)
setColWidths(wb_deseq, "DESeq2_spxA2KO_vs_WT",
             cols = 1:ncol(res_a2ko_df), widths = "auto")

addWorksheet(wb_deseq, "Concordance_NanoString_DESeq2")
writeData(wb_deseq, "Concordance_NanoString_DESeq2",
          concordance_df[order(concordance_df$Label), ],
          headerStyle = header_style)
freezePane(wb_deseq, "Concordance_NanoString_DESeq2", firstRow = TRUE)

# Prism-ready concordance scatter data (DESeq2 logFC on x-axis)
for (strain in c("spxA1KO", "spxA2KO")) {
  x_col  <- paste0("DESeq2_log2FC_", strain)
  y_col  <- paste0("NS_logFC_", strain, "_vs_WT")
  sig_col <- paste0("NS_Sig_", strain)
  prism_df <- concordance_df[
    !is.na(concordance_df[[x_col]]) &
    !is.na(concordance_df[[y_col]]),
    c("Label", y_col, x_col, sig_col,
      paste0("DESeq2_Sig_", strain))]
  colnames(prism_df) <- c("Gene", "NanoString_logFC",
                            "DESeq2_logFC",
                            "Sig_NanoString", "Sig_DESeq2")
  prism_df <- prism_df[order(prism_df$Gene), ]
  sheet_nm <- paste0("Prism_", strain, "_concordance")
  addWorksheet(wb_deseq, sheet_nm)
  writeData(wb_deseq, sheet_nm, prism_df, headerStyle = header_style)
  setColWidths(wb_deseq, sheet_nm, cols = 1:5, widths = "auto")
}

saveWorkbook(wb_deseq,
  file.path(outdir, "DESeq2_RNAseq_results.xlsx"),
  overwrite = TRUE)

message("Output files saved.")
message("  ST_Validation_NanoString_RNAseq.xlsx")
message("  DESeq2_RNAseq_results.xlsx")

# =============================================================================
# MODULE 2.10 — UPDATE WORKSPACE
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
  locus_key, tpm_wide, ns_rnaseq_map,
  concordance_df, crosswalk,
  c9q_to_name, spym3_map, name_to_c9q,
  dds, res_a1ko_df, res_a2ko_df,
  deseq_a1ko_mapped, deseq_a2ko_mapped,
  r_wt, r2_wt, r_deseq_a1ko, r_deseq_a2ko,
  ns_de_a1ko, ns_de_a2ko,
  both_sig_a1ko, both_sig_a2ko,
  strain_colors, outdir, data_dir,
  file = file.path(outdir, "NanoString_workspace_noMo.RData")
)

message("Workspace updated: NanoString_workspace_noMo.RData")
message("\nModule 2 complete.")
message("\n=== FINAL VALIDATION METRICS ===")
print(val_metrics)

# =============================================================================
# SESSION INFO
# =============================================================================
cat("\n=== Session Info ===\n")
writeLines(capture.output(sessionInfo()),
           file.path(OUT_DIR, "RNAseq_sessionInfo.txt"))
sessionInfo()
