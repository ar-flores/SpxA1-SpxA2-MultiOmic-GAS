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
# ChIP-exo Analysis — CovR Differential Binding
# =============================================================================
# Description:
#   DiffBind/DESeq2 differential CovR binding analysis, composite visualization,
#   FIMO motif enrichment, priority target annotation, and figure generation.
#
#   NOTE ON SpxA2OE STRAIN:
#   The "SpxA2OE" condition uses the LiaSQ146A (liaS-Q146A) isogenic mutant,
#   which encodes a constitutively active LiaS sensor kinase that drives high
#   SpxA2 expression independent of AMP induction. This strain is used as a
#   SpxA2 overexpression surrogate throughout the analysis. BioSample accessions:
#   SAMN60523488 (rep1), SAMN60523489 (rep2), SAMN60523490 (rep3).
#
# Input files:
#   BAM files in ChIPexo/bam/ (9 samples)
#   Peak files in ChIPexo/peaks/ (narrowPeak format, from chipexo_pipeline.sh)
#   Composite .out files in ChIPexo/output_summits/ (from ScriptManager)
#   FIMO results in ChIPexo/motif_analysis/
#
# Output files (ChIPexo/results/):
#   plots/    — composite, priority targets, motif composition, QC figures
#   tables/   — ST1-ST5 supplemental tables, Prism-ready data
#   chipexo_analysis_final.RData
#   sessionInfo.txt
# =============================================================================

library(DiffBind)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(tidyr)

# Set working directory
# =============================================================================
# SECTION 0: PATHS — set BASE_DIR to repository root before running
# =============================================================================
BASE_DIR <- "."
CHIP_DIR <- file.path(BASE_DIR, "ChIPexo")
BAM_DIR  <- file.path(CHIP_DIR, "bam")
PEAK_DIR <- file.path(CHIP_DIR, "peaks")
OUT_DIR  <- file.path(CHIP_DIR, "results")
MOTIF_DIR<- file.path(CHIP_DIR, "motif_analysis")
SUMM_DIR <- file.path(CHIP_DIR, "output_summits")
dir.create(file.path(OUT_DIR, "plots"),  showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(OUT_DIR, "tables"), showWarnings=FALSE, recursive=TRUE)

# Create output directories
dir.create("diffbind_results", showWarnings=FALSE)
dir.create(file.path(OUT_DIR, "plots", showWarnings=FALSE)
dir.create(file.path(OUT_DIR, "tables", showWarnings=FALSE)

# =============================================================================
# SECTION 1: DIFFBIND - LOAD DATA
# =============================================================================

# Create sample sheet
# NOTE: Each condition uses the same merged peak file for all replicates
# This is intentional - peaks called on merged BAMs for maximum sensitivity
samplesheet <- data.frame(
    SampleID = c("WT_rep1","WT_rep2","WT_rep3",
                 "SpxA2KO_rep1","SpxA2KO_rep2","SpxA2KO_rep3",
                 "SpxA2OE_rep1","SpxA2OE_rep2","SpxA2OE_rep3"),
    Condition = c("WT","WT","WT",
                  "SpxA2KO","SpxA2KO","SpxA2KO",
                  "SpxA2OE","SpxA2OE","SpxA2OE"),
    Replicate = c(1,2,3,1,2,3,1,2,3),
    bamReads = c(
        file.path(BAM_DIR, "WT_rep1.bam",
        file.path(BAM_DIR, "WT_rep2.bam",
        file.path(BAM_DIR, "WT_rep3.bam",
        file.path(BAM_DIR, "SpxA2KO_rep1.bam",
        file.path(BAM_DIR, "SpxA2KO_rep2.bam",
        file.path(BAM_DIR, "SpxA2KO_rep3.bam",
        file.path(BAM_DIR, "SpxA2OE_rep1.bam",
        file.path(BAM_DIR, "SpxA2OE_rep2.bam",
        file.path(BAM_DIR, "SpxA2OE_rep3.bam"),
    Peaks = c(
        file.path(PEAK_DIR, "WT_diffbind_peaks.narrowPeak",
        file.path(PEAK_DIR, "WT_diffbind_peaks.narrowPeak",
        file.path(PEAK_DIR, "WT_diffbind_peaks.narrowPeak",
        file.path(PEAK_DIR, "SpxA2KO_diffbind_peaks.narrowPeak",
        file.path(PEAK_DIR, "SpxA2KO_diffbind_peaks.narrowPeak",
        file.path(PEAK_DIR, "SpxA2KO_diffbind_peaks.narrowPeak",
        file.path(PEAK_DIR, "SpxA2OE_diffbind_peaks.narrowPeak",
        file.path(PEAK_DIR, "SpxA2OE_diffbind_peaks.narrowPeak",
        file.path(PEAK_DIR, "SpxA2OE_diffbind_peaks.narrowPeak"),
    PeakCaller = rep("narrowPeak", 9)
)

# Load into DiffBind
dba <- dba(sampleSheet=samplesheet)
print(dba)
# Expected: 9 Samples, 701 sites in matrix

# =============================================================================
# SECTION 2: DIFFBIND - COUNT READS
# =============================================================================

# Count reads at consensus peaks
# summits=75 creates 150bp windows centered on peak summits
dba <- dba.count(dba,
                 summits=75,
                 bParallel=TRUE)
print(dba)
# Expected FRiP: 0.21-0.24 across all samples (high quality)

# =============================================================================
# SECTION 3: DIFFBIND - NORMALIZE
# =============================================================================

# Library-size normalization using DESeq2
# NOTE: No input control available - ChIP-exo strand specificity
# provides inherent background discrimination (Rhee & Pugh, 2011)
dba <- dba.normalize(dba,
                     method=DBA_DESEQ2,
                     normalize=DBA_NORM_LIB,
                     library=DBA_LIBSIZE_FULL)

# Check normalization factors
dba.normalize(dba, bRetrieve=TRUE)

# =============================================================================
# SECTION 4: DIFFBIND - DIFFERENTIAL BINDING ANALYSIS
# =============================================================================

# Set contrasts with WT as reference
dba <- dba.contrast(dba,
                    reorderMeta=list(Condition="WT"))
print(dba)
# Expected: 3 contrasts (WT vs KO, WT vs OE, OE vs KO)

# Run DESeq2 differential binding analysis
dba <- dba.analyze(dba,
                   method=DBA_DESEQ2,
                   bParallel=TRUE)

# Show summary
dba.show(dba, bContrasts=TRUE)
# Expected: WT vs KO = 439, WT vs OE = 76, OE vs KO = 378

# Extract results (FDR < 0.05, no fold change filter)
res_WT_KO <- dba.report(dba, contrast=1, th=0.05, bUsePval=FALSE)
res_WT_OE <- dba.report(dba, contrast=2, th=0.05, bUsePval=FALSE)
res_OE_KO <- dba.report(dba, contrast=3, th=0.05, bUsePval=FALSE)

cat("WT vs KO differential sites:", length(res_WT_KO), "\n")
cat("WT vs OE differential sites:", length(res_WT_OE), "\n")
cat("OE vs KO differential sites:", length(res_OE_KO), "\n")

# Convert to dataframes and add direction
df_WT_KO <- as.data.frame(res_WT_KO)
df_WT_OE <- as.data.frame(res_WT_OE)
df_OE_KO <- as.data.frame(res_OE_KO)

df_WT_KO$direction <- ifelse(df_WT_KO$Fold > 0, "More_in_WT", "More_in_KO")
df_WT_OE$direction <- ifelse(df_WT_OE$Fold > 0, "More_in_WT", "More_in_OE")

cat("WT vs KO directions:\n")
print(table(df_WT_KO$direction))
cat("\nWT vs OE directions:\n")
print(table(df_WT_OE$direction))

# =============================================================================
# SECTION 5: PRIORITY TARGET ANNOTATION
# =============================================================================

# Known CovR priority targets with updated NZ_CP067090.1 coordinates
# Note: Locus tags updated from Horstmann 2023 SpyM3 nomenclature
priority_coords <- data.frame(
    gene = c("covR","ska","hasA","scpC","mga","nra","slo",
             "ideS","sagA","codY","ahpC","brnQ","emm","scl1","scl2"),
    start = c(283636,1677322,1843888,334612,1711522,115763,153255,
              685851,573409,1497347,1755232,275221,1709593,1681671,795619),
    end   = c(284323,1678645,1845148,339559,1713133,117299,154971,
              686871,573571,1498130,1755793,276592,1710973,1682517,797445)
)

# Annotation function
annotate_peaks <- function(df, priority_coords, window=300){
    df$gene <- NA
    for(i in 1:nrow(priority_coords)){
        gene   <- priority_coords$gene[i]
        gstart <- priority_coords$start[i]
        gend   <- priority_coords$end[i]
        idx <- which(df$start <= (gend + window) &
                     df$end   >= (gstart - window))
        if(length(idx) > 0) df$gene[idx] <- gene
    }
    return(df)
}

# Annotate results
df_WT_KO_ann <- annotate_peaks(df_WT_KO, priority_coords)
df_WT_OE_ann <- annotate_peaks(df_WT_OE, priority_coords)

# Print priority target results
cat("Priority CovR targets - WT vs KO:\n")
print(df_WT_KO_ann[!is.na(df_WT_KO_ann$gene),
                    c("seqnames","start","end","Fold","FDR","direction","gene")])

cat("\nPriority CovR targets - WT vs OE:\n")
print(df_WT_OE_ann[!is.na(df_WT_OE_ann$gene),
                    c("seqnames","start","end","Fold","FDR","direction","gene")])

# =============================================================================
# SECTION 6: SAVE RESULTS
# =============================================================================

# Save annotated results
write.csv(df_WT_KO_ann,
          file.path(OUT_DIR, "WT_vs_SpxA2KO_annotated.csv",
          row.names=TRUE)
write.csv(df_WT_OE_ann,
          file.path(OUT_DIR, "WT_vs_SpxA2OE_annotated.csv",
          row.names=TRUE)
write.csv(as.data.frame(res_OE_KO),
          file.path(OUT_DIR, "SpxA2OE_vs_SpxA2KO.csv",
          row.names=TRUE)

# Export peak BED files for motif analysis
# KO gained (more binding in KO than WT)
df_KO_gained <- df_WT_KO[df_WT_KO$direction=="More_in_KO",]
df_WT_gained <- df_WT_KO[df_WT_KO$direction=="More_in_WT",]

write.table(data.frame(
    chr=df_KO_gained$seqnames, start=df_KO_gained$start,
    end=df_KO_gained$end,
    name=paste0("peak_", 1:nrow(df_KO_gained)), score=0, strand="."),
    file.path(MOTIF_DIR, "KO_gained_peaks.bed",
    sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

write.table(data.frame(
    chr=df_WT_gained$seqnames, start=df_WT_gained$start,
    end=df_WT_gained$end,
    name=paste0("peak_", 1:nrow(df_WT_gained)), score=0, strand="."),
    file.path(MOTIF_DIR, "WT_gained_peaks.bed",
    sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

cat("Results saved\n")

# =============================================================================
# SECTION 7: QC FIGURES
# =============================================================================

# PCA plot
png(file.path(OUT_DIR, "plots/PCA_plot.png", width=800, height=600)
dba.plotPCA(dba, attributes=DBA_CONDITION, label=DBA_ID)
dev.off()

# Correlation heatmap
png(file.path(OUT_DIR, "plots/correlation_heatmap.png", width=900, height=800)
dba.plotHeatmap(dba, attributes=DBA_ID, colScheme="Greens")
dev.off()

# MA plots
png(file.path(OUT_DIR, "plots/MA_WT_vs_KO.png", width=800, height=600)
dba.plotMA(dba, contrast=1)
dev.off()

png(file.path(OUT_DIR, "plots/MA_WT_vs_OE.png", width=800, height=600)
dba.plotMA(dba, contrast=2)
dev.off()

cat("QC figures saved\n")

# =============================================================================
# SECTION 8: PRIORITY TARGETS BARPLOT (Figure 2)
# =============================================================================

priority_summary <- df_WT_KO_ann %>%
    filter(!is.na(gene)) %>%
    group_by(gene) %>%
    slice_min(FDR) %>%
    ungroup() %>%
    arrange(Fold)

png(file.path(OUT_DIR, "plots/priority_targets_WT_vs_KO.png",
    width=1000, height=600, res=150)

ggplot(priority_summary,
       aes(x=reorder(gene, Fold), y=Fold, fill=direction)) +
    geom_bar(stat="identity") +
    geom_hline(yintercept=0, linetype="dashed") +
    scale_fill_manual(
        values=c("More_in_KO"="#E41A1C", "More_in_WT"="#377EB8"),
        labels=c("More_in_KO"="Higher in KO", "More_in_WT"="Higher in WT")) +
    labs(title="CovR Binding at Priority Targets: WT vs SpxA2 KO",
         x="Gene",
         y="Log2 Fold Change (WT - SpxA2KO)",
         fill="Direction") +
    theme_bw() +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=12),
          axis.text.y=element_text(size=12),
          axis.title=element_text(size=14),
          plot.title=element_text(size=16, hjust=0.5),
          legend.text=element_text(size=10)) +
    geom_text(aes(label=ifelse(FDR<0.001, "***",
                        ifelse(FDR<0.01, "**",
                        ifelse(FDR<0.05, "*", "")))),
              vjust=ifelse(priority_summary$Fold>0, -0.5, 1.5),
              size=5)

dev.off()
cat("Priority targets figure saved\n")

# =============================================================================
# SECTION 9: MOTIF COMPOSITION FIGURE (Figure 3)
# =============================================================================
# NOTE: Run AFTER bash script Section 6 (FIMO analysis)

# Load FIMO results
ko_dimer   <- read.table(file.path(MOTIF_DIR, "fimo_KO_final/fimo.tsv",
                          header=TRUE, comment.char="#", sep="\t")
ko_monomer <- read.table(file.path(MOTIF_DIR, "fimo_KO_monomer/fimo.tsv",
                          header=TRUE, comment.char="#", sep="\t")
wt_dimer   <- read.table(file.path(MOTIF_DIR, "fimo_WT_final/fimo.tsv",
                          header=TRUE, comment.char="#", sep="\t")
wt_monomer <- read.table(file.path(MOTIF_DIR, "fimo_WT_monomer/fimo.tsv",
                          header=TRUE, comment.char="#", sep="\t")

# Get unique sequence names with hits
ko_dimer_seqs   <- unique(ko_dimer$sequence_name)
ko_monomer_seqs <- unique(ko_monomer$sequence_name)
wt_dimer_seqs   <- unique(wt_dimer$sequence_name)
wt_monomer_seqs <- unique(wt_monomer$sequence_name)

# Calculate overlap categories
ko_both        <- sum(ko_dimer_seqs %in% ko_monomer_seqs)
ko_dimer_only  <- sum(!ko_dimer_seqs %in% ko_monomer_seqs)
ko_monomer_only<- sum(!ko_monomer_seqs %in% ko_dimer_seqs)
ko_neither     <- 198 - length(union(ko_dimer_seqs, ko_monomer_seqs))

wt_both        <- sum(wt_dimer_seqs %in% wt_monomer_seqs)
wt_dimer_only  <- sum(!wt_dimer_seqs %in% wt_monomer_seqs)
wt_monomer_only<- sum(!wt_monomer_seqs %in% wt_dimer_seqs)
wt_neither     <- 241 - length(union(wt_dimer_seqs, wt_monomer_seqs))

cat("KO gained motif composition:\n")
cat("  Dimer only:", ko_dimer_only, "\n")
cat("  Both:", ko_both, "\n")
cat("  Monomer only:", ko_monomer_only, "\n")
cat("  Neither:", ko_neither, "\n")

cat("\nWT gained motif composition:\n")
cat("  Dimer only:", wt_dimer_only, "\n")
cat("  Both:", wt_both, "\n")
cat("  Monomer only:", wt_monomer_only, "\n")
cat("  Neither:", wt_neither, "\n")

# Statistical tests
# Dimer motif enrichment
dimer_table <- matrix(c(ko_dimer_only + ko_both, ko_neither + ko_monomer_only,
                         wt_dimer_only + wt_both, wt_neither + wt_monomer_only),
                       nrow=2,
                       dimnames=list(c("KO_gained","WT_gained"),
                                     c("Has_dimer","No_dimer")))
fisher_dimer <- fisher.test(dimer_table)
cat("\nDimer motif enrichment (KO vs WT gained):\n")
cat("KO:", round((ko_dimer_only+ko_both)/198*100,1), "%\n")
cat("WT:", round((wt_dimer_only+wt_both)/241*100,1), "%\n")
cat("p-value:", fisher_dimer$p.value, "\n")
cat("Odds ratio:", fisher_dimer$estimate, "\n")

# Monomer-only enrichment
monomer_only_table <- matrix(c(ko_monomer_only, 198-ko_monomer_only,
                                wt_monomer_only, 241-wt_monomer_only),
                              nrow=2,
                              dimnames=list(c("KO_gained","WT_gained"),
                                            c("Monomer_only","Other")))
fisher_monomer <- fisher.test(monomer_only_table)
cat("\nMonomer-only enrichment (WT vs KO gained):\n")
cat("KO:", round(ko_monomer_only/198*100,1), "%\n")
cat("WT:", round(wt_monomer_only/241*100,1), "%\n")
cat("p-value:", fisher_monomer$p.value, "\n")
cat("Odds ratio:", fisher_monomer$estimate, "\n")

# Motif composition stacked bar plot
motif_composition <- data.frame(
    Class = c(rep("KO gained\n(n=198)", 4), rep("WT gained\n(n=241)", 4)),
    Category = rep(c("Dimer only","Both motifs","Monomer only","Neither"), 2),
    Count   = c(ko_dimer_only, ko_both, ko_monomer_only, ko_neither,
                wt_dimer_only, wt_both, wt_monomer_only, wt_neither),
    Percent = c(round(ko_dimer_only/198*100,1), round(ko_both/198*100,1),
                round(ko_monomer_only/198*100,1), round(ko_neither/198*100,1),
                round(wt_dimer_only/241*100,1), round(wt_both/241*100,1),
                round(wt_monomer_only/241*100,1), round(wt_neither/241*100,1))
)

motif_composition$Category <- factor(
    motif_composition$Category,
    levels=c("Neither","Monomer only","Both motifs","Dimer only"))

png(file.path(OUT_DIR, "plots/motif_composition.png",
    width=800, height=600, res=150)

ggplot(motif_composition,
       aes(x=Class, y=Percent, fill=Category)) +
    geom_bar(stat="identity", width=0.6) +
    scale_fill_manual(values=c(
        "Dimer only"   ="#2166AC",
        "Both motifs"  ="#92C5DE",
        "Monomer only" ="#F4A582",
        "Neither"      ="#BABABA")) +
    labs(title="CovR Binding Motif Composition\nby Peak Class",
         x="", y="Percentage of peaks (%)", fill="Motif") +
    theme_bw() +
    theme(axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title=element_text(size=14),
          plot.title=element_text(size=14, hjust=0.5),
          legend.text=element_text(size=10)) +
    geom_text(aes(label=paste0(Percent,"%")),
              position=position_stack(vjust=0.5),
              size=3.5, color="black")

dev.off()
cat("Motif composition figure saved\n")

# =============================================================================
# SECTION 10: COMPOSITE PLOT (Figure 1)
# =============================================================================
# NOTE: Run AFTER bash script Section 5 (ScriptManager TagPileup)

read_composite_wide <- function(file, condition){
    df <- read.table(file, header=TRUE, sep="\t",
                     row.names=1, check.names=FALSE)
    positions <- as.numeric(colnames(df))
    result <- data.frame(
        position  = positions,
        sense     = as.numeric(df[1,]),
        anti      = as.numeric(df[2,]),
        condition = condition
    )
    return(result)
}

wt <- read_composite_wide(file.path(SUMM_DIR, "WT_merged_composite.out",    "WT")
ko <- read_composite_wide(file.path(SUMM_DIR, "SpxA2KO_composite.out", "SpxA2KO")
oe <- read_composite_wide(file.path(SUMM_DIR, "SpxA2OE_composite.out", "SpxA2OE")

all_data <- rbind(wt, ko, oe)
all_data$condition <- factor(all_data$condition,
                              levels=c("WT","SpxA2KO","SpxA2OE"))

all_data_long <- all_data %>%
    pivot_longer(cols=c("sense","anti"),
                 names_to="strand", values_to="count")

png(file.path(OUT_DIR, "plots/composite_combined.png",
    width=1000, height=700, res=150)

ggplot(all_data_long,
       aes(x=position, y=count, color=condition, linetype=strand)) +
    geom_line(size=0.8) +
    scale_color_manual(values=c("WT"="#2166AC",
                                 "SpxA2KO"="#D6604D",
                                 "SpxA2OE"="#4DAC26")) +
    scale_linetype_manual(values=c("sense"="solid", "anti"="dashed"),
                           labels=c("sense"="Sense (+)", "anti"="Anti (-)")) +
    labs(title="CovR ChIP-exo Composite Plot",
         subtitle="Averaged across 362 WT peak summits \u00b1 100bp",
         x="Distance from Peak Summit (bp)",
         y="Average Tag Count",
         color="Condition", linetype="Strand") +
    theme_bw() +
    theme(axis.text=element_text(size=11),
          axis.title=element_text(size=13),
          plot.title=element_text(size=14, hjust=0.5),
          plot.subtitle=element_text(size=11, hjust=0.5),
          legend.text=element_text(size=11)) +
    geom_vline(xintercept=0, linetype="dotted", color="gray50", alpha=0.5)

dev.off()
cat("Composite figure saved\n")

# =============================================================================
# SECTION 11: EXPORT TABLES FOR PRISM
# =============================================================================

# Priority targets table
priority_export <- priority_summary %>%
    select(gene, Fold, FDR, direction) %>%
    mutate(
        neglog10_FDR   = -log10(FDR),
        significance   = ifelse(FDR<0.001, "***",
                         ifelse(FDR<0.01,  "**",
                         ifelse(FDR<0.05,  "*", "ns")))
    ) %>%
    arrange(Fold)

write.csv(priority_export,
          file.path(OUT_DIR, "tables/priority_targets_table.csv",
          row.names=FALSE)

# Motif composition table
write.csv(motif_composition,
          file.path(OUT_DIR, "tables/motif_composition_table.csv",
          row.names=FALSE)

# Composite data for Prism
prism_composite <- data.frame(
    position        = wt$position,
    WT_sense        = wt$sense,     WT_anti        = wt$anti,
    SpxA2KO_sense   = ko$sense,     SpxA2KO_anti   = ko$anti,
    SpxA2OE_sense   = oe$sense,     SpxA2OE_anti   = oe$anti
)
write.csv(prism_composite,
          file.path(OUT_DIR, "tables/composite_data_prism.csv",
          row.names=FALSE)

# ST1: All 700 consensus peaks
all_peaks    <- dba.peakset(dba, bRetrieve=TRUE)
all_peaks_df <- as.data.frame(all_peaks)
write.csv(all_peaks_df,
          file.path(OUT_DIR, "tables/ST1_all_700_peaks.csv",
          row.names=FALSE)

# ST4: FIMO hits
make_hits_df <- function(df, motif, peak_class){
    out <- df[, c("sequence_name","start","stop","strand",
                  "score","p.value","q.value","matched_sequence")]
    out$motif      <- motif
    out$peak_class <- peak_class
    return(out)
}

ST4 <- rbind(
    make_hits_df(ko_dimer,   "dimer",   "KO_gained"),
    make_hits_df(wt_dimer,   "dimer",   "WT_gained"),
    make_hits_df(ko_monomer, "monomer", "KO_gained"),
    make_hits_df(wt_monomer, "monomer", "WT_gained")
)
write.csv(ST4,
          file.path(OUT_DIR, "tables/ST4_FIMO_hits.csv",
          row.names=FALSE)

cat("All tables exported\n")

# =============================================================================
# SECTION 12: SAVE SESSION
# =============================================================================

save.image(file.path(OUT_DIR, "chipexo_analysis_final.RData")
writeLines(capture.output(sessionInfo()),
           file.path(OUT_DIR, "sessionInfo.txt")

cat("\n=== Analysis complete ===\n")
cat("Figures saved to: diffbind_results/plots/\n")
cat("Tables saved to:  diffbind_results/tables/\n")
cat("Full results in:  diffbind_results/\n")
