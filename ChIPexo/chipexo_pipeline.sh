#!/bin/bash
# =============================================================================
# SpxA2-CovR ChIP-exo Analysis Pipeline
# =============================================================================
# Study: SpxA2 modulates CovR DNA occupancy in Streptococcus pyogenes MGAS10870
# Strain: emm3 GAS MGAS10870 (NZ_CP067090.1)
# Conditions: WT, SpxA2KO (deltaSpxA2), SpxA2OE
# Replicates: 3 per condition (9 total)
# Author: Anthony R. Flores MD MPH PhD
# Institution: Vanderbilt University Medical Center
# Contact:     flores-lab.org
# Journal:     mBio (in submission)
# GEO:         PRJNA1472884
# GitHub:      https://github.com/flores-lab/SpxA1-SpxA2-MultiOmic-GAS
# =============================================================================

# =============================================================================
# SECTION 0: ENVIRONMENT SETUP
# =============================================================================

# Software requirements:
# - conda environment "chipexo": macs3, deeptools, meme-suite, bedtools, samtools, biopython
# - R v4.5.3 with DiffBind v3.20, DESeq2 v1.50.2, ggplot2, dplyr, tidyr, pheatmap
# - ScriptManager v0.14 (https://github.com/CEGRcode/scriptmanager)
# - Java (for ScriptManager)

# Activate conda environment
conda activate chipexo

# Set base directory
BASEDIR="${BASE_DIR:-$(pwd)/ChIPexo}"  # override with BASE_DIR env variable if set
BAMDIR=${BASEDIR}/841-Flores
PEAKDIR=${BASEDIR}/peaks
OUTPUTDIR=${BASEDIR}/output_summits
MOTIFDIR=${BASEDIR}/motif_analysis
SCRIPTMANAGER="${SCRIPTMANAGER_PATH:-$HOME/gensoft/ScriptManager-v0.14.jar}"

# =============================================================================
# SECTION 1: REFERENCE GENOME PREPARATION
# =============================================================================

# Generate BED file from NCBI GenBank file (NZ_CP067090.1)
# Handles reannotated RefSeq locus tags (C9Q_RS prefix)
python3 << 'EOF'
from Bio import SeqIO

with open("CP067090.gbk") as in_handle, \
     open("CP067090_final.bed", "w") as out_handle:
    for record in SeqIO.parse(in_handle, "genbank"):
        for feature in record.features:
            if feature.type in ("gene", "pseudogene"):
                chrom = "NZ_CP067090.1"
                start = int(feature.location.start)
                end = int(feature.location.end)
                strand = "+" if feature.location.strand == 1 else "-"
                # Gene name first, fall back to RS locus tag
                gene_name = feature.qualifiers.get("gene", [None])[0]
                locus_tags = feature.qualifiers.get("locus_tag", [])
                rs_tag = next((t for t in locus_tags if "RS" in t), None)
                old_tag = next((t for t in locus_tags if "RS" not in t), None)
                name = gene_name or rs_tag or old_tag or "unknown"
                out_handle.write(f"{chrom}\t{start}\t{end}\t{name}\t0\t{strand}\n")
EOF

# Manually add scl2 (pseudogene not captured by standard extraction)
echo -e "NZ_CP067090.1\t795619\t797445\tscl2\t0\t+" >> CP067090_final.bed

# Generate FASTA from GenBank
python3 << 'EOF'
from Bio import SeqIO
with open("CP067090.gbk") as in_handle, \
     open("NZ_CP067090.1.fasta", "w") as out_handle:
    for record in SeqIO.parse(in_handle, "genbank"):
        record.id = "NZ_CP067090.1"
        record.description = "Streptococcus pyogenes MGAS10870"
        SeqIO.write(record, out_handle, "fasta")
EOF

echo "Reference preparation complete"

# =============================================================================
# SECTION 2: BAM FILE RENAMING
# =============================================================================
# Cornell core naming convention: XXXXX_CovR_custom_MGAS10870_-_THY_-_BX_FilteredBAM.bam

cd ${BAMDIR}

# WT replicates (43721-43723)
for i in 1 2 3; do
    num=$((43720 + i))
    mv "${num}_CovR_custom_MGAS10870_-_THY_-_BX_FilteredBAM.bam" "WT_rep${i}.bam"
    mv "${num}_CovR_custom_MGAS10870_-_THY_-_BX_FilteredBAM.bam.bai" "WT_rep${i}.bam.bai"
done

# SpxA2KO replicates (43724-43726)
for i in 1 2 3; do
    num=$((43723 + i))
    mv "${num}_CovR_custom_MGAS10870_-_THY_-_BX_FilteredBAM.bam" "SpxA2KO_rep${i}.bam"
    mv "${num}_CovR_custom_MGAS10870_-_THY_-_BX_FilteredBAM.bam.bai" "SpxA2KO_rep${i}.bam.bai"
done

# NOTE: SpxA2OE condition uses the LiaSQ146A (liaS-Q146A) isogenic mutant,
# which drives high SpxA2 expression as a constitutively active LiaS sensor
# kinase. This is used as a SpxA2 overexpression surrogate throughout.
# BioSample accessions: SAMN60523488-490 (reps 1-3)
# SpxA2OE replicates (43727-43729; LiaSQ146A strain)
for i in 1 2 3; do
    num=$((43726 + i))
    mv "${num}_CovR_custom_MGAS10870_-_THY_-_BX_FilteredBAM.bam" "SpxA2OE_rep${i}.bam"
    mv "${num}_CovR_custom_MGAS10870_-_THY_-_BX_FilteredBAM.bam.bai" "SpxA2OE_rep${i}.bam.bai"
done

echo "BAM renaming complete"

# =============================================================================
# SECTION 3: BAM MERGING AND INDEXING
# =============================================================================

cd ${BAMDIR}

# Merge replicates per condition
for strain in WT SpxA2KO SpxA2OE; do
    echo "Merging ${strain}..."
    samtools merge -f -@ 54 \
        ${strain}_merged.bam \
        ${strain}_rep1.bam \
        ${strain}_rep2.bam \
        ${strain}_rep3.bam &
done
wait

# Index merged BAMs
for strain in WT SpxA2KO SpxA2OE; do
    samtools index -@ 54 ${strain}_merged.bam &
done
wait

# QC check
for strain in WT SpxA2KO SpxA2OE; do
    echo "=== ${strain}_merged ==="
    samtools flagstat ${strain}_merged.bam
done

echo "BAM merging complete"

# =============================================================================
# SECTION 4: PEAK CALLING (MACS3)
# =============================================================================

mkdir -p ${PEAKDIR}

# Call peaks on merged BAMs
# Note: MACS3 installed in chipexo conda environment
# Note: --keep-dup all because Cornell core already removed duplicates
for strain in WT SpxA2KO SpxA2OE; do
    echo "Calling peaks for ${strain}..."
    macs3 callpeak \
        -t ${BAMDIR}/${strain}_merged.bam \
        -f BAM \
        -g 1840000 \
        -n ${strain}_merged_v2 \
        --outdir ${PEAKDIR} \
        --nomodel \
        --extsize 50 \
        --shift -25 \
        -p 1e-10 \
        --min-length 50 \
        --keep-dup all \
        2>&1 | tee ${PEAKDIR}/${strain}_macs3.log &
done
wait

# Filter peaks by score
# score > 5000 for visualization (362/425/355 peaks per condition)
# score > 3000 for DiffBind input (541/632/547 peaks per condition)
for strain in WT SpxA2KO SpxA2OE; do
    # Visualization peaks
    awk '$5 > 5000' \
        ${PEAKDIR}/${strain}_merged_v2_peaks.narrowPeak \
        > ${PEAKDIR}/${strain}_filtered_peaks.narrowPeak

    # DiffBind input peaks
    awk '$5 > 3000' \
        ${PEAKDIR}/${strain}_merged_v2_peaks.narrowPeak \
        > ${PEAKDIR}/${strain}_diffbind_peaks.narrowPeak

    echo -n "${strain} filtered peaks: "
    wc -l ${PEAKDIR}/${strain}_filtered_peaks.narrowPeak
done

echo "Peak calling complete"

# =============================================================================
# SECTION 5: CHIP-EXO COMPOSITE VISUALIZATION (ScriptManager)
# =============================================================================

mkdir -p ${OUTPUTDIR}

# Create summit-centered BED from WT filtered peaks
# Column 10 of narrowPeak is summit offset from start
awk '{print $1"\t"($2+$10-100)"\t"($2+$10+100)"\t"$4"\t"$5"\t+"}' \
    ${PEAKDIR}/WT_filtered_peaks.narrowPeak \
    > ${MOTIFDIR}/WT_summits_200bp.bed

# Run TagPileup for each condition using WT summit coordinates
for strain in WT SpxA2KO SpxA2OE; do
    # Use merged BAM naming convention
    if [ "$strain" == "WT" ]; then
        BAM="${BAMDIR}/WT_merged.bam"
    else
        BAM="${BAMDIR}/${strain}_merged.bam"
    fi

    java -Xmx50g -jar ${SCRIPTMANAGER} read-analysis tag-pileup \
        --cpu=54 \
        ${MOTIFDIR}/WT_summits_200bp.bed \
        ${BAM} \
        -o ${OUTPUTDIR}/${strain}_composite.out \
        -M ${OUTPUTDIR}/${strain}_matrix &
done
wait

echo "TagPileup complete"
echo "Load composite.out files in R for visualization (see chipexo_analysis.R)"

# =============================================================================
# SECTION 6: MOTIF ANALYSIS (MEME/FIMO)
# =============================================================================

mkdir -p ${MOTIFDIR}

# --- 6A: Extract sequences for MEME/FIMO ---
# NOTE: Run AFTER DiffBind analysis (Section 7) to get KO/WT gained peak BED files

# Extract sequences for differentially bound peak classes
# (KO_gained_peaks.bed and WT_gained_peaks.bed generated by DiffBind R script)
for class in KO_gained WT_gained; do
    bedtools getfasta \
        -fi ${BASEDIR}/NZ_CP067090.1.fasta \
        -bed ${MOTIFDIR}/${class}_peaks.bed \
        -nameOnly \
        -fo ${MOTIFDIR}/${class}_peaks_named.fasta

    # Add condition prefix to sequence names
    prefix=$(echo $class | cut -d'_' -f1)
    sed "s/>peak_/>${prefix}_peak_/" \
        ${MOTIFDIR}/${class}_peaks_named.fasta > \
        ${MOTIFDIR}/${class}_peaks_final.fasta
done

# --- 6B: De novo motif discovery (MEME) ---
# KO gained peaks - expect dimer motif (~20-25bp AT-rich)
meme ${MOTIFDIR}/KO_gained_peaks_final.fasta \
    -dna \
    -oc ${MOTIFDIR}/meme_KO_gained \
    -nmotifs 3 \
    -minw 10 \
    -maxw 25 \
    -mod zoops \
    -revcomp \
    -p 8 \
    2>&1 | tee ${MOTIFDIR}/meme_KO_gained.log &

# WT gained peaks - expect monomer motif (~6bp ATTARA)
meme ${MOTIFDIR}/WT_gained_peaks_final.fasta \
    -dna \
    -oc ${MOTIFDIR}/meme_WT_gained \
    -nmotifs 3 \
    -minw 4 \
    -maxw 15 \
    -mod zoops \
    -revcomp \
    -p 8 \
    2>&1 | tee ${MOTIFDIR}/meme_WT_gained.log &
wait

# --- 6C: Known motif scanning (FIMO) ---
# Create known CovR motif reference file (from Horstmann et al. 2023, emm3)
cat > ${MOTIFDIR}/CovR_known_motifs.txt << 'MOTIFEOF'
MEME version 4

ALPHABET= ACGT

strands: + -

Background letter frequencies
A 0.329 C 0.171 G 0.171 T 0.329

MOTIF CovR_dimer
letter-probability matrix: alength= 4 w= 19
 0.25  0.00  0.00  0.75
 0.50  0.00  0.00  0.50
 0.75  0.00  0.00  0.25
 0.50  0.00  0.00  0.50
 0.50  0.00  0.00  0.50
 1.00  0.00  0.00  0.00
 0.50  0.00  0.00  0.50
 1.00  0.00  0.00  0.00
 1.00  0.00  0.00  0.00
 0.75  0.00  0.00  0.25
 1.00  0.00  0.00  0.00
 0.75  0.25  0.00  0.00
 0.75  0.00  0.00  0.25
 0.50  0.00  0.00  0.50
 0.75  0.00  0.00  0.25
 0.25  0.00  0.25  0.50
 1.00  0.00  0.00  0.00
 0.50  0.00  0.00  0.50
 1.00  0.00  0.00  0.00
MOTIFEOF

# Monomer motif file (separate for permissive scanning)
cat > ${MOTIFDIR}/CovR_monomer_only.txt << 'MOTIFEOF'
MEME version 4

ALPHABET= ACGT

strands: + -

Background letter frequencies
A 0.329 C 0.171 G 0.171 T 0.329

MOTIF CovR_monomer
letter-probability matrix: alength= 4 w= 6
 0.95  0.02  0.02  0.02
 0.50  0.00  0.00  0.50
 0.50  0.00  0.00  0.50
 0.95  0.02  0.02  0.02
 0.50  0.00  0.50  0.00
 0.95  0.02  0.02  0.02
MOTIFEOF

# Scan KO gained peaks for dimer motif (p < 0.001)
fimo \
    --oc ${MOTIFDIR}/fimo_KO_final \
    --thresh 0.001 \
    ${MOTIFDIR}/CovR_known_motifs.txt \
    ${MOTIFDIR}/KO_gained_peaks_final.fasta

# Scan WT gained peaks for dimer motif
fimo \
    --oc ${MOTIFDIR}/fimo_WT_final \
    --thresh 0.001 \
    ${MOTIFDIR}/CovR_known_motifs.txt \
    ${MOTIFDIR}/WT_gained_peaks_final.fasta

# Scan both for monomer motif (p < 0.001)
fimo \
    --oc ${MOTIFDIR}/fimo_KO_monomer \
    --thresh 0.001 \
    ${MOTIFDIR}/CovR_monomer_only.txt \
    ${MOTIFDIR}/KO_gained_peaks_final.fasta

fimo \
    --oc ${MOTIFDIR}/fimo_WT_monomer \
    --thresh 0.001 \
    ${MOTIFDIR}/CovR_monomer_only.txt \
    ${MOTIFDIR}/WT_gained_peaks_final.fasta

echo "Motif analysis complete"

# =============================================================================
# SECTION 7: HORSTMANN 2023 VALIDATION
# =============================================================================

mkdir -p ${BASEDIR}/validation/horstmann2023

# Create BED file of Horstmann 2023 priority CovR sites
# Note: Locus tags updated for NZ_CP067090.1 reannotation
cat > ${BASEDIR}/validation/horstmann2023/horstmann_sites.bed << 'EOF'
NZ_CP067090.1	231	1587	dnaA	0	+
NZ_CP067090.1	115763	117299	nra	0	-
NZ_CP067090.1	151386	152742	nga	0	+
NZ_CP067090.1	283636	284323	covR	0	+
NZ_CP067090.1	1705714	1709224	scpA/B	0	-
NZ_CP067090.1	573409	573571	sagA	0	+
NZ_CP067090.1	1362260	1363178	trxB	0	-
NZ_CP067090.1	1497347	1498130	codY	0	-
NZ_CP067090.1	1677322	1678645	ska	0	+
NZ_CP067090.1	1711522	1713133	mga	0	-
NZ_CP067090.1	1843888	1845148	hasA	0	+
NZ_CP067090.1	903685	905248	guaA	0	+
NZ_CP067090.1	1034409	1034796	grab	0	-
NZ_CP067090.1	965999	966773	cfa	0	-
NZ_CP067090.1	1529639	1532639	endoS	0	-
NZ_CP067090.1	1755232	1755793	ahpC	0	+
NZ_CP067090.1	334612	339559	scpC	0	+
NZ_CP067090.1	685851	686871	ideS	0	-
NZ_CP067090.1	1681671	1682517	scl1	0	-
NZ_CP067090.1	795619	797445	scl2	0	+
NZ_CP067090.1	282835	283369	dahA	0	+
NZ_CP067090.1	275221	276592	brnQ	0	-
NZ_CP067090.1	1709593	1710973	emm	0	-
NZ_CP067090.1	349411	349663	spyA	0	+
NZ_CP067090.1	153255	154971	slo	0	+
EOF

# Sort BED file
sort -k1,1 -k2,2n \
    ${BASEDIR}/validation/horstmann2023/horstmann_sites.bed > \
    ${BASEDIR}/validation/horstmann2023/horstmann_sites_sorted.bed

# Intersect with WT peaks using 500bp window
bedtools slop \
    -i ${BASEDIR}/validation/horstmann2023/horstmann_sites_sorted.bed \
    -g <(echo "NZ_CP067090.1	1863912") \
    -b 500 | \
bedtools intersect \
    -a - \
    -b ${PEAKDIR}/WT_filtered_peaks.narrowPeak \
    -wa -u | \
    awk '{print $4}' > \
    ${BASEDIR}/validation/horstmann2023/found_with_500bp_window.txt

echo "Horstmann 2023 overlap:"
echo "Found: $(wc -l < ${BASEDIR}/validation/horstmann2023/found_with_500bp_window.txt)/25 sites"

echo "Pipeline complete - proceed to chipexo_analysis.R for DiffBind and visualization"
