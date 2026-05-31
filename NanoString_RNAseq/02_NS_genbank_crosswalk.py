#!/usr/bin/env python3
"""
=============================================================================
NanoString Probe → MGAS10870 Locus Tag Crosswalk
=============================================================================
SpxA1 and SpxA2 Function as a Stoichiometry-Dependent Regulatory Rheostat
Governing Virulence Gene Expression in Group A Streptococcus

Flores Streptococcal Laboratory
Division of Pediatric Infectious Diseases
Vanderbilt University Medical Center

Authors:  Anthony R. Flores, Matthew A. Sanson, Luis A. Vega [et al.]
Contact:  flores-lab.org
Journal:  mBio (in submission)
GEO:      PRJNA1472884

Description:
    Builds a mapping from NanoString nCounter probe IDs (M5005_Spy_ / SpyM3_
    format) to MGAS10870 C9Q_ locus identifiers for cross-platform validation
    of NanoString and RNA-seq data.

    Four-tier matching strategy (in order):
      1. SpyM3_ probes  — direct SpyM3 tag lookup via MGAS10870 locus key
      2. Exact protein  — identical protein sequence (MGAS5005 vs MGAS10870)
      3. Gene name      — gene name match in MGAS10870 GenBank annotation
      4. k-mer          — Jaccard similarity of 10-mer sets (threshold >= 0.5)

    Complement-strand probes (M5005_Spy_XXXX c suffix) resolved by stripping
    trailing 'c' before MGAS5005 lookup.

    Five probes assigned manually after lab verification:
      sclA, sclB, fasC, emm3.0, ralp4

Usage:
    python3 02_NS_genbank_crosswalk.py \
        --mgas5005   /path/to/M1_MGAS5005.gb \
        --mgas10870  /path/to/M3_CP067090.gb \
        --locus_key  /path/to/MGAS10870_locus-tag_key.xlsx \
        --output     /path/to/NS_crosswalk_final.csv

Requirements:
    Python >= 3.11
    openpyxl (pip install openpyxl)

Input files (available from NCBI):
    M1_MGAS5005.gb     — MGAS5005 GenBank (CP000017)
    M3_CP067090.gb     — MGAS10870 GenBank (NZ_CP067090.1)
    MGAS10870_locus-tag_key.xlsx — C9Q_ to SpyM3 locus tag mapping

Output:
    NS_crosswalk_final.csv — probe_id, C9Q_tag, SpyM3_tag, match_type,
                              kmer_identity columns
=============================================================================
"""


import re
import csv
import argparse
import sys
from pathlib import Path

try:
    import openpyxl
except ImportError:
    sys.exit("ERROR: openpyxl not found. Run: pip install openpyxl")


# ---------------------------------------------------------------------------
# GENBANK PARSER
# ---------------------------------------------------------------------------

def parse_genbank_proteins(filepath):
    """
    Parse CDS features from a GenBank flat file.

    Returns a list of dicts with keys:
        locus_tag, new_locus_tag, gene, product, protein_seq
    For MGAS5005, locus_tag is set to old_locus_tag (M5005_Spy_XXXX format).
    For MGAS10870, locus_tag is the current C9Q_ tag.
    """
    records = []
    with open(filepath, encoding="utf-8", errors="replace") as f:
        content = f.read()

    cds_pattern = re.compile(
        r'     CDS\s+(?:complement\()?\d+\.\.\d+\)?'
        r'\n(.*?)(?=\n     \w|\nORIGIN)',
        re.DOTALL
    )

    for match in cds_pattern.finditer(content):
        body = match.group(1)

        old_lt  = re.search(r'/old_locus_tag="([^"]+)"', body)
        new_lt  = re.search(r'/locus_tag="([^"]+)"',     body)
        gene    = re.search(r'/gene="([^"]+)"',           body)
        product = re.search(r'/product="([^"]+)"',        body)
        transl  = re.search(r'/translation="([^"]+)"',    body, re.DOTALL)

        locus     = old_lt.group(1) if old_lt else (
                    new_lt.group(1) if new_lt else None)
        new_locus = new_lt.group(1) if new_lt else None

        seq = re.sub(r'\s+', '', transl.group(1)) if transl else None

        records.append({
            'locus_tag'     : locus,
            'new_locus_tag' : new_locus,
            'gene'          : gene.group(1)    if gene    else None,
            'product'       : product.group(1) if product else None,
            'protein_seq'   : seq,
        })

    return records


# ---------------------------------------------------------------------------
# K-MER SIMILARITY
# ---------------------------------------------------------------------------

def kmer_identity(seq1, seq2, k=10):
    """
    Jaccard similarity of k-mer sets between two protein sequences.
    Returns float in [0, 1]; 0.0 if either sequence is empty or too short.
    """
    if not seq1 or not seq2 or min(len(seq1), len(seq2)) < k:
        return 0.0
    kmers1 = set(seq1[i:i+k] for i in range(len(seq1) - k + 1))
    kmers2 = set(seq2[i:i+k] for i in range(len(seq2) - k + 1))
    union  = kmers1 | kmers2
    return len(kmers1 & kmers2) / len(union) if union else 0.0


# ---------------------------------------------------------------------------
# NANOSTRING PROBE LIST
# ---------------------------------------------------------------------------

# All 119 probes in the SpxA1/SpxA2 codeset (target genes + endogenous controls)
NS_PROBES = [
    "M5005_Spy_0013", "M5005_Spy_0017", "M5005_Spy_0018",
    "M5005_Spy_0034", "M5005_Spy_0036", "M5005_Spy_0040",
    "M5005_Spy_0077", "M5005_Spy_0106c","M5005_Spy_0107",
    "M5005_Spy_0109", "M5005_Spy_0114", "M5005_Spy_0115c",
    "M5005_Spy_0124", "M5005_Spy_0139", "M5005_Spy_0142",
    "M5005_Spy_0149", "M5005_Spy_0161", "M5005_Spy_0182",
    "M5005_Spy_0186c","M5005_Spy_0195", "M5005_Spy_0205",
    "M5005_Spy_0274c","M5005_Spy_0282", "M5005_Spy_0317",
    "M5005_Spy_0340", "M5005_Spy_0341", "M5005_Spy_0351c",
    "M5005_Spy_0356c","M5005_Spy_0367c","M5005_Spy_0368",
    "M5005_Spy_0407", "M5005_Spy_0413", "M5005_Spy_0424",
    "M5005_Spy_0435", "M5005_Spy_0473c","M5005_Spy_0474",
    "M5005_Spy_0476", "M5005_Spy_0478c","M5005_Spy_0484",
    "M5005_Spy_0508", "M5005_Spy_0517c","M5005_Spy_0543",
    "M5005_Spy_0556", "M5005_Spy_0559c","M5005_Spy_0561",
    "M5005_Spy_0562", "M5005_Spy_0563", "M5005_Spy_0571",
    "M5005_Spy_0598c","M5005_Spy_0651", "M5005_Spy_0668c",
    "M5005_Spy_0680", "M5005_Spy_0691", "M5005_Spy_0701c",
    "M5005_Spy_0713", "M5005_Spy_0777", "M5005_Spy_0805",
    "M5005_Spy_0831c","M5005_Spy_0874", "M5005_Spy_0948c",
    "M5005_Spy_0959c","M5005_Spy_0981c","M5005_Spy_0996",
    "M5005_Spy_1073c","M5005_Spy_1106c","M5005_Spy_1145c",
    "M5005_Spy_1169", "M5005_Spy_1237c","M5005_Spy_1275c",
    "M5005_Spy_1291c","M5005_Spy_1307c","M5005_Spy_1334c",
    "M5005_Spy_1395c","M5005_Spy_1402", "M5005_Spy_1407",
    "M5005_Spy_1415c","M5005_Spy_1474c","M5005_Spy_1479",
    "M5005_Spy_1512c","M5005_Spy_1528c","M5005_Spy_1531c",
    "M5005_Spy_1540c","M5005_Spy_1557", "M5005_Spy_1571c",
    "M5005_Spy_1576", "M5005_Spy_1611c","M5005_Spy_1612c",
    "M5005_Spy_1625c","M5005_Spy_1636c","M5005_Spy_1680c",
    "M5005_Spy_1684", "M5005_Spy_1687c","M5005_Spy_1702",
    "M5005_Spy_1711c","M5005_Spy_1715c","M5005_Spy_1718c",
    "M5005_Spy_1719c","M5005_Spy_1720c","M5005_Spy_1723c",
    "M5005_Spy_1725c","M5005_Spy_1732c","M5005_Spy_1735c",
    "M5005_Spy_1737", "M5005_Spy_1738c","M5005_Spy_1751c",
    "M5005_Spy_1782c","M5005_Spy_1798c","M5005_Spy_1850c",
    "M5005_Spy_1851", "M5005_Spy_1857c","M5005_Spy_1865",
    # emm3-specific SpyM3 probes
    "SpyM3_0097", "SpyM3_0098", "SpyM3_0100", "SpyM3_0104",
    "SpyM3_0131", "SpyM3_0307", "SpyM3_0738", "SpyM3_1409",
]

# Manual overrides verified by laboratory curator
# Applied after automated matching to correct diverged/mismatched probes.
MANUAL_OVERRIDES = {
    # Probe           C9Q_tag      Gene name   SpyM3 tag    Match type
    "M5005_Spy_1687c": {
        'C9Q_tag'       : 'C9Q_08620',
        'C9Q_gene_name' : 'sclA',
        'SpyM3_tag'     : '',
        'Match_type'    : 'manual_verified',
        'Notes'         : ('sclA (scl2) — collagen-like repeat protein; '
                           'diverged repeat region; verified C9Q_08620'),
    },
    "M5005_Spy_0777": {
        'C9Q_tag'       : 'C9Q_04075',
        'C9Q_gene_name' : 'sclB',
        'SpyM3_tag'     : 'SpyM3_0737',
        'Match_type'    : 'manual_verified',
        'Notes'         : ('sclB (scl1) — collagen-like repeat protein; '
                           'verified C9Q_04075 / SpyM3_0737'),
    },
    "M5005_Spy_0205": {
        'C9Q_tag'       : 'C9Q_01125',
        'C9Q_gene_name' : 'fasC',
        'SpyM3_tag'     : '',
        'Match_type'    : 'manual_verified',
        'Notes'         : 'fasC — verified C9Q_01125',
    },
    "M5005_Spy_1719c": {
        'C9Q_tag'       : 'C9Q_08745',
        'C9Q_gene_name' : 'emm3.0',
        'SpyM3_tag'     : '',
        'Match_type'    : 'manual_emm_crosshybrid',
        'Notes'         : ('emm1.0 probe cross-hybridizes with emm3.0; '
                           'sequences too diverged for automated match; '
                           'verified C9Q_08745'),
    },
    "M5005_Spy_0186c": {
        'C9Q_tag'       : '',
        'C9Q_gene_name' : '',
        'SpyM3_tag'     : '',
        'Match_type'    : 'absent_M10870_verified',
        'Notes'         : ('ralp4 — confirmed absent from MGAS10870; '
                           'excluded from RNA-seq validation'),
    },
}

# Probes known to target M1-specific genes absent from emm3
M1_SPECIFIC = {
    "M5005_Spy_0106c",  # rofA — M1 FCT region
    "M5005_Spy_0107",   # cpa  — M1 FCT region
    "M5005_Spy_0109",   # prtF — M1 FCT region
    "M5005_Spy_0114",   # srtB — M1 FCT region
    "M5005_Spy_0356c",  # speJ — M1 phage exotoxin
    "M5005_Spy_0561",   # epf  — M1 FCT region variant
    "M5005_Spy_0805",   # srtK — absent in emm3
    "M5005_Spy_1169",   # spd3 — phage DNase
    "M5005_Spy_1415c",  # sdaD2 — phage streptodornase
    "M5005_Spy_1702",   # smeZ — M1 phage exotoxin
    "M5005_Spy_1718c",  # sic1 — M1-specific complement inhibitor
}


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Build NanoString → MGAS10870 locus tag crosswalk")
    parser.add_argument("--mgas5005",  required=True,
                        help="MGAS5005 GenBank file (M1_MGAS5005.gb)")
    parser.add_argument("--mgas10870", required=True,
                        help="MGAS10870 GenBank file (M3_CP067090.gb)")
    parser.add_argument("--locus_key", required=True,
                        help="MGAS10870_locus-tag_key.xlsx")
    parser.add_argument("--output",    required=True,
                        help="Output CSV path (NS_crosswalk_final.csv)")
    args = parser.parse_args()

    # ── Parse GenBank files ────────────────────────────────────────────────
    print(f"Parsing MGAS5005: {args.mgas5005}")
    m5005  = parse_genbank_proteins(args.mgas5005)
    print(f"  {len(m5005)} CDS features parsed")

    print(f"Parsing MGAS10870: {args.mgas10870}")
    m10870 = parse_genbank_proteins(args.mgas10870)
    print(f"  {len(m10870)} CDS features parsed")

    # ── Build MGAS5005 lookup ──────────────────────────────────────────────
    # Index by both old_locus_tag (M5005_Spy_XXXX) and new_locus_tag.
    # MGAS5005 GenBank uses old_locus_tag WITHOUT the 'c' suffix;
    # our NanoString annotation appends 'c' to complement-strand genes.
    # Strip trailing 'c' when looking up.
    m5005_by_locus = {}
    for r in m5005:
        for lt in (r['locus_tag'], r['new_locus_tag']):
            if lt:
                m5005_by_locus[lt] = r

    def lookup_m5005(probe_id):
        if probe_id in m5005_by_locus:
            return m5005_by_locus[probe_id]
        stripped = re.sub(r'c$', '', probe_id)
        return m5005_by_locus.get(stripped)

    # ── Build MGAS10870 lookups ────────────────────────────────────────────
    m10870_by_seq   = {}
    m10870_by_locus = {}
    m10870_by_gene  = {}
    for r in m10870:
        if r['protein_seq']:
            m10870_by_seq.setdefault(r['protein_seq'], []).append(r)
        if r['locus_tag']:
            m10870_by_locus[r['locus_tag']] = r
        if r['gene']:
            m10870_by_gene.setdefault(r['gene'], []).append(r)

    # ── Load locus key ─────────────────────────────────────────────────────
    print(f"Loading locus key: {args.locus_key}")
    wb = openpyxl.load_workbook(args.locus_key, read_only=True)
    ws = wb.active
    spym3_map   = {}   # C9Q_ → SpyM3 tag
    c9q_to_name = {}   # C9Q_ → gene name
    name_to_c9q = {}   # gene name → [C9Q_ tags]

    for row in ws.iter_rows(min_row=2, values_only=True):
        name, c9q, spym3 = str(row[0] or ''), str(row[1] or ''), str(row[2] or '')
        if c9q and spym3 and spym3 != '#N/A':
            spym3_map[c9q]  = spym3
        if c9q and name:
            c9q_to_name[c9q] = name
            name_to_c9q.setdefault(name, []).append(c9q)

    spym3_to_c9q = {v: k for k, v in spym3_map.items()}
    print(f"  {len(spym3_map)} SpyM3 mappings loaded")

    # ── Build crosswalk ────────────────────────────────────────────────────
    crosswalk = []

    for probe in NS_PROBES:
        row = {
            'NanoString_ProbeID' : probe,
            'M5005_gene'         : '',
            'M5005_product'      : '',
            'C9Q_tag'            : '',
            'C9Q_gene_name'      : '',
            'SpyM3_tag'          : '',
            'Match_type'         : '',
            'Kmer_identity'      : '',
            'Notes'              : '',
        }

        # ── Tier 0: SpyM3_ probes — direct lookup ─────────────────────────
        if probe.startswith('SpyM3_'):
            row['SpyM3_tag']  = probe
            row['Match_type'] = 'SpyM3_direct'
            c9q = spym3_to_c9q.get(probe, '')
            row['C9Q_tag']       = c9q
            row['C9Q_gene_name'] = c9q_to_name.get(c9q, '')
            crosswalk.append(row)
            continue

        # ── Tier 1: Exact protein sequence match ──────────────────────────
        m5_rec = lookup_m5005(probe)
        if not m5_rec:
            row['Match_type'] = 'not_in_MGAS5005'
            crosswalk.append(row)
            continue

        row['M5005_gene']    = m5_rec.get('gene')    or ''
        row['M5005_product'] = m5_rec.get('product') or ''
        seq5 = m5_rec.get('protein_seq')

        if not seq5:
            row['Match_type'] = 'no_protein_seq'
            crosswalk.append(row)
            continue

        if seq5 in m10870_by_seq:
            matches = m10870_by_seq[seq5]
            if len(matches) == 1:
                m3 = matches[0]
                row['C9Q_tag']       = m3['locus_tag']
                row['C9Q_gene_name'] = m3.get('gene') or \
                                        c9q_to_name.get(m3['locus_tag'], '')
                row['SpyM3_tag']     = spym3_map.get(m3['locus_tag'], '')
                row['Match_type']    = 'exact_protein'
                row['Kmer_identity'] = '1.000'
            else:
                # Multiple exact matches — try gene name to resolve
                gene_match = [m for m in matches
                              if m.get('gene') == row['M5005_gene']]
                if len(gene_match) == 1:
                    m3 = gene_match[0]
                    row['C9Q_tag']       = m3['locus_tag']
                    row['C9Q_gene_name'] = m3.get('gene') or ''
                    row['SpyM3_tag']     = spym3_map.get(m3['locus_tag'], '')
                    row['Match_type']    = 'exact_protein_gene_resolved'
                    row['Kmer_identity'] = '1.000'
                else:
                    row['C9Q_tag']   = '|'.join(m['locus_tag']
                                                  for m in matches)
                    row['SpyM3_tag'] = '|'.join(
                        spym3_map.get(m['locus_tag'], '')
                        for m in matches if spym3_map.get(m['locus_tag']))
                    row['Match_type'] = 'paralog_exact'
                    row['Notes']      = f"{len(matches)} identical copies"
            crosswalk.append(row)
            continue

        # ── Tier 2: Gene name match in MGAS10870 ──────────────────────────
        gene5 = row['M5005_gene']
        if gene5 and gene5 in m10870_by_gene:
            gmatches = m10870_by_gene[gene5]
            if len(gmatches) == 1:
                m3  = gmatches[0]
                sim = kmer_identity(seq5, m3['protein_seq'] or '')
                row['C9Q_tag']       = m3['locus_tag']
                row['C9Q_gene_name'] = m3.get('gene') or ''
                row['SpyM3_tag']     = spym3_map.get(m3['locus_tag'], '')
                row['Match_type']    = 'gene_name_match'
                row['Kmer_identity'] = f"{sim:.3f}"
                crosswalk.append(row)
                continue

        # ── Tier 3: Locus key name lookup ─────────────────────────────────
        if gene5 and gene5 in name_to_c9q:
            c9qs = name_to_c9q[gene5]
            if len(c9qs) == 1:
                c9q = c9qs[0]
                row['C9Q_tag']       = c9q
                row['C9Q_gene_name'] = c9q_to_name.get(c9q, '')
                row['SpyM3_tag']     = spym3_map.get(c9q, '')
                row['Match_type']    = 'locus_key_name_match'
                if c9q in m10870_by_locus:
                    m3  = m10870_by_locus[c9q]
                    sim = kmer_identity(seq5, m3.get('protein_seq') or '')
                    row['Kmer_identity'] = f"{sim:.3f}"
                crosswalk.append(row)
                continue

        # ── Tier 4: k-mer similarity scan ─────────────────────────────────
        best_sim  = 0.0
        best_recs = []
        for m3_seq, m3_recs in m10870_by_seq.items():
            if not m3_seq:
                continue
            sim = kmer_identity(seq5, m3_seq)
            if sim > best_sim:
                best_sim  = sim
                best_recs = m3_recs
            elif sim == best_sim and sim > 0:
                best_recs.extend(m3_recs)

        if best_sim >= 0.5 and len(best_recs) == 1:
            m3 = best_recs[0]
            row['C9Q_tag']       = m3['locus_tag']
            row['C9Q_gene_name'] = m3.get('gene') or \
                                    c9q_to_name.get(m3['locus_tag'], '')
            row['SpyM3_tag']     = spym3_map.get(m3['locus_tag'], '')
            row['Match_type']    = 'kmer_similarity'
            row['Kmer_identity'] = f"{best_sim:.3f}"
        elif best_sim >= 0.5 and len(best_recs) > 1:
            row['C9Q_tag']       = '|'.join(m['locus_tag']
                                              for m in best_recs)
            row['SpyM3_tag']     = '|'.join(
                spym3_map.get(m['locus_tag'], '')
                for m in best_recs if spym3_map.get(m['locus_tag']))
            row['Match_type']    = 'kmer_paralog'
            row['Kmer_identity'] = f"{best_sim:.3f}"
            row['Notes']         = f"{len(best_recs)} tied best matches"
        else:
            m1_note = ' (M1-specific, expected absent in emm3)' \
                      if probe in M1_SPECIFIC else ''
            row['Match_type'] = 'no_match'
            row['Notes']      = (f"Best kmer identity: {best_sim:.3f}"
                                  f"{m1_note}")

        crosswalk.append(row)

    # ── Apply manual overrides ─────────────────────────────────────────────
    for r in crosswalk:
        pid = r['NanoString_ProbeID']
        if pid in MANUAL_OVERRIDES:
            r.update(MANUAL_OVERRIDES[pid])

    # ── Compute Has_C9Q flag ───────────────────────────────────────────────
    for r in crosswalk:
        r['Has_C9Q'] = (
            bool(r['C9Q_tag']) and
            '|' not in r['C9Q_tag'] and
            r['Match_type'] not in ('absent_M10870_verified',
                                     'not_in_MGAS5005', 'no_match')
        )

    # ── Summary ───────────────────────────────────────────────────────────
    by_type = {}
    for r in crosswalk:
        by_type[r['Match_type']] = by_type.get(r['Match_type'], 0) + 1

    print("\n=== CROSSWALK SUMMARY ===")
    for k, v in sorted(by_type.items()):
        print(f"  {k:<40s}: {v}")
    print(f"\nTotal probes:        {len(crosswalk)}")
    print(f"Has C9Q_ assigned:   {sum(r['Has_C9Q'] for r in crosswalk)}")
    print(f"Excluded from RNA-seq validation: "
          f"{sum(not r['Has_C9Q'] for r in crosswalk)}")

    # ── Write CSV ─────────────────────────────────────────────────────────
    fieldnames = ['NanoString_ProbeID', 'M5005_gene', 'M5005_product',
                  'C9Q_tag', 'C9Q_gene_name', 'SpyM3_tag',
                  'Match_type', 'Kmer_identity', 'Notes', 'Has_C9Q']

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames,
                                 extrasaction='ignore')
        writer.writeheader()
        writer.writerows(crosswalk)

    print(f"\nCrosswalk saved: {output_path}")


if __name__ == '__main__':
    main()
