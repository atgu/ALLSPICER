#!/usr/bin/env python3
"""Pull gnomAD v4.1.1 exome annotations for pLoF variants in APOB and TET2, to
stratify the ALLSPICE Task-2 (indel-removal) pLoF signal by:
  - annotation confidence / LoF-escape : LOFTEE lof, lof_flags, lof_filter
  - consequence class                  : frameshift vs stop_gained vs splice
  - NMD escape / isoform               : exon rank/total, MANE, CDS/protein pos
  - technical artifact                 : filters, region_flags, AS_VQSLOD, pab_max
  - frequency sanity                   : adj AF/AN, grpmax AF

Runs in a Hail/GCS environment. Exports a flat TSV keyed by locus+alleles that
joins to the ALLSPICER pLoF variant table (locus 'chrN:pos', alleles ["ref","alt"]).
"""
import os
import sys

# This directory contains a sibling `pickle.py` (and `utils.py`); Python prepends the
# script's own directory to sys.path, so numpy's `import pickle` would resolve to that
# local file and crash. Drop the script dir from sys.path before importing hail/numpy.
_self_dir = os.path.dirname(os.path.abspath(__file__))
sys.path = [p for p in sys.path if os.path.abspath(p or os.getcwd()) != _self_dir]

import hail as hl
# Same worst-consequence-per-gene logic the project uses to build the pLoF set
# (utils/annotations.py -> create_gene_map_ht). Requires the `gnomad` package.
from gnomad.utils.vep import process_consequences

# --- config ---
GNOMAD_HT = "gs://gcp-public-data--gnomad/release/4.1.1/ht/exomes/gnomad.exomes.v4.1.1.sites.ht"
OUTPUT    = "gs://aou_wlu/allspice/apob_tet2_plof_gnomad_v4.1.1.tsv.bgz"  # <- set to a writable bucket
GENES     = hl.set(["APOB", "TET2"])

# GRCh38 gene loci (padded) -- interval-filter first so we scan ~kb, not the whole exome.
INTERVALS = ["chr2:21000000-21050000",      # APOB
             "chr4:105140000-105280000"]    # TET2

# Consequence sets from utils/annotations.py (used for the project annotation label).
PLOF_CSQ = hl.set(["transcript_ablation", "splice_acceptor_variant",
                   "splice_donor_variant", "stop_gained", "frameshift_variant"])
MISSENSE_CSQ = hl.set(["stop_lost", "start_lost", "transcript_amplification",
                       "inframe_insertion", "inframe_deletion", "missense_variant"])
SYNONYMOUS_CSQ = hl.set(["stop_retained_variant", "synonymous_variant"])

ht = hl.read_table(GNOMAD_HT)
ht = hl.filter_intervals(ht, [hl.parse_locus_interval(i, "GRCh38") for i in INTERVALS])

# Pick the worst consequence per gene on the CANONICAL transcript exactly as the project
# does (utils/annotations.py: filter transcript_consequences to ENSG, then
# process_consequences -> worst_csq_by_gene_canonical), so the pLoF set here is built
# identically to the one we join against.
ht = ht.annotate(vep=ht.vep.annotate(
    transcript_consequences=ht.vep.transcript_consequences.filter(
        lambda x: x.gene_id.startswith("ENSG"))))
ht = process_consequences(ht, has_polyphen=False)  # v4+ dropped PolyPhen from vep struct
ht = ht.explode(ht.vep.worst_csq_by_gene_canonical)
ht = ht.annotate(tc=ht.vep.worst_csq_by_gene_canonical)
ht = ht.filter(GENES.contains(ht.tc.gene_symbol) & ht.tc.gene_id.startswith("ENSG"))

# Project annotation label (annotation_case_builder): LOFTEE HC -> pLoF, LC -> LC,
# else by most_severe_consequence. The pLoF set we join to is exactly lof == 'HC'.
ht = ht.annotate(annotation=(hl.case(missing_false=True)
    .when(ht.tc.lof == "HC", "pLoF")
    .when(ht.tc.lof == "LC", "LC")
    .when(MISSENSE_CSQ.contains(ht.tc.most_severe_consequence), "missense")
    .when(SYNONYMOUS_CSQ.contains(ht.tc.most_severe_consequence), "synonymous")
    .or_missing()))
# Keep pLoF-relevant variants (HC/LC or a pLoF consequence); the downstream join
# restricts to the HC pLoF set.
ht = ht.filter((ht.tc.lof == "HC") | (ht.tc.lof == "LC") |
               ht.tc.consequence_terms.any(lambda c: PLOF_CSQ.contains(c)))

# exon is "rank/total" (e.g. "5/29"); parse for the last-exon NMD-escape heuristic.
exon_parts = hl.or_missing(hl.is_defined(ht.tc.exon) & ht.tc.exon.contains("/"),
                           ht.tc.exon.split("/"))
exon_rank  = hl.or_missing(hl.is_defined(exon_parts), hl.int(exon_parts[0]))
exon_total = hl.or_missing(hl.is_defined(exon_parts), hl.int(exon_parts[1]))

ref, alt = ht.alleles[0], ht.alleles[1]
out = ht.select(
    # --- join key / identity ---
    variant_id = ht.locus.contig + ":" + hl.str(ht.locus.position) + ":" + ref + ":" + alt,
    gene_symbol = ht.tc.gene_symbol,
    annotation = ht.annotation,                       # project label: pLoF (HC) / LC / ...
    # --- consequence class ---
    most_severe_consequence = ht.vep.most_severe_consequence,
    consequence_terms = hl.delimit(ht.tc.consequence_terms, "|"),
    variant_class = ht.vep.variant_class,
    # --- LoF confidence / escape (LOFTEE) ---
    lof = ht.tc.lof, lof_flags = ht.tc.lof_flags, lof_filter = ht.tc.lof_filter,
    # --- transcript / isoform / NMD-escape ---
    transcript_id = ht.tc.transcript_id,
    mane_select = ht.tc.mane_select,
    canonical = ht.tc.canonical == 1,
    exon = ht.tc.exon, intron = ht.tc.intron,
    exon_rank = exon_rank, exon_total = exon_total,
    in_last_exon = hl.or_missing(hl.is_defined(exon_rank), exon_rank == exon_total),
    cds_start = ht.tc.cds_start, cds_end = ht.tc.cds_end,
    protein_start = ht.tc.protein_start, protein_end = ht.tc.protein_end,
    hgvsc = ht.tc.hgvsc, hgvsp = ht.tc.hgvsp,
    # --- variant type ---
    allele_type = ht.allele_info.allele_type,         # snv / ins / del / mixed
    was_mixed = ht.allele_info.was_mixed,
    is_snv = (hl.len(ref) == 1) & (hl.len(alt) == 1),
    indel_len = hl.abs(hl.len(alt) - hl.len(ref)),
    is_frameshift = ht.tc.consequence_terms.contains("frameshift_variant"),
    # --- technical QC ---
    filters = hl.delimit(hl.array(ht.filters), "|"),
    pass_qc = hl.len(ht.filters) == 0,
    AS_VQSLOD = ht.info.AS_VQSLOD,
    AS_pab_max = ht.info.AS_pab_max,
    monoallelic = ht.info.monoallelic,
    lcr = ht.region_flags.lcr,
    segdup = ht.region_flags.segdup,
    fail_interval_qc = ht.region_flags.fail_interval_qc,
    # --- frequency sanity (freq[0] = overall adj) ---
    AC_adj = ht.freq[0].AC, AF_adj = ht.freq[0].AF, AN_adj = ht.freq[0].AN,
    nhomalt_adj = ht.freq[0].homozygote_count,
    grpmax_AF = ht.grpmax.gnomad.AF, grpmax_anc = ht.grpmax.gnomad.gen_anc,
    # --- splice support (for splice pLoF) ---
    spliceai_ds_max = ht.in_silico_predictors.spliceai_ds_max,
    pangolin_largest_ds = ht.in_silico_predictors.pangolin_largest_ds,
    cadd_phred = ht.in_silico_predictors.cadd.phred,
)

print(f"APOB+TET2 pLoF variants: {out.count()}")
out.group_by(out.gene_symbol, out.lof, out.allele_type).aggregate(n=hl.agg.count()).show(50)
out.export(OUTPUT)
print(f"wrote {OUTPUT}")
