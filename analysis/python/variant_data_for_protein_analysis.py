from ukb_exomes import *
import hail as hl
from gnomad.utils.vep import process_consequences
from gnomad.resources.grch38.reference_data import clinvar
from gnomad.resources.resource_utils import DataException
from gnomad.utils.vep import vep_or_lookup_vep

my_bucket = 'gs://ukbb-exome-pharma-wlu'

def subset_gene_variant_pheno_ht(
        gene_symbol: str,
        phenos: set,
        tranche=CURRENT_TRANCHE
):
    print(tranche)
    var = hl.read_matrix_table(get_results_mt_path('variant', tranche=tranche))
    sub = var.filter_rows((var.gene == gene_symbol) & (hl.is_defined(var.annotation)))
    sub = sub.filter_cols(hl.literal(phenos).contains(sub.phenocode))

    vep = hl.read_table(var_annotations_ht_path('vep', *TRANCHE_DATA[tranche]))
    vep = process_consequences(vep)
    vep = vep.explode(vep.vep.worst_csq_by_gene_canonical)
    vep = vep.explode(vep.vep.transcript_consequences)
    annotation = annotation_case_builder(vep.vep.worst_csq_by_gene_canonical)

    vep = vep.annotate(annotation=hl.if_else(hl.literal({'missense', 'LC'}).contains(annotation), 'missense|LC', annotation))
    vep = vep.filter(hl.is_defined(vep.annotation) & (vep.annotation != 'non-coding'))

    vep = vep.annotate(hgvsp = vep.vep.transcript_consequences.hgvsp,
                       gene_id = vep.vep.worst_csq_by_gene_canonical.gene_id,
                       gene_symbol = vep.vep.worst_csq_by_gene_canonical.gene_symbol,
                       annotation = vep.annotation)
    vep = vep.filter((vep.gene_symbol == gene_symbol))

    sub = sub.drop('AF', 'AC')
    sub = sub.annotate_rows(
        hgvsp=vep[sub.row_key]['hgvsp'],
        AC=sub.call_stats.AC,
        AF=sub.call_stats.AF)
    sub = sub.entries()
    sub = sub.select('AC', 'AF', 'gene', 'annotation', 'BETA', 'Pvalue', 'hgvsp')

    return sub

gene_list = hl.import_table(f'{my_bucket}/wrap_up_results/ukb_pleiotropy_sig_genes_1e_5_for_siwei.tsv', delimiter='\t')
gene_list = gene_list.key_by('gene')
phenos = gene_list.value.collect()
genes = gene_list.gene.collect()

info = {}
genes = ['ALPL', 'ALB']
for i in range(len(genes)):
    if genes[i] in info.keys():
        if type(info[genes[i]]) == list:
            info[genes[i]].append(phenos[i])
        else:
            info[genes[i]] = [info[genes[i]], phenos[i]]
    else:
        info[genes[i]]=phenos[i]

for i in range(len(info)):
    print(list(info.keys())[i])
    print(list(info.values())[i])
    if not hl.hadoop_exists(f'gs://ukbb-exome-pharma-wlu/wrap_up_results/for_siwei/{list(info.keys())[i]}_var_500k_for_siwei.txt'):
        ht = subset_gene_variant_pheno_ht(list(info.keys())[i], list(info.values())[i])
        ht = ht.filter((ht.AF <= 1e-4) & (ht.annotation == 'missense'))
        ht = ht.annotate(
            chromosome = ht.locus.contig.replace('chr', ''),
            start = ht.locus.position,
            end = ht.locus.position + ht.alleles[1].length() - ht.alleles[0].length(),
            allele = hl.if_else(ht.alleles[0].length() == ht.alleles[1].length(), ht.alleles[0] + '/' + ht.alleles[1],
                                hl.if_else(ht.alleles[0].length() > ht.alleles[1].length(), ht.alleles[0][1:] + '/-', '-/' + ht.alleles[1][1:])
                                )
        )
        ht.export(f'gs://ukbb-exome-pharma-wlu/wrap_up_results/for_siwei/{list(info.keys())[i]}_var_500k_for_siwei.txt')
    # break