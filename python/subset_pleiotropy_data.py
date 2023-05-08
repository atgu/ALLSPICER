import hail as hl
import logging
import argparse
from ukb_exomes import *
from ukbb_common import *

CURRENT_TRANCHE = "300k"
logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("subset_pleiotropy_data")
logger.setLevel(logging.INFO)

def main(args):
    print(CURRENT_TRANCHE)
    logger.info(f"Current tranche: {CURRENT_TRANCHE} data... ")
    if args.random_phenos:
        gene = hl.read_matrix_table(get_results_mt_path("gene", random_phenos=True))
        var = hl.read_matrix_table(get_results_mt_path("variant", random_phenos=True))

        gene = gene.annotate_rows(n_var=hl.agg.collect_as_set(gene.total_variants))
        gene = gene.explode_rows(gene.n_var)
        gene = gene.drop("total_variants")
        gene = gene.annotate_rows(total_variants=gene.n_var)
        gene = gene.select(P_VALUE_FIELDS[args.test_type], "n_cases")
        var = var.select("AC", "AF", "gene", "annotation", "BETA", "Pvalue", "n_cases")

    if args.real_phenos:
        if args.filter_rows:
            gene = load_final_sumstats_table(result_type="gene", extension="mt", tranche=CURRENT_TRANCHE)
            gene.count()
            gene = gene.filter_rows(
                gene[f"keep_gene_{args.test_type}"]
                & gene['keep_gene_expected_ac' if CURRENT_TRANCHE=='500k' else 'keep_gene_caf']
                & gene.keep_gene_coverage
                & gene.keep_gene_n_var
            )
            var = load_final_sumstats_table(result_type="variant", extension="mt", tranche=CURRENT_TRANCHE)
            var.count()
            var = annotate_additional_info_mt(var, "variant", tranche=CURRENT_TRANCHE)
            var = var.filter_rows(var.keep_var_annt)
        else:
            gene = hl.read_matrix_table(get_results_mt_path("gene"))
            var = hl.read_matrix_table(get_results_mt_path("variant"))
        if args.filter_cols:
            gene = gene.filter_cols(gene[f"keep_pheno_{args.test_type}"])
            var = var.filter_cols(var[f"keep_pheno_{args.test_type}"])
        if args.phenocode_lst is not None:
            gene = gene.filter_cols(
                hl.literal(args.phenocode_lst).contains(gene.phenocode)
            )
            var = var.filter_cols(
                hl.literal(args.phenocode_lst).contains(var.phenocode)
            )
        if args.irnt:
            gene = gene.filter_cols(hl.literal({"irnt"}).contains(gene.modifier))
            var = var.filter_cols(hl.literal({"irnt"}).contains(var.modifier))
        if args.icd:
            gene = gene.annotate_cols(trait_type2 = hl.if_else(
                    gene.trait_type == "icd_first_occurrence", "icd10", gene.trait_type
                ))
            var = var.annotate_cols(trait_type2 = hl.if_else(
                    var.trait_type == "icd_first_occurrence", "icd10", var.trait_type
                ))
            gene = gene.filter_cols(gene.trait_type2 == "icd10")
            var = var.filter_cols(var.trait_type2 == "icd10")
        if args.continuous:
            gene = gene.filter_cols(gene.trait_type == "continuous")
            var = var.filter_cols(var.trait_type == "continuous")
        if args.gene_lst is not None:
            gene = gene.filter_rows(
                hl.literal(args.gene_lst).contains(gene.gene_symbol)
            )
            var = var.filter_rows(hl.literal(args.gene_lst).contains(var.gene))
        if CURRENT_TRANCHE == "500k":
            gene = gene.filter_rows(gene.annotation != "pLoF|missense|LC")
            var = var.drop("AC", "AF")
            var = var.annotate_rows(**var.call_stats)

        gene = gene.annotate_rows(
            **{
                f"sig_{args.test_type}": hl.agg.count_where(
                    gene[P_VALUE_FIELDS[args.test_type]]
                    < EMPIRICAL_P_THRESHOLDS[args.test_type]
                )
            }
        )
        gene = gene.filter_rows((gene[f"sig_{args.test_type}"] > 1))
        gene = gene.entries()
        var = var.entries()

        gene = gene.select(
            P_VALUE_FIELDS[args.test_type],
            f"{'n_cases_defined' if CURRENT_TRANCHE=='500k' else 'n_cases'}",
            "description",
        )
        var = var.select(
            "AC",
            "AF",
            "gene",
            "annotation",
            "BETA",
            "Pvalue",
            "n_cases",
            "description",
            "pathogenicity",
            "polyphen2",
        )

        var = var.filter(
            hl.literal(
                list(dict(gene.aggregate(hl.agg.counter(gene.gene_symbol))).keys())
            ).contains(var.gene)
            & (var.AF < args.af_upper)
        )
        var = var.annotate_globals(
            af_upper=args.af_upper
        )

    gene.write(
        f"gs://ukbb-exome-pharma-wlu/pleiotropy_oct2022/{CURRENT_TRANCHE}/corr_testing_{args.test_type}_gene_{CURRENT_TRANCHE}.ht",
        overwrite=args.overwrite,
    )
    var.write(
        f"gs://ukbb-exome-pharma-wlu/pleiotropy_oct2022/{CURRENT_TRANCHE}/corr_testing_{args.test_type}_var_{CURRENT_TRANCHE}_new_AF.ht",
        overwrite=args.overwrite,
    )
    gene = hl.read_table(
        f"gs://ukbb-exome-pharma-wlu/pleiotropy_oct2022/{CURRENT_TRANCHE}/corr_testing_{args.test_type}_gene_{CURRENT_TRANCHE}.ht"
    )
    gene.export(
        f"gs://ukbb-exome-pharma-wlu/pleiotropy_oct2022/{CURRENT_TRANCHE}/corr_testing_{args.test_type}_gene_{CURRENT_TRANCHE}.txt.bgz"
    )
    var = hl.read_table(
        f"gs://ukbb-exome-pharma-wlu/pleiotropy_oct2022/{CURRENT_TRANCHE}/corr_testing_{args.test_type}_var_{CURRENT_TRANCHE}_new_AF.ht"
    )
    var.export(
        f"gs://ukbb-exome-pharma-wlu/pleiotropy_oct2022/{CURRENT_TRANCHE}/corr_testing_{args.test_type}_var_{CURRENT_TRANCHE}_new_AF.txt.bgz"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--random_phenos", help="use random phenotypes", action="store_true"
    )
    parser.add_argument(
        "--real_phenos", help="use real phenotypes", action="store_true"
    )
    parser.add_argument(
        "--filter_rows",
        help="whether to filter to high quality set of gene& variant",
        action="store_true",
    )
    parser.add_argument(
        "--filter_cols",
        help="whether to filter to high quality set of phenotypes",
        action="store_true",
    )
    parser.add_argument(
        "--phenocode_lst", type=list, help="subset to a list of phenotypes"
    )
    parser.add_argument("--gene_lst", type=list, help="subset to a list of genes")
    parser.add_argument("--irnt", help="use irnt phenotypes", action="store_true")
    parser.add_argument("--icd", help="use icd phenotypes", action="store_true")
    parser.add_argument(
        "--continuous", help="use continuous phenotypes", action="store_true"
    )
    parser.add_argument(
        "--test_type",
        type=str,
        help="Test results to apply lambda filters on: skato OR burden",
    )
    parser.add_argument("--overwrite", help="Overwrite everything", action="store_true")
    parser.add_argument(
        "--af_upper",
        type=float,
        help="Keep genes/variants with higher cumulative allele frequency",
    )
    parser.add_argument(
        "--p_upper",
        type=float,
        help="Keep genes/variants with higher cumulative allele frequency",
    )

    args = parser.parse_args()
    print(args)

main(args)
# hailctl dataproc start wlu --packages ukbb_common,"git+https://github.com/broadinstitute/gnomad_methods.git@main" --requester-pays-allow-all -w 10 --max-idle 10m
# hailctl dataproc submit wlu subset_pleiotropy_data.py --real_phenos --filter_rows --test_type burden --filter_cols --af_upper 0.0001
