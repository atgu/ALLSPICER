import hail as hl
import logging
import argparse
from ukb_exomes import *
from ukbb_common import *

CURRENT_TRANCHE = "500k"
logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("subset_pleiotropy_data")
logger.setLevel(logging.INFO)

def main(args):
    print(CURRENT_TRANCHE)
    logger.info(f"Current tranche: {CURRENT_TRANCHE} data... ")
    final_pheno_ht = hl.read_table('gs://ukbb-exome-pharma-wlu/wrap_up_results/final_599_phenotypes.ht')

    gene = load_final_sumstats_table(result_type="gene", extension="mt", tranche=CURRENT_TRANCHE)
    print(f'Number of genes pre-filter: {gene.count()}')
    gene = gene.filter_rows(
        gene[f"keep_gene_burden"]
        & gene['keep_gene_expected_ac' if CURRENT_TRANCHE=='500k' else 'keep_gene_caf']
        & gene.keep_gene_coverage
        & gene.keep_gene_n_var
    )
    gene = gene.filter_cols(hl.is_defined(final_pheno_ht[gene.col_key]))
    print(f'Number of genes post-filter: {gene.count()}')

    var = load_final_sumstats_table(result_type="variant", extension="mt", tranche=CURRENT_TRANCHE)
    print(f'Number of variants pre-filter: {var.count()}')
    var = var.filter_rows(var.keep_var_annt)
    var = var.filter_cols(hl.is_defined(final_pheno_ht[var.col_key]))
    print(f'Number of variants post-filter: {var.count()}')


    if args.continuous:
        gene = gene.filter_cols(gene.trait_type == "continuous")
        var = var.filter_cols(var.trait_type == "continuous")

    if CURRENT_TRANCHE == "500k":
        gene = gene.filter_rows(gene.annotation != "pLoF|missense|LC")
        var = var.drop("AC", "AF")
        var = var.annotate_rows(**var.call_stats)

    gene = gene.annotate_rows(
        **{
            f"sig_burden": hl.agg.count_where(
                gene.Pvalue_Burden < 2.5e-6
            )
        }
    )
    gene = gene.filter_rows((gene[f"sig_burden"] > 1))
    gene = gene.entries()
    gene.describe()
    var = var.entries()
    var.describe()

    gene = gene.select(
        'Pvalue_Burden',
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
    if args.p_upper is not None:
        var = var.filter(
            var.Pvalue < args.p_upper
        )
        var = var.annotate_globals(
            p_upper=args.p_upper
        )

    gene.write(
        f"gs://ukbb-exome-pharma-wlu/wrap_up_results/corr_testing_burden_gene_{CURRENT_TRANCHE}.ht",
        overwrite=args.overwrite,
    )
    var.write(
        f"gs://ukbb-exome-pharma-wlu/wrap_up_results/corr_testing_burden_var_{args.af_upper*100}{'' if args.p_upper is None else f'_p_{args.p_upper}'}_{CURRENT_TRANCHE}.ht",
        overwrite=args.overwrite,
    )
    gene = hl.read_table(
        f"gs://ukbb-exome-pharma-wlu/wrap_up_results/corr_testing_burden_gene_{CURRENT_TRANCHE}.ht"
    )
    gene.export(
        f"gs://ukbb-exome-pharma-wlu/wrap_up_results/corr_testing_burden_gene_{CURRENT_TRANCHE}.txt.bgz"
    )
    var = hl.read_table(
        f"gs://ukbb-exome-pharma-wlu/wrap_up_results/corr_testing_burden_var_{args.af_upper*100}{'' if args.p_upper is None else f'_p_{args.p_upper}'}_{CURRENT_TRANCHE}.ht"
    )
    var.export(
        f"gs://ukbb-exome-pharma-wlu/wrap_up_results/corr_testing_burden_var_{args.af_upper*100}{'' if args.p_upper is None else f'_p_{args.p_upper}'}_{CURRENT_TRANCHE}.txt.bgz"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--continuous", help="use continuous phenotypes", action="store_true"
    )
    parser.add_argument("--overwrite", help="Overwrite everything", action="store_true")
    parser.add_argument(
        "--af_upper",
        type=float,
        help="Keep variants with higher cumulative allele frequency",
    )
    parser.add_argument(
        "--p_upper",
        type=float,
        help="Keep single variants with significant association",
        default=None
    )
    args = parser.parse_args()
    print(args)

main(args)
# hailctl dataproc start wlu --packages ukbb_common,"git+https://github.com/broadinstitute/gnomad_methods.git@main" --requester-pays-allow-all -w 10 --max-idle 10m
# hailctl dataproc submit wlu subset_pleiotropy_data.py --real_phenos --filter_rows --test_type burden --filter_cols --af_upper 0.0001
