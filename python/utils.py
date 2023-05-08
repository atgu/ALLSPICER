from ukb_exomes import *
import hail as hl

RBC = {'30010', '30020', '30030', '30040', '30050', '30060', '30070'}
WBC = {'30000', '30120', '30130', '30140', '30150', '30160', '30180', '30190', '30200', '30210', '30220'}
PLT = {'30080', '30090', '30100', '30110'}
RET = {'30250', '30240', '30260', '30270', '30280', '30290', '30300'}

BLOOD = set.union(RBC, WBC, PLT, RET)
BONE = {'30610', '30680', '30890'}
CARDIOVASCULAR = {'30630', '30640', '30710', '30690', '30760', '30780', '30790', '30870'}
DIABETE = {'30740', '30750'}
HORMONE = {'30770', '30830', '30850'}
LIVER = {'30620', '30600', '30650', '30660', '30730', '30840'}
RENAL = {'30700', '30510', '30720', '30500', '30810', '30520', '30530', '30860', '30880', '30670'}

P_VALUE_FIELDS = {'skato': 'Pvalue', 'skat': 'Pvalue_SKAT', 'burden': 'Pvalue_Burden', 'variant':'Pvalue'}
my_bucket = 'gs://ukbb-exome-pharma-wlu'


def subset_gene_variant_pheno_ht(
        gene_symbol: str,
        phenos: set,
        tranche=CURRENT_TRANCHE
):

    var = hl.read_matrix_table(get_results_mt_path('variant', tranche=tranche))
    sub = var.filter_rows((var.gene == gene_symbol) & (hl.is_defined(var.annotation)))
    sub = sub.filter_cols(hl.literal(phenos).contains(sub.phenocode))

    vep = hl.read_table(var_annotations_ht_path('vep', *TRANCHE_DATA[tranche]))
    vep = process_consequences(vep)
    vep = vep.explode(vep.vep.worst_csq_by_gene_canonical)
    vep = vep.explode(vep.vep.transcript_consequences)
    annotation = annotation_case_builder(vep.vep.worst_csq_by_gene_canonical)

    vep = vep.annotate(
        annotation=hl.if_else(hl.literal({'missense', 'LC'}).contains(annotation), 'missense|LC', annotation))
    vep = vep.filter(hl.is_defined(vep.annotation) & (vep.annotation != 'non-coding'))

    vep = vep.annotate(hgvsp=vep.vep.transcript_consequences.hgvsp,
                       gene_id=vep.vep.worst_csq_by_gene_canonical.gene_id,
                       gene_symbol=vep.vep.worst_csq_by_gene_canonical.gene_symbol,
                       annotation=vep.annotation)
    vep = vep.filter((vep.gene_symbol == gene_symbol))
    sub = sub.drop('AF', 'AC')
    sub = sub.annotate_rows(hgvsp=vep[sub.row_key]['hgvsp'],
                            **sub.call_stats)
    sub = sub.entries()
    sub = sub.select('AC', 'AF', 'gene', 'annotation', 'BETA', 'Pvalue', 'hgvsp')

    return sub

def filtered_mt(test_type, result_type, message=True):
    # pheno_filter = ['keep_pheno_unrelated', f'keep_pheno_{test_type}']
    pheno_filter = [f'keep_pheno_{test_type}']
    if result_type == 'gene':
        row_filter = [f'keep_gene_{test_type}', 'keep_gene_coverage', 'keep_gene_expected_ac', 'keep_gene_n_var']
    else:
        row_filter = ['keep_var_annt', 'keep_var_expected_ac']
    test_type = 'variant' if result_type=='variant' else test_type
    if message:
        print(f'Type of test: {test_type.upper()}')
    mt = load_final_sumstats_table(result_type, extension = 'mt', tranche = CURRENT_TRANCHE)
    if message:
        print(f'[PRE-FILTER] number of {result_type}s: {mt.count_rows()}, number of phenotypes: {mt.count_cols()} \n')
    for f in row_filter:
        mt = mt.filter_rows(mt[f])
    for f in pheno_filter:
        mt = mt.filter_cols(mt[f])
    mt = mt.filter_entries(mt.keep_entry_expected_ac)
    if message:
        print(f'[POST-FILTER] number of {result_type}s: {mt.count_rows()}, number of phenotypes: {mt.count_cols()} \n')
    return mt

def raw_pleiotropy_count(test_type, result_type):
    mt = filtered_mt(test_type, result_type)
    test_type = 'variant' if result_type=='variant' else test_type
    mt = mt.annotate_rows(sig_cnt = hl.agg.count_where(mt[P_VALUE_FIELDS[test_type]]<EMPIRICAL_P_THRESHOLDS[test_type]))
    print(f'Percentage of {result_type}s with 1+ signals: {mt.aggregate_rows(hl.agg.fraction(mt.sig_cnt>1))}\n')
    print(f'Number of {result_type}s with 1+ signals: {mt.aggregate_rows(hl.agg.count_where(mt.sig_cnt>1))}\n')
    print(f'Percentage of {result_type}s with 1+ signals stratified by annotation:{mt.aggregate_rows(hl.agg.group_by(mt.annotation, hl.agg.fraction(mt.sig_cnt>1)))}')
    print(f'Number of {result_type}s with 1+ signals stratified by annotation:{mt.aggregate_rows(hl.agg.group_by(mt.annotation, hl.agg.count_where(mt.sig_cnt>1)))}')

def annotate_icd_chapters_mt(mt: hl.matrixtable):
        mt = mt.annotate_cols(
            pheno_group=hl.case()
                .when(mt.trait_type == "icd10", mt.phenocode[0])
                .when(mt.trait_type == "icd_first_occurrence", mt.description.split("first", 2)[0][5])
                .or_missing(),
            icd10_index=hl.int32(
                hl.case()
                    .when(mt.trait_type == "icd10", mt.phenocode[1:3])
                    .when(mt.trait_type == "icd_first_occurrence", mt.description.split("first", 2)[0][6:8])
                    .or_missing()
            )
        )
        mt = mt.annotate_cols(pheno_group=hl.case()
                              .when((mt.pheno_group == 'B'), 'A')
                              .when((mt.pheno_group == 'D') & (mt.icd10_index < 50), 'C')
                              .when((mt.pheno_group == 'H') & (mt.icd10_index < 60), 'H1')
                              .when((mt.pheno_group == 'H') & (mt.icd10_index >= 60), 'H2')
                              .default(mt.pheno_group))
        return mt


def generate_domain_data(test_type, result_type, domain_type):
    mt = filtered_mt(test_type, result_type, message=False)
    print(mt.count())
    mt = mt.drop('n_cases', 'n_controls', 'heritability', 'saige_version', 'inv_normalized',
                 'n_cases_females', 'n_cases_males', 'description_more', 'coding_description')
    if domain_type == 'icd':
        mt = mt.filter_cols(hl.literal({'icd10', 'icd_first_occurrence'}).contains(mt.trait_type))
        mt = annotate_icd_chapters_mt(mt)
    elif domain_type == 'biomarker':
        mt = mt.filter_cols(mt.phenocode.startswith('30') & (hl.len(mt.phenocode) == 5))
        mt = mt.annotate_cols(pheno_group=hl.case()
                              .when(hl.literal(BONE).contains(mt.phenocode), 'bone and joint')
                              .when(hl.literal(CARDIOVASCULAR).contains(mt.phenocode), 'cardiovascular')
                              .when(hl.literal(DIABETE).contains(mt.phenocode), 'diabetes')
                              .when(hl.literal(HORMONE).contains(mt.phenocode), 'hormone')
                              .when(hl.literal(LIVER).contains(mt.phenocode), 'liver')
                              .when(hl.literal(RENAL).contains(mt.phenocode), 'renal')
                              .when(hl.literal(BLOOD).contains(mt.phenocode), 'blood')
                              .default('other'))
    elif domain_type == 'both':
        mt = mt.filter_cols(hl.literal({'icd10', 'icd_first_occurrence'}).contains(mt.trait_type) |
                            (mt.phenocode.startswith('30') & (hl.len(mt.phenocode) == 5)))
        mt = annotate_icd_chapters_mt(mt)
        mt = mt.annotate_cols(pheno_group=hl.case()
                              .when(hl.literal(BONE).contains(mt.phenocode), 'bone and joint')
                              .when(hl.literal(CARDIOVASCULAR).contains(mt.phenocode), 'cardiovascular')
                              .when(hl.literal(DIABETE).contains(mt.phenocode), 'diabetes')
                              .when(hl.literal(HORMONE).contains(mt.phenocode), 'hormone')
                              .when(hl.literal(LIVER).contains(mt.phenocode), 'liver')
                              .when(hl.literal(RENAL).contains(mt.phenocode), 'renal')
                              .when(hl.literal(BLOOD).contains(mt.phenocode), 'blood')
                              .default(mt.pheno_group))
    elif domain_type == 'blood':
        mt = mt.filter_cols(hl.literal(BLOOD).contains(mt.phenocode))
        mt = mt.annotate_cols(pheno_group=hl.case()
                              .when(hl.literal(RBC).contains(mt.phenocode), 'Red blood cells')
                              .when(hl.literal(WBC).contains(mt.phenocode), 'White blood cells')
                              .when(hl.literal(PLT).contains(mt.phenocode), 'Platelet')
                              .when(hl.literal(RET).contains(mt.phenocode), 'Reticulocyte')
                              .default('other'))

    else:
        mt = mt.annotate_cols(pheno_group=hl.case()
                              .when(mt.category.matches('Mental health'), 'Mental Health')
                              .when(mt.category.matches('Medications'), 'Treatment/medications')
                              .when((mt.category.matches('Medical conditions') & (hl.len(mt.phenocode) == 5)) | (
            mt.category.matches('Medical information')), 'Diseases')
                              # .when(mt.category.matches('Touchscreen'), 'Touchscreen')
                              .when(hl.literal({'icd10', 'icd_first_occurrence'}).contains(mt.trait_type) |
                                    (mt.phenocode == 'COVID19') |
                                    ((mt.modifier == 'custom') & (mt.trait_type == 'categorical')), 'Diseases')
                              .when(mt.category.matches('Physical measures'), 'Physical measures')
                              .when(mt.category.matches('Imaging'), 'Imaging')
                              .when(mt.category.matches('inpatient'), 'Hospital inpatient operations')
                              .when(mt.description.matches('Operation code'), 'Verbal interview operations')
                              # .when(mt.category.matches('Diet by 24-hour recall'), 'Diet')
                              .when((mt.category.matches('Biological samples')), 'Biomarkers')
                              .or_missing())
    test_type = 'variant' if result_type == 'variant' else test_type
    mt = mt.annotate_entries(pheno_sig=hl.if_else(mt[P_VALUE_FIELDS[test_type]] <
                                                  EMPIRICAL_P_THRESHOLDS[test_type], 1, 0))
    mt = mt.annotate_entries(pheno_sig=hl.if_else(hl.is_missing(mt.pheno_sig), 0, mt.pheno_sig))
    print(f'number of phenotypes per group:{mt.aggregate_cols(hl.agg.counter(mt.pheno_group))}')
    sub = mt.group_cols_by(mt.pheno_group).aggregate(pheno_group_sig_cnt=hl.agg.sum(mt.pheno_sig))
    sub = sub.annotate_entries(pheno_group_sig=sub.pheno_group_sig_cnt > 0)
    return sub


def domain_pleiotropy_count(test_type, result_type, domain_type):
    data = generate_domain_data(test_type, result_type, domain_type)
    data = data.annotate_rows(sig_cnt=hl.agg.count_where(data.pheno_group_sig))
    print(
        f'Number of {result_type}s associated with 1+ domain: {data.aggregate_rows(hl.agg.count_where(data.sig_cnt > 1))}\n')
    print(
        f'Number of {result_type}s associated with 1+ domain stratified by annotation:{data.aggregate_rows(hl.agg.group_by(data.annotation, hl.agg.count_where(data.sig_cnt > 1)))}')

def write_domain_level_data(test_type = 'skato', result_type = 'gene', group_type = 'biomarker'):
    sub = generate_domain_data(test_type, result_type, group_type)
    sub = sub.select_rows('interval', 'total_variants')
    sub = sub.entries()
    sub.export(f'{my_bucket}/pleiotropy/{CURRENT_TRANCHE}/pleiotropy_{group_type}_domain_level_{result_type}_sig_{test_type}_{CURRENT_TRANCHE}.txt.bgz')

def get_combined_group_analysis_mt(test_type):
    mt = get_qc_result_mt(result_type='gene', test_type=test_type, tranche = '500k')
    mt = mt.filter_rows(mt.annotation != 'synonymous')
    gene = mt.group_rows_by('gene_id', 'gene_symbol').aggregate(
    **{
        f"{annotation}_sig_{test}":
        hl.agg.count_where((mt.annotation == annotation) & (mt[P_VALUE_FIELDS[test]]< EMPIRICAL_P_THRESHOLDS[test]))
        for annotation in ['pLoF', 'missense|LC', 'pLoF|missense|LC']
        for test in ['skato', 'burden']
    }
    )
    gene = gene.annotate_cols(trait_type2 = hl.if_else(hl.literal({'icd10', 'icd_first_occurrence'}).contains(gene.trait_type), 'icd10', gene.trait_type))
    gene = gene.annotate_rows(
    **{
        f"combined_sig_{test}_{trait}":
        hl.agg.count_where((gene[f'pLoF|missense|LC_sig_{test}'] >0) &
                           (gene[f'pLoF_sig_{test}'] == 0) &
                           (gene[f'missense|LC_sig_{test}'] == 0) &
                           (gene.trait_type2 == trait)
                          )
        for test in ['skato', 'burden']
        for trait in TRAIT_TYPES
    },
    **{
        f"combined_sig_{test}":
        hl.agg.count_where((gene[f'pLoF|missense|LC_sig_{test}'] >0) &
                           (gene[f'pLoF_sig_{test}'] == 0) &
                           (gene[f'missense|LC_sig_{test}'] == 0)
                          )
        for test in ['skato', 'burden']
    },
    )
    gene = gene.annotate_cols(
    **{
        f"combined_gene_sig_{test}":
        hl.agg.count_where((gene[f'pLoF|missense|LC_sig_{test}'] >0) &
                           (gene[f'pLoF_sig_{test}'] == 0) &
                           (gene[f'missense|LC_sig_{test}'] == 0)
                          )
        for test in ['skato', 'burden']
    },
    )
    return gene