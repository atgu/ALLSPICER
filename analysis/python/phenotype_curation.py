from ukb_exomes import *
from utils import *
import hail as hl
import argparse
from gnomad.utils.vep import process_consequences
from gnomad.resources.grch38.reference_data import clinvar
from gnomad.resources.resource_utils import DataException
from gnomad.utils.vep import vep_or_lookup_vep
import subprocess
import re
import pickle


hl.init(tmp_dir = f'{my_bucket}/tmp/')


def get_pleiotropy_gene_mt(result_type:str = 'gene'):
    ht = hl.read_table(correlated_pheno_ht)
    mt = load_final_sumstats_table(result_type=result_type, extension="mt")
    mt = mt.annotate_cols(
        keep_pheno_unrelated=hl.is_missing(
            ht.key_by(
                **{field: ht.node[field] for field in mt.col_key}
            )[mt.col_key]
        ),
    )
    return mt


def get_filtered_pleiotropy_mt(individual_level, result_type:str='gene', test_type:str='burden', message=True):
    if individual_level:
        pheno_filter = ['keep_pheno_unrelated', f'keep_pheno_{test_type}']
    else:
        pheno_filter = [f'keep_pheno_{test_type}']
    if result_type == 'gene':
        row_filter = [f'keep_gene_{test_type}', 'keep_gene_coverage', 'keep_gene_expected_ac', 'keep_gene_n_var']
    else:
        row_filter = ['keep_var_annt', 'keep_var_expected_ac']
        # row_filter = ['keep_var_annt']
    test_type = 'variant' if result_type=='variant' else test_type
    if message:
        print(f'Type of test: {test_type.upper()}')
        print(CURRENT_TRANCHE)
    mt = get_pleiotropy_gene_mt(result_type=result_type)
    if message:
        print(f'[PRE-FILTER] number of {result_type}s: {mt.count_rows()}, number of phenotypes: {mt.count_cols()} \n')
    for f in row_filter:
        mt = mt.filter_rows(mt[f])
    for f in pheno_filter:
        mt = mt.filter_cols(mt[f])
    # mt = mt.filter_entries(mt.keep_entry_expected_ac)
    if message:
        print(f'[POST-FILTER] number of {result_type}s: {mt.count_rows()}, number of phenotypes: {mt.count_cols()} \n')
    return mt


def format_phenotypes(mt):
    # Edit phenotypes with missing categories
    mt = mt.annotate_cols(
        category=hl.if_else((mt.modifier == 'custom') & (hl.is_missing(mt.category)), 'custom', mt.category),
        description=hl.if_else((mt.modifier == 'custom') & (hl.is_missing(mt.description)), 'custom', mt.description))
    mt = mt.annotate_cols(category=hl.if_else(mt.phenocode == 'COVID19',
                                              'Health-related outcomes > First occurrences > Respiratory system disorders',
                                              mt.category))
    mt = mt.annotate_cols(
        category=hl.if_else(mt.category == '', 'UK Biobank Assessment Centre > Physical measures', mt.category))
    # Edit self-reported disease phenotype
    mt = mt.annotate_cols(note=hl.if_else((mt.description == 'Cancer code, self-reported') |
                                          (mt.description == 'Non-cancer illness code, self-reported') |
                                          (mt.category.matches(
                                              'UK Biobank Assessment Centre > Touchscreen > | Online follow-up |UK Biobank Assessment Centre > Verbal interview >')),
                                          'self-reported', ''),
                          description=hl.if_else((mt.description == 'Cancer code, self-reported') | (
                                      mt.description == 'Non-cancer illness code, self-reported') | (
                                                             mt.description == 'Mental health problems ever diagnosed by a professional'),
                                                 mt.coding_description, mt.description))

    # annotate_pheno_groups
    mt = mt.annotate_cols(pheno_group=hl.case()
                          .when(mt.category.matches('Diet by 24-hour recall'), 'Diet')
                          .when(mt.category.matches('Cognitive function'), 'Cognitive function')
                          .when(mt.phenocode == 'Fluid_intelligence_score_custom', 'Cognitive function')
                          .when(
        (mt.phenocode == '20544') | (mt.phenocode == '20126') | (mt.phenocode == 'Depressive_symptoms_custom'),
        'Diseases')
                          .when(mt.category.matches('Mental health'), 'Mental health')
                          .when(
        hl.literal({"Depressive_symptoms_custom", "insomnia_1200_quant"}).contains(mt.phenocode), 'Mental health')
                          .when((mt.category.matches('Lifestyle and environment')), 'Lifestyle')
                          .when(mt.category.matches('Medication'), 'Treatment/medications')
                          .when(
        (mt.category.matches('Medical conditions|Medical information|Health and medical history')), 'Diseases')
                          .when(hl.literal({'icd10', 'icd_first_occurrence'}).contains(mt.trait_type) |
                                (mt.phenocode == 'COVID19') |
                                ((mt.modifier == 'custom') & (mt.trait_type == 'categorical')), 'Diseases')
                          .when(mt.category.matches('Physical measures'), 'Physical measures')
                          .when(mt.phenocode == '20022', 'Physical measures')
                          .when(mt.category.matches('Brain MRI'), 'Brain imaging')
                          .when((mt.category.matches('Biological samples')), 'Biomarkers')
                          .when(mt.modifier == 'custom', 'Eye measures')
                          .when(mt.category.matches('Imaging'), 'Physical measures')
                          .or_missing())
    mt = mt.annotate_cols(pheno_group=hl.case()
                          .when((mt.category.matches(
        'UK Biobank Assessment Centre > Touchscreen > Early life factors')) | mt.description.matches('Age | age ') |
                                hl.literal(
                                    {'120', '189', '630', '2724', '22200', 'Touchscreen_duration_custom', '12651',
                                     '2217'}).contains(mt.phenocode) |
                                ((mt.phenocode == '6152') & (mt.coding == '100')), 'Others')
                          .when((mt.category.matches('inpatient') |
                                 mt.description.matches('Operation code')), 'Operations')
                          .default(mt.pheno_group))

    return mt


icd_names = {'A': 'Infectious', 'B': 'Infectious', 'C': 'Neoplasms', 'D': 'Blood/immune', 'E': 'Endocrine/metabolic',
             'F': 'Mental/behavioral', 'G': 'Nervous', 'H1': 'Eye',
             'H2': 'Ear', 'I': 'Circulatory', 'J': 'Respiratory', 'K': 'Digestive', 'L': 'Skin/subcutaneous',
             'M': 'Musculoskeletal', 'N': 'Genitourinary', 'O': 'Pregnancy',
             'P': 'Perinatal', 'Q': 'Congenital', 'R': 'Symptoms', 'S': 'Injury/poison', 'T': 'Injury/poison',
             'V': 'External causes', 'Y': 'External causes', 'Z': 'Health Factors'}


def annotate_pheno_group_mt(mt: hl.matrixtable):
    # Annotate ICD chapters
    mt = mt.annotate_cols(
        disease_group=hl.case()
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
    mt = mt.annotate_cols(disease_group=hl.case()
                          .when((mt.disease_group == 'B'), 'A')
                          .when((mt.disease_group == 'D') & (mt.icd10_index < 50), 'C')
                          .when((mt.disease_group == 'H') & (mt.icd10_index < 60), 'H1')
                          .when((mt.disease_group == 'H') & (mt.icd10_index >= 60), 'H2')
                          .default(mt.disease_group))
    # Group other diseases
    mt = mt.annotate_cols(description=hl.if_else(mt.description == 'custom', mt.phenocode, mt.description))
    mt = mt.annotate_cols(description=mt.description.replace('_', ' '))
    mt = mt.annotate_cols(disease_group=hl.if_else((mt.pheno_group == 'Diseases') & hl.is_missing(mt.disease_group),
                                                   hl.case()
                                                   .when(mt.description == 'eye trauma', 'S')
                                                   .when(mt.description_more.contains('Code for cancer'), 'C')
                                                   .default(mt.disease_group), mt.disease_group))
    mt = mt.annotate_cols(
        disease_group=hl.if_else((mt.pheno_group == 'Diseases') & hl.is_missing(mt.disease_group), hl.case()
                                 .when(mt.description.matches('tuberculosis|sepsis|chickenpox|shingles'), 'A')
                                 .when(mt.description.matches('pituitary adenoma/tumour'), 'C')
                                 .when(mt.description.matches('Cancer|cancer'), 'C')
                                 .when(mt.description.matches('anaemia|haematological'), 'D')
                                 .when(mt.description.matches('Diabetes|diabetes|thyroiditis|cholesterol|goitre'), 'E')
                                 .when(
            mt.description.matches('Depression|depression|anxiety|Schizophrenia|alcohol|Bipolar|Depressive'), 'F')
                                 .when(mt.phenocode == '20544', 'F')
                                 .when(mt.description.matches(
            'Alzheimer|alzheimer|sleep|nerve|Epilepsy|Parkinsons| Sclerosis| sclerosis|ALS|insomnia|neuralgia|dementia|carpal tunnel|fatigue'),
                                       'G')
                                 .when(mt.description.matches('Eye|eye|retinal|cataract|iritis|Glaucoma|macular|sicca'),
                                       'H1')
                                 .when(mt.description.matches('vestibular|labyrinthitis|vertigo|otosclerosis|tiniitis'),
                                       'H2')
                                 .when(hl.literal({'4803'}).contains(mt.phenocode), 'H2')
                                 .when(mt.description.matches(
            'Heart|heart|vascular|hypertension|haemorrhage|claudication|stroke|Stroke|hemorrhage|vein|aneurysm|VTE|blood|CAD|Atrial|raynaud|svt'),
                                       'I')
                                 .when(mt.description.matches(
            'respiratory|asthma|bronchiectasis|throat|pneumothorax|emphysema|COPD|rhinitis|nasal polyp|bronchitis|dust'),
                                       'J')
                                 .when(mt.category.matches('Hepatobiliary'), 'K')
                                 .when(mt.description.matches(
            'Crohn|colitis|Colitis|stomach|oesophageal|liver|oesophagus|gastritis|gall|cholecystitis|IBD|Celiac|Bowel|haemorrhoids|dyspepsia|constipation|hernia|hepatitis'),
                                       'K')
                                 .when(
            mt.description.matches('Atopic Dermatitis|acne|Alopecia|Psoriasis|skin|urticaria|cellulitis'), 'L')
                                 .when(mt.description.matches(
            'Chest pain|pleurisy|abdominal|urinary frequency|breakdown|Weigh|headaches|jaundice|breath|Wheeze'), 'R')
                                 .when(mt.description.matches(
            'pain|Ankylosing spondylitis|bone|back|Lupus|lupus|joint|muscle|RA|sciatica|disc|fibromyalgia|elbow|tendonitis'),
                                       'M')
                                 .when(mt.description.matches('renal|kidney|ureteric|cystitis|UTI'), 'N')
                                 .when(mt.description.matches('injury|fracture'), 'S')
                                 .when(mt.description.matches('burns|anaphylaxis'), 'T')
                                 .when(mt.description.matches('suicide'), 'Y')
                                 .when(mt.description.matches('stress|Falls|COVID-19|food|drug'), 'Z')
                                 .default(mt.disease_group), mt.disease_group))
    mt = mt.annotate_cols(disease_group=hl.or_missing(mt.pheno_group == 'Diseases', mt.disease_group))
    # Annotate disease name
    mt = mt.annotate_cols(disease_name=hl.literal(icd_names).get(mt.disease_group))

    # Annotate biomarker group
    mt = mt.annotate_cols(biomarker_group=hl.case()
                          .when(hl.literal(BONE).contains(mt.phenocode), 'bone and joint')
                          .when(hl.literal(CARDIOVASCULAR).contains(mt.phenocode), 'cardiovascular')
                          .when(hl.literal(DIABETE).contains(mt.phenocode), 'diabetes')
                          .when(hl.literal(HORMONE).contains(mt.phenocode), 'hormone')
                          .when(hl.literal(LIVER).contains(mt.phenocode), 'liver')
                          .when(hl.literal(RENAL).contains(mt.phenocode), 'renal')
                          .when(hl.literal(BLOOD).contains(mt.phenocode), 'blood')
                          .or_missing())
    # Annotate treatment description
    mt = mt.annotate_cols(
        description=hl.if_else(mt.pheno_group == 'Treatment/medications', mt.coding_description, mt.description))

    return mt

def get_phenotype_curated_mt(pheno_groups_to_keep: list=None, result_type: str = 'gene', test_type: str='burden'):
    mt = get_filtered_pleiotropy_mt(individual_level = True, result_type=result_type, test_type=test_type, message=False)
    mt = format_phenotypes(mt)
    mt = annotate_pheno_group_mt(mt)
    if pheno_groups_to_keep is not None:
        mt = mt.filter_cols(hl.literal(pheno_groups_to_keep).contains(mt.pheno_group))
    return mt

def generate_domain_data(pheno_groups_to_keep, test_type, result_type, domain_type):
    mt = get_phenotype_curated_mt(pheno_groups_to_keep, result_type, test_type)
    mt = mt.drop('n_controls', 'heritability', 'saige_version', 'inv_normalized',
                 'n_cases_females', 'n_cases_males', 'description_more', 'coding_description')
    GROUP_NAME = 'pheno_group' if domain_type == 'pheno' else ('disease_group' if domain_type == 'icd' else 'biomarker_group')
    test_type = 'variant' if result_type=='variant' else test_type
    mt = mt.annotate_entries(pheno_sig=hl.if_else(mt[P_VALUE_FIELDS[test_type]] <
                                                  PHEWAS_P_CUTOFF, 1, 0))
    mt = mt.annotate_entries(pheno_sig=hl.if_else(hl.is_missing(mt.pheno_sig), 0, mt.pheno_sig))
    print(f'number of phenotypes per group:{mt.aggregate_cols(hl.agg.counter(mt[GROUP_NAME]))}')
    sub = mt.group_cols_by(mt[GROUP_NAME]).aggregate(pheno_group_sig_cnt = hl.agg.sum(mt.pheno_sig))
    sub = sub.annotate_entries(pheno_group_sig = sub.pheno_group_sig_cnt > 0)
    return sub


def write_domain_level_data(pheno_groups_to_keep, test_type = 'skato', result_type = 'gene', group_type = 'biomarker'):
    sub = generate_domain_data(pheno_groups_to_keep, test_type, result_type, group_type)
    sub = sub.select_rows('total_variants', 'CAF', 'mean_coverage')
    sub = sub.entries()
    sub.export(f'{my_bucket}/pleiotropy_{group_type}_domain_level_{result_type}_sig_{test_type}_{CURRENT_TRANCHE}.txt.bgz')


def main(args):
    ## Check number of phenotypes per pheno_group
    pheno_groups_to_keep = None
    mt = get_phenotype_curated_mt(pheno_groups_to_keep=pheno_groups_to_keep, result_type='gene', test_type='burden')
    print(mt.aggregate_cols(hl.agg.counter(mt.pheno_group)))

    ## Compute p_phewas
    pheno_groups_to_keep = {'Biomarkers', 'Brain imaging', 'Cognitive function', 'Diseases', 'Physical measures'}
    mt = get_phenotype_curated_mt(pheno_groups_to_keep=pheno_groups_to_keep, result_type='gene', test_type='burden')
    n_pheno = mt.count_cols()
    PHEWAS_P_CUTOFF = 0.05 / n_pheno
    print(f'PHEWAS_P_CUTOFF: {PHEWAS_P_CUTOFF}')

    ## Write Domain-level data
    write_domain_level_data(
        pheno_groups_to_keep={'Biomarkers', 'Brain imaging', 'Cognitive function', 'Diseases', 'Physical measures'},
        test_type='burden', result_type='gene', group_type='icd')
    write_domain_level_data(
        pheno_groups_to_keep={'Biomarkers', 'Brain imaging', 'Cognitive function', 'Diseases', 'Physical measures'},
        test_type='burden', result_type='gene', group_type='pheno')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    args = parser.parse_args()
    main(args)
