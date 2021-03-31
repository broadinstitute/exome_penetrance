import hail as hl


PLINK_DATA_SOURCES = ["UKBB_regeneron", "UKBB_gatk"]
VCF_DATA_SOURCES = ["52K", "UKBB_regeneron_vcf", "UKBB_gatk_vcf"]
DATA_SOURCES = ['UKBB_regeneron',"UKBB_gatk", "52K"]


gene_set = {"ABCC8", "AKT2", "APOA1", "APOA5", "APOB", "APPL1", "BLK", "CEL", "CETP",
            "GATA4", "GATA6", "GCK", "HNF1A", "HNF1B", "HNF4A", "INS", "KCNJ11",
            "KLF11", "LDLR", "LMNA", "LPL", "MC4R", "NEUROD1", "PAX4", "PCSK9", "PDX1",
            "PLIN1", "PPARG"}

# Removed CEL and PLIN1 because they seem to be specific framshift variants not necessarily LoF mechanism
lof_gene_set = {"APOA5", "APOB", "CEL", "CETP", "GCK", "HNF1A", "HNF1B", "HNF4A", "INS", "LDLR",
                "LMNA", "LPL", "MC4R", "PCSK9", "PDX1", "PPARG"}



path_to_52K_samples_to_filter = "gs://gnomad-julia/penetrance/52K/samples_to_keep.txt"

tissues = ['Spleen','Brain_FrontalCortex_BA9_','SmallIntestine_TerminalIleum','Skin_SunExposed_Lowerleg_',
           'Artery_Coronary','Brain_Hippocampus','Esophagus_Muscularis','Brain_Nucleusaccumbens_basalganglia_',
           'Artery_Tibial','Brain_Hypothalamus','Adipose_Visceral_Omentum_','Nerve_Tibial',
           'Brain_CerebellarHemisphere','Breast_MammaryTissue','Liver','Skin_NotSunExposed_Suprapubic_',
           'AdrenalGland','Pancreas','Lung','Pituitary','Muscle_Skeletal','Colon_Transverse','Artery_Aorta',
           'Heart_AtrialAppendage','Adipose_Subcutaneous','Esophagus_Mucosa','Heart_LeftVentricle','Brain_Cerebellum',
           'Brain_Cortex','Thyroid','Stomach','WholeBlood','Brain_Anteriorcingulatecortex_BA24_',
           'Brain_Putamen_basalganglia_','Brain_Caudate_basalganglia_','Colon_Sigmoid',
           'Esophagus_GastroesophagealJunction','Brain_Amygdala']


v7_tissues_to_drop = ["Bladder", "Brain_Spinalcord_cervicalc_1_", "Brain_Substantianigra",
                      "Cervix_Ectocervix","Cervix_Endocervix", "FallopianTube", "Kidney_Cortex",
                      "MinorSalivaryGland", "Uterus", "Ovary","Testis", "Vagina",
                      "Cells_EBV_transformedlymphocytes", "Cells_Transformedfibroblasts", "Prostate"]
v8_tissues_to_drop = ["Bladder", "Brain_Spinalcord_cervicalc_1_", "Brain_Substantianigra",
                      "Cervix_Ectocervix","Cervix_Endocervix", "FallopianTube", "Kidney_Cortex",
                      "MinorSalivaryGland", "Uterus", "Ovary","Testis", "Vagina",
                      "Cells_EBV_transformedlymphocytes","Cells_Culturedfibroblasts", "Prostate"]

def get_tx_annotation_kt_path(data_source: str) -> str:
    if data_source == "52K":
        return "gs://gnomad-public/papers/2019-tx-annotation/data/GTEx.v7.gene_expression_per_gene_per_tissue.120518.kt"
    else:
        return "gs://gnomad-public/papers/2019-tx-annotation/data/GTEx.v8.gene_expression_per_gene_per_tissue.042319.kt"


def get_tx_annotation_path(data_source: str) -> str:
    if data_source == "52K":
        return "gs://gnomad-public/papers/2019-tx-annotation/data/GTEx.V7.tx_medians.110818.ht"
    else:
        return "gs://gnomad-public/papers/2019-tx-annotation/data/GTEx.V8.tx_medians.042319.ht"

def get_tx_annotation_tissues_to_drop(data_source: str) -> str:
    if data_source == "52K":
        return v7_tissues_to_drop
    else:
        return v8_tissues_to_drop

def data_prefix(data_source: int) -> str:
    if data_source not in DATA_SOURCES:
        raise DataException("This data_source is currently not present")
    return f'gs://gnomad-julia/penetrance/{data_source}'


def get_penetrance_vcf_data(data_source):
    if data_source != "52K":
        raise DataException("This data_source doesn't have vcf data")
    mt = hl.import_vcf("gs://ukbb_josep_t2d/55k.clean.all.vcf.bgz")
    return mt


def get_ukbb_plink_data(data_source) -> str:
    if data_source == "UKBB_regeneron":
        mt = hl.import_plink(bed="gs://fc-fdd512d3-61cc-4e5e-8701-c33342a9feb4/wave01/plink/ukb_evc_chr1_v1.bed",
                            bim="gs://fc-fdd512d3-61cc-4e5e-8701-c33342a9feb4/wave01/plink/ukb_spb_exm_chrall_v1.bim",
                            fam=f'{data_prefix(data_source)}/ukb27892_evc_chr1_v1_s49959.fam',
                            reference_genome="GRCh38",
                            skip_invalid_loci=True)
    elif data_source == "UKBB_gatk":
        mt = hl.import_plink(bed="gs://fc-72d33328-e60d-4e5a-96e2-03fe2a0c8ae8/wave01/plink/ukb_efe_chr1_v1.bed",
                            bim="gs://fc-72d33328-e60d-4e5a-96e2-03fe2a0c8ae8/wave01/plink/ukb_fe_exm_chrall_v1.bim",
                             fam=f'{data_prefix(data_source)}/ukb27892_efe_chr1_v1_s49959.fam',
                             reference_genome="GRCh38",
                             skip_invalid_loci=True)
    else:
        raise DataException("This data_source doesn't have plink data")
    return mt


def get_gencode_annotation_path(data_source: str) -> str:
    if data_source == "52K":
        return "gs://gnomad-zach/data/external/gencode/gencode_genes.json.gz"
    else:
        return "gs://gnomad-julia/penetrance/resources/gencode.v29.annotation.json"

def get_vep_config_path(data_source: str) -> str:
    if data_source == "52K":
        return "gs://hail-common/vep/vep/vep85-loftee-gcloud.json"
    else:
        return "gs://hail-common/vep/vep/vep95-GRCh38-loftee-gcloud.json"


def get_penetrance_path(data_source: str, data_type: str) -> str:
    """
    Wrapper function to get paths to penetrance data
    :param str data_source: UKBB_regeneron, or UKBB_gatk
    :param str data_type: raw, split, vep, qc, lof, lof_pext
    :return: Path to chosen MT
    :rtype: str
    """
    if data_source not in DATA_SOURCES:
        raise DataException("This data_source is currently not present")
    data_type = f'.{data_type}' if data_type != "raw" else ''
    mt_or_ht = 'ht' if data_type in ['vep_csq','HNF1A_vep_csq'] else 'mt'
    return f'{data_prefix(data_source)}/{data_source}{data_type}.{mt_or_ht}'


def get_path_to_filtered_variant_file(data_source: str) -> str:
    return f'{data_prefix(data_source)}/{data_source}.filtered_variants.gz'

def get_path_to_lof_vcf(data_source: str) -> str:
    return f'{data_prefix(data_source)}/{data_source}.lof.vcf.bgz'

def get_path_to_qc_vcf(data_source: str) -> str:
    return f'{data_prefix(data_source)}/{data_source}.qc.vcf.bgz'

def get_path_to_pext_worstcsq_tsv(data_source: str) -> str:
    return f'{data_prefix(data_source)}/{data_source}.lof.pext.worstcsq.tsv.bgz'



ukbb_calling_intervals_path = 'gs://gnomad-julia/penetrance/resources/ukbb_exome_calling.interval_list'
lcr_intervals_path = 'gs://gnomad-julia/penetrance/resources/LCRFromHengH38_chr1-22_XY.txt'

# Sample QC files
def sample_qc_prefix(data_source: str) -> str:
    if data_source not in DATA_SOURCES:
        raise DataException("This data_source is currently not present")
    return f'{data_prefix(data_source)}/sample_qc'



def sex_ht_path(data_source: str) -> str:
    return f'{sample_qc_prefix(data_source)}/sex_check/sex.ht'

def qc_mt_path(data_source: str, ld_pruned: bool = False) -> str:
    """
    Returns path of MatrixTable for sample QC purposes
    :param bool ld_pruned: Should the qc matrix be LD pruned
    :return: Path MatrixTable for sample QC purposes
    :rtype: str
    """
    ld_pruned = '.pruned' if ld_pruned else ''
    return f'{sample_qc_prefix(data_source)}/high_callrate_common_biallelic_snps{ld_pruned}.mt'

def qc_ht_path(data_source: str) -> str:
    return f'{sample_qc_prefix(data_source)}/high_callrate_common_biallelic_snps.ht'

def raw_qc_ht_path(data_source: str) -> str:
    return f'{sample_qc_prefix(data_source)}/sample_qc.ht'

def callrate_mt_path(data_source: str) -> str:
    return f'{sample_qc_prefix(data_source)}/platform_pca/callrate.mt'

def platform_pca_scores_ht_path(data_source: str) -> str:
    return f'{sample_qc_prefix(data_source)}/platform_pca/platform_pca_scores.ht'

def platform_pca_loadings_ht_path(data_source: str) -> str:
    return f'{sample_qc_prefix(data_source)}/platform_pca/platform_pca_loadings.ht'

def platform_pca_results_ht_path(data_source: str) -> str:
    return f'{sample_qc_prefix(data_source)}/platform_pca/platform_pca_results.ht'

def relatedness_pca_scores_ht_path(data_source: str) -> str:
    return f'{sample_qc_prefix(data_source)}/relatedness/pruned.pca_scores.ht'

def relatedness_ht_path(data_source: str) -> str:
    return f'{sample_qc_prefix(data_source)}/relatedness/relatedness.ht'

def duplicates_ht_path(data_source: str, dup_sets: bool = False) -> str:
    dup_sets = f'_sets' if dup_sets else ''
    return f'{sample_qc_prefix(data_source)}/relatedness/duplicate{dup_sets}.ht'

def inferred_ped_path(data_source: str) -> str:
    return f'{sample_qc_prefix(data_source)}/relatedness/ped.txt'

def related_drop_path(data_source: str) -> str:
    return f'{sample_qc_prefix(data_source)}/relatedness/related_samples_to_drop.ht'

def ancestry_pca_scores_ht_path(data_source: str, population: str = None) -> str:
    pop = f'.{population}' if population else ''
    return f'{sample_qc_prefix(data_source)}/population_pca/unrelated.pca_scores{pop}.ht'

def ancestry_pca_loadings_ht_path(data_source: str, population: str = None) -> str:
    pop = f'.{population}' if population else ''
    return f'{sample_qc_prefix(data_source)}/population_pca/unrelated.pca_loadings{pop}.ht'

def ancestry_pc_project_scores_ht_path(data_source: str, data_type: str = None) -> str:
    """
    Returns path of Table for scores and pop assignments from pc_project on gnomAD PCs
    :param str data_type: either None for UKBB only or joint for merged UKBB and gnomAD
    :return: Path to Table
    :rtype: str
    """
    data_type = f'.{data_type}' if data_type else ''
    return f'{sample_qc_prefix(data_source)}/population_pca/pc_project_scores_pop_assign{data_type}.ht'

def platform_pop_outlier_ht_path(data_source: str, variant_class_prefix: str) -> str:
    return f'{sample_qc_prefix(data_source)}/outlier_detection_{variant_class_prefix}.ht'

def qc_temp_data_prefix(data_source: str):
    return f'{sample_qc_prefix(data_source)}/temp/'


def meta_path(data_source: str) -> str:
    return f'{sample_qc_prefix(data_source)}/meta.ht'


class DataException(Exception):
    pass