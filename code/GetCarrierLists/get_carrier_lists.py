from gnomad_hail import *
import hail as hl
import argparse
import json
from tx_annotation import *
from penetrance.resources import *


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("ukbb_penetrance")
logger.setLevel(logging.INFO)


def recode_interval(i, buffer=10):
    chrom, window = i.split(':')
    start, end = map(int, window.split('-'))
    start -= buffer
    end += buffer
    if start < 1:
        start = 1
    return "{0}:{1}-{2}".format(chrom, start, end)


def main(args):
    data_source = args.data_source

    if data_source == "52K":
        reference = 'GRCh37'
    else:
        reference = 'GRCh38'

    hl.init(log='/ukb_penetrance.log', default_reference=reference, tmp_dir='hdfs:///write_vcf.tmp/')

    if not args.skip_import_data:
        logger.info(f'Importing {data_source} into hail Matrix Table...')
        if data_source in PLINK_DATA_SOURCES:
            mt = get_ukbb_plink_data(data_source)
        if data_source in VCF_DATA_SOURCES:
            mt = get_penetrance_vcf_data(data_source)
            if data_source == "52K":
                logger.info('Filtering forbidden samples')
                filter_ht = hl.import_table(path_to_52K_samples_to_filter, impute=True).key_by("s")
                mt = mt.filter_cols(filter_ht[mt.col_key].is_forbidden, keep=False)
        mt = mt.checkpoint(get_penetrance_path(data_source, "raw"), overwrite=args.overwrite)
        variants, samples = mt.count()
        logger.info(f'{variants} variants and {samples} samples found in {data_source} data')


    if not args.skip_filter_intervals:
        logger.info(f'Import gencode data...')
        gencode = json.load(hl.utils.hadoop_open(get_gencode_annotation_path(data_source)))
        assert all(gene in gencode for gene in gene_set)
        intervals = [hl.parse_locus_interval(recode_interval(gencode[gene]['interval'], buffer=100000)) for gene in
                     gene_set]
        mt = hl.read_matrix_table(get_penetrance_path(data_source, "raw"))
        logger.info(f'filtering to intervals of interest...')
        mt = hl.filter_intervals(mt, intervals)

        logger.info(f'Splitting multi-allelic variants...')
        if data_source in PLINK_DATA_SOURCES:
            mt = hl.split_multi(mt)
        if data_source in VCF_DATA_SOURCES:
            mt = hl.split_multi_hts(mt)
        mt.write(get_penetrance_path(data_source, "split"), overwrite=args.overwrite)


    if not args.skip_genotype_filter:
        if data_source in PLINK_DATA_SOURCES:
            raise DataException("Genotype filtering not available for PLINK dataset")
        mt = hl.read_matrix_table(get_penetrance_path(data_source, "split"))
        ab = mt.AD[1] / hl.sum(mt.AD)
        filter_condition = ((mt.GQ < 20) | (mt.DP < 10) | (mt.GT.is_het() & (ab < 0.25)))
        mt = mt.annotate_entries(GT=hl.cond(filter_condition, hl.null(hl.tcall), mt.GT))
        mt.write(get_penetrance_path(data_source, "gt_filter"), overwrite=args.overwrite)


    if not args.skip_variant_qc:
        logger.info(f'Running variant qc...')
        if data_source in PLINK_DATA_SOURCES:
            mt = hl.read_matrix_table(get_penetrance_path(data_source, "split"))
        else:
            mt = hl.read_matrix_table(get_penetrance_path(data_source, "gt_filter"))
        mt = hl.variant_qc(mt)
        mt = mt.filter_rows(mt.variant_qc.AC[1] > 0)
        mt = mt.checkpoint(get_penetrance_path(data_source, "qc"), overwrite=args.overwrite)

        logger.info(f'Exporting qc vcf...')
        mt = mt.annotate_rows(info=hl.struct(
                AC=mt.variant_qc.AC[1], AN=mt.variant_qc.AN, AF=mt.variant_qc.AF[1],
                n_hom=mt.variant_qc.homozygote_count[1]))
        hl.export_vcf(mt, get_path_to_qc_vcf(data_source))


    if not args.skip_vep:
        # Note: Need to start cluster with --vep --vep-reference
        logger.info(f'Running VEP on the filtered and split MT...')
        mt = hl.read_matrix_table(get_penetrance_path(data_source, "qc"))
        vep_ht = mt.rows().select()
        vep_ht = hl.vep(vep_ht, get_vep_config_path(data_source), csq=True)
        vep_ht.write(get_penetrance_path(data_source, "vep_csq"), overwrite=args.overwrite)
        mt = hl.vep(mt, get_vep_config_path(data_source))
        mt.write(get_penetrance_path(data_source, "vep"), overwrite=args.overwrite)


    if not args.skip_filter_csq:
        logger.info(f'Filtering to variants with a consequence on any of the genes of interest...')
        mt = hl.read_matrix_table(get_penetrance_path(data_source, "vep"))
        mt = mt.annotate_globals(gene_set=gene_set)
        mt = mt.annotate_rows(vep_transcript_consequences=mt.vep.transcript_consequences.filter(
            lambda t: mt.gene_set.contains(t.gene_symbol)))
        mt = mt.filter_rows(mt.vep_transcript_consequences.length() >= 1, keep=True)
        mt = mt.drop('gene_set')
        mt.write(get_penetrance_path(data_source, "filtered"), overwrite=args.overwrite)


    if not args.skip_export_filtered:
        logger.info(f'Exporting list of all variants with consequence on genes of interest...')
        mt = hl.read_matrix_table(get_penetrance_path(data_source, "filtered"))
        ht = mt.rows()
        ht = ht.annotate(variant=hl.str(ht.locus) + ':' + hl.delimit(ht.alleles, ':'))
        ht = ht.select('variant', ht.variant_qc.n_called, ht.variant_qc.n_not_called, ht.variant_qc.AC,
                       ht.variant_qc.AF, ht.variant_qc.n_non_ref, ht.variant_qc.n_het)
        ht.export(get_path_to_filtered_variant_file(data_source))


    if not args.skip_filter_to_lof:
        logger.info(f'Filtering variants to those with a high-confidence loss-of-function annotation...')
        mt = hl.read_matrix_table(get_penetrance_path(data_source, "filtered"))
        mt = mt.annotate_globals(gene_set=lof_gene_set)
        mt = mt.annotate_rows(vep_transcript_consequences_HC=mt.vep_transcript_consequences.filter(
            lambda t: ((t.lof == 'HC') & mt.gene_set.contains(t.gene_symbol))))
        mt = mt.filter_rows(mt.vep_transcript_consequences_HC.length() >= 1, keep=True)
        mt.write(get_penetrance_path(data_source, "lof"), overwrite=args.overwrite)


    if not args.skip_tx_annotation:
        logger.info(f'Adding transcript annotation...')
        mt, gtex = read_tx_annotation_tables(get_penetrance_path(data_source, "lof"),
                                             get_tx_annotation_path(data_source), "mt")
        mt = tx_annotate_mt(mt, gtex, tissues_to_filter=get_tx_annotation_tissues_to_drop(data_source),
                            gene_maximums_kt_path=get_tx_annotation_kt_path(data_source),
                            tx_annotation_type="proportion", filter_to_csqs=all_coding_csqs,
                            filter_to_genes=lof_gene_set, gene_column_in_mt="gene_symbol")

        logger.info(f'Pulling out pext values for the worst consequence annotation...')
        ht_worst_csq = pull_out_worst_from_tx_annotate(mt.filter_rows(~hl.is_missing(mt.tx_annotation))).rows()
        ht_worst_csq.export(get_path_to_pext_worstcsq_tsv(data_source))
        print(ht_worst_csq.describe())
        ht_worst_csq = ht_worst_csq.annotate(worst_csq_pext=hl.delimit([t + ":" + hl.format('%.5e', ht_worst_csq[t]) for t in tissues], ","))
        print(ht_worst_csq.describe())
        print(ht_worst_csq.worst_csq_pext.show())

        mt = mt.annotate_rows(worst_csq_mean_proportion=ht_worst_csq[mt.row_key].mean_proportion,
                              worst_csq_pext=ht_worst_csq[mt.row_key].worst_csq_pext,
                              lof_flag=ht_worst_csq[mt.row_key].lof_flag)
        mt.write(get_penetrance_path(data_source, "lof_pext"), overwrite=args.overwrite)


    if not args.skip_liftover:
        if data_source == "52K":
            raise DataException("52K data on GRCh37 build so don't need liftover use --skip_liftover")
        logger.info('Lifting over data from GRCh38 to GRCh37...')
        mt = hl.read_matrix_table(get_penetrance_path(data_source, "lof_pext"))
        rg37 = hl.get_reference('GRCh37')
        rg38 = hl.get_reference('GRCh38')
        rg38.add_liftover('gs://hail-common/references/grch38_to_grch37.over.chain.gz', rg37)

        mt = mt.annotate_rows(new_locus=hl.liftover(mt.locus, 'GRCh37', include_strand=True))
        mt = mt.annotate_rows(liftover_37=hl.str(mt.new_locus.result) + ":" + hl.delimit(mt.alleles, ':'))
        mt = mt.write(get_penetrance_path(data_source, "lof_pext_liftover"), overwrite=args.overwrite)


    if not args.skip_export_lof_vcf:
        logger.info(f'Exporting LoF vcf...')
        vep_ht = hl.read_table(get_penetrance_path(data_source, "vep_csq"))
        if data_source == "52K":
            mt = hl.read_matrix_table(get_penetrance_path(data_source, "lof_pext"))
            mt = mt.annotate_rows(info=hl.struct(
                AC=mt.variant_qc.AC[1], AN=mt.variant_qc.AN, AF=mt.variant_qc.AF[1],
                n_hom=mt.variant_qc.homozygote_count[1], CSQ=vep_ht[mt.row_key].vep, worst_csq_pext=mt.worst_csq_pext,
                worst_csq_mean_proportion=mt.worst_csq_mean_proportion, lof_flag=mt.lof_flag))
        else:
            mt = hl.read_matrix_table(get_penetrance_path(data_source, "lof_pext_liftover"))
            mt = mt.annotate_rows(info=hl.struct(
                AC=mt.variant_qc.AC[1], AN=mt.variant_qc.AN, AF=mt.variant_qc.AF[1], n_called=mt.variant_qc.n_called,
                n_not_called=mt.variant_qc.n_not_called, n_non_ref=mt.variant_qc.n_non_ref,n_het=mt.variant_qc.n_het,
                n_hom=mt.variant_qc.homozygote_count[1],
                CSQ=vep_ht[mt.row_key].vep, liftover_37=mt.liftover_37, worst_csq_pext=mt.worst_csq_pext,
                worst_csq_mean_proportion=mt.worst_csq_mean_proportion, lof_flag=mt.lof_flag))
        hl.export_vcf(mt, get_path_to_lof_vcf(data_source), metadata={
            'info': {'CSQ': {'Description': vep_ht.vep_csq_header.collect()[0]}}
        })


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--overwrite', help='Overwrite all data from this subset (default: False)',
                        action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('-d', '--data_source', help='Data source', choices=['UKBB_regeneron', 'UKBB_gatk', '52K'],
                        required=True)

    parser.add_argument('--skip_import_data', help='Skip import plink data', action='store_true')
    parser.add_argument('--skip_filter_intervals', help='Skip filter to genes of interest', action='store_true')
    parser.add_argument('--skip_genotype_filter', help='Skip filtering by GQ, DP, and AB. Must skip for PLINK data sets',
                        action='store_true')
    parser.add_argument('--skip_variant_qc', help='Skip adding variant qc to the MT', action='store_true')
    parser.add_argument('--skip_vep', help='Skip running vep on intervals of interest', action='store_true')
    parser.add_argument('--skip_filter_csq', help='Skip filtering to variants with csq on genes', action='store_true')
    parser.add_argument('--skip_export_filtered', help='Skip exporting variants with csq on genes', action='store_true')
    parser.add_argument('--skip_filter_to_lof', help='Skip filtering to lof variants', action='store_true')
    parser.add_argument('--skip_tx_annotation', help='Skip transcript annotation', action='store_true')
    parser.add_argument('--skip_liftover', help='Skip adding liftover', action='store_true')
    parser.add_argument('--skip_export_lof_vcf', help='Skip exporting lof vcf', action='store_true')


    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)