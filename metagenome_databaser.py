'''
    Front end for metagenome databasing tool.
    Only used to parse user commands, then pass the result on to the appropriate pipeline.
'''
import sys, argparse

from scripts.DatabaseManipulator import DatabaseManipulator
from scripts.DatabaseExportFactory import DatabaseExportFactory as export_factory

def main():

    '''
        Define the version of the database software. This is passed into the database connection function.
            If creation, the database is stamped with the software version.
            If loading, the database is tested to ensure it is compatible with the current software version.
    '''
    __version_info__ = ('0', '1', '0')
    version_str = '.'.join(__version_info__)

    ''' Set up the argument/subargument parsers '''
    parser = argparse.ArgumentParser( description='' )
    parser.add_argument('-v', '--version', action='version', version='%(prog)s version {}'.format( version_str ) )

    ''' Begin loading sub-commands '''
    subparser = parser.add_subparsers(dest='command')
    populate_create(subparser)
    populate_add(subparser)
    populate_export(subparser)
    populate_view(subparser)

    # Remove/delete subparser
    #parser_remove = subparser.add_parser('remove', help='Remove entries from an existing database')
    #parser_remove.add_argument('-d', '--db', help='Database name', required=True)
    #parser_remove.add_argument('--contig', help='Contigs or scaffolds to remove. If removed, will also remove any associated genes, coverage values, or bin associations', required=False)
    #parser_remove.add_argument('--gene', help='Genes to remove. If removed, will also remove any associated annotations or transcription records', required=False)
    #parser_remove.add_argument('--annotation', help='Annotations to remove', required=False)
    #parser_remove.add_argument('--coverage', help='Coverage records to remove', required=False)
    #parser_remove.add_argument('--transcript', help='Transcription records to remove', required=False)
    #parser_remove.add_argument('--bin', help='Bin associations to remove. Removing bins does not affect the underlying contigs/scaffolds', required=False)

    args = parser.parse_args()

    ''' Process flow control
        For options other than create, ensure software and database are compatible before proceeding '''    
    if args.command == 'create':
        create_database(args.db, version_str)

    elif args.command == 'add':
        verify_versions(args.db, version_str)
        add_features(args.db, version_str, vars(args) )

    elif args.command == 'export':
        verify_versions(args.db, version_str)
        export_data(args.db, version_str, args.output, vars(args) )

    elif args.command == 'view':
        verify_versions(args.db, version_str)
        summarise_database(args.db, version_str)

###############################################################################
#region Subparser
#   Functions to tidy away the creation of subparsers for the user interface

def populate_create(subparser):
    parser_create = subparser.add_parser('create', help='Create a new database file')
    parser_create.add_argument('-d', '--db', help='Database name', required=True)

def populate_add(subparser):
    parser_add = subparser.add_parser('add', help='Add entries to an existing database')
    parser_add.add_argument('-d', '--db', help='Database name', required=True)
    parser_add.add_argument('--contig', help='Assembled contig or scaffold file for basing the database', required=False)
    parser_add.add_argument('--prodigal_aa', help='Genes (amino acids) predicted using prodigal', required=False)
    parser_add.add_argument('--prodigal_nt', help='Genes (nucleotide) predicted using prodigal', required=False)
    parser_add.add_argument('--rrna_ssu', help='Small subunit (16S, 18S) rRNA gene fasta file, predicted by MeTaxa2. Automatically creates annotations.', required=False)
    parser_add.add_argument('--rrna_lsu', help='Large subunit (23S, 28S) rRNA gene fasta file, predicted by MeTaxa2. Automatically creates annotations.', required=False)
    parser_add.add_argument('--trna', help='Transfer RNA (tRNA) gene fasta file, predicted by Aragorn. Automatically creates annotations.', required=False)
    parser_add.add_argument('--trna_contigs', help='The contigs file used in Aragorn prediction. This is required for adding tRNA records.', required=False)
    parser_add.add_argument('--annotations', help='Annotation table for genes present in the database', required=False)
    parser_add.add_argument('--coverage', help='Coverage table for contigs present in the database', required=False)
    parser_add.add_argument('--transcript', help='Transcription/expression table for genes present in the database', required=False)
    parser_add.add_argument('--bin_file', help='Table of bin descriptions to add to database', required=False)
    parser_add.add_argument('--bin_dir', help='Directory of bins to add to database', required=False)

def populate_export(subparser):
    parser_export = subparser.add_parser('export', help='Retrieve entries from an existing database.')
    parser_export.add_argument('-d', '--db', help='Database name', required=True)
    parser_export.add_argument('-o', '--output', help='Target file handle for export. Feature-specific suffixes will be appended to this value', required=True)
    parser_export.add_argument('--seqs_as_tables', help='Export contigs or genes as tables instead of fasta files (Default: False)', action='store_true', required=False)
    parser_export.add_argument('--all_records', help='Export all annotation records. (Default: Current version only)', action='store_true', required=False)
    parser_export.add_argument('--contig', help='Export all contigs', action='store_true', required=False)
    parser_export.add_argument('--gene', help='Export all genes', action='store_true', required=False)
    parser_export.add_argument('--gene_by_method', help='Export all genes predicted using a particular tool', required=False)
    parser_export.add_argument('--annotation', help='Export annotations associated with a set of genes', action='store_true', required=False)
    parser_export.add_argument('--coverage', help='Export all contig coverage values in the database', action='store_true', required=False)
    parser_export.add_argument('--transcript', help='Export all gene expression values in the database', action='store_true', required=False)
    parser_export.add_argument('--bin', help='Export contigs, genes, and current annotations associated with the current bin set', action='store_true', required=False)
    parser_export.add_argument('--bin_all', help='Export contigs, genes, and current annotations associated with all bins', action='store_true', required=False)
    parser_export.add_argument('--bin_by_version', help='Export contigs, genes, and current annotations associated with bins of a user-specifed binning version', required=False)

def populate_view(subparser):
    parser_view = subparser.add_parser('view', help='Report the current contents of the database to stdout')
    parser_view.add_argument('-d', '--db', help='Database name', required=True)

#endregion
###############################################################################
#region User input flow control
#   Parse the user parameters into commands to the DatabaseManipulator
def create_database(db_name, db_version):

        ''' Open a connection and create the blank tables '''
        db_connection = DatabaseManipulator(db_name, db_version)
        status, msg = db_connection.create_blank_database()

        if status != 1:
            print(msg)

def verify_versions(db_name, db_version):
    db_connection = DatabaseManipulator(db_name, db_version)
    db_connection.validate_database_version()

def add_features(db_name, db_version, add_params):
    
    ''' Open a connection then attempt to add each feature provided.
        The insert operations are ordered so if dependent features are added at the same time as their dependencies, there will be no error. '''
    try:

        db_connection = DatabaseManipulator(db_name, db_version)
        db_connection.open_connection()

        if add_params['contig']:
            db_connection.add_contigs(add_params['contig'])

        if add_params['prodigal_aa'] or add_params['prodigal_nt']: db_connection.add_genes_prodigal(aa_file=add_params['prodigal_aa'], nt_file=add_params['prodigal_nt'])

        if add_params['rrna_ssu']: db_connection.add_genes_metaxa(add_params['rrna_ssu'], 'SSU')

        if add_params['rrna_lsu']: db_connection.add_genes_metaxa(add_params['rrna_lsu'], 'LSU')

        if add_params['trna'] and add_params['trna_contigs']:
            db_connection.add_genes_aragorn(add_params['trna_contigs'], add_params['trna'])
        elif add_params['trna'] or add_params['trna_contigs']:
            print('Missing contig file required for tRNA annotation, please check input options.')

        if add_params['annotations']: db_connection.add_annotations_by_table(add_params['annotations'])

        if add_params['coverage']: db_connection.add_coverage_table(add_params['coverage'])

        if add_params['transcript']: db_connection.add_transcript_table(add_params['transcript'])

        if add_params['bin_file'] and add_params['bin_dir']:
            db_connection.add_bin_table(add_params['bin_file'], add_params['bin_dir'])
        elif add_params['bin_file'] or add_params['bin_dir']:
            print('Missing required input parameters for bins, please check input options.')

        db_connection.conn.commit()
        db_connection.close_connection()

    except Exception as e:
        db_connection.conn.rollback()
        print( 'Error encountered: {}'.format(e))

def export_data(db_name, db_version, output_prefix, export_params):

    try:

        db_connection = DatabaseManipulator(db_name, db_version)
        db_connection.open_connection()

        if export_params['contig']:
            contig_df = db_connection.get_contigs()
            export_factory.contig_export(contig_df, output_prefix, export_params['seqs_as_tables'])

        if export_params['gene']:
            gene_df = db_connection.get_genes()
            export_factory.gene_export(gene_df, output_prefix, export_params['seqs_as_tables'])

        if export_params['gene_by_method']:
            gene_df = db_connection.get_genes(prediction_method=export_params['gene_by_method'])
            export_factory.gene_export(gene_df, '{}.{}'.format(output_prefix, export_params['gene_by_method'].lower()), export_params['seqs_as_tables'])

        if export_params['annotation']:
            annotation_df = db_connection.get_annotations(old_annotations=True) if export_params['all_records'] else db_connection.get_annotations()
            export_factory.annotation_export(annotation_df, output_prefix)

        if export_params['coverage']:
            coverage_df = db_connection.get_coverage()
            export_factory.coverage_export(coverage_df, output_prefix)

        if export_params['transcript']:
            transcript_df = db_connection.get_transcript()
            export_factory.transcript_export(transcript_df, output_prefix)

        if export_params['bin']:
            bin_df = db_connection.get_current_bins()
            export_factory.bin_export(bin_df, db_connection, output_prefix)

        if export_params['bin_all']:
            bin_df = db_connection.get_all_bins()
            export_factory.bin_export(bin_df, db_connection, output_prefix, record_version='all')

        if export_params['bin_by_version']:
            bin_df = db_connection.get_bins_by_version(export_params['bin_by_version'])
            export_factory.bin_export(bin_df, db_connection, output_prefix, record_version=export_params['bin_by_version'])

        db_connection.close_connection()

    except Exception as e:
        print( 'Error encountered: {}'.format(e))

def summarise_database(db_name, db_version):

    try:
        db_connection = DatabaseManipulator(db_name, db_version)
        db_connection.open_connection()

        export_factory.summarise_database(db_name, db_connection)

        db_connection.close_connection()

    except Exception as e:
        print( 'Error encountered: {}'.format(e))

#endregion
###############################################################################
if __name__ == '__main__':
    main()