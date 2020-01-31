'''
    Front end for metagenome databasing tool.

    Only used to parse user commands, then pass the result on to the appropriate pipeline.

    python metagenome_databaser.py create -d asd.db -c tests/contigs.fna
'''
import sys, argparse

from scripts.DatabaseManipulator import DatabaseManipulator

def main():

    ''' Set up the argument/subargument parsers '''
    parser = argparse.ArgumentParser( description='' )
    subparser = parser.add_subparsers(dest='command')

    # Creation subparser
    parser_create = subparser.add_parser('create', help='Create a new database file')
    parser_create.add_argument('-d', '--db', help='Database name', required=True)

    # Add subparser
    parser_add = subparser.add_parser('add', help='Add entries to an existing database')
    parser_add.add_argument('-d', '--db', help='Database name', required=True)
    parser_add.add_argument('--contig', help='Assembled contig or scaffold file for basing the database', required=False)
    parser_add.add_argument('--prodigal_aa', help='Genes (amino acids) predicted using prodigal', required=False)
    parser_add.add_argument('--prodigal_nt', help='Genes (nucleotide) predicted using prodigal', required=False)
    parser_add.add_argument('--rrna_ssu', help='Small subunit (16S, 18S) rRNA gene fasta file, predicted by MeTaxa2. Automatically creates annotations.', required=False)
    parser_add.add_argument('--rrna_lsu', help='Large subunit (23S, 28S) rRNA gene fasta file, predicted by MeTaxa2. Automatically creates annotations.', required=False)
    parser_add.add_argument('--trna', help='Transfer RNA (tRNA) gene fasta file, predicted by Aragorn. Automatically creates annotations.', required=False)
    parser_add.add_argument('--trna_meta', help='Metadata file to pair with Aragorn prediction file. This is required for adding tRNA records.', required=False)
    parser_add.add_argument('--annotations', help='Annotation table for genes present in the database', required=False)
    parser_add.add_argument('--coverage', help='Coverage table for contigs present in the database', required=False)
    parser_add.add_argument('--transcript', help='Transcription/expression table for genes present in the database', required=False)
    parser_add.add_argument('--bin_file', help='Table of bin descriptions to add to database', required=False)
    parser_add.add_argument('--bin_dir', help='Directory of bins to add to database', required=False)

    # Export subparser
    parser_export = subparser.add_parser('export', help='Retrieve entries from an existing database.')
    parser_export.add_argument('-d', '--db', help='Database name', required=True)
    parser_export.add_argument('-o', '--output', help='Target file name for export. Feature-specific suffixes will be appended to this value', required=True)
    parser_export.add_argument('--seqs_as_tables', help='Export contigs or genes as tables instead of fasta files (Default: False)', action='store_true', required=False)
    parser_export.add_argument('--contig', help='Export all contigs', action='store_true', required=False) #get_contigs
    parser_export.add_argument('--gene', help='Export all genes', action='store_true', required=False) #get_contigs
    parser_export.add_argument('--annotation', help='Export annotations associated with a set of genes. Requires a text file of gene names', required=False) #get_contigs
    parser_export.add_argument('--coverage', help='Export all contig coverage values in the database', action='store_true', required=False) #get_contigs
    parser_export.add_argument('--coverage_by_sample', help='Export all contig coverage values found in a set of samples. Requires a text file of sample names', required=False) #get_contigs
    parser_export.add_argument('--transcript', help='Export all gene expression values in the database', action='store_true', required=False) #get_contigs
    parser_export.add_argument('--transcript_by_sample', help='Export all gene expression values found in a set of samples. Requires a text file of sample names', required=False) #get_contigs
    parser_export.add_argument('--bin', help='Export contigs, genes, and annotations associated with a bin. Requires a text file of bin names', required=False) #get_contigs

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

    ''' Process flow control '''
    if args.command == 'create':
        create_database(args.db)

    elif args.command == 'add':
        add_features(args.db, vars(args) )

    elif args.command == 'export':
        pass

###############################################################################
#
# Parse the user parameters into commands to the DatabaseManipulator
#

def create_database(db_name):

        ''' Open a connection and create the blank tables '''
        db_connection = DatabaseManipulator(db_name)

        status, msg = db_connection.create_blank_database()

        if status != 1:
            print(msg)

        db_connection.close_connection()

def add_features(db_name, add_params):
    
    ''' Open a connection then attempt to add each feature provided.
        The insert operations are ordered so if dependent features are added at the same time as their dependencies, there will be no error. '''
    try:

        db_connection = DatabaseManipulator(db_name)

        if add_params['contig']:
            db_connection.add_contigs(add_params['contig'])

        if add_params['prodigal_aa'] or add_params['prodigal_nt']: db_connection.add_genes_prodigal(aa_file=add_params['prodigal_aa'], nt_file=add_params['prodigal_nt'])

        if add_params['rrna_ssu']: db_connection.add_genes_metaxa(wargs['rrna_ssu'], 'SSU')

        if add_params['rrna_lsu']: db_connection.add_genes_metaxa(add_params['rrna_lsu'], 'LSU')

        if add_params['trna'] and add_params['trna_meta']:
            db_connection.create_aragorn_trna(add_params['trna'], add_params['trna_meta'])
        elif add_params['trna'] or add_params['trna_meta']:
            print('Missing required input parameters for tRNA annotation, please check input options.')

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

        print(e)

'''
--seqs_as_tables', help='Export contigs or genes as tables instead of fasta files (Default: False)', action='store_true', required=False)
    parser_export.add_argument('--contig', help='Export all contigs', action='store_true', required=False) #get_contigs
    parser_export.add_argument('--gene', help='Export all genes', action='store_true', required=False) #get_contigs
    parser_export.add_argument('--annotation', help='Export annotations associated with a set of genes. Requires a text file of gene names', required=False) #get_contigs
    parser_export.add_argument('--coverage', help='Export all contig coverage values in the database', action='store_true', required=False) #get_contigs
    parser_export.add_argument('--coverage_by_sample', help='Export all contig coverage values found in a set of samples. Requires a text file of sample names', required=False) #get_contigs
    parser_export.add_argument('--transcript', help='Export all gene expression values in the database', action='store_true', required=False) #get_contigs
    parser_export.add_argument('--transcript_by_sample', help='Export all gene expression values found in a set of samples. Requires a text file of sample names', required=False) #get_contigs
    parser_export.add_argument('--bin', help='Expo
'''
def export_data(db_name, output_prefix, export_contigs=None, export_genes=None, export_annotations=None, export_coverage=None, export_cov_by_sample=None, export_transcript=None, export_transcript_by_sample=None):


    return

###############################################################################
if __name__ == '__main__':
    main()