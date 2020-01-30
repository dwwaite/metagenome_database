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
    parser_create = subparser.add_parser('create', help='Create the basic database file and populate with contigs')
    parser_create.add_argument('-d', '--db', help='Database name', required=True)

    # Add subparser
    parser_add = subparser.add_parser('add', help='Add entries to an existing database')
    parser_add.add_argument('-d', '--db', help='Database name', required=True)
    parser_add.add_argument('-c', '--contig', help='Assembled contig or scaffold file for basing the database', required=False)
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

    # Remove/delete subparser

    #args = parser.parse_args()
    args = parser.parse_args()

    ''' Process flow control '''
    if args.command == 'create':
        create_database(args.db)

    elif args.command == 'add':
        add_features(args.db,
                     contig_file=args.contig,
                     aa_file=args.prodigal_aa,
                     nt_file=args.prodigal_nt,
                     ssu_file=args.rrna_ssu,
                     lsu_file=args.rrna_lsu,
                     trna_file=args.trna,
                     trna_meta_file=args.trna_meta,
                     ann_file=args.annotations,
                     coverage_file=args.coverage,
                     transcript_file=args.transcript,
                     bin_file=args.bin_file,
                     bin_dir=args.bin_dir)

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

def add_features(db_name, contig_file=None, aa_file=None, nt_file=None, ssu_file=None, lsu_file=None, trna_file=None, trna_meta_file=None, ann_file=None, coverage_file=None, transcript_file=None, bin_file=None, bin_dir=None):

    ''' Open a connection then attempt to add each feature provided.
        The insert operations are ordered so if dependent features are added at the same time as their dependencies, there will be no error. '''
    try:

        db_connection = DatabaseManipulator(db_name)

        if contig_file:
            db_connection.add_contigs(contig_file)

        if aa_file or nt_file:
            db_connection.add_genes_prodigal(aa_file=aa_file, nt_file=nt_file)

        if ssu_file:
            db_connection.add_genes_metaxa(ssu_file, 'SSU')

        if lsu_file:
            db_connection.add_genes_metaxa(lsu_file, 'LSU')

        if trna_file and trna_meta_file:
            db_connection.create_aragorn_trna(trna_file, trna_meta_file)

        if ann_file:
            db_connection.add_annotations_by_table(ann_file)

        if coverage_file:
            db_connection.add_coverage_table(coverage_file)

        if transcript_file:
            db_connection.add_transcript_table(transcript_file)

        if bin_file and bin_dir:
            db_connection.add_bin_table(bin_file, bin_dir)

        db_connection.conn.commit()
        db_connection.close_connection()

    except Exception as e:

        print(e)

###############################################################################
if __name__ == '__main__':
    main()