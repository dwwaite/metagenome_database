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
    parser_create.add_argument('-c', '--contig', help='Assembled contig or scaffold file for basing the database', required=True)

    # Add subparser
    parser_add = subparser.add_parser('add', help='Add entries to an existing database')
    parser_add.add_argument('-d', '--db', help='Database name', required=True)
    parser_add.add_argument('--prodigal_aa', help='Genes (amino acids) predicted using prodigal', required=False)
    parser_add.add_argument('--prodigal_nt', help='Genes (nucleotide) predicted using prodigal', required=False)
    parser_add.add_argument('--rrna_ssu', help='Small subunit (16S, 18S) rRNA gene fasta file, predicted by MeTaxa2', required=False)
    parser_add.add_argument('--rrna_lsu', help='Large subunit (23S, 28S) rRNA gene fasta file, predicted by MeTaxa2', required=False)
    parser_add.add_argument('--trna', help='Transfer RNA (tRNA) gene fasta file, predicted by Aragorn', required=False)

    # Remove/delete subparser

    #args = parser.parse_args()
    args = parser.parse_args()

    ''' Process flow control '''
    if args.command == 'create':
        create_database(args.db, args)

###############################################################################
#
# Parse the user parameters into commands to the DatabaseManipulator
#

def create_database(db_name, arguments):

        ''' Open a connection and create the blank tables '''
        db_connection = DatabaseManipulator(db_name)

        status, msg = db_connection.create_blank_database()

        if status < 0:
            print(msg)
            sys.exit()

        ''' No errors reported, load in the contigs '''
        db_connection.add_contigs(arguments.contig)

def add_features():
    pass
    """
        ''' With the base created, add in any additional features provided '''
        if arguments.prodigal_aa or arguments.prodigal_nt:
            db_connection.add_genes_prodigal(aa_file=arguments.prodigal_aa, nt_file=arguments.prodigal_nt)

        if arguments.rrna_ssu:
            db_connection.add_genes_metaxa(arguments.rrna_ssu, 'SSU')

        if arguments.rrna_lsu:
            db_connection.add_genes_metaxa(arguments.rrna_lsu, 'LSU')

        if arguments.trna:
            pass
    """

###############################################################################
if __name__ == '__main__':
    main()