import os
import re
from itertools import compress
import pandas as pd

class DatabasePreprocessFactory():

    @staticmethod
    def blast_to_df(blast_file, database_name, annot_method, output_name):

        blast_df = pd.read_csv(blast_file, sep='\t')
        blast_df.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']

        ''' Discard columns that will not be relevant to the database '''
        blast_df.drop(['mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'bitscore'], axis=1, inplace=True)

        ''' Update the columns of the BLAST table to match the expected format '''
        blast_df.rename( columns={ 'qseqid': 'gene', 'sseqid': 'hit', 'pident': 'identity', 'length': 'coverage' }, inplace=True)

        ''' Populate descriptive columns '''
        blast_df['database'] = database_name
        blast_df['method'] = annot_method
        blast_df['description'] = ''

        ''' Reorder columns and export '''
        blast_df = blast_df[ ['gene', 'method', 'database', 'hit', 'identity', 'coverage', 'evalue', 'description'] ]
        blast_df.to_csv('{}.blast_table.txt'.format(output_name), sep='\t', index=False)

    @staticmethod
    def _hmm_line_to_dict(hmm_line):

        line_entries = [ entry for entry in hmm_line.split(' ') if entry ]

        ''' Cast the filtered line entries to a dict. Value mappings are:
                target name => gene
                query name => description
                accession => hit
                E-value (domain) => evalue
        '''
        col_mask = [1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        return { k: v for k, v in zip(['gene', 'description', 'hit', 'evalue'], compress(line_entries, col_mask)) }

    @staticmethod
    def hmm_to_df(hmm_file, database_name, annot_method, output_name):

        hmm_content = [ line.strip() for line in open(hmm_file, 'r') if not re.search('^#', line) ]
        hmm_buffer = [ DatabasePreprocessFactory._hmm_line_to_dict(hmm_entry) for hmm_entry in hmm_content ]

        hmm_df = pd.DataFrame(hmm_buffer)

        ''' Populate empty columns '''
        hmm_df['database'] = database_name
        hmm_df['method'] = annot_method
        hmm_df['identity'] = -1
        hmm_df['coverage'] = -1

        ''' Reorder the columns and export '''
        hmm_df = hmm_df[ ['gene', 'method', 'database', 'hit', 'identity', 'coverage', 'evalue', 'description'] ]
        hmm_df.to_csv('{}.hmmer_table.txt'.format(output_name), sep='\t', index=False)

    @staticmethod
    def _append_taxonomies(checkm_df, gtdbtk_file):

        ''' Load in the GTDB-TK taxonomy file '''
        gtdbtk_df = pd.read_csv(gtdbtk_file, sep='\t')
        gtdbtk_dict = { u: c for u, c in zip(gtdbtk_df.user_genome, gtdbtk_df.classification) }

        ''' If any of the CheckM bins is in the file, update the taxonomy '''
        checkm_df.Taxonomy = [ gtdbtk_dict[b] if b in gtdbtk_dict else t for b, t in zip(checkm_df.Bin, checkm_df.Taxonomy) ]

    @staticmethod
    def gtdb_to_df(checkm_file, output_name, bact_file=None, arch_file=None):

        ''' Pull in the CheckM file, and reduce/rename down to the values of interest '''
        checkm_df = pd.read_csv(checkm_file, sep='\t')

        checkm_df.drop(['Marker lineage', '# genomes', '# markers', '# marker sets', '0', '1', '2', '3', '4', '5+'], axis=1, inplace=True)
        checkm_df.rename(columns={'Bin Id': 'Bin', 'Strain heterogeneity': 'Heterogeneity'}, inplace=True)

        ''' Add in a blank Taxonomy column, then populate it if possible '''
        checkm_df['Taxonomy'] = ''

        if bact_file:
            DatabasePreprocessFactory._append_taxonomies(checkm_df, bact_file)

        if arch_file:
            DatabasePreprocessFactory._append_taxonomies(checkm_df, arch_file)

        checkm_df.to_csv( '{}.bin_import.txt'.format(output_name), sep='\t', index=False )