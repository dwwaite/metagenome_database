import os
import pandas as pd

class DatabaseExportFactory():

    #region Internal functions

    @staticmethod
    def _as_df(df, output_name, identifier):
        df.to_csv('{}.{}.txt'.format(output_name, identifier), sep='\t', index=False)

    @staticmethod
    def _as_fasta(df, name_column, seq_column, output_name, ext):

        output_writer = open('{}.{}'.format(output_name, ext), 'w')
        for _, row in df.iterrows():
            _ = output_writer.write( '>{}\n{}\n'.format(row[name_column], row[seq_column]) )

        output_writer.close()

    #endregion

    #region Exposed functions

    @staticmethod
    def contig_export(contigs_df, output_name, seqs_as_tables=None):
        if seqs_as_tables:
            DatabaseExportFactory._as_df(contigs_df, output_name, 'contigs')
        else:
            DatabaseExportFactory._as_fasta(contigs_df, 'contig_name', 'sequence', output_name, 'contigs.fna')

    @staticmethod
    def gene_export(genes_df, output_name, seqs_as_tables=None):
        if seqs_as_tables:
            DatabaseExportFactory._as_df(genes_df, output_name, 'genes')
        else:
            DatabaseExportFactory._as_fasta(genes_df, 'gene_name', 'sequence_aa', output_name, 'genes.faa')
            DatabaseExportFactory._as_fasta(genes_df, 'gene_name', 'sequence_nt', output_name, 'genes.fna')

    @staticmethod
    def annotation_export(annotations_df, output_name):
        DatabaseExportFactory._as_df(annotations_df, output_name, 'annotations')

    @staticmethod
    def coverage_export(coverage_df, output_name):
        DatabaseExportFactory._as_df(coverage_df, output_name, 'coverage')

    @staticmethod
    def transcript_export(transcript_df, output_name):
        DatabaseExportFactory._as_df(transcript_df, output_name, 'transcript')

    @staticmethod
    def bin_export(bin_df, dbm, output_name, record_version=None):

        ''' Attempt to create the directory for export. Return if it already exists. '''
        output_dir = '{}.bin_export.{}'.format(output_name, record_version) if record_version else '{}.bin_export'.format(output_name)

        if os.path.isdir(output_dir):
            return

        os.mkdir(output_dir)
        output_prefix = os.path.join(output_dir, output_name)

        ''' Write the metadata '''
        DatabaseExportFactory._as_df(bin_df, output_prefix, 'bin_metadata')

        ''' Read content from the database, ready for exporting '''
        master_contig_df = dbm.get_contigs()
        master_gene_df = dbm.get_genes()
        master_annotation_df = dbm.get_annotations()

        ''' For each bin, export data '''
        for _, row in bin_df.iterrows():

            ''' Contigs '''
            contig_list = dbm.get_contigs_in_bin( row['bin_name'] )
            bin_contig_df = master_contig_df[ master_contig_df.contig_name.isin(contig_list) ]
            DatabaseExportFactory._as_fasta(bin_contig_df, 'contig_name', 'sequence', output_prefix, '{}.fna'.format(row['bin_name']) )

            ''' Genes '''
            bin_gene_df = master_gene_df[ master_gene_df.contig_name.isin(contig_list) ]
            DatabaseExportFactory._as_fasta(bin_gene_df, 'gene_name', 'sequence_aa', output_prefix, '{}.genes.faa'.format(row['bin_name']) )
            DatabaseExportFactory._as_fasta(bin_gene_df, 'gene_name', 'sequence_nt', output_prefix, '{}.genes.fna'.format(row['bin_name']) )

            ''' Current annotations '''
            bin_annotation_df = master_annotation_df[ master_annotation_df.gene_name.isin( bin_gene_df.gene_name ) ]
            DatabaseExportFactory._as_df(bin_annotation_df, output_prefix, '{}.annotations'.format(row['bin_name']) )

    #endregion