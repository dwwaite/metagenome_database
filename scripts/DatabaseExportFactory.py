import os
import pandas as pd
import numpy as np

class DatabaseExportFactory():

    #region Internal functions

    @staticmethod
    def _as_df(df, output_name, identifier):
        df.to_csv('{}.{}.txt'.format(output_name, identifier), sep='    ', index=False)

    @staticmethod
    def _as_fasta(df, name_column, seq_column, output_name, ext):

        output_writer = open('{}.{}'.format(output_name, ext), 'w')
        for _, row in df.iterrows():

            if row[seq_column] != '':
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

    @staticmethod
    def summarise_database(db_name, dbm):

        ''' Report statistics for database '''
        print('')
        print('')
        print( 'Generating report for database: {}'.format(db_name) )
        user, version = dbm.read_version
        print( '    Built under software version {}'.format(version) )
        print( '    Built by user {}'.format(user) )

        ''' Start with contigs '''
        contig_df = dbm.get_contigs()

        print('')
        print('--== CONTIGS ==--')
        print( 'The database contains a total of {} contigs'.format( contig_df.shape[0] ) )
        print( '    Smallest contig length: {}'.format( min(contig_df.length) ) )
        print( '    Median contig length: {}'.format( np.median(contig_df.length) ) )
        print( '    Longest contig length: {}'.format( max(contig_df.length) ) )

        ''' Move to genes '''
        gene_df = dbm.get_genes()
        print('')
        print('--== GENES ==--')
        print( 'The database contains a total of {} genes, using {} methods.'.format( gene_df.shape[0],
                                                                                                len( gene_df.prediction_tool.unique() ) ) )

        for pt, df in gene_df.groupby('prediction_tool'):
            print( 'Genes predicted by {}: {}'.format(pt, df.shape[0]) )

        ''' Now annotations '''
        annotation_df = dbm.get_annotations(old_annotations=True)
        print('')
        print('--== ANNOTATIONS ==--')
        print( 'The database contains a total of {} annotations, using {} methods.'.format( annotation_df.shape[0],
                                                                                                len( annotation_df.method.unique() ) ) )
        print( 'There are currently {} annotations marked as current.'.format( annotation_df.iscurrent.sum() ) )

        for m, df in annotation_df.groupby('method'):

            if m in set(['MeTaxa2', 'Aragorn']):
                print( '    Genes annotated using {}: {} ({} total)'.format(m, df.iscurrent.sum(), df.shape[0]) )

            else:
                for db, _ in df.groupby('annotation_db'):
                    print( '    Genes annotated using {} ({}): {} ({} total)'.format(m, db, _.iscurrent.sum(), _.shape[0]) )

        ''' Coverage... '''
        coverage_df = dbm.get_coverage()
        print('')
        print('--== COVERAGE ==--')
        print( 'A total of {} reads have been mapped to {} contigs, across {} samples.'.format(coverage_df.reads_mapped.sum(),
                                                                                              len( coverage_df.contig_name.unique() ),
                                                                                              len( coverage_df['sample'].unique() ) ) )

        ''' Transcription... '''
        transcript_df = dbm.get_transcript()
        print('')
        print('--== TRANSCRIPTION ==--')
        print( 'A total of {} reads have been mapped to {} genes, across {} samples.'.format(transcript_df.reads_mapped.sum(),
                                                                                            len( transcript_df.gene_name.unique() ),
                                                                                            len( transcript_df['sample'].unique() ) ) )
        
        ''' Binning... '''
        all_bins_df = dbm.get_all_bins()
        current_bins_df = dbm.get_current_bins()

        print('')
        print('--== BINNING ==--')
        print( 'A total of {} bins have been stored in this database, with {} marked as current.'.format(all_bins_df.shape[0], current_bins_df.shape[0]) )


        n_submissions = len( all_bins_df.version_name.unique() )
        print('')
        print( 'A total of {} bin submission{} {} been made to this database.'.format(n_submissions,
                                                                                      '' if n_submissions == 1 else 's',
                                                                                      'has' if n_submissions == 1 else 'have' ) )

        print( 'In the current bin set ({})...'.format( current_bins_df.version_name[0] ) )
        print( '    The average completeness is {}% (current), or {}% across all bins.'.format( np.median(current_bins_df.completeness),
                                                                                              np.median(all_bins_df.completeness) ) )
        print( '    The average contamination is {}% (current), or {}% across all bins'.format( np.median(current_bins_df.contamination),
                                                                                            np.median(all_bins_df.contamination) ) )

        print('')
        print('Accross every binning version:')
        for v, df in all_bins_df.groupby('version_name'):
            print( '    Submission \'{}\': {} bins, average completeness is {}%.'.format(v, df.shape[0], np.median(df.completeness)) )

    #endregion