'''
    Automated testing of the creation, insertion, update, and delete functionality of the DatabaseManipulator class

    Logic flows in the order of dependencies, rather than organising by database operations

    For example:
        Stage 0 - Basic data handling functions
        Stage 1 - Database creation and establishing/closing a connection
        Stage 2 - Creation and manipulation of contigs
        Stage 3 - Creation and manipulation of genes, which require contigs as FK
        Stage 4 - Creation and manipulation of annotations, which require genes as FK
        ...and so on.
'''
import sys, os, io, time
import unittest
import sqlite3
import pandas as pd
from datetime import datetime

class TestDatabaseManipulator(unittest.TestCase):

    def setUp(self):

        ''' Create an in-memory database for testing purposes '''
        __version_info__ = ('0', '1', '0')

        self.db_m = DatabaseManipulator(':memory:', __version_info__)
        self.db_m.create_blank_database()

        #cursor = self.db_m.conn.cursor()
        #cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
        #print(cursor.fetchall())

        self.temp_file_buffer = []
        self.temp_folder_buffer = []

    def tearDown(self):

        for temp_file in self.temp_file_buffer:
            os.remove(temp_file)

        for temp_folder in self.temp_folder_buffer:
            os.rmdir(temp_folder)

    # region File manipulation for test purposes

    def spawn_persistant_database(self, db_version):

        db_name = 'mock.testing.hdb'
        self.db_m = DatabaseManipulator(db_name, db_version)

        self.db_m.create_blank_database()
        self.db_m.close_connection()

        self.temp_file_buffer.append(db_name)
        return db_name

    def make_directory(self, d_name):
        os.mkdir(d_name)
        self.temp_folder_buffer.append(d_name)

    def push_to_fasta(self, f_name, seq_names, seq_values, seq_meta=None):

        mock_fna = open(f_name, 'w')
        self.temp_file_buffer.append(f_name)

        if seq_meta:
            for k, m, s in zip(seq_names, seq_meta, seq_values):

                mock_fna.write( '>{} {}\n{}\n'.format(k, m, s) )

        else:
            for k, s in zip(seq_names, seq_values):

                            mock_fna.write( '>{}\n{}\n'.format(k, s) )

        mock_fna.close()

    def append_to_fasta(self, f_name, seq_name, seq_value):

        writer = open(f_name, 'a')
        writer.write( '>{}\n{}\n'.format(seq_name, seq_value) )
        writer.close()

    def push_to_table(self, f_name, df):

        df.to_csv(f_name, index=False, sep='\t')

        self.temp_file_buffer.append(f_name)

    # endregion

    # region Stage 0 - Row casting functions

    def test_gene_row_to_dict(self):

        mock_row = ('abc_1', 'abc', 'ATCGAT', 'ID', 6, 2, 1, 6, 1)
        gene_dict = self.db_m._gene_row_to_dict(mock_row)

        for exp_key, val in zip(['gene_name', 'contig_name', 'prediction_tool', 'sequence_nt', 'sequence_aa', 'length_nt', 'length_aa', 'start_pos', 'stop_pos', 'orientation'], mock_row):

            self.assertTrue( exp_key in gene_dict )
            self.assertTrue( gene_dict[exp_key] == val )

    def test_contig_row_to_dict(self):

        mock_row = ('abc', 'ATCGATCG', 8, 0.5)
        contig_dict = self.db_m._contig_row_to_dict(mock_row)

        for exp_key, val in zip(['contig_name', 'sequence', 'length', 'gc'], mock_row):

            self.assertTrue( exp_key in contig_dict )
            self.assertTrue( contig_dict[exp_key] == val )

    def test_annotation_row_to_dict(self):

        mock_row = ('abc', 'blast', 'nothing', '12345', 70.5, 30.2, 0.001, 'fake annotation', True)
        annot_dict = self.db_m._annotation_row_to_dict(mock_row)

        for exp_key, val in zip(['gene_name', 'method', 'annotation_db', 'accession', 'identity', 'coverage', 'evalue', 'description', 'iscurrent'], mock_row):

            self.assertTrue( exp_key in annot_dict )
            self.assertTrue( annot_dict[exp_key] == val )

    def test_coverage_row_to_dict(self):

        mock_row = ('abs', 'sample1', 123)
        coverage_dict = self.db_m._coverage_row_to_dict(mock_row)

        for exp_key, val in zip(['contig_name', 'sample', 'reads_mapped'], mock_row):

            self.assertTrue( exp_key in coverage_dict )
            self.assertTrue( coverage_dict[exp_key] == val )

    def test_transcript_row_to_dict(self):

        mock_row = ('abs', 'sample1', 123)
        coverage_dict = self.db_m._transcript_row_to_dict(mock_row)

        for exp_key, val in zip(['gene_name', 'sample', 'reads_mapped'], mock_row):

            self.assertTrue( exp_key in coverage_dict )
            self.assertTrue( coverage_dict[exp_key] == val )

    def test_bin_row_to_dict(self):

        mock_row = ( 'bin1', 'a;b;c', 1.00, 2.00, 3.00, 'v1', datetime.now().isoformat() )
        bin_dict = self.db_m._bin_row_to_dict(mock_row)

        for exp_key, val in zip(['bin_name', 'taxonomy', 'completeness', 'contamination', 'strainheterogeneity', 'version_name', 'date_created'], mock_row):

            self.assertTrue( exp_key in bin_dict )
            self.assertTrue( bin_dict[exp_key] == val )

    # endregion

    # region Stage 1 - Create and evaluate database

    def test_open_connection(self):

        ''' There does not appear to be a way to actually test a connection is open.
            Most commonly recommended way is to simply execute a query, and in the absence of error, assume the connection is valid. '''
        self.db_m.open_connection()
        cursor = self.db_m.conn.cursor()

        try:

            cursor.execute(''' SELECT count(name) FROM sqlite_master WHERE type='table' AND name='contig' ''')
            self.assertEqual(1, 1)

        except:
            self.assertEqual(0, 1)

    def test_db_creation(self):

        version_created = ('0', '1', '0')
        _ = self.spawn_persistant_database(version_created)

        self.db_m.open_connection()

        cursor = self.db_m.conn.cursor()
        cursor.execute(''' SELECT count(name) FROM sqlite_master WHERE type='table' ''' )

        for tbl in ['tbl_version', 'contig', 'gene', 'annotation', 'coverage', 'transcript', 'bin']:
            cursor.execute(''' SELECT count(name) FROM sqlite_master WHERE type='table' AND name=? ''', (tbl,) )
            self.assertEqual(cursor.fetchone()[0], 1 )

        self.db_m.close_connection()

    def test_db_creation_alreadyexists(self):

        version_created = ('0', '1', '0')
        database_name = self.spawn_persistant_database(version_created)

        ''' Test '''
        second_db = DatabaseManipulator(database_name, version_created)
        _, err = second_db.create_blank_database()

        ''' Evaluate '''
        self.assertIn( 'Database file {} already exists. Did you mean to invoke the update command?'.format(database_name), err )

    def test_database_version_error(self):

        version_created = ('0', '1', '0')
        version_attached = ('0', '2', '0')

        database_name = self.spawn_persistant_database(version_created)
        db_m = DatabaseManipulator(database_name, version_attached)

        ''' Evaluate '''
        with self.assertRaises(Exception) as cm:
            db_m.validate_database_version()

        self.assertIn( 'Software version is {}, but database was created under version {}. Unable to proceed.'.format('.'.join(version_attached),
                                                                                                                      '.'.join(version_created)), str(cm.exception) )

    @unittest.skip('Database structure is not explicitly tested in working code yet')
    def test_database_exists_structure_error(self):

        ''' Open database connection, create cursor, drop table, test creation '''
        self.db_m.create_blank_database()
        self.db_m.open_connection()
        cursor = self.db_m.conn.cursor()

        cursor.execute( "DROP TABLE bin" )

        self.assertEqual( self.db_m.database_exists(), -1 )

    def test_close_connection(self):

        self.db_m.open_connection()
        cursor = self.db_m.conn.cursor()

        self.db_m.close_connection()

        try:

            cursor.execute(''' SELECT count(name) FROM sqlite_master WHERE type='table' AND name='contig' ''')
            self.assertEqual(0, 1, 'Close connection to database')

        except:

            self.assertEqual(1, 1, 'Close connection to database')

    # endregion

    # region Stage 2 - Contig functions

    def create_contigs_file(self, f_name):

        exp_keys = ['contig1', 'contig2', 'contig3']
        exp_seq = ['AAATTTCCCGGG', 'AAATTTCCC', 'AAATTT']

        self.push_to_fasta(f_name, exp_keys, exp_seq)

        return exp_keys, exp_seq

    def test_obtain_contig_set(self):

        ''' Set up files for test '''
        file_name = 'mock.contigs.fasta'
        exp_contigs, exp_seqs = self.create_contigs_file(file_name)

        ''' Test '''
        self.db_m.create_blank_database()
        self.db_m.add_contigs(file_name)

        contig_set = self.db_m._obtain_contig_set()
        self.assertEqual( len( set(exp_contigs) & contig_set), len(exp_contigs) )

    def test_validate_contigs(self):

        ''' Set up files for test '''
        file_name = 'mock.contigs.fasta'
        exp_contigs, exp_seqs = self.create_contigs_file(file_name)

        exp_results = { c: Contig(c, s) for c, s in zip(exp_contigs, exp_seqs) }

        self.db_m.create_blank_database()

        ''' Evaluate '''
        result_buffer = self.db_m._validate_contigs(file_name)

        for obs_contig in result_buffer:

            exp_contig = exp_results[ obs_contig.name ]

            self.assertIn( obs_contig.name, exp_contig.name )
            self.assertEqual( obs_contig.sequence, exp_contig.sequence )  
            self.assertEqual( obs_contig.length, exp_contig.length )
            self.assertEqual( obs_contig.gc, exp_contig.gc )

    def test_add_contig_alreadyexists(self):

        ''' Set up files for test '''
        file_name = 'mock.contigs.fasta'
        exp_contigs, exp_seqs = self.create_contigs_file(file_name)

        self.db_m.create_blank_database()
        self.db_m.add_contigs(file_name)

        ''' Evaluate - try to add the contigs a second time '''
        with self.assertRaises(Exception) as cm:
            self.db_m.add_contigs(file_name)

        self.assertIn( 'Duplicate contig detected', str(cm.exception) )

        contig_df = self.db_m.get_contigs()
        self.assertEqual(contig_df.shape[0], len(exp_contigs) )

    def test_add_contigs(self):

        ''' Set up files for test '''
        file_name = 'mock.contigs.fasta'
        exp_contigs, _ = self.create_contigs_file(file_name)

        self.db_m.create_blank_database()
        self.db_m.add_contigs(file_name)

        ''' Evaluate '''
        contig_df = self.db_m.get_contigs()
        self.assertEqual(contig_df.shape[0], len(exp_contigs) )

    def test_get_contigs(self):

        ''' Set up files for test '''
        file_name = 'mock.contigs.fasta'
        exp_contigs, exp_seqs = self.create_contigs_file(file_name)

        exp_dict = { c: Contig(c, s) for c, s in zip(exp_contigs, exp_seqs) }

        self.db_m.create_blank_database()
        self.db_m.add_contigs(file_name)

        ''' Evaluate '''
        contig_df = self.db_m.get_contigs()
        for _, row in contig_df.iterrows():

            self.assertIn( row['contig_name'], exp_dict )
            obs_contig = exp_dict[ row['contig_name'] ]

            self.assertEqual( row['sequence'], obs_contig.sequence )
            self.assertEqual( row['length'], obs_contig.length )
            self.assertEqual( row['gc'], obs_contig.gc )

    def test_get_contig_by_name_exists(self):

        ''' Set up files for test '''
        file_name = 'mock.contigs.fasta'
        exp_contigs, exp_seqs = self.create_contigs_file(file_name)

        self.db_m.create_blank_database()
        self.db_m.add_contigs(file_name)

        ''' Evaluate '''
        contig_dict = self.db_m._get_contig_by_name( exp_contigs[0] )

        self.assertIsNotNone(contig_dict)
        self.assertEqual( contig_dict['contig_name'], exp_contigs[0] )

    def test_get_contig_by_name_notexists(self):

        ''' Set up files for test '''
        file_name = 'mock.contigs.fasta'
        exp_contigs, exp_seqs = self.create_contigs_file(file_name)

        self.db_m.create_blank_database()
        self.db_m.add_contigs(file_name)

        ''' Evaluate '''
        contig_dict = self.db_m._get_contig_by_name('xxxxx')

        self.assertIsNone(contig_dict)

    # endregion

    # region Stage 3.0 - Gene functions - generic

    def create_database_with_contigs(self):
        
        file_name = 'mock.contigs.fasta'
        exp_contigs, exp_seqs = self.create_contigs_file(file_name)

        self.db_m.create_blank_database()
        self.db_m.add_contigs(file_name)

        return exp_contigs

    def create_mock_gene(self, contig_name):
        mock_gene = Gene(gene_name='contig1_1', contig_name=contig_name, prediction_tool='test_spawn',
                         sequence_nt='ATCGAT', sequence_aa='CC',
                         start=1, stop=6, orientation=1)
        return mock_gene

    def test_obtain_gene_set(self):

        ''' Test '''
        contigs = self.create_database_with_contigs()
        first_gene = self.create_mock_gene( contigs[0] )
        second_gene = self.create_mock_gene( contigs[0] )
        second_gene._gene_name = '{}_v2'.format( first_gene.name )

        ''' Evaluate '''
        self.db_m._add_genes( [first_gene, second_gene] )
        gene_set = self.db_m._obtain_gene_set()

        self.assertEqual( set([ first_gene.name, second_gene.name ]), gene_set )

    def test_validate_genes(self):

        ''' Test '''
        contigs = self.create_database_with_contigs()
        mock_gene = self.create_mock_gene( contigs[0] )

        ''' Evaluate '''
        success_state = self.db_m._validate_genes( [mock_gene] )
        self.assertTrue( success_state )

    def test_validate_genes_nocontig(self):

        ''' Test '''
        _ = self.create_database_with_contigs()

        mock_gene = self.create_mock_gene( 'xxx' )

        ''' Evaluate '''
        with self.assertRaises(Exception) as cm:
            self.db_m._validate_genes( [mock_gene] )

        self.assertIn( 'Gene {} is not linked to an existing contig. Aborting...'.format(mock_gene.name), str(cm.exception) )

    def test_validate_genes_badint(self):

        ''' Test '''
        contigs = self.create_database_with_contigs()

        mock_gene = self.create_mock_gene( contigs[0] )
        mock_gene._start = 'a'

        ''' Evaluate '''
        with self.assertRaises(Exception) as cm:
            self.db_m._validate_genes( [mock_gene] )

        self.assertIn( 'Error converting values to int in gene {}. Aborting...'.format(mock_gene.name), str(cm.exception) )

    def test_validate_genes_duplicategene_exists(self):

        ''' Test '''
        contigs = self.create_database_with_contigs()
        mock_gene = self.create_mock_gene( contigs[0] )
        self.db_m._add_genes( [mock_gene] )

        ''' Evaluate - Fail is tested by adding the gene a second time '''
        with self.assertRaises(Exception) as cm:
            self.db_m._validate_genes( [mock_gene] )

        self.assertIn( 'Gene {} is already in the database. Aborting...'.format(mock_gene.name), str(cm.exception) )

    def test_validate_genes_duplicategene_adding(self):

        ''' Test '''
        _ = self.db_m.create_blank_database()
        mock_gene = self.create_mock_gene( 'xxx' )

        ''' Evaluate - Fail is tested by adding the gene a second time '''
        with self.assertRaises(Exception) as cm:
            self.db_m._validate_genes( [mock_gene, mock_gene] )

        self.assertIn( 'There are duplicate genes in the input file. Aborting...', str(cm.exception) )

    def test_add_genes(self):

        ''' Test '''
        contigs = self.create_database_with_contigs()
        mock_gene = self.create_mock_gene( contigs[0] )
        self.db_m._add_genes( [mock_gene] )

        ''' Evaluate '''
        gene_df = self.db_m.get_genes()
        self.assertEqual( gene_df.shape[0], 1)

    def test_get_genes(self):

        ''' Test '''
        contigs = self.create_database_with_contigs()
        mock_gene = self.create_mock_gene( contigs[0] )
        self.db_m._add_genes( [mock_gene] )

        ''' Evaluate '''
        gene_df = self.db_m.get_genes()
        self.assertEqual( gene_df.shape[0], 1 )

        for _, row in gene_df.iterrows():

            self.assertEqual( row['gene_name'], mock_gene.name )
            self.assertEqual( row['contig_name'], mock_gene.contig_name )
            self.assertEqual( row['prediction_tool'], mock_gene.prediction_tool )
            self.assertEqual( row['sequence_nt'], mock_gene.sequence_nt )
            self.assertEqual( row['sequence_aa'], mock_gene.sequence_aa )
            self.assertEqual( row['length_nt'], mock_gene.length_nt )
            self.assertEqual( row['length_aa'], mock_gene.length_aa )
            self.assertEqual( row['start_pos'], mock_gene.start )
            self.assertEqual( row['stop_pos'], mock_gene.stop )
            self.assertEqual( row['orientation'], mock_gene.orientation )

    def test_get_genes_bypred(self):

        ''' Test '''
        contigs = self.create_database_with_contigs()
        first_gene = self.create_mock_gene( contigs[0] )
        second_gene = self.create_mock_gene( contigs[0] )
        second_gene._gene_name = '{}_v2'.format( first_gene.name )

        alt_method = '{}_v2'.format(first_gene.prediction_tool)
        second_gene._prediction_tool = alt_method

        self.db_m._add_genes( [ first_gene, second_gene ] )
    
        ''' Evaluate '''
        gene_df = self.db_m.get_genes(prediction_method=alt_method)
        self.assertEqual( gene_df.shape[0], 1)

        for _, row in gene_df.iterrows():
            self.assertEqual( row['prediction_tool'], alt_method )

    def test_get_gene_by_name(self):

        ''' Test '''
        contigs = self.create_database_with_contigs()
        first_gene = self.create_mock_gene( contigs[0] )
        second_gene = self.create_mock_gene( contigs[0] )
        second_gene._gene_name = '{}_v2'.format(first_gene.name)
    
        ''' Evaluate '''
        self.db_m._add_genes( [ first_gene, second_gene ] )
        gene_dict = self.db_m.get_gene_by_name( first_gene.name )
        
        self.assertIsNotNone( gene_dict )
        self.assertEqual( gene_dict['gene_name'], first_gene.name )

    def test_get_gene_by_name_notexists(self):

        ''' Test '''
        contigs = self.create_database_with_contigs()
        mock_gene = self.create_mock_gene( contigs[0] )
        self.db_m._add_genes( [mock_gene] )

        ''' Evaluate '''
        gene_dict = self.db_m.get_gene_by_name('xxxxx')
        self.assertIsNone(gene_dict)

    def test_get_genes_by_contig(self):

        ''' Set up files for test '''
        contigs = self.create_database_with_contigs()
        first_gene = self.create_mock_gene( contigs[0] )
        second_gene = self.create_mock_gene( contigs[1] )
        second_gene._gene_name = 'asdasdasd'

        self.db_m._add_genes( [ first_gene, second_gene ] )

        ''' Evaluate '''
        gene_df = self.db_m.get_genes_by_contig( contigs[0] )
        self.assertEqual( gene_df.shape[0], 1 )

    def test_get_genes_by_contig_notexists(self):

        ''' Test '''
        contigs = self.create_database_with_contigs()
        mock_gene = self.create_mock_gene( contigs[0] )

        self.db_m._add_genes( [mock_gene] ) 

        ''' Evaluate '''
        gene_df = self.db_m.get_genes_by_contig('xxx')
        self.assertEqual( gene_df.shape[0], 0 )

    # endregion

    # region Stage 3.1 - Gene functions - Prodigal

    def create_genes_prod_nt(self, f_name):

        exp_keys = ['contig1_1', 'contig1_2', 'contig2_1', 'contig3_1']
        exp_meta = ['# 2 # 25 # 1 # ID=1_1;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc_cont=0.50',
                    '# 4 # 46 # 1 # ID=1_1;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc_cont=0.50',
                    '# 5 # 56 # -1 # ID=1_1;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc_cont=0.50',
                    '# 25 # 46 # 1 # ID=1_1;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc_cont=0.50']
        exp_seq = ['TCGATCGATCGATCGATCGATCG',
                   'GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT',
                   'TCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG',
                   'CGATCGATCGATCGATCGATC']

        return exp_keys, exp_meta, exp_seq

    def create_genes_prod_aa(self, f_name):

        exp_keys = ['contig1_1', 'contig1_2', 'contig2_1', 'contig3_1']
        exp_meta = ['# 2 # 25 # 1 # ID=1_1;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc_cont=0.50',
                    '# 4 # 46 # 1 # ID=1_1;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc_cont=0.50',
                    '# 5 # 56 # -1 # ID=1_1;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc_cont=0.50',
                    '# 25 # 46 # 1 # ID=1_1;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc_cont=0.50']
        exp_seq = ['SIDRSID*', 'DRSIDRSIDRSIDR*', 'SIDRSIDRSIDRSIDRS', 'RSIDRSI']

        return exp_keys, exp_meta, exp_seq

    def test_extract_prodigal_contig_simple(self):

        contig_name = self.db_m._extract_prodigal_contig('abc_1')
        self.assertEqual(contig_name, 'abc')

    def test_extract_prodigal_contig_complex(self):

        contig_name = self.db_m._extract_prodigal_contig('ab_c_1')
        self.assertEqual(contig_name, 'ab_c')

    def test_extract_prodigal_metadata(self):

        start, stop, orient = self.db_m._extract_prodigal_metadata('# 2 # 25 # -1 # ID=1_1;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc_cont=0.50')
        self.assertEqual(start, '2')
        self.assertEqual(stop, '25')
        self.assertEqual(orient, '-1')

    def test_add_genes_prodigal_aa_only(self):

        ''' Set up files for test '''
        _ = self.create_database_with_contigs()

        aa_file_name = 'mock.prod.faa'
        aa_keys, aa_meta, aa_seq = self.create_genes_prod_aa(aa_file_name)
        self.push_to_fasta(aa_file_name, aa_keys, aa_seq, seq_meta=aa_meta)

        self.db_m.add_genes_prodigal(aa_file=aa_file_name)
        exp_dict = { k: (m, aa) for k, m, aa in zip(aa_keys, aa_meta, aa_seq) }

        ''' Evaluate '''        
        genes_df = self.db_m.get_genes()
        self.assertEqual( genes_df.shape[0], len(aa_keys) )

        for _, row in genes_df.iterrows():

            self.assertTrue( row['gene_name'] in exp_dict )

            ''' Bundle the expected values for comparison'''
            exp_meta, exp_aa = exp_dict[ row['gene_name'] ]
            exp_start, exp_stop, exp_orient = self.db_m._extract_prodigal_metadata(exp_meta)

            for e, o in zip( (exp_aa, '', exp_start, exp_stop, exp_orient), ['sequence_aa', 'sequence_nt', 'start_pos', 'stop_pos', 'orientation'] ):
                self.assertEqual( e, str( row[o] ) )

    def test_add_genes_prodigal_nt_only(self):

        ''' Set up files for test '''
        _ = self.create_database_with_contigs()

        nt_file_name = 'mock.prod.fna'
        nt_keys, nt_meta, nt_seq = self.create_genes_prod_nt(nt_file_name)
        self.push_to_fasta(nt_file_name, nt_keys, nt_seq, seq_meta=nt_meta)

        self.db_m.add_genes_prodigal(nt_file=nt_file_name)  
        exp_dict = { k: (m, nt) for k, m, nt in zip(nt_keys, nt_meta, nt_seq) }
        
        ''' Evaluate '''
        gene_df = self.db_m.get_genes()
        self.assertEqual( gene_df.shape[0], len(nt_keys) )

        for _, row in gene_df.iterrows():

            self.assertTrue( row['gene_name'] in exp_dict )

            ''' Bundle the expected values for comparison'''
            exp_meta, exp_nt = exp_dict[ row['gene_name'] ]
            exp_start, exp_stop, exp_orient = self.db_m._extract_prodigal_metadata(exp_meta)

            for e, o in zip( ('', exp_nt, exp_start, exp_stop, exp_orient), ['sequence_aa', 'sequence_nt', 'start_pos', 'stop_pos', 'orientation'] ):
                self.assertEqual( e, str( row[o] ) )

    def test_add_genes_prodigal_both(self):

        ''' Set up files for test '''
        _ = self.create_database_with_contigs()

        aa_file_name = 'mock.prod.faa'
        aa_keys, aa_meta, aa_seq = self.create_genes_prod_aa(aa_file_name)
        self.push_to_fasta(aa_file_name, aa_keys, aa_seq, seq_meta=aa_meta)

        nt_file_name = 'mock.prod.fna'
        nt_keys, nt_meta, nt_seq = self.create_genes_prod_nt(nt_file_name)
        self.push_to_fasta(nt_file_name, nt_keys, nt_seq, seq_meta=nt_meta)

        self.db_m.add_genes_prodigal(aa_file=aa_file_name, nt_file=nt_file_name)
        exp_dict = { k: (m, aa, nt) for k, m, aa, nt in zip(aa_keys, aa_meta, aa_seq, nt_seq) }

        ''' Evaluate '''
        gene_df = self.db_m.get_genes()
        self.assertEqual( gene_df.shape[0], len(exp_dict) )

        for _, row in gene_df.iterrows():

            self.assertTrue( row['gene_name'] in exp_dict )

            ''' Bundle the expected values for comparison'''
            exp_meta, exp_aa, exp_nt = exp_dict[ row['gene_name'] ]
            exp_start, exp_stop, exp_orient = self.db_m._extract_prodigal_metadata(exp_meta)

            for e, o in zip( (exp_aa, exp_nt, exp_start, exp_stop, exp_orient), ['sequence_aa', 'sequence_nt', 'start_pos', 'stop_pos', 'orientation'] ):
                self.assertEqual( e, str( row[o] ) )

    def test_add_genes_prodigal_failedcoord(self):

        ''' Set up files for test '''
        _ = self.create_database_with_contigs()

        aa_file_name = 'mock.prod.faa'
        aa_keys, aa_meta, aa_seq = self.create_genes_prod_aa(aa_file_name)

        aa_meta[0] = aa_meta[0].replace('# 2 # 25 #', '# abc # 25 #')
        self.push_to_fasta(aa_file_name, aa_keys, aa_seq, seq_meta=aa_meta)

        ''' Test '''
        with self.assertRaises(Exception) as cm:
            self.db_m.add_genes_prodigal(aa_file=aa_file_name)

        self.assertIn( 'Error converting values to int in gene {}. Aborting...'.format(aa_keys[0]), str(cm.exception) )

        genes_df = self.db_m.get_genes()
        self.assertEqual( genes_df.shape[0], 0 )

    # endregion

    # region Stage 3.2 - Gene functions - MeTaxa

    def create_metaxa_ssu(self, f_name):

        exp_keys = ['contig1|B', 'contig2|E']
        exp_meta = ['Predicted Bacterial SSU rRNA (34 bp) From domain V1l to V9r on main strand',
                    'Predicted Eukaryotic SSU rRNA (23 bp) From domain V1l to V9r on complementary strand']
        exp_seq = ['CGATCGATCGATCGATCGATCGATCGATCGATCG', 'TCGATCGATCGATCGATCGATCG']

        self.push_to_fasta(f_name, exp_keys, exp_seq, seq_meta=exp_meta)
        return exp_keys, exp_meta, exp_seq

    def create_metaxa_key(self, gene_string):
        *prefix, rna_tax, rna_type = gene_string.split('_')
        prefix = '_'.join(prefix)
        return '{}|{}'.format(prefix, rna_tax)
    
    def test_metaxa_extract_names(self):

        ''' String pulled from an example of the *.extraction.fasta file:
                NODE_10065_length_8699_cov_9.406272|E Predicted Eukaryotic LSU rRNA (3505 bp) From domain C01 to C16 on main strand '''

        test_str = 'NODE_10065_length_8699_cov_9.406272|E'
        exp_contig = 'NODE_10065_length_8699_cov_9.406272'
        exp_gene = 'NODE_10065_length_8699_cov_9.406272_E_mock'

        obs_contig, obs_gene = self.db_m._metaxa_extract_names(test_str, 'mock')

        self.assertEqual( exp_contig, obs_contig )
        self.assertEqual( exp_gene, obs_gene )

    def test_metaxa_extract_orientation(self):

        ''' String pulled from an example of the *.extraction.fasta file:
                NODE_10065_length_8699_cov_9.406272|E Predicted Eukaryotic LSU rRNA (3505 bp) From domain C01 to C16 on main strand '''
        test_str_forward = 'Predicted Eukaryotic LSU rRNA (3505 bp) From domain C01 to C16 on main strand'
        test_str_reverse = 'Predicted Eukaryotic LSU rRNA (3505 bp) From domain C01 to C16 on complementary strand'

        obs_forward = self.db_m._metaxa_extract_orientation(test_str_forward)
        obs_reverse = self.db_m._metaxa_extract_orientation(test_str_reverse)

        self.assertEqual( 1, obs_forward )
        self.assertEqual( -1, obs_reverse )

    def test_metaxa_locate_coords_forward(self):

        ''' Set up files for test '''
        contigs = self.create_database_with_contigs()

        ''' Evaluate '''
        obs_start, obs_stop = self.db_m._metaxa_locate_coords('contig1', 'AAATTTCCC', 1)

        self.assertEqual( 1, obs_start )
        self.assertEqual( 9, obs_stop )

    def test_metaxa_locate_coords_reverse(self):

        ''' Set up files for test '''
        contigs = self.create_database_with_contigs()

        ''' Evaluate '''
        obs_start, obs_stop = self.db_m._metaxa_locate_coords('contig1', 'CGGGAAA', -1)

        self.assertEqual( 10, obs_start )
        self.assertEqual( 4, obs_stop )

    def test_add_genes_metaxa(self):

        ''' Set up files for test '''
        contigs = self.create_database_with_contigs()

        file_name = 'mock.ssu.fna'
        exp_keys, *_ = self.create_metaxa_ssu(file_name)

        ''' Evaluate '''
        self.db_m.add_genes_metaxa(file_name, 'SSU')
        
        gene_df = self.db_m.get_genes()
        self.assertEqual( gene_df.shape[0], len(exp_keys) )

        annot_df = self.db_m.get_annotations()
        self.assertEqual( annot_df.shape[0], len(exp_keys) )

    # endregion

    # region Stage 3.3 - Gene functions - Aragorn

    def create_aragorn_trna(self, f_name):

        exp_keys = ['1-1', '1-2', '2-1']
        exp_meta = ['tRNA-Arg(gcg) [3,16]', 'tRNA-Pro(aaa) [28,32]', 'tRNA-His(ggg) [18,28]']
        exp_seq = ['CGATCGATCGATCG', 'GATCG', 'CGATCGATCG']

        self.push_to_fasta(f_name, exp_keys, exp_seq, seq_meta=exp_meta)
        return exp_keys, exp_meta, exp_seq

    def create_aragorn_meta(self, contig_file, gene_names, metadata_values):

        contig_sequence = self.db_m._read_contig_sequence(contig_file, False)

        for i, (gene_name, metadata) in enumerate( zip(gene_names, metadata_values) ):

            contig_i, gene_i = self.db_m._aragorn_split_indicies(gene_name)
            contig_name = contig_sequence[ int(contig_i) - 1 ]

            trna_type = self.db_m._aragorn_extract_type(metadata)
            gene_names[i] = '{}_{}_{}'.format(contig_name, gene_i, trna_type)

        return gene_names

    def test_aragorn_split_indicies(self):

        mock_string = '1-2'
        contig_i, gene_i = self.db_m._aragorn_split_indicies(mock_string)

        self.assertEqual( contig_i, 1 )
        self.assertEqual( gene_i, 2 )

    def test_aragorn_extract_type(self):

        mock_str = 'tRNA-Arg(gcg) [58280,58368]'
        exp_trna_type = 'tRNA-Arg-gcg'
        obs_trna_type = self.db_m._aragorn_extract_type(mock_str)

        self.assertEqual(exp_trna_type, obs_trna_type)

    def test_aragorn_extract_coords(self):

        mock_str = 'tRNA-Arg(gcg) [58280,58368]'
        obs_start, obs_stop = self.db_m._aragorn_extract_coords(mock_str)

        self.assertEqual('58280', obs_start)
        self.assertEqual('58368', obs_stop)

    def test_aragorn_extract_names(self):

        exp_sequence = [ 'a', 'b' ]
        exp_metadata = 'tRNA-Arg(gcg) [58280,58368]'

        obs_contig, obs_gene = self.db_m._aragorn_extract_names(exp_sequence, exp_metadata, 1, 1)
        self.assertEqual( obs_contig, exp_sequence[0] )
        self.assertEqual( obs_gene, 'a_1_tRNA-Arg-gcg' )
    
    def test_add_genes_aragorn(self):

        ''' Test '''
        contigs = self.create_database_with_contigs()

        aragorn_file_name = 'mock.trna.fna'
        exp_keys, exp_meta, exp_seq = self.create_aragorn_trna(aragorn_file_name)

        obs_keys = self.create_aragorn_meta('mock.contigs.fasta', exp_keys, exp_meta)
        exp_dict = { k: (m, seq) for k, m, seq in zip(obs_keys, exp_meta, exp_seq) }

        ''' Evaluate '''
        #trna_input_str = '{},{}'.format('mock.contigs.fasta', aragorn_file_name)
        self.db_m.add_genes_aragorn('mock.contigs.fasta', aragorn_file_name)

        gene_df = self.db_m.get_genes()
        self.assertEqual( gene_df.shape[0], len(exp_keys) )

    # endregion

    # region Stage 4 - Annotation functions

    def create_database_for_gene(self):

        contigs = self.create_database_with_contigs()

        mock_gene = self.create_mock_gene( contigs[0] )
        self.db_m._add_genes( [mock_gene] )

        return contigs, mock_gene

    def create_mock_annotation(self, base_gene):

        mock_annotation = Annotation(gene_name=base_gene.name, method='test_method', annotation_database='test_db', accession='DW12345',
                                     identity=70.5, coverage=30.2, evalue=0.001, description='Random text goes here')
        return mock_annotation

    def test_annotations_are_current_true(self):
        self.assertTrue( self.db_m._annotations_are_current(True) )

    def test_annotations_are_current_false(self):
        self.assertFalse( self.db_m._annotations_are_current(False) )

    def test_annotations_are_current_none(self):
        self.assertTrue( self.db_m._annotations_are_current(None) )

    def test_validate_annotations(self):

        ''' Test '''
        contigs, mock_gene = self.create_database_for_gene()
        mock_annotation = self.create_mock_annotation(mock_gene)

        ''' Evaluate '''
        success_state = self.db_m._validate_annotations( [mock_annotation] )
        self.assertTrue(success_state)
    
    def test_validate_annotations_badfloat(self):

        ''' Test '''
        contigs, mock_gene = self.create_database_for_gene()
        mock_annotation = self.create_mock_annotation(mock_gene)
        mock_annotation._coverage = 'xxx'

        ''' Evaluate '''
        with self.assertRaises(Exception) as cm:
            self.db_m._validate_annotations( [mock_annotation] )

        self.assertIn( 'Error converting values to float in gene {}. Aborting...'.format(mock_annotation.gene_name), str(cm.exception) )

    def test_validate_annotations_missinggene(self):

        ''' Test '''
        _ = self.create_database_for_gene()
        mock_annotation = Annotation(gene_name='xxx', method='test_method', annotation_database='test_db', accession='DW12345',
                                     identity=70.5, coverage=30.2, evalue=0.001, description='Random text goes here')

        ''' Evaluate '''
        with self.assertRaises(Exception) as cm:
            self.db_m._validate_annotations( [mock_annotation] )

        self.assertIn( 'Annotation entry (method {}, database {}) is not matched to a valid gene (xxx). Aborting...'.format(mock_annotation.method, mock_annotation.annotation_database), str(cm.exception) )

    def test_add_annotation(self):

        ''' Test '''
        contigs, mock_gene = self.create_database_for_gene()
        mock_annotation = self.create_mock_annotation(mock_gene)

        ''' Evaluate '''
        self.db_m._add_annotations( [mock_annotation] )
        annot_df = self.db_m.get_annotations()
        self.assertEqual(annot_df.shape[0], 1)

    def test_get_annotations(self):

        ''' Test '''
        contigs, mock_gene = self.create_database_for_gene()
        mock_annotation = self.create_mock_annotation(mock_gene)

        self.db_m._add_annotations( [mock_annotation] )

        ''' Evaluate '''        
        annot_df = self.db_m.get_annotations()
        self.assertEqual(annot_df.shape[0], 1)

        for _, row in annot_df.iterrows():
            self.assertEqual( row['method'], mock_annotation.method )
            self.assertEqual( row['annotation_db'], mock_annotation.annotation_database )
            self.assertEqual( row['accession'], mock_annotation.accession )
            self.assertEqual( row['identity'], mock_annotation.identity )
            self.assertEqual( row['coverage'], mock_annotation.coverage )
            self.assertEqual( row['evalue'], mock_annotation.evalue )
            self.assertEqual( row['description'], mock_annotation.description )
            self.assertEqual( row['iscurrent'], mock_annotation.is_current )

    def test_get_annotations_by_gene(self):

        ''' Test '''
        contigs, first_gene = self.create_database_for_gene()
        second_gene = self.create_mock_gene( contigs[1] )

        first_gene._gene_name = 'mock_1'
        first_annotation = self.create_mock_annotation(first_gene)
        second_gene._gene_name = 'mock_2'
        second_annotation = self.create_mock_annotation(second_gene)

        self.db_m._add_genes( [first_gene, second_gene] )
        self.db_m._add_annotations( [ first_annotation, second_annotation ] )

        ''' Evaluate '''
        annot_df = self.db_m.get_annotations()
        self.assertEqual( annot_df.shape[0], 2 )

        annot_df = self.db_m.get_annotations_by_gene(first_gene.name)
        self.assertEqual( annot_df.shape[0], 1 )

    def test_get_annotation_by_gene_method(self):

        ''' Test
            Create 2 genes, the first gene linked to 2 annotations and the second linked to 1 annotatation
            On search, would expect the following:
                If get_annotations() -> 3 hits
                If get_annotations_by_gene() -> 2 hits
                If get_annotation_by_gene_method -> 1 hit '''
        contigs, first_gene = self.create_database_for_gene()
        second_gene = self.create_mock_gene( contigs[1] )

        first_gene._gene_name = 'mock_1'
        first_annotation = self.create_mock_annotation(first_gene)

        second_annotation = self.create_mock_annotation(first_gene)
        second_annotation._method = 'xxx'

        second_gene._gene_name = 'mock_2'
        third_annotation = self.create_mock_annotation(second_gene)
        third_annotation._method = 'xxx'

        self.db_m._add_genes( [first_gene, second_gene] )
        self.db_m._add_annotations( [ first_annotation, second_annotation, third_annotation ] )

        ''' Evaluate '''
        annot_df = self.db_m.get_annotations()
        self.assertEqual( annot_df.shape[0], 3 )

        annot_df = self.db_m.get_annotation_by_gene_method(first_gene.name, first_annotation.method)
        self.assertEqual( annot_df.shape[0], 1 )

    def test_set_annotation_historic(self):

        ''' Test '''
        contigs, mock_gene = self.create_database_for_gene()
        mock_annotation = self.create_mock_annotation(mock_gene)

        self.db_m._add_annotations( [mock_annotation] )

        ''' Evaluate '''
        annot_df = self.db_m.get_annotations()
        self.assertEqual( annot_df.shape[0], 1 )

        self.db_m._set_annotations_historic( [mock_annotation] )

        annot_df = self.db_m.get_annotations()
        self.assertEqual( annot_df.shape[0], 0 )

    def test_add_annotations_by_table(self):

        ''' Test '''
        contigs, mock_gene = self.create_database_for_gene()
        mock_annotation = self.create_mock_annotation(mock_gene)

        fake_df_buffer = [ { 'gene': mock_annotation.gene_name, 'method': mock_annotation.method, 'database': mock_annotation.annotation_database,
                             'hit': mock_annotation.accession, 'identity': mock_annotation.identity, 'coverage': mock_annotation.coverage, 'evalue': mock_annotation.evalue,
                             'description': mock_annotation.description } ]

        annot_file_name = 'mock.annot.txt'
        self.push_to_table(annot_file_name, pd.DataFrame(fake_df_buffer))

        ''' Evaluate '''
        self.db_m.add_annotations_by_table(annot_file_name)

        annot_df = self.db_m.get_annotations()
        self.assertEqual(annot_df.shape[0], len(fake_df_buffer) )

        for i, row in annot_df.iterrows():
            self.assertEqual( row['gene_name'], mock_annotation.gene_name )
            self.assertEqual( row['method'], mock_annotation.method )
            self.assertEqual( row['annotation_db'], mock_annotation.annotation_database )
            self.assertEqual( row['accession'], mock_annotation.accession )
            self.assertEqual( row['identity'], mock_annotation.identity )
            self.assertEqual( row['coverage'], mock_annotation.coverage )
            self.assertEqual( row['evalue'], mock_annotation.evalue )
            self.assertEqual( row['description'], mock_annotation.description )
            self.assertEqual( row['iscurrent'], mock_annotation.is_current )

    def test_add_annotations_by_table_failedvalue(self):

        ''' Test '''
        contigs, mock_gene = self.create_database_for_gene()
        mock_annotation = self.create_mock_annotation(mock_gene)

        fake_df_buffer = [ { 'gene': mock_annotation.gene_name, 'method': mock_annotation.method, 'database': mock_annotation.annotation_database,
                             'hit': mock_annotation.accession, 'identity': 'xxx', 'coverage': mock_annotation.coverage, 'evalue': mock_annotation.evalue,
                             'description': mock_annotation.description } ]

        annot_file_name = 'mock.annot.txt'
        self.push_to_table(annot_file_name, pd.DataFrame(fake_df_buffer))

        ''' Evaluate '''
        with self.assertRaises(Exception) as cm:
            self.db_m.add_annotations_by_table(annot_file_name)

        self.assertIn( 'Aborting', str(cm.exception) )

    # endregion

    # region Stage 5 - Coverage functions

    def create_mock_mapping(self, contig_name):
        return Mapping(target_name=contig_name, sample_name='asd', reads_mapped=12345)

    def test_validate_mapping_contig(self):

        ''' Test '''
        contigs = self.create_database_with_contigs()
        mock_mapping = self.create_mock_mapping( contigs[0] )

        ''' Evaluate '''
        success_state = self.db_m._validate_mapping( [mock_mapping], target='contig' )
        self.assertTrue(success_state)

    def test_validate_mapping_contig_badint(self):

        ''' Test '''
        contigs = self.create_database_with_contigs()
        mock_mapping = self.create_mock_mapping( contigs[0] )
        mock_mapping._reads_mapped = 'x'

        ''' Evaluate '''
        with self.assertRaises(Exception) as cm:
            self.db_m._validate_mapping( [mock_mapping], target='contig' )

        self.assertIn( 'Error converting contig mapping value to int in {}. Aborting...'.format(mock_mapping.target_name), str(cm.exception) )

    def test_validate_mapping_contig_missingcontig(self):

        ''' Test '''
        self.db_m.create_blank_database()
        mock_mapping = self.create_mock_mapping( 'xxx' )

        ''' Evaluate '''
        with self.assertRaises(Exception) as cm:
            self.db_m._validate_mapping( [mock_mapping], target='contig' )

        self.assertIn( 'Record {} is not linked to an existing contig. Aborting...'.format(mock_mapping.target_name), str(cm.exception) )

    def test_add_coverage_table(self):

        ''' Test '''
        contigs = self.create_database_with_contigs()
        first_mapping = self.create_mock_mapping( contigs[0] )
        second_mapping = self.create_mock_mapping( contigs[1] )

        coverage_buffer = [ { 'Contig': first_mapping.target_name, first_mapping.sample_name: first_mapping.reads_mapped },
                            { 'Contig': first_mapping.target_name, first_mapping.sample_name: first_mapping.reads_mapped }  ]
        coverage_df = pd.DataFrame(coverage_buffer)
        self.push_to_table('mock.coverage.txt', coverage_df)

        ''' Evaluate '''
        self.db_m.add_coverage_table('mock.coverage.txt')

        coverage_df = self.db_m.get_coverage()
        self.assertEqual( coverage_df.shape[0], len(coverage_buffer) )

    def test_add_coverage_table_nocontig(self):

        ''' Test '''
        contigs = self.create_database_with_contigs()
        mock_mapping = self.create_mock_mapping( contigs[0] )

        coverage_buffer = [ { 'Contig': mock_mapping.target_name, mock_mapping.sample_name: mock_mapping.reads_mapped },
                            { 'Contig': 'xxx', mock_mapping.sample_name: mock_mapping.reads_mapped }  ]
        coverage_df = pd.DataFrame(coverage_buffer)
        self.push_to_table('mock.coverage.txt', coverage_df)

        ''' Evaluate '''
        with self.assertRaises(Exception) as cm:
            self.db_m.add_coverage_table('mock.coverage.txt')

        self.assertIn( 'Record xxx is not linked to an existing contig. Aborting...', str(cm.exception) )

        coverage_df = self.db_m.get_coverage()
        self.assertEqual( coverage_df.shape[0], 0 )

    def test_get_coverage(self):

        ''' Test '''
        contigs = self.create_database_with_contigs()
        first_mapping = self.create_mock_mapping( contigs[0] )
        second_mapping = self.create_mock_mapping( contigs[1] )

        coverage_buffer = [ { 'Contig': first_mapping.target_name, first_mapping.sample_name: first_mapping.reads_mapped },
                            { 'Contig': second_mapping.target_name, second_mapping.sample_name: second_mapping.reads_mapped }  ]
        coverage_df = pd.DataFrame(coverage_buffer)
        self.push_to_table('mock.coverage.txt', coverage_df)

        self.db_m.add_coverage_table('mock.coverage.txt')

        ''' Evaluate '''
        coverage_df = self.db_m.get_coverage()
        self.assertEqual( coverage_df.shape[0], len(coverage_buffer) )

        for i, row in coverage_df.iterrows():

            exp_vals = coverage_buffer[i]
            self.assertEqual( row['contig_name'], exp_vals['Contig'] )
            self.assertEqual( row['reads_mapped'], exp_vals[ row['sample'] ] )
            self.assertIn( row['sample'], exp_vals )
 
    # endregion

    # region Stage 6 - Transcript functions

    def test_validate_mapping_gene(self):

        ''' Test '''
        contigs, mock_gene = self.create_database_for_gene()
        mock_mapping = self.create_mock_mapping( mock_gene.name )

        ''' Evaluate '''
        success_state = self.db_m._validate_mapping( [mock_mapping], target='gene' )
        self.assertTrue(success_state)

    def test_validate_mapping_gene_badint(self):

        ''' Test '''
        contigs, mock_gene = self.create_database_for_gene()
        mock_mapping = self.create_mock_mapping( mock_gene.name )
        mock_mapping._reads_mapped = 'x'

        ''' Evaluate '''
        with self.assertRaises(Exception) as cm:
            self.db_m._validate_mapping( [mock_mapping], target='gene' )

        self.assertIn( 'Error converting gene mapping value to int in {}. Aborting...'.format(mock_mapping.target_name), str(cm.exception) )
 
    def test_validate_mapping_gene_missinggene(self):

        ''' Test '''
        self.db_m.create_blank_database()
        mock_mapping = self.create_mock_mapping( 'xxx' )

        ''' Evaluate '''
        with self.assertRaises(Exception) as cm:
            self.db_m._validate_mapping( [mock_mapping], target='gene' )

        self.assertIn( 'Record {} is not linked to an existing gene. Aborting...'.format(mock_mapping.target_name), str(cm.exception) )

    def test_add_transcript_table(self):

        ''' Test '''
        contigs, mock_gene = self.create_database_for_gene()
        mock_mapping = self.create_mock_mapping( mock_gene.name )

        transcript_buffer = [ { 'Gene': mock_mapping.target_name, mock_mapping.sample_name: mock_mapping.reads_mapped } ]
        transcript_df = pd.DataFrame(transcript_buffer)
        self.push_to_table('mock.transcript.txt', transcript_df)

        ''' Evaluate '''
        self.db_m.add_transcript_table('mock.transcript.txt')

        transcript_df = self.db_m.get_transcript()
        self.assertEqual( transcript_df.shape[0], len(transcript_buffer) )

    def test_add_transcript_table_nogene(self):

        ''' Test '''
        contigs, mock_gene = self.create_database_for_gene()
        mock_mapping = self.create_mock_mapping( 'xxx' )

        transcript_buffer = [ { 'Gene': mock_mapping.target_name, mock_mapping.sample_name: mock_mapping.reads_mapped } ]
        transcript_df = pd.DataFrame(transcript_buffer)
        self.push_to_table('mock.transcript.txt', transcript_df)

        ''' Evaluate '''
        with self.assertRaises(Exception) as cm:
            self.db_m.add_transcript_table('mock.transcript.txt')

        self.assertIn( 'Record xxx is not linked to an existing gene. Aborting...', str(cm.exception) )

        transcript_df = self.db_m.get_transcript()
        self.assertEqual( transcript_df.shape[0], 0 )

    def test_get_transcript(self):

        ''' Test '''
        contigs, mock_gene = self.create_database_for_gene()
        mock_mapping = self.create_mock_mapping( mock_gene.name )

        transcript_buffer = [ { 'Gene': mock_mapping.target_name, mock_mapping.sample_name: mock_mapping.reads_mapped } ]
        transcript_df = pd.DataFrame(transcript_buffer)
        self.push_to_table('mock.transcript.txt', transcript_df)

        self.db_m.add_transcript_table('mock.transcript.txt')

        ''' Evaluate '''
        transcript_df = self.db_m.get_transcript()
        self.assertEqual( transcript_df.shape[0], len(transcript_buffer) )

        for i, row in transcript_df.iterrows():

            exp_vals = transcript_buffer[0]
            self.assertEqual( row['gene_name'], exp_vals['Gene'] )
            self.assertEqual( row['reads_mapped'], exp_vals[ row['sample'] ] )
            self.assertIn( row['sample'], exp_vals )

    # endregion

    # region Stage 7 - Bin functions

    def create_mock_bin(self):

        mock_bin = Bin(bin_name= 'bin1', taxonomy='a;b;c',
                       completeness=95.12, contamination=2.97, strain_heterogeneity=12.25,
                       version_name='v1', date_created=datetime.now().isoformat())
        return mock_bin

    def create_dummy_bintable(self, f_name, bin_list):

        buffer = [ { 'Bin': b.name, 'Taxonomy': b.taxonomy, 'Completeness': b.completeness, 'Contamination': b.contamination,
                     'Heterogeneity': b.strain_heterogeneity } for b in bin_list ]

        self.push_to_table(f_name, pd.DataFrame(buffer) )

    def test_obtain_bin_set(self):

        ''' Test '''
        first_bin = self.create_mock_bin()
        second_bin = self.create_mock_bin()
        second_bin._bin_name = 'bin2'

        self.db_m.create_blank_database()
        self.db_m._add_bins( [first_bin, second_bin] )

        ''' Evaluate '''
        bin_set = self.db_m._obtain_bin_set()

        exp_bins = set( [ first_bin.name, second_bin.name ] )
        self.assertEqual( len( exp_bins & bin_set), len(exp_bins) )

    def test_validate_bins(self):

        ''' Test '''
        self.db_m.create_blank_database()
        mock_bin = self.create_mock_bin()

        ''' Evaluate '''
        success_state = self.db_m._validate_bins( [mock_bin] )
        self.assertTrue(success_state)

    def test_validate_bins_badfloat(self):

        ''' Test '''
        self.db_m.create_blank_database()
        mock_bin = self.create_mock_bin()
        mock_bin._completeness = 'x'

        ''' Evaluate '''
        with self.assertRaises(Exception) as cm:
            self.db_m._validate_bins( [mock_bin] )

        self.assertIn( 'Error converting values to float in bin {}. Aborting...'.format(mock_bin.name), str(cm.exception) )

    def test_validate_bins_duplicateentry_db(self):

        ''' Test '''
        self.db_m.create_blank_database()
        mock_bin = self.create_mock_bin()

        self.db_m._add_bins( [mock_bin] )

        ''' Evaluate '''
        with self.assertRaises(Exception) as cm:
            self.db_m._validate_bins( [mock_bin] )

        self.assertIn( 'Duplicate entry for bin {} detected. Aborting...'.format(mock_bin.name), str(cm.exception) )

    def test_validate_bins_duplicateentry_file(self):

        ''' Test '''
        self.db_m.create_blank_database()
        mock_bin = self.create_mock_bin()

        ''' Evaluate '''
        with self.assertRaises(Exception) as cm:
            self.db_m._validate_bins( [mock_bin, mock_bin] )

        self.assertIn( 'There are duplicate bins in the input file. Aborting...', str(cm.exception) )

    def test_index_contig_folder(self):

        ''' Test '''
        self.make_directory('mock_bin_files/')
        exp_contigs, exp_seqs = self.create_contigs_file('mock.contigs.fna')

        self.push_to_fasta('mock_bin_files/mock.bin1.fasta', exp_contigs, exp_seqs)
        self.push_to_fasta('mock_bin_files/mock.bin2.fasta', exp_contigs, exp_seqs)

        ''' Evaluate '''
        contig_binning_dict = self.db_m._index_contig_folder('mock_bin_files/')

        self.assertIn( 'mock.bin1', contig_binning_dict)
        self.assertIn( 'mock.bin2', contig_binning_dict)

        self.assertEqual( exp_contigs, contig_binning_dict['mock.bin1'] )
        self.assertEqual( exp_contigs, contig_binning_dict['mock.bin2'] )

    def test_validate_contig_binning(self):

        ''' Test '''
        contigs = self.create_database_with_contigs()

        mock_bin = self.create_mock_bin()
        mock_dict = { mock_bin.name: contigs }

        ''' Evaluate '''
        success_state = self.db_m._validate_contig_binning( [mock_bin], mock_dict )
        self.assertTrue(success_state)

    def test_validate_contig_binning_missingcontig(self):

        ''' Test '''
        contigs = self.create_database_with_contigs()
        contigs.append('xxx')

        mock_bin = self.create_mock_bin()
        mock_dict = { mock_bin.name: contigs }

        ''' Evaluate '''
        with self.assertRaises(Exception) as cm:
            self.db_m._validate_contig_binning( [mock_bin], mock_dict )

        self.assertIn( 'The contig xxx in bin {} does not exist in the database. Aborting...'.format(mock_bin.name), str(cm.exception) )

    def test_add_bin_table(self):

        ''' Test '''
        contigs = self.create_database_with_contigs()
        mock_bin = self.create_mock_bin()

        file_name = 'mock.binning_table.txt'
        self.create_dummy_bintable(file_name, [mock_bin])

        dir_name = 'mock_bin_files'
        self.make_directory(dir_name)
        self.push_to_fasta('{}/{}.fasta'.format(dir_name, mock_bin.name), contigs, ['aaa', 'ttt', 'ccc'])

        ''' Evaluate '''
        self.db_m.add_bin_table(file_name, dir_name)
        bin_df = self.db_m.get_all_bins()

        self.assertEqual( bin_df.shape[0], 1 )

        for _, row in bin_df.iterrows():

            self.assertEqual(row['bin_name'], mock_bin.name)
            self.assertEqual(row['taxonomy'], mock_bin.taxonomy)
            self.assertEqual(row['completeness'], mock_bin.completeness)
            self.assertEqual(row['contamination'], mock_bin.contamination)
            self.assertEqual(row['strainheterogeneity'], mock_bin.strain_heterogeneity)
            self.assertEqual(row['version_name'], 'mock.binning_table')

    def test_add_bin_table_failedvalue(self):

        ''' Test '''
        contigs = self.create_database_with_contigs()
        mock_bin = self.create_mock_bin()

        file_name = 'mock.binning_table.txt'
        self.create_dummy_bintable(file_name, [mock_bin])

        dir_name = 'mock_bin_files'
        self.make_directory(dir_name)
        self.push_to_fasta('{}/{}.fasta'.format(dir_name, mock_bin.name), contigs, ['aaa', 'ttt', 'ccc'])

        ''' Set the completeness to a string value '''
        df = pd.read_csv(file_name, sep='\t')
        df.iloc[ 0 , 1 ] = 'a'
        df.to_csv(file_name, sep='\t', index=False)

        ''' Evaluate '''
        with self.assertRaises(Exception) as cm:
            self.db_m.add_bin_table(file_name, dir_name)

        self.assertIn( 'Error converting values to float in bin', str(cm.exception) )

        bin_df = self.db_m.get_all_bins()
        self.assertEqual( bin_df.shape[0], 0)

    def test_get_all_bins(self):

        ''' Test '''
        contigs = self.create_database_with_contigs()
        first_bin = self.create_mock_bin()

        second_bin = self.create_mock_bin()
        second_bin._bin_name = 'bin2'

        self.db_m._add_bins( [first_bin, second_bin] )

        ''' Evaluate '''
        bin_df = self.db_m.get_all_bins()
        self.assertEqual( bin_df.shape[0], 2 )

    def test_get_bins_by_version(self):

        ''' Test '''
        contigs = self.create_database_with_contigs()
        first_bin = self.create_mock_bin()

        second_bin = self.create_mock_bin()
        second_bin._bin_name = 'bin2'
        second_bin._version_name = 'V2'
        
        self.db_m._add_bins( [first_bin, second_bin] )

        ''' Evaluate '''
        bin_df = self.db_m.get_bins_by_version(first_bin.version_name)
        self.assertEqual( bin_df.shape[0], 1 )

        for _, row in bin_df.iterrows():
            self.assertEqual( row['version_name'], first_bin.version_name )

    def test_get_current_bins(self):

        ''' Test '''
        contigs = self.create_database_with_contigs()
        first_bin = self.create_mock_bin()

        ''' Induce a short sleep, to ensure the second entry is measurably different from te first one '''
        time.sleep(2)
        second_bin = self.create_mock_bin()
        second_bin._bin_name = 'bin2'
        
        self.db_m._add_bins( [first_bin, second_bin] )

        ''' Evaluate '''
        bin_df = self.db_m.get_current_bins()
        self.assertEqual( bin_df.shape[0], 1 )

        for _, row in bin_df.iterrows():
            self.assertEqual( row['date_created'], second_bin.date_created )

    def test_get_contigs_in_bin(self):

        ''' Test '''
        contigs = self.create_database_with_contigs()
        mock_bin = self.create_mock_bin()

        file_name = 'mock.binning_table.txt'
        self.create_dummy_bintable(file_name, [mock_bin])

        self.make_directory('mock_bin_files/')
        self.push_to_fasta('mock_bin_files/{}.fasta'.format(mock_bin.name), contigs[0:2], ['aaa', 'ttt'])

        self.db_m.add_bin_table(file_name, 'mock_bin_files/')

        ''' Evaluate '''
        binned_contigs = self.db_m.get_contigs_in_bin(mock_bin.name)

        self.assertEqual( len(binned_contigs), 2 )
        self.assertEqual( binned_contigs, contigs[0:2] )

    # endregion

    # region Auxillary functions

    def test_pull_fasta_metadata(self):

        contig_str = 'ABC DEF GHI'
        seq_name, metadata = self.db_m._pull_fasta_metadata(contig_str)

        self.assertEqual( seq_name, 'ABC' )
        self.assertEqual( metadata, 'DEF GHI' )
    
    def test_pull_fasta_metadata_nospace(self):

        contig_str = 'ABC'
        seq_name, metadata = self.db_m._pull_fasta_metadata(contig_str)

        self.assertEqual( seq_name, 'ABC' )
        self.assertEqual( metadata, '' )

    def test_parse_fasta_nometadata(self):

        ''' Test '''
        file_name = 'mock.contigs.fasta'
        exp_names, _, exp_seqs = self.create_genes_prod_nt(file_name)
        self.push_to_fasta(file_name, exp_names, exp_seqs)

        ''' Evaluate '''
        d = self.db_m._parse_fasta(file_name, True)
        for k, s in zip(exp_names, exp_seqs):

            self.assertTrue( k in d )
            self.assertEqual( d[k][0], s )

    def test_parse_fasta_metadata(self):

        ''' Test '''
        file_name = 'mock.contigs.fasta'
        exp_names, exp_meta, exp_seqs = self.create_genes_prod_nt(file_name)
        self.push_to_fasta(file_name, exp_names, exp_seqs, seq_meta=exp_meta)

        ''' Evaluate '''
        d = self.db_m._parse_fasta(file_name, True)
        for k, m, s in zip(exp_names, exp_meta, exp_seqs):

            self.assertTrue( k in d )

            seq, meta = d[k]
            self.assertEqual( seq, s )
            self.assertEqual( meta, m )

    @unittest.skip('Function moved to DataBaseObjectFactory class, unittests to be moved.')
    def test_calc_gc(self):

        self.assertEqual( self.db_m._calc_gc('ATCGATCG'), 50.0 )

    def test_reverse_complement(self):

        seq = 'ATCGATNG'
        rev = 'CNATCGAT'
        self.assertEqual( self.db_m._reverse_complement(seq), rev )

    def test_read_contig_sequence_keepmetadata(self):

        ''' Test '''
        file_name = 'mock.contigs.fasta'
        exp_names, exp_meta, exp_seqs = self.create_genes_prod_nt(file_name)
        self.push_to_fasta(file_name, exp_names, exp_seqs, seq_meta=exp_meta)

        ''' Evaluate '''
        contig_sequence = self.db_m._read_contig_sequence(file_name, True)
        for exp_contig, exp_met, (obs_contig, obs_meta) in zip(exp_names, exp_meta, contig_sequence):

            self.assertTrue( exp_contig == obs_contig )
            self.assertTrue( exp_met == obs_meta )

    def test_read_contig_sequence_dropmetadata(self):

        ''' Test '''
        file_name = 'mock.contigs.fasta'
        exp_names, exp_meta, exp_seqs = self.create_genes_prod_nt(file_name)
        self.push_to_fasta(file_name, exp_names, exp_seqs, seq_meta=exp_meta)

        ''' Evaluate '''
        contig_sequence = self.db_m._read_contig_sequence(file_name, False)
        for exp_contig, obs_contig in zip(exp_names, contig_sequence):

            self.assertTrue( exp_contig == obs_contig )

    def test_validate_columns(self):

        ''' Test '''
        exp_columns = ['a', 'b', 'c']
        fake_df_buffer = [ {'a': 1, 'b': 2, 'c': 3} ]

        fake_df = pd.DataFrame(fake_df_buffer)

        ''' Evaluate '''
        result, _ = self.db_m._validate_columns(fake_df, exp_columns)
        self.assertTrue(result)

    def test_validate_columns_missingcol(self):

        ''' Test '''
        exp_columns = ['a', 'b', 'c', 'd']
        fake_df_buffer = [ {'a': 1, 'b': 2, 'c': 3} ]

        fake_df = pd.DataFrame(fake_df_buffer)

        ''' Evaluate '''
        result, msg = self.db_m._validate_columns(fake_df, exp_columns)
        self.assertFalse(result)
        self.assertEqual(msg, 'd')

    def test_cast_as_ints(self):

        exp_ints = ( '1', '2', '3' )

        success_state, (obs_1, obs_2, obs_3) = self.db_m._cast_as_ints( exp_ints )

        self.assertTrue( success_state )
        self.assertEqual( int( exp_ints[0] ), obs_1 )
        self.assertEqual( int( exp_ints[1] ), obs_2 )
        self.assertEqual( int( exp_ints[2] ), obs_3 )

    def test_cast_as_ints_fail(self):

        exp_ints = ( '1', '2', 'a' )

        success_state, *obs_vals = self.db_m._cast_as_ints( exp_ints )

        self.assertFalse( success_state )

    def test_cast_as_floats(self):

        exp_floats = ( '1.0', '.2', '3' )

        success_state, (obs_1, obs_2, obs_3) = self.db_m._cast_as_floats( exp_floats )

        self.assertTrue( success_state )
        self.assertEqual( float( exp_floats[0] ), obs_1 )
        self.assertEqual( float( exp_floats[1] ), obs_2 )
        self.assertEqual( float( exp_floats[2] ), obs_3 )

    def test_cast_as_floats_fail(self):

        exp_floats = ( '1.0', '.02', 'a' )

        success_state, *obs_vals = self.db_m._cast_as_floats( exp_floats )

        self.assertFalse( success_state )

    # endregion

if __name__ == '__main__':

    ''' Import the parent path, so that we can import the scripts folder '''
    sys.path.insert(0, '..')
    from scripts.DatabaseManipulator import DatabaseManipulator
    from scripts.DatabaseObjectFactory import Contig, Gene, Annotation, Mapping, Bin

    unittest.main()