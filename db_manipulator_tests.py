'''
    Automated testing of the creation, insertion, update, and delete functionality of the DatabaseManipulator class
'''
import unittest, os

from scripts.DatabaseManipulator import DatabaseManipulator

class TestDatabaseManipulator(unittest.TestCase):

    def setUp(self):

        # Create an in-memory database for testing purposes
        self.db_m = DatabaseManipulator(':memory:')

    # region Test file creation

    def create_contigs_noheader(self, f_name):
        mock_fna = open(f_name, 'w')
        mock_fna.write( '>A\nATCG\n>B\nGCTA\n' )
        mock_fna.close()

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

        self.push_to_fasta(f_name, exp_keys, exp_meta, exp_seq)
        return exp_keys, exp_meta, exp_seq

    def create_genes_prod_aa(self, f_name):

        exp_keys = ['contig1_1', 'contig1_2', 'contig2_1', 'contig3_1']
        exp_meta = ['# 2 # 25 # 1 # ID=1_1;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc_cont=0.50',
                    '# 4 # 46 # 1 # ID=1_1;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc_cont=0.50',
                    '# 5 # 56 # -1 # ID=1_1;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc_cont=0.50',
                    '# 25 # 46 # 1 # ID=1_1;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc_cont=0.50']
        exp_seq = ['SIDRSID*', 'DRSIDRSIDRSIDR*', 'SIDRSIDRSIDRSIDRS', 'RSIDRSI']

        self.push_to_fasta(f_name, exp_keys, exp_meta, exp_seq)
        return exp_keys, exp_meta, exp_seq

    def create_metaxa_ssu(self, f_name):

        exp_keys = ['contig1|B', 'contig2|E']
        exp_meta = ['Predicted Bacterial SSU rRNA (34 bp) From domain V1l to V9r on main strand',
                    'Predicted Eukaryotic SSU rRNA (23 bp) From domain V1l to V9r on complementary strand']
        exp_seq = ['CGATCGATCGATCGATCGATCGATCGATCGATCG', 'TCGATCGATCGATCGATCGATCG']

        self.push_to_fasta(f_name, exp_keys, exp_meta, exp_seq)
        return exp_keys, exp_meta, exp_seq


    def push_to_fasta(self, f_name, seq_names, seq_meta, seq_values):

        mock_fna = open(f_name, 'w')
        for k, m, s in zip(seq_names, seq_meta, seq_values):

            mock_fna.write( '>{} {}\n{}\n'.format(k, m, s) )

        mock_fna.close()

    def remove_temp_file(self, f_name):

        os.remove(f_name)

    # endregion

    # region Create and evaluate database

    def test_open_connection(self):

        # There does not appear to be a way to actually test a connection is open.
        # Most recommended way is to simply execute a query, and in the absence of error, assume the connection is valid.
        cursor = self.db_m.conn.cursor()

        try:

            cursor.execute(''' SELECT count(name) FROM sqlite_master WHERE type='table' AND name='contig' ''')
            self.assertEqual(1, 1, 'Open connection to database')

        except:

            self.assertEqual(0, 1, 'Open connection to database')

    def test_db_creation(self):

        self.db_m.create_blank_database()

        cursor = self.db_m.conn.cursor()

        for tbl in ['contig', 'gene', 'annotation', 'coverage', 'bin']:

            cursor.execute(''' SELECT count(name) FROM sqlite_master WHERE type='table' AND name=? ''', (tbl,) )

            self.assertEqual(cursor.fetchone()[0], 1, 'Missing table: {}'.format(tbl) )

    def test_db_creation_alreadyexists(self):

        self.db_m.create_blank_database()
        i, _ = self.db_m.create_blank_database()

        self.assertEqual(i, -1, 'Database already exists' )

    def test_database_exists_notexists(self):

        ''' As above, but without invoking the create_blank_database() function '''
        self.assertEqual( self.db_m.database_exists(), 0 )

    def test_database_exists_othererror(self):

        ''' Open database connection, create cursor, drop table, test creation '''
        self.db_m.create_blank_database()
        cursor = self.db_m.conn.cursor()

        cursor.execute( "DROP TABLE bin" )

        self.assertEqual( self.db_m.database_exists(), -1 )

    def test_close_connection(self):

        cursor = self.db_m.conn.cursor()

        self.db_m.close_connection()

        try:

            cursor.execute(''' SELECT count(name) FROM sqlite_master WHERE type='table' AND name='contig' ''')
            self.assertEqual(0, 1, 'Close connection to database')

        except:

            self.assertEqual(1, 1, 'Close connection to database')

    # endregion

    # region Insert functions

    def test_add_contigs(self):

        ''' Set up files for test '''
        file_name = 'test_parse_fasta.fasta'
        self.create_contigs_noheader(file_name)

        ''' Test '''
        self.db_m.create_blank_database()
        self.db_m.add_contigs(file_name)

        ''' Evaluate '''
        c = self.db_m.conn.cursor()
        c.execute( ''' SELECT * FROM contig ''' )
        rows = c.fetchall()

        self.assertEqual( len(rows), 2 )
        self.remove_temp_file(file_name)

    def test_add_genes_prodigal_both(self):

        ''' Set up files for test '''
        contigs_file_name = 'test_parse_fasta.fasta'
        self.create_contigs_noheader(contigs_file_name)
        aa_file_name = 'mock.prod.faa'
        aa_keys, aa_meta, aa_seq = self.create_genes_prod_aa(aa_file_name)
        nt_file_name = 'mock.prod.fna'
        *_, nt_seq = self.create_genes_prod_nt(nt_file_name)

        exp_dict = { k: (m, aa, nt) for k, m, aa, nt in zip(aa_keys, aa_meta, aa_seq, nt_seq) }

        ''' Test '''
        self.db_m.create_blank_database()
        self.db_m.add_contigs(contigs_file_name)
        self.db_m.add_genes_prodigal(aa_file=aa_file_name, nt_file=nt_file_name)

        ''' Evaluate '''
        c = self.db_m.conn.cursor()
        c.execute("SELECT * FROM gene;")

        for row in c.fetchall():

            ''' Bundle the database row for comparison '''
            gene_dict = self.db_m._gene_row_to_dict(row)

            self.assertTrue( gene_dict['gene_name'] in exp_dict )

            ''' Bundle the expected values for comparison'''
            exp_meta, exp_aa, exp_nt = exp_dict[ gene_dict['gene_name'] ]
            exp_start, exp_stop, exp_orient = self.db_m._extract_prodigal_metadata(exp_meta)

            for e, o in zip( (exp_aa, exp_nt, exp_start, exp_stop, exp_orient),
                             ['sequence_aa', 'sequence_nt', 'start_pos', 'stop_pos', 'orientation'] ):

                self.assertEqual(e, gene_dict[o], 'Gene {}, column {}'.format(gene_dict['gene_name'], o))

        self.remove_temp_file(contigs_file_name)
        self.remove_temp_file(aa_file_name)
        self.remove_temp_file(nt_file_name)

    def test_add_genes_prodigal_aa_only(self):

        ''' Set up files for test '''
        contigs_file_name = 'test_parse_fasta.fasta'
        self.create_contigs_noheader(contigs_file_name)
        aa_file_name = 'mock.prod.faa'
        aa_keys, aa_meta, aa_seq = self.create_genes_prod_aa(aa_file_name)

        exp_dict = { k: (m, aa) for k, m, aa in zip(aa_keys, aa_meta, aa_seq) }

        ''' Test '''
        self.db_m.create_blank_database()
        self.db_m.add_contigs(contigs_file_name)
        self.db_m.add_genes_prodigal(aa_file=aa_file_name)

        ''' Evaluate '''
        c = self.db_m.conn.cursor()
        c.execute("SELECT * FROM gene;")

        for row in c.fetchall():

            ''' Bundle the database row for comparison '''
            gene_dict = self.db_m._gene_row_to_dict(row)

            self.assertTrue( gene_dict['gene_name'] in exp_dict )

            ''' Bundle the expected values for comparison'''
            exp_meta, exp_aa = exp_dict[ gene_dict['gene_name'] ]
            exp_start, exp_stop, exp_orient = self.db_m._extract_prodigal_metadata(exp_meta)

            for e, o in zip( (exp_aa, '', exp_start, exp_stop, exp_orient),
                             ['sequence_aa', 'sequence_nt', 'start_pos', 'stop_pos', 'orientation'] ):

                self.assertEqual(e, gene_dict[o], 'Gene {}, column {}'.format(gene_dict['gene_name'], o))

        self.remove_temp_file(contigs_file_name)
        self.remove_temp_file(aa_file_name)

    def test_add_genes_prodigal_nt_only(self):

        ''' Set up files for test '''
        contigs_file_name = 'test_parse_fasta.fasta'
        self.create_contigs_noheader(contigs_file_name)
        nt_file_name = 'mock.prod.fna'
        nt_keys, nt_meta, nt_seq = self.create_genes_prod_nt(nt_file_name)

        exp_dict = { k: (m, nt) for k, m, nt in zip(nt_keys, nt_meta, nt_seq) }

        ''' Test '''
        self.db_m.create_blank_database()
        self.db_m.add_contigs(contigs_file_name)
        self.db_m.add_genes_prodigal(nt_file=nt_file_name)

        ''' Evaluate '''
        c = self.db_m.conn.cursor()
        c.execute("SELECT * FROM gene;")

        for row in c.fetchall():

            ''' Bundle the database row for comparison '''
            gene_dict = self.db_m._gene_row_to_dict(row)

            self.assertTrue( gene_dict['gene_name'] in exp_dict )

            ''' Bundle the expected values for comparison'''
            exp_meta, exp_nt = exp_dict[ gene_dict['gene_name'] ]
            exp_start, exp_stop, exp_orient = self.db_m._extract_prodigal_metadata(exp_meta)

            for e, o in zip( ('', exp_nt, exp_start, exp_stop, exp_orient),
                             ['sequence_aa', 'sequence_nt', 'start_pos', 'stop_pos', 'orientation'] ):

                self.assertEqual(e, gene_dict[o], 'Gene {}, column {}'.format(gene_dict['gene_name'], o))

        self.remove_temp_file(contigs_file_name)
        self.remove_temp_file(nt_file_name)

    def test_add_genes_metaxa(self):

        ''' Set up files for test '''
        contigs_file_name = 'test_parse_fasta.fasta'
        self.create_contigs_noheader(contigs_file_name)
        file_name = 'mock.ssu.fna'
        exp_keys, exp_meta, exp_seq = self.create_metaxa_ssu(file_name)

        ''' Test '''
        self.db_m.create_blank_database()
        self.db_m.add_contigs(contigs_file_name)
        self.db_m.add_genes_metaxa(file_name, 'SSU')

        ''' Evaluate '''
        c = self.db_m.conn.cursor()
        c.execute("SELECT * FROM gene;")

        for row in c.fetchall():

            # FUNCTION-IZE THE GENE MAPPING SYSTEM, BECAUSE IT's GOING TO GET A LOT OF USE IN HERE...
    
    # endregion

    # region Aux functions

    def test_parse_fasta_nometadata(self):

        ''' Prepare - File contents and existance are specific to this test, so isolated here '''
        file_name = 'test_parse_fasta.fasta'
        self.create_contigs_noheader(file_name)

        ''' Execute '''
        d = self.db_m._parse_fasta(file_name, False)

        ''' Evaluate '''
        for k, v in zip(['A', 'B'], ['ATCG', 'GCTA']):

            self.assertTrue( k in d )
            self.assertTrue( d[k] == v )

        self.remove_temp_file(file_name)

    def test_parse_fasta_metadata(self):

        file_name = 'test_parse_fasta.fasta'
        exp_keys, exp_meta, exp_seq = self.create_genes_prod_nt(file_name)

        ''' Execute '''
        d = self.db_m._parse_fasta(file_name, True)

        ''' Evaluate '''
        for k, m, s in zip(exp_keys, exp_meta, exp_seq):

            self.assertTrue( k in d )

            seq, meta = d[k]
            self.assertTrue( seq == s )
            self.assertTrue( meta == m )

        self.remove_temp_file(file_name)

    def test_calc_gc(self):

        self.assertEqual( self.db_m._calc_gc('ATCGATCG'), 50.0, '__calc_gc function' )

    def test_extract_prodigal_contig_simple(self):

        contig_name = self.db_m._extract_prodigal_contig('abc_1')
        self.assertEqual(contig_name, 'abc')

    def test_extract_prodigal_contig_complex(self):

        contig_name = self.db_m._extract_prodigal_contig('ab_c_1')
        self.assertEqual(contig_name, 'ab_c')

    def test_extract_prodigal_metadata(self):

        start, stop, orient = self.db_m._extract_prodigal_metadata('# 2 # 25 # -1 # ID=1_1;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc_cont=0.50')
        self.assertEqual(start, 2)
        self.assertEqual(stop, 25)
        self.assertEqual(orient, -1)

    def test_reverse_complement(self):

        seq = 'ATCGATNG'
        rev = 'CNATCGAT'
        self.assertEqual( self.db_m._reverse_complement(seq), rev )

    # endregion

    # region Row casting functions

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

    # endregion

    # region Update statements

    # endregion

    # region Delete statements

    # endregion

if __name__ == '__main__':
    unittest.main()