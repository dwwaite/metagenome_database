'''
    Backend for processing the database

    Useful query scripts

import sqlite3
con = sqlite3.connect('asd.db')
cursor = con.cursor()

# GENE BASED
cursor.execute("SELECT * FROM gene;")

for row in cursor.fetchall():
    print(row)
    _ = cursor.execute("SELECT * FROM contig WHERE contigname=?;", (row[1],))
    cursor.fetchall()
    print()

# CONTIG BASED
cursor.execute("SELECT * FROM contig;")
for row in cursor.fetchall():
    print(row)
    _ = cursor.execute("SELECT * FROM gene WHERE contigname=?;", (row[0],))
    cursor.fetchall()
    print()
'''
import sqlite3, sys
from sqlite3 import Error
from collections import Counter

# region Front end handling and option validation

class DatabaseManipulator():

    def __init__(self, database_name):
        self.database_name = database_name

        self.open_connection()

    def create_blank_database(self):

        '''
            Attempt to create a blank database schema in the connection.
            If the database is already formatted, abort.

            Using the self.conn object to create cursors, all SQL statements are created within this function
        '''

        if self.database_exists():

            return -1, 'Database file {} is already populated. Did you mean to invoke the update command?'.format(self.database_name)

        try:

            ''' Create a cursor and populate the expected tables '''
            c = self.conn.cursor()

            c.execute('''CREATE TABLE contig (contigname TEXT PRIMARY KEY,
                                              sequence TEXT,
                                              length INTEGER,
                                              gc FLOAT);''')

            c.execute('''CREATE TABLE gene (genename TEXT PRIMARY KEY,
                                            contigname TEXT,
                                            sequence_nt TEXT,
                                            sequence_aa TEXT,
                                            length_nt INTEGER,
                                            length_aa INTEGER,
                                            start INTEGER,
                                            stop INTEGER,
                                            orientation INTEGER,
                                            FOREIGN KEY(contigname) REFERENCES contig(contigname));''')

            c.execute('''CREATE TABLE annotation (annotationid INTEGER PRIMARY KEY,
                                                  genename TEXT,
                                                  method TEXT,
                                                  annotation_schema TEXT,
                                                  accession TEXT,
                                                  description TEXT,
                                                  iscurrent BOOLEAN,
                                                  FOREIGN KEY(genename) REFERENCES gene(genename));''')

            c.execute('''CREATE TABLE coverage (coverageid INTEGER PRIMARY KEY,
                                               contigname TEXT,
                                               sample TEXT,
                                               readsmapped INTEGER,
                                               FOREIGN KEY(contigname) REFERENCES contig(contigname));''')

            c.execute('''CREATE TABLE bin (binid INTEGER PRIMARY KEY,
                                           contigname TEXT,
                                           binname TEXT,
                                           taxonomy TEXT,
                                           version INTEGER,
                                           FOREIGN KEY(contigname) REFERENCES contig(contigname));''')

            self.conn.commit()

            return 1, None

        except Error as e:

            return -1, e

    #region Connection handling

    def open_connection(self):

        try:
    
            self.conn = sqlite3.connect( self.database_name )

        except Error as e:
            print(e)

    def database_exists(self):

        ''' Count each of the expected tables, then evaluate if the number of successful observation is correct '''
        obs_tables = 0
        tables_defined = ['contig', 'gene', 'annotation', 'coverage', 'bin']
        exp_tables = len( tables_defined )

        c = self.conn.cursor()
    
        for tbl in tables_defined:

            c.execute(''' SELECT count(name) FROM sqlite_master WHERE type='table' AND name=? ''', (tbl,) )

            if c.fetchone()[0]: obs_tables += 1

        ''' Evaluate the outcome
            1. Obs == Exp, database exists, return True
            2. Obs == 0, database does not exist, return False
            3. Obs > 0, Obs != Exp, tables exist, but not as many as expected. Something has gone wrong.
        '''

        if obs_tables == exp_tables:

            return 1

        elif obs_tables == 0:

            return 0

        else:
       
            return -1
    
    def close_connection(self):

        self.conn.close()

    # endregion

    # region Record entry

    def add_contigs(self, contig_file):

        c = self.conn.cursor()

        for contig, sequence in self._parse_fasta(contig_file, False).items():

            ''' Polish the contig name '''
            if ' ' in contig:
                contig = contig.split(' ')[0]

            seq_len = len(sequence)
            seq_gc = self._calc_gc(sequence)

            c.execute('INSERT INTO contig(contigname, sequence, length, gc) VALUES(?, ?, ?, ?)', (contig, sequence, seq_len, seq_gc) )

        self.conn.commit()

    def _add_gene(self, gene_name, contig_name, sequence_nt, sequence_aa, length_nt, length_aa, start_pos, stop_pos, orientation):

        c = self.conn.cursor()

        try:

            c.execute( ''' INSERT INTO gene(genename, contigname, sequence_nt, sequence_aa, length_nt, length_aa, start, stop, orientation) VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?) ''',
                       (gene_name, contig_name, sequence_nt, sequence_aa, length_nt, length_aa, start_pos, stop_pos, orientation) )
            self.conn.commit()

        except sqlite3.Error as e:

            print(e)
            print( 'Did not add entry {}'.format(gene_name) )

    def _extract_prodigal_contig(self, gene_name):
        contig_pieces = gene_name.split('_')[0:-1]
        return '_'.join(contig_pieces)

    def _extract_prodigal_metadata(self, metadata):
        start, stop, orient, *remainder = metadata.split('#')[1:]
        return int(start), int(stop), int(orient)

    def add_genes_prodigal(self, aa_file=None, nt_file=None):

        ''' To begin, I need to line up amino acid and nucleotide predictions if both are provided
            Achieve this by creating empty dictionaries, then replacing them with real ones if requested.
            Build a list of all keys, then iterate through and extract sequences and required '''
        aa_map = {}
        if aa_file:
            aa_map = self._parse_fasta(aa_file, True)

        nt_map = {}
        if nt_file:
            nt_map = self._parse_fasta(nt_file, True)

        all_keys = set( aa_map.keys() ) | set( nt_map.keys() )

        for key in all_keys:

            aa = ''; nt = ''; meta = ''
            if key in aa_map:
                aa, meta = aa_map[key]
            if key in nt_map:
                nt, meta = nt_map[key]

            contig_name = self._extract_prodigal_contig(key)

            start, stop, orientation = self._extract_prodigal_metadata(meta)
            self._add_gene(key, contig_name, nt, aa, len(nt), len(aa), start, stop, orientation)

    def  add_genes_metaxa():

        return 0

    def add_genes_aragorn():
    
        return 0

    # endregion

    # region Record removal

    # endregion

    # region Retrieval

    '''
        To begin, created a range of functions that simply take database rows and cast them into dicts. This isn't groupbreaking, but helps to keep
        track of what variables I have at other parts of the program, and in unit testing.
    '''

    def _gene_row_to_dict(self, _row):
        return { n: r for n, r in zip( ['gene_name', 'contig_name', 'sequence_nt', 'sequence_aa', 'length_nt', 'length_aa', 'start_pos', 'stop_pos', 'orientation'],
                                       _row ) }

    def _contig_row_to_dict(self, _row):
        return { n: r for n, r in zip( ['contig_name', 'sequence', 'length', 'gc'], _row ) }

    # endregion

    # region File handling and sequence manupulation

    def _parse_fasta(self, f_name, keep_metadata):

        fna_dict = {}

        ''' Not really necessary to handle this way, but unittests flag open file reader warnings if I don't explicitly close '''
        fna_reader = open(f_name, 'r')
        fna_content = fna_reader.read()
        fna_reader.close()

        for entry in fna_content.split('>')[1:]:

            header, *seq = entry.split('\n')
            seq = ''.join(seq)

            if keep_metadata:

                seq_name, *metadata = header.split(' ')
                metadata = ' '.join(metadata)

                fna_dict[seq_name] = (seq, metadata)

            else:

                fna_dict[header] = seq

        return fna_dict

    def _calc_gc(self, sequence):

        nuclotide_counter = Counter(sequence)

        gc_frac = float( nuclotide_counter['G'] + nuclotide_counter['C'] ) / sum( nuclotide_counter.values() )
        return gc_frac * 100

    # endregion
