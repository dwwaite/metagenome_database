import sqlite3, sys, os, glob
import pandas as pd
import dateutil.parser as date_parser

from datetime import datetime
from sqlite3 import Error

from scripts.DatabaseObjectFactory import Contig, Gene, Annotation, Mapping, Bin

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
                                            prediction_tool TEXT,
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
                                                  annotation_db TEXT,
                                                  accession TEXT,
                                                  description TEXT,
                                                  identity FLOAT,
                                                  coverage FLOAT,
                                                  evalue FLOAT,
                                                  iscurrent BOOLEAN,
                                                  FOREIGN KEY(genename) REFERENCES gene(genename));''')

            c.execute('''CREATE TABLE coverage (coverageid INTEGER PRIMARY KEY,
                                               contigname TEXT,
                                               sample TEXT,
                                               readsmapped INTEGER,
                                               FOREIGN KEY(contigname) REFERENCES contig(contigname));''')

            c.execute('''CREATE TABLE transcript (transcriptid INTEGER PRIMARY KEY,
                                               genename TEXT,
                                               sample TEXT,
                                               readsmapped INTEGER,
                                               FOREIGN KEY(genename) REFERENCES gene(genename));''')

            c.execute('''CREATE TABLE bin (binname TEXT PRIMARY KEY,
                                           taxonomy TEXT,
                                           completeness FLOAT,
                                           contamination FLOAT,
                                           strainheterogeneity FLOAT,
                                           version TEXT,
                                           datecreated TEXT );''')

            c.execute('''CREATE TABLE bin_contig (bcid INTEGER PRIMARY KEY,
                                                  contigname TEXT,
                                                  binname TEXT,
                                                  FOREIGN KEY(contigname) REFERENCES contig(contigname),
                                                  FOREIGN KEY(binname) REFERENCES bin(binname));''')

            self.conn.commit()

            return 1, None

        except Error as e:

            return -1, e

    #region Connection handling

    def open_connection(self):

        try:
    
            self.conn = sqlite3.connect( self.database_name )
            self.conn.executescript( ''' pragma foreign_keys=on; ''' )

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

    # region Record entry - contigs

    def add_contigs(self, contig_file):

        entry_valid, contig_data = self._validate_contigs(contig_file)

        if entry_valid:

            c = self.conn.cursor()
            c.executemany('INSERT INTO contig(contigname, sequence, length, gc) VALUES(?, ?, ?, ?)',
                          [ (c.name, c.sequence, c.length, c.gc) for c in contig_data ] )

        else:
            raise Exception('Duplicate contig detected ({}). Aborting...'.format(contig_data) )

    def _obtain_contig_set(self):

        c = self.conn.cursor()
        c.execute( '''SELECT contigname FROM contig''' )
        return set( [ row[0] for row in c.fetchall() ] )

    def _validate_contigs(self, contig_file):

        contig_buffer = []
        current_contig_set = self._obtain_contig_set()

        success_state = True

        for name, sequence in self._parse_fasta(contig_file, False).items():

            ''' Polish the contig name '''
            if ' ' in name:
                name = name.split(' ')[0]

            if name in current_contig_set:
                success_state = False

            contig_buffer.append( Contig(name, sequence) )

        return success_state, contig_buffer

    # endregion

    # region Record entry - genes (general purpose)

    def _obtain_gene_set(self):

        c = self.conn.cursor()
        c.execute( '''SELECT genename FROM gene''' )
        return set( [ row[0] for row in c.fetchall() ] )

    def _validate_genes(self, gene_list):

        current_contigs = self._obtain_contig_set()
        current_genes = self._obtain_gene_set()

        exp_genes = set( [ gene.name for gene in gene_list ] )
        if len(exp_genes) < len(gene_list):
            raise Exception( 'There are duplicate genes in the input file. Aborting...' )

        for gene in gene_list:

            ints_are_valid, _ = self._cast_as_ints( (gene.length_nt, gene.length_aa, gene.start, gene.stop, gene.orientation) )
            if not ints_are_valid:
                raise Exception( 'Error converting values to int in gene {}. Aborting...'.format(gene.name) )

            if not gene.contig_name in current_contigs:
                raise Exception( 'Gene {} is not linked to an existing contig. Aborting...'.format(gene.name) )

            if gene.name in current_genes:
                raise Exception( 'Gene {} is already in the database. Aborting...'.format(gene.name) )

        return True

    def _add_genes(self, gene_list):

        c = self.conn.cursor()

        try:

            c.executemany( ''' INSERT INTO gene(genename, contigname, prediction_tool, sequence_nt, sequence_aa, length_nt, length_aa, start, stop, orientation) VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?) ''',
                           [ (g.name, g.contig_name, g.prediction_tool, g.sequence_nt, g.sequence_aa, g.length_nt, g.length_aa, g.start, g.stop, g.orientation) for g in gene_list ])
            self.conn.commit()

        except sqlite3.Error as e:
            print(e)

    # endregion
    
    # region Record entry - genes (Prodigal files)

    def _extract_prodigal_contig(self, gene_name):
        contig_pieces = gene_name.split('_')[0:-1]
        return '_'.join(contig_pieces)

    def _extract_prodigal_metadata(self, metadata):
        start, stop, orient, *remainder = metadata.split('#')[1:]
        return start.strip(), stop.strip(), orient.strip()

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

        ''' Once the gene sequences are mapped, iterate the all_keys set and create gene tuples for each entry,
            filling in blanks as needed. '''
        gene_buffer = []

        for key in all_keys:

            aa = ''; nt = ''; meta = ''
            if key in aa_map: aa, meta = aa_map[key]
            if key in nt_map: nt, meta = nt_map[key]

            contig_name = self._extract_prodigal_contig(key)
            start, stop, orientation = self._extract_prodigal_metadata(meta)

            prod_gene = Gene(gene_name=key, contig_name=contig_name, prediction_tool='Prodigal',
                             sequence_nt=nt, sequence_aa=aa,
                             start=start, stop=stop, orientation=orientation)
            gene_buffer.append( prod_gene )

        ''' Validate the genes, if they pass, add them.  '''
        genes_are_valid = self._validate_genes( gene_buffer )

        if genes_are_valid:
            self._add_genes( gene_buffer )


    # endregion

    # region Record entry - genes (MeTaxa files)

    def _metaxa_extract_names(self, contig_string, gene_type):

        ''' MeTaxa2 uses a | to denote the start of metadata, so my standard split won't work.
                However, it only reports a single character before the first space, which is a single letter marking which taxa the prediction occured from '''

        *contig_name, taxa_flag = contig_string.split('|')
        contig_name = '|'.join(contig_name)

        gene_name = '{}_{}_{}'.format(contig_name, taxa_flag, gene_type)
        return contig_name, gene_name

    def _metaxa_extract_orientation(self, metadata):
        return 1 if metadata.split(' ')[-2] == 'main' else -1

    def _metaxa_locate_coords(self, parent_contig_name, rrna_sequence, orientation):

        parent_contig = self._get_contig_by_name(parent_contig_name)

        if parent_contig and orientation == 1:

                start = parent_contig['sequence'].find(rrna_sequence)
                stop = start + len(rrna_sequence)
                return start + 1, stop

        elif parent_contig and orientation == -1:

                rev_parent = self._reverse_complement( parent_contig['sequence'] )
                start = len(rev_parent) - rev_parent.find(rrna_sequence)
                stop = start - len(rrna_sequence)
                return start, stop + 1

        return -1, -1

    def add_genes_metaxa(self, rrna_file, gene_type):

        ''' To begin, parse the file into a list of gene tuples. Create a matching annotation to go with each rRNA gene '''

        rna_predictions = self._parse_fasta(rrna_file, True)

        gene_buffer = []
        annotation_buffer = []

        for contig, (sequence, metadata) in rna_predictions.items():

            contig_name, gene_name = self._metaxa_extract_names(contig, gene_type)
            
            ''' Find the start/stop/orientation values for the sequence '''
            orientation = self._metaxa_extract_orientation(metadata)
            start, stop = self._metaxa_locate_coords(contig_name, sequence.upper(), orientation)

            metaxa_gene = Gene(gene_name=gene_name, contig_name=contig_name, prediction_tool='MeTaxa2',
                               sequence_nt=sequence, sequence_aa='',
                               start=start, stop=stop, orientation=orientation)

            metaxa_annotation = Annotation(gene_name=gene_name, method='MeTaxa2', annotation_database=gene_type, accession='-',
                                            identity=-1.0, coverage=-1.0, evalue=-1.0, description=metadata)

            gene_buffer.append( metaxa_gene )
            annotation_buffer.append( metaxa_annotation )

        ''' Validate the genes, if they pass, add them. '''
        genes_are_valid = self._validate_genes( gene_buffer )

        if genes_are_valid:

            self._add_genes( gene_buffer )
            self._add_annotations( annotation_buffer )

    # endregion

    # region Record entry - genes (Aragorn files)

    def _aragorn_extract_names(self, trna_sequence, metadata, contig_i, gene_i):

        contig_name = trna_sequence[ contig_i - 1 ]
        trna_type = self._aragorn_extract_type(metadata)
        gene_name = '{}_{}_{}'.format(contig_name, gene_i, trna_type)
        
        return contig_name, gene_name

    def _aragorn_split_indicies(self, contig_str):

        ''' Contig_str is a string of the form 1-1 '''
        contig_i, gene_i = contig_str.split('-')
        return int(contig_i), int(gene_i)

    def _aragorn_extract_type(self, trna_metadata):

        ''' Metadata is a string of the form "tRNA-Arg(gcg) [58280,58368]" '''
        trna_type = trna_metadata.split(' ')[0]
        return trna_type.replace('(', '-').replace(')', '')

    def _aragorn_extract_coords(self, trna_metadata):

        ''' Metadata is a string of the form "tRNA-Arg(gcg) [58280,58368]" '''
        trna_coords = trna_metadata.split(' ')[-1]
        start_str, stop_str = trna_coords.split(',')

        return start_str.replace('[', ''), stop_str.replace(']', '')

    def _aragorn_determine_orientation(self, start_str, stop_str):

        ints_are_valid, (start, stop) = self._cast_as_ints( (start_str, stop_str) )

        if ints_are_valid:
                return 1 if stop > start else -1
        else:
            return 'not numeric'

    def add_genes_aragorn(self, trna_contigs, trna_predictions):
    
        ''' Aragorn is a bit tricky, predictions are named according to the order of input contigs, but are relabeled using
            [contig#]-[gene#] pairs.
            This means I need an ordered copy of the contigs file used during prediction '''
        trna_sequence = self._read_contig_sequence(trna_contigs, False)
        trna_predictions = self._parse_fasta(trna_predictions, True)

        ''' Create a buffer for gene tuples, then iterate the files and extract necessary information '''
        gene_buffer = []
        annotation_buffer = []

        for aragorn_contig, (sequence, metadata) in trna_predictions.items():

            ''' Translate the uninformative Aragorn gene name into the more useful format [Contig]_[Gene number]_[tRNA type] '''
            contig_i, gene_i = self._aragorn_split_indicies(aragorn_contig)

            contig_name, gene_name = self._aragorn_extract_names(trna_sequence, metadata, contig_i, gene_i)
            
            ''' Compute the coordinates and orientation of the tRNA prediction '''
            start, stop = self._aragorn_extract_coords(metadata)
            orientation = self._aragorn_determine_orientation(start, stop)

            aragorn_gene = Gene(gene_name=gene_name, contig_name=contig_name, prediction_tool='Aragorn',
                                sequence_nt=sequence, sequence_aa='',
                                start=start, stop=stop, orientation=orientation)

            aragorn_annotation = Annotation(gene_name=gene_name, method='Aragorn', annotation_database='tRNA prediction', accession='-',
                                            identity=-1.0, coverage=-1.0, evalue=-1.0, description=metadata)

            gene_buffer.append( aragorn_gene )
            annotation_buffer.append( aragorn_annotation )

        ''' Validate the genes, if they pass, add them. '''
        genes_are_valid = self._validate_genes( gene_buffer )

        if genes_are_valid:

            self._add_genes( gene_buffer )
            self._add_annotations( annotation_buffer )

    # endregion

    # region Record entry - annotations

    def _set_annotations_historic(self, annotatation_list):

        c = self.conn.cursor()

        for annotation in annotatation_list:
            c.execute('''UPDATE annotation SET iscurrent = 0 WHERE genename = ? AND method = ? AND annotation_db = ?''',
                      (annotation.gene_name, annotation.method, annotation.annotation_database) )

        self.conn.commit()

    def _validate_annotations(self, annotatation_list):

        current_genes = self._obtain_gene_set()

        for annotation in annotatation_list:

            floats_are_valid, _ = self._cast_as_floats( (annotation.identity, annotation.coverage, annotation.evalue) )
            if not floats_are_valid:
                raise Exception( 'Error converting values to float in gene {}. Aborting...'.format(annotation.gene_name) )

            if not annotation.gene_name in current_genes:
                raise Exception( 'Annotation entry (method {}, database {}) is not matched to a valid gene ({}). Aborting...'.format(annotation.method, annotation.annotation_database, annotation.gene_name) )

        return True

    def _add_annotations(self, annotatation_list):

        ''' Check the database for an annotation with the same gene/method/db '''
        self._set_annotations_historic(annotatation_list)
            
        c = self.conn.cursor()

        try:

            c.executemany( ''' INSERT INTO annotation(genename, method, annotation_db, accession, identity, coverage, evalue, description, iscurrent) VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?) ''',
                           [ (a.gene_name, a.method, a.annotation_database, a.accession, a.identity, a.coverage, a.evalue, a.description, a.is_current) for a in annotatation_list ] )

        except sqlite3.Error as e:
            print(e)

    def add_annotations_by_table(self, annotation_table):

        ''' Verify the necessary columns are present '''
        anntotation_df = pd.read_csv(annotation_table, sep='\t')
        table_is_valid, col_err = self._validate_columns(anntotation_df, ['gene', 'method', 'database', 'hit', 'identity', 'coverage', 'evalue', 'description'])

        ''' If there is a problem with the table, report it and abort the process. '''
        if not table_is_valid:
            raise Exception( 'Unable to find column {} in table {}. Aborting...'.format(col_err, annotation_table) )

        ''' Create annotation tuples from the DataFrame rows '''
        annotation_buffer = []
        for _, row in anntotation_df.iterrows():

            annotation = Annotation(gene_name=row['gene'], method=row['method'], annotation_database=row['database'], accession=row['hit'],
                                    identity=row['identity'], coverage=row['coverage'], evalue=row['evalue'], description=row['description'])
            annotation_buffer.append( annotation )

        annotations_are_valid = self._validate_annotations( annotation_buffer )

        if annotations_are_valid:
            self._add_annotations( annotation_buffer )

    # endregion
    
    # region Record entry - coverage

    def _validate_mapping(self, mapping_list, target):

        if target == 'contig':
            current_targets = self._obtain_contig_set()
        elif target == 'gene':
            current_targets = self._obtain_gene_set()

        for mapping_record in mapping_list:

            ints_are_valid, _ = self._cast_as_ints( (mapping_record.reads_mapped, ) )
            if not ints_are_valid:
                raise Exception( 'Error converting {} mapping value to int in {}. Aborting...'.format(target, mapping_record.target_name) )

            if not mapping_record.target_name in current_targets:
                raise Exception( 'Record {} is not linked to an existing {}. Aborting...'.format(mapping_record.target_name, target) )

        return True

    def _add_coverages(self, mapping_list):

        c = self.conn.cursor()

        try:
            c.executemany('''INSERT INTO coverage(contigname, sample, readsmapped) VALUES(?, ?, ?)''',
                          [ (m.target_name, m.sample_name, m.reads_mapped) for m in mapping_list ] )
            self.conn.commit()

        except sqlite3.Error as e:
            print(e)

    def add_coverage_table(self, coverage_table):

        ''' There should be some validation of the table parameters prior to getting here.
            Read assumption is that the first column is contig name, all subsequence columns are sample observations '''

        coverage_df = pd.read_csv(coverage_table, sep='\t')

        contig_col = coverage_df.columns[0]
        sample_names = coverage_df.columns[1:]

        mapping_list = []

        for _, row in coverage_df.iterrows():

            for col_name in sample_names:

                mapping_record = Mapping(target_name=row[contig_col], sample_name=col_name, reads_mapped=row[col_name])
                mapping_list.append( mapping_record )

        coverages_are_valid = self._validate_mapping( mapping_list, target='contig' )
        
        if coverages_are_valid:
            self._add_coverages( mapping_list )

    # endregion

    # region Record entry - transcription

    def _add_transcript(self, mapping_list):

        c = self.conn.cursor()

        try:
            c.executemany('''INSERT INTO transcript(genename, sample, readsmapped) VALUES(?, ?, ?)''',
                          [ (m.target_name, m.sample_name, m.reads_mapped) for m in mapping_list ] )
            self.conn.commit()

        except sqlite3.Error as e:
            print(e)

    def add_transcript_table(self, transcript_table):

        ''' There should be some validation of the table parameters prior to getting here.
            Read assumption is that the first column is gene name, all subsequence columns are sample observations '''
        transcript_df = pd.read_csv(transcript_table, sep='\t')

        gene_col = transcript_df.columns[0]
        sample_names = transcript_df.columns[1:]

        mapping_list = []

        for _, row in transcript_df.iterrows():

            for gene_name in sample_names:

                mapping_record = Mapping(target_name=row[gene_col], sample_name=gene_name, reads_mapped=row[gene_name])
                mapping_list.append( mapping_record )

        transcripts_are_valid = self._validate_mapping( mapping_list, target='gene' )
        
        if transcripts_are_valid:
            self._add_transcript( mapping_list )

    # endregion

    # region Record entry - bins and contig associations

    def _obtain_bin_set(self):

        c = self.conn.cursor()
        c.execute( '''SELECT binname FROM bin''' )
        return set( [ row[0] for row in c.fetchall() ] )

    def _validate_bins(self, bin_list):

        current_bins = self._obtain_bin_set()

        ''' Check for duplicates in the input file '''
        exp_bins = set( [ b.name for b in bin_list ] )
        if len(exp_bins) < len(bin_list):
            raise Exception( 'There are duplicate bins in the input file. Aborting...' )

        for b in bin_list:

            floats_are_valid, _ = self._cast_as_floats( (b.completeness, b.contamination, b.strain_heterogeneity) )
            if not floats_are_valid:
                raise Exception( 'Error converting values to float in bin {}. Aborting...'.format(b.name) )
            
            if b.name in current_bins:
                raise Exception( 'Duplicate entry for bin {} detected. Aborting...'.format(b.name) )

        return True

    def _validate_contig_binning(self, bin_list, contig_binning_dict):

        current_contigs = self._obtain_contig_set()

        for b in bin_list:

            contig_selection = set( contig_binning_dict[ b.name ] )

            invalid_contigs = contig_selection - current_contigs
            if len(invalid_contigs) > 0:
                raise Exception( 'The contig {} in bin {} does not exist in the database. Aborting...'.format( ','.join(invalid_contigs), b.name) )

        return True

    def _index_contig_folder(self, contig_folder):

        contig_binning_dict = {}

        for contig_file in glob.glob( contig_folder + '/*'):

            bin_file = os.path.split(contig_file)[-1]
            bin_name = os.path.splitext(bin_file)[0]

            contig_binning_dict[bin_name] = self._read_contig_sequence(contig_file, False)

        return contig_binning_dict

    def add_bin_table(self, bin_table, contig_folder):

        '''
            There should be some validation of the table parameters prior to getting here.
            Format assumption is that the table contains the following columns:
                1. Bin name (Bin)
                2. Bin taxonomy (Taxonomy)
                3. Completeness estimator, CheckM or similar (Completeness)
                4. Contamination estimator, CheckM or similar (Contamination)
                5. Strain heterogeneity estimator, CheckM or similar (Heterogeneity)

            Version is taken as the bin_table file name with the extension removed.
            Date created is calcualted by a system call at the time of execution.
        '''
        bin_df = pd.read_csv(bin_table, sep='\t')

        binning_version = os.path.split( bin_table )[-1]
        binning_version = os.path.splitext(binning_version)[0]

        bin_timestamp = datetime.now().isoformat()

        bin_buffer = []
        for _, row in bin_df.iterrows():

            current_bin = Bin(bin_name=row['Bin'], taxonomy=row['Taxonomy'],
                              completeness=row['Completeness'], contamination=row['Contamination'], strain_heterogeneity=row['Heterogeneity'],
                              version_name=binning_version, date_created=bin_timestamp)
            bin_buffer.append( current_bin )

        contig_binning_dict = self._index_contig_folder(contig_folder)

        ''' Validate '''
        success_state_bin = self._validate_bins(bin_buffer)
        success_state_contig = self._validate_contig_binning(bin_buffer, contig_binning_dict)
        
        if success_state_bin and success_state_contig:
            self._add_bins(bin_buffer)
            self._add_contigs_to_bins(contig_binning_dict)
    
    def _add_bins(self, bin_buffer):

        c = self.conn.cursor()

        try:
            c.executemany('''INSERT INTO bin(binname, taxonomy, completeness, contamination, strainheterogeneity, version, datecreated) VALUES(?, ?, ?, ?, ?, ?, ?)''',
                          [ (b.name, b.taxonomy, b.completeness, b.contamination, b.strain_heterogeneity, b.version_name, b.date_created) for b in bin_buffer ] )
            self.conn.commit()

        except sqlite3.Error as e:
            print(e)

    def _add_contigs_to_bins(self, contig_binning_dict):

        c = self.conn.cursor()

        ''' Unpack the dictionary into a list of tuples '''
        contig_associations = []
        for bin_name, contig_set in contig_binning_dict.items():
            contig_associations.extend( [ (bin_name, contig) for contig in contig_set ] )

        try:
            c.executemany('''INSERT INTO bin_contig(binname, contigname) VALUES(?, ?)''', contig_associations )
            self.conn.commit()

        except sqlite3.Error as e:
            print(e)

    # endregion

    # region Record retrieval - contigs

    def _contig_row_to_dict(self, _row):
        return { n: r for n, r in zip( ['contig_name', 'sequence', 'length', 'gc'], _row ) }

    def get_contigs(self):

        ''' Returns a DataFrame of contigs '''

        ''' Create cursor, execute command and retrieve rows '''
        c = self.conn.cursor()
        c.execute('''SELECT contigname, sequence, length, gc FROM contig''')

        frame_buffer = [ self._contig_row_to_dict(row) for row in c.fetchall() ]

        return pd.DataFrame(frame_buffer)

    def _get_contig_by_name(self, contig_name):

            ''' Returns a dict of contig data for a single contig in the database, or None if the entry is not observed '''

            ''' Create cursor, execute command and retrieve the row '''
            c = self.conn.cursor()
            c.execute('''SELECT contigname, sequence, length, gc FROM contig WHERE contigname = ?''', (contig_name,) )

            ''' Parse the results '''
            result = c.fetchone()

            if result:
                return self._contig_row_to_dict(result)

            else:
                return None

    # endregion

    # region Record retrieval - genes

    def _gene_row_to_dict(self, _row):
        return { n: r for n, r in zip( ['gene_name', 'contig_name', 'prediction_tool', 'sequence_nt', 'sequence_aa', 'length_nt', 'length_aa', 'start_pos', 'stop_pos', 'orientation'],
                                       _row ) }

    def get_genes(self, prediction_method=None):

        ''' Returns a DataFrame of all genes, optionally filtered by a prediction method '''

        ''' Create cursor, execute command and retrieve rows '''
        c = self.conn.cursor()
        
        if prediction_method:

            c.execute('''SELECT genename, contigname, prediction_tool, sequence_nt, sequence_aa, length_nt, length_aa, start, stop, orientation FROM gene WHERE prediction_tool = ?''', (prediction_method,) )
        
        else: 

            c.execute('''SELECT genename, contigname, prediction_tool, sequence_nt, sequence_aa, length_nt, length_aa, start, stop, orientation FROM gene''')

        frame_buffer = [ self._gene_row_to_dict(row) for row in c.fetchall() ]

        return pd.DataFrame(frame_buffer)

    def get_gene_by_name(self, gene_name):

        ''' Returns a dict of contig data for a single contig in the database, or None if the entry is not observed '''

        ''' Create cursor, execute command and retrieve the row '''
        c = self.conn.cursor()
        c.execute('''SELECT genename, contigname, prediction_tool, sequence_nt, sequence_aa, length_nt, length_aa, start, stop, orientation FROM gene WHERE genename = ?''', (gene_name,) )

        result = c.fetchone()

        if result:
            return self._gene_row_to_dict(result)

        else:
            return None

    def get_genes_by_contig(self, contig_name):

        ''' Returns a DataFrame of genes predicted on a user-specified  contig in the database '''

        ''' Create cursor, execute command and retrieve the row '''
        c = self.conn.cursor()
        c.execute('''SELECT genename, contigname, prediction_tool, sequence_nt, sequence_aa, length_nt, length_aa, start, stop, orientation FROM gene WHERE contigname = ?''', (contig_name,) )

        frame_buffer = [ self._gene_row_to_dict(row) for row in c.fetchall() ]
        return pd.DataFrame(frame_buffer)

    # endregion

    # region Record retrieval - annotations

    def _annotation_row_to_dict(self, _row):
        return { n: r for n, r in zip( ['gene_name', 'method', 'annotation_db', 'accession', 'identity', 'coverage', 'evalue', 'description', 'iscurrent'], _row ) }

    def _annotations_are_current(self, old_annotations):
        if old_annotations is None:
            return True
        
        else:
            return old_annotations

    def get_annotations(self, old_annotations=None):

        ''' Returns a DataFrame of annotation data. By default, only current annotations are returned '''

        ''' Create cursor, execute command and retrieve rows '''
        c = self.conn.cursor()

        current_gene_only = self._annotations_are_current(old_annotations)
        if current_gene_only:
            c.execute('''SELECT genename, method, annotation_db, accession, identity, coverage, evalue, description, iscurrent FROM annotation WHERE iscurrent = ?''', (current_gene_only,) )

        else:
            c.execute('''SELECT genename, method, annotation_db, accession, identity, coverage, evalue, description, iscurrent FROM annotation''')

        frame_buffer = [ self._annotation_row_to_dict(row) for row in c.fetchall() ]
        return pd.DataFrame(frame_buffer)

    def get_annotations_by_gene(self, gene_name, old_annotations=None):

        ''' Returns a DataFrame of annotation data for a single gene in the database. By default, only the current annotation is returned '''

        ''' Create cursor, execute command and retrieve the row '''
        c = self.conn.cursor()

        current_gene_only = self._annotations_are_current(old_annotations)
        if current_gene_only:
            c.execute('''SELECT genename, method, annotation_db, accession, description, iscurrent FROM annotation WHERE genename = ? AND iscurrent = ?''', (gene_name, current_gene_only) )

        else:
            c.execute('''SELECT genename, method, annotation_db, accession, description, iscurrent FROM annotation WHERE genename = ?''', (gene_name,) )

        frame_buffer = [ self._annotation_row_to_dict(row) for row in c.fetchall() ]
        return pd.DataFrame(frame_buffer)

    def get_annotation_by_gene_method(self, gene_name, method, old_annotations=None):

        ''' Returns a DataFrame of annotation data for a single gene/method in the database. By default, only the current annotation is returned '''
        gene_df = self.get_annotations_by_gene(gene_name, old_annotations)
        return gene_df[ gene_df.method == method ]

    # endregion

    # region Record retrieval - coverage

    def _coverage_row_to_dict(self, _row):
        return { n: r for n, r in zip( ['contig_name', 'sample', 'reads_mapped'], _row ) }

    def get_coverage(self):

        ''' Returns a DataFrame of read mapping values '''

        ''' Create cursor, execute command and retrieve rows '''
        c = self.conn.cursor()
        c.execute('''SELECT contigname, sample, readsmapped FROM coverage''')

        frame_buffer = [ self._coverage_row_to_dict(row) for row in c.fetchall() ]
        return pd.DataFrame(frame_buffer)

    # endregion

    # region Record retrieval - transcript

    def _transcript_row_to_dict(self, _row):
        return { n: r for n, r in zip( ['gene_name', 'sample', 'reads_mapped'], _row ) }

    def get_transcript(self):

        ''' Returns a DataFrame of read mapping values '''

        ''' Create cursor, execute command and retrieve rows '''
        c = self.conn.cursor()
        c.execute('''SELECT genename, sample, readsmapped FROM transcript''')

        frame_buffer = [ self._transcript_row_to_dict(row) for row in c.fetchall() ]
        return pd.DataFrame(frame_buffer)

    # endregion

    # region Record retrieval - bins and contig associations

    def _bin_row_to_dict(self, _row):
        return { n: r for n, r in zip( ['bin_name', 'taxonomy', 'completeness', 'contamination', 'strainheterogeneity', 'version_name', 'date_created'], _row ) }

    def get_all_bins(self):

        ''' Returns a DataFrame of bin information - returns ALL bin versions '''

        ''' Create cursor, execute command and retrieve rows '''
        c = self.conn.cursor()
        c.execute('''SELECT binname, taxonomy, completeness, contamination, strainheterogeneity, version, datecreated FROM bin''')

        frame_buffer = [ self._bin_row_to_dict(row) for row in c.fetchall() ]
        return pd.DataFrame(frame_buffer)

    def get_bins_by_version(self, version_name):

        all_bins_df = self.get_all_bins()
        return all_bins_df[ all_bins_df.version_name == version_name ]

    def get_current_bins(self):

        all_bins_df = self.get_all_bins()

        observed_dates = [ date_parser.parse(d) for d in all_bins_df.date_created.unique() ]
        observed_dates.sort()

        most_recent_date = observed_dates[-1].isoformat()

        return all_bins_df[ all_bins_df.date_created == most_recent_date ]

    def _get_contig_bin_df(self):

        ''' Only used for internal functionality, not exposed to user '''

        ''' Create cursor, execute command and retrieve rows '''
        c = self.conn.cursor()
        c.execute('''SELECT contigname, binname FROM bin_contig''' )

        rows = [ row for row in c.fetchall() ]
        return pd.DataFrame.from_records(rows, columns=['contig_name', 'bin_name'])

    def get_contigs_in_bin(self, bin_name):

        ''' Returns a list of contigs associated with a particular bin. Returns None if no contigs are found. '''

        ''' Create cursor, execute command and retrieve rows '''
        c = self.conn.cursor()
        c.execute('''SELECT contigname FROM bin_contig WHERE binname = ?''', (bin_name,) )

        associated_contigs = [ row[0] for row in c.fetchall() ]

        if len(associated_contigs) > 0:
            return associated_contigs

        else:
            return None

    # endregion

    # region Input validation

    def _cast_as_ints(self, exp_ints):

        exp_ints = list(exp_ints)

        for i, exp_int in enumerate(exp_ints):

            try:

                exp_ints[i] = int(exp_int)

            except ValueError as ve:

                return False, exp_ints

        return True, tuple(exp_ints)

    def _cast_as_floats(self, exp_floats):

        exp_floats = list(exp_floats)

        for i, exp_float in enumerate(exp_floats):

            try:

                exp_floats[i] = float(exp_float)

            except ValueError as ve:

                return False, exp_floats

        return True, tuple(exp_floats)

    def _validate_columns(self, anntotation_df, exp_columns):

        try:

            for exp_column in exp_columns:

                assert( exp_column in anntotation_df.columns ), exp_column
            
            return True, ''

        except AssertionError as ae:

            return False, str(ae)
    # endregion

    # region File handling and sequence manupulation

    def _pull_fasta_metadata(self, read_header):

        if ' ' in read_header:
            seq_name, *metadata = read_header.split(' ')
            metadata = ' '.join(metadata)

        else:
            seq_name = read_header
            metadata = ''
        
        return seq_name, metadata

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
                seq_name, metadata = self._pull_fasta_metadata(header)

                fna_dict[seq_name] = (seq, metadata)

            else:
                fna_dict[header] = seq

        return fna_dict

    def _reverse_complement(self, sequence):

        backmap = { 'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N' }

        rev_comp = [''] * len(sequence)

        for i, nt in enumerate( list(sequence) ):
            rev_comp[i] = backmap[nt]

        return ''.join( rev_comp[::-1] )

    def _read_contig_sequence(self, contig_file, keep_metadata):

        ''' Not really necessary to handle this way, but unittests flag open file reader warnings if I don't explicitly close '''
        fna_reader = open(contig_file, 'r')
        fna_content = fna_reader.read()
        fna_reader.close()

        contigs_pass = [ x.strip().replace('>', '') for x in fna_content.split('\n') if '>' in x ]

        if keep_metadata:
            return [ self._pull_fasta_metadata(contig_header) for contig_header in contigs_pass ]

        else:
            return [ self._pull_fasta_metadata(contig_header)[0] for contig_header in contigs_pass ]

    # endregion

    # region Generation of user-specified data contexts

    '''
        This area is a bit inconsistent in implementation. In some cases, I use inner joins to perform cross-table queries, 
            and in others individual tables are pulled using existing functions then joined using pandas.

        This is a trade off between the complexity of writing the statements and how easily I can rely on the simple get() functions to acheive the desired effect.
    '''
    def present_coverage_by_contigs(self):

        ''' Return a DataFrame of reads mapped per bin, spread over sample sites '''

        cov_df = self.get_coverage()

        spread_df = cov_df.pivot(index='contig_name',columns='sample',values='reads_mapped').fillna(0)
        spread_df.reset_index(inplace=True)

        contig_df = self.get_contigs()
        contig_df.set_index('contig_name')

        #spread_df[ 'Contig_length' ] = 


        return spread_df

    def present_coverage_by_bins(self, current=False):

        ''' Return a DataFrame of reads mapped per bin, spread over sample sites. Optionally filter so that only the most current
            binning iteration is reported. '''
        return 0

    # endregion