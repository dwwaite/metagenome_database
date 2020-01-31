'''
    Automated testing of the interation between the front- and back-end of the database tool.
'''
import sys
import argparse
import os
import io
import unittest
import sqlite3
import time
from functools import reduce
import pandas as pd

class TestFrontend(unittest.TestCase):

    # region Overhead functions
 
    def setUp(self):
        self.temp_file_buffer = []

    def tearDown(self):
        ''' Add a small delay, to make sure the database connection has closed before removing the file '''
        time.sleep(0.5)
        _ = [ os.remove(temp_file) for temp_file in self.temp_file_buffer if os.path.exists(temp_file) ]               

    def start_logging_stdout(self):
        self.print_capture = io.StringIO()
        sys.stdout = self.print_capture

    def stop_logging_stdout(self):
        sys.stdout = sys.__stdout__
        return self.print_capture.getvalue()

    # endregion

    #region Creation

    def create_database(self):

        database_name = 'mock.hdb'
        fe_create_database(database_name)
        self.temp_file_buffer.append(database_name)
        return database_name

    def test_create_database(self):

        ''' Attempt to create a fresh database, populating with the contents of tests/contigs.fna '''

        ''' Test '''
        database_name = self.create_database()

        ''' Evaluate '''
        self.assertTrue( os.path.exists(database_name) )

    def test_create_database_already_exists(self):

        ''' Make sure that if a database already exists, the correct error is reported.
            The actual error code is masked away by the fe_create_database() function, but it passses up the error
                message spawned by the underlying DatabaseManipulator, so just catching and reading this message is sufficient. '''

        ''' Test '''
        database_name = self.create_database()

        self.start_logging_stdout()
        _ = self.create_database()

        ''' Evaluate '''
        msg = self.stop_logging_stdout()
        self.assertIn( 'Database file {} is already populated. Did you mean to invoke the update command?'.format(database_name), msg )

    #endregion

    #region Add
    
    def create_blank_params(self):
        return { 'contig': None, 'prodigal_aa': None, 'prodigal_nt': None,
                 'rrna_ssu': None, 'rrna_lsu': None, 'trna': None, 'trna_meta': None,
                 'annotations': None, 'coverage': None, 'transcript': None, 'bin_file': None, 'bin_dir': None }

    def test_add_features(self):

        ''' Won't test this exhaustively, since individual unit tests already ensure these functions work as expected.
            Just to a simple pair of contig/gene additions to ensure the dependency is satisfied. '''

        ''' Set up the database and params dict '''
        database_name = self.create_database()
        params = self.create_blank_params()
        params['contig'] = 'contigs.fna'
        params['prodigal_aa'] = 'genes_prod_aa.faa'

        fe_add_features(database_name, params )

        ''' Evaluate '''
        db_m = DatabaseManipulator(database_name)

        contig_df = db_m.get_contigs()
        self.assertEqual( contig_df.shape[0], 3 )

        gene_df = db_m.get_genes()
        self.assertEqual( gene_df.shape[0], 4 )

    def test_add_features_failed_dependencies(self):

        ''' Set up the database and params dict '''
        database_name = self.create_database()
        params = self.create_blank_params()
        params['prodigal_aa'] = 'genes_prod_aa.faa'

        ''' Evaluate '''
        self.start_logging_stdout()
        fe_add_features(database_name, params)

        msg = self.stop_logging_stdout()
        gene_names = ['contig1_1', 'contig1_2', 'contig2_1', 'contig3_1']

        err_mask = [ 'Gene {} is not linked to an existing contig. Aborting...'.format(g) in msg for g in gene_names ]
        self.assertTrue( reduce(lambda a, b : a or b, err_mask) )

        db_m = DatabaseManipulator(database_name)
        gene_df = db_m.get_genes()
        self.assertEqual( gene_df.shape[0], 0 )

    def test_add_trna_missing_param(self):

        ''' Set up the database and params dict '''
        database_name = self.create_database()
        params = self.create_blank_params()
        params['contig'] = 'contigs.fna'
        params['trna'] = 'genes_trna.fna'

        ''' Evaluate '''
        self.start_logging_stdout()
        fe_add_features(database_name, params )

        msg = self.stop_logging_stdout()
        self.assertIn( 'Missing required input parameters for tRNA annotation, please check input options.', msg )

    #endregion

if __name__ == '__main__':

    ''' Import the parent path, so that we can import the metagenome_databaser functions.
        When importing front-end functions, alias them with the fe_ prefix, so that it's easy to spot them in code '''
    sys.path.insert(0, '..')
    from metagenome_databaser import create_database as fe_create_database
    from metagenome_databaser import add_features as fe_add_features
    from scripts.DatabaseManipulator import DatabaseManipulator

    unittest.main()