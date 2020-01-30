'''
    Automated testing of the interation between the front- and back-end of the database tool.

    Useful snippet for checking database programmatically
    
import sqlite3
db = sqlite3.connect('mock.hdb')
cursor = db.cursor()
cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
print(cursor.fetchall())

cursor = db.cursor()
cursor.execute("SELECT * FROM contig;")
print(cursor.fetchall())

'''
import sys, argparse, os, io
import unittest
import sqlite3
import pandas as pd

class TestFrontend(unittest.TestCase):

    # region Overhead functions
 
    def setUp(self):
        self.temp_file_buffer = []

    def tearDown(self):
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

        ''' Tidy up '''
        os.remove(database_name)

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

    """
        db_m = DatabaseManipulator(database_name)
        contig_df = db_m.get_contigs()
    """

    #endregion

if __name__ == '__main__':

    ''' Import the parent path, so that we can import the metagenome_databaser functions.
        When importing front-end functions, alias them with the fe_ prefix, so that it's easy to spot them in code '''
    sys.path.insert(0, '..')
    from metagenome_databaser import create_database as fe_create_database
    from scripts.DatabaseManipulator import DatabaseManipulator

    unittest.main()