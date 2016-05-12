from __future__ import unicode_literals, division, print_function
import unittest
import shutil
import os
import tempfile
from pymatgen.io.abinit.scheduler_error_parsers import get_parser, ALL_PARSERS, AbstractErrorParser
from pymatgen.io.abinit.scheduler_error_parsers import SlurmErrorParser
from pymatgen.io.abinit.scheduler_error_parsers import SubmitError, FullQueueError


__author__ = "Michiel van Setten"
__copyright__ = " "
__version__ = "0.9"
__maintainer__ = "Michiel van Setten"
__email__ = "mjvansetten@gmail.com"
__date__ = "Mar 24, 2014"


class QueueErrorParseAPITest(unittest.TestCase):

    def test_get_parser(self):
        """
        Testing the get parser API
        """
        for (scheduler) in ALL_PARSERS:
            err_file = 'dummy'
            parser = get_parser(scheduler, err_file=err_file)
            self.assertIsInstance(parser.files, dict)
            self.assertIsInstance(parser.errors, list)
            self.assertIsInstance(parser.error_definitions, dict)
            self.assertEqual(parser.parse(), len(parser.error_definitions))

    def test_parser_abstract_methods(self):
        """
        Testing the methods of the abstract error parser class
        """

        with self.assertRaises(TypeError):
            AbstractErrorParser(err_file='dummy')

        class DummyErrorParser(AbstractErrorParser):

            @property
            def error_definitions(self):
                return dict()

        parser = DummyErrorParser(err_file='dummy')
        lines = ['A 12 B', 'A message B']
        meta_filer = {'broken_limit': [r"A (\d+) B", 1]}

        self.assertEqual(parser.parse(), 0)
        self.assertIsInstance(parser.error_definitions, dict)
        self.assertEqual(parser.extract_metadata(lines=lines, meta_filter=meta_filer), {u'broken_limit': [u'12']})
        self.assertEqual(parser.parse_single({}), (False, None, None))


class QueueErrorParseTestSLURM(unittest.TestCase):

    def setUp(self):
        pass

    def test_SubmitError(self):
        """
        Testing the SLURM error parsing of submit que errors
        """
        test_file = tempfile.mktemp()
        err_file = tempfile.mktemp()
        with open(err_file, 'w') as f:
            f.write("")
        string = "before\nBatch job submission failed\nafter"
        with open(test_file, 'w') as f:
            f.write(string)
        parser = SlurmErrorParser(err_file=err_file, batch_err_file=test_file)
        parser.parse()
        self.assertEqual(len(parser.errors), 1)
        self.assertIsInstance(parser.errors[0], SubmitError)
        os.remove(test_file)
        os.remove(err_file)

    def test_FullQueueError(self):
        """
        Testing the SLURM error parsing of full que errors
        """
        test_file = tempfile.mktemp()
        err_file = tempfile.mktemp()
        with open(err_file, 'w') as f:
            f.write("")
        string = "before\nsbatch: error: Batch job submission failed: Job violates accounting/QOS policy\nafter"
        with open(test_file, 'w') as f:
            f.write(string)
        parser = SlurmErrorParser(err_file=err_file, batch_err_file=test_file)
        parser.parse()
        self.assertEqual(len(parser.errors), 2)
        self.assertIsInstance(parser.errors[1], FullQueueError)
        os.remove(test_file)
        os.remove(err_file)

    def test_MemoryCancelError(self):
        """
        Testing the SLURM error parsing of submit que errors
        """

        definition = {
                           'err': {
                               'string': "Exceeded job memory limit",
                               'meta_filter': {}
                           }
                       },


    def test_TimeCancelError(self):
        """
        Testing the SLURM error parsing of submit que errors
        """

        definition = {
                         'err': {
                             'string': "DUE TO TIME LIMIT",
                             'meta_filter': {
                                 'time_of_cancel': [r"JOB (\d+) CANCELLED AT (\S*) DUE TO TIME LIMIT", 1]
                             }
                         }
                     },

    def NodeFailureError(self):
        """
        Testing the SLURM error parsing of submit que errors
        """

        definition = {
                          'run_err': {
                              'string': "can't open /dev/ipath, network down",
                              'meta_filter': {
                                  'nodes': [r"node(\d+)\.(\d+)can't open (\S*), network down \(err=26\)", 1]
                              }
                          }
                      },


if __name__ == '__main__':
    unittest.main()
