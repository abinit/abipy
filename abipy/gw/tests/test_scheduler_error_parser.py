from __future__ import unicode_literals, division, print_function
import unittest
import os
import tempfile
from pymatgen.io.abinit.scheduler_error_parsers import get_parser, ALL_PARSERS, AbstractErrorParser
from pymatgen.io.abinit.scheduler_error_parsers import SlurmErrorParser, PBSErrorParser
from pymatgen.io.abinit.scheduler_error_parsers import SubmitError, FullQueueError, MemoryCancelError, TimeCancelError
from pymatgen.io.abinit.scheduler_error_parsers import NodeFailureError

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

        #with self.assertRaises(TypeError):
        #    AbstractErrorParser(err_file='dummy')

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
        string = "before\nBatch job submission failed\nafter"
        test_file = tempfile.mktemp()
        err_file = tempfile.mktemp()
        with open(err_file, 'w') as f:
            f.write("")
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
        string = "before\nsbatch: error: Batch job submission failed: Job violates accounting/QOS policy\nafter"
        test_file = tempfile.mktemp()
        err_file = tempfile.mktemp()
        with open(err_file, 'w') as f:
            f.write("")
        with open(test_file, 'w') as f:
            f.write(string)
        parser = SlurmErrorParser(err_file=err_file, batch_err_file=test_file)
        parser.parse()
        self.assertEqual(len(parser.errors), 2)
        print(parser.errors)
        self.assertTrue(isinstance(parser.errors[0], FullQueueError) or isinstance(parser.errors[1], FullQueueError))
        os.remove(test_file)
        os.remove(err_file)

    def test_MemoryCancelError(self):
        """
        Testing the SLURM error parsing of MemoryCancelError errors
        """
        string = "Exceeded job memory limit"
        err_file = tempfile.mktemp()
        with open(err_file, 'w') as f:
            f.write(string)
        parser = SlurmErrorParser(err_file=err_file)
        parser.parse()
        self.assertEqual(len(parser.errors), 1)
        self.assertIsInstance(parser.errors[0], MemoryCancelError)
        self.assertEqual(parser.errors[0].meta_data, {})
        os.remove(err_file)

    def test_TimeCancelError(self):
        """
        Testing the SLURM error parsing of TimeCancelError que errors
        """
        # 'err'
        string = "JOB 999 CANCELLED AT 888 DUE TO TIME LIMIT"
        err_file = tempfile.mktemp()
        with open(err_file, 'w') as f:
            f.write(string)
        parser = SlurmErrorParser(err_file=err_file)
        parser.parse()
        self.assertEqual(len(parser.errors), 1)
        self.assertIsInstance(parser.errors[0], TimeCancelError)
        print(parser.errors[0])
        self.assertEqual(parser.errors[0].meta_data['time_of_cancel'][0], u'888')
        os.remove(err_file)

    def test_NodeFailureError(self):
        """
        Testing the SLURM error parsing of NodeFailureError que errors
        """
        string = "node123.456can't open /dev/ipath, network down (err=26)"
        err_file = tempfile.mktemp()
        test_file = tempfile.mktemp()
        with open(err_file, 'w') as f:
            f.write("")
        with open(test_file, 'w') as f:
            f.write(string)
        parser = SlurmErrorParser(err_file=err_file, run_err_file=test_file)
        parser.parse()
        self.assertEqual(len(parser.errors), 1)
        self.assertIsInstance(parser.errors[0], NodeFailureError)
        self.assertEqual(parser.errors[0].meta_data['nodes'], [u'123'])
        os.remove(err_file)
        os.remove(test_file)


class QueueErrorParseTestPBS(unittest.TestCase):

    def test_TimeCancelError(self):
        """
        Testing the PBS error parsing of TimeCancelError errors
        """
        string = "before\n=>> PBS: job killed: walltime 1001 exceeded limit 1000\nafter"
        test_file = tempfile.mktemp()
        err_file = tempfile.mktemp()
        with open(err_file, 'w') as f:
            f.write("")
        with open(test_file, 'w') as f:
            f.write(string)
        parser = PBSErrorParser(err_file=err_file, out_file=test_file)
        parser.parse()
        self.assertEqual(len(parser.errors), 1)
        print(parser.errors)
        self.assertIsInstance(parser.errors[0], TimeCancelError)
        os.remove(test_file)
        os.remove(err_file)

    def test_MemoryCancelError(self):
        """
        Testing the PBS error parsing of MemoryCancelError errors
        """
        string = "before\n(.*)job killed: vmem 4001kb exceeded limit 4000kb\nafter"
        test_file = tempfile.mktemp()
        err_file = tempfile.mktemp()
        with open(err_file, 'w') as f:
            f.write("")
        with open(test_file, 'w') as f:
            f.write(string)
        parser = PBSErrorParser(err_file=err_file, out_file=test_file)
        parser.parse()
        self.assertEqual(len(parser.errors), 1)
        print(parser.errors)
        self.assertIsInstance(parser.errors[0], MemoryCancelError)
        os.remove(test_file)
        os.remove(err_file)


if __name__ == '__main__':
    unittest.main()
