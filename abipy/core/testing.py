"""
Common test support for all abipy test scripts.

This single module should provide all the common functionality for abipy tests
in a single location, so that test scripts can just import it and work right away.
"""
from __future__ import print_function, division

from pymatgen.util.io_utils import which
from pymatgen.util.testing import PymatgenTest

__all__ = [
    "AbipyTest",
    "AbipyFileTest",
]

class AbipyTest(PymatgenTest):
    """Extend TestCase with functions from numpy.testing.utils that support ndarrays."""

    @staticmethod
    def which(program):
        """Returns full path to a executable. None if not found or not executable."""
        return which(program)

    def serialize_with_pickle(self, objects, protocols=None):
        """
        Test whether the object can be serialized and deserialized with pickle.

        Args:
            objects:
                Object or list of objects. The object must define the __eq__ operator.
            protocols:
                List of pickle protocols to test.

        Returns:
            List of objects deserialized with the specified protocols.
        """
        import tempfile
        import cPickle as pickle

        # Build a list even when we receive a single object.
        if not isinstance(objects, (list, tuple)):
            objects = [objects]

        # By default, all pickle protocols are tested.
        if protocols is None:
            protocols = set([0, 1, 2] + [pickle.HIGHEST_PROTOCOL])

        # This list will contains the deserialized object for each protocol.
        objects_by_protocol = []

        for protocol in protocols:
            # Serialize and deserialize the object.
            mode = "w" if protocol == 0 else "wb"
            fd, tmpfile = tempfile.mkstemp(text="b" not in mode)

            with open(tmpfile, mode) as fh:
                pickle.dump(objects, fh, protocol=protocol)

            with open(tmpfile, "r") as fh:
                new_objects = pickle.load(fh)

            # Save the deserialized objects and test for equality.
            objects_by_protocol.append(new_objects)

            for old_obj, new_obj in zip(objects, new_objects):
                self.assert_equal(old_obj, new_obj)
                #self.assertTrue(str(old_obj) == str(new_obj))

        return objects_by_protocol


class AbipyFileTest(AbipyTest):
    """
    Test class for files with a __str__ attribute.
    At setup, must set the 'file' attribute of the AbipyFileTest.
    """
    file = None

    @staticmethod
    def normalize(string):
        string = string.replace('\t', '  ')
        string = string.replace("$", " CASH ")
        string = string.replace("(", " LP ")
        string = string.replace(")", " RP ")
        string = string.replace("*", " STAR ")
        string = string.strip()
        string = '\n'.join([line.strip() for line in string.splitlines()])
        return string

    def assertContains(self, expression):
        """
        Assert that the string representation of the file contains 'expression'
        'expression' is trimmed of leading new line.
        Each line of 'expression' is trimmed of blank spaces.
        Empty lines are ignored.
        """
        expression = self.normalize(expression)
        ref = self.normalize(str(self.file))
        return  self.assertRegexpMatches(ref, expression)

    def assertContainsNot(self, expression):
        """
        Assert that the string representation of the file does not contain
        'expression'.
        """
        expression = self.normalize(expression)
        ref = self.normalize(str(self.file))
        return  self.assertNotRegexpMatches(ref, expression)

    def assertEmpty(self):
        """Assert the string representation is empty."""
        s = str(self.file).strip()
        self.assertFalse(bool(s))

    def assertOrder(self, expression1, expression2):
        """
        Assert that the string representation of the file
        contains 'expression1' before 'expression2'.
        """
        expression1 = self.normalize(expression1)
        expression2 = self.normalize(expression2)
        ref = self.normalize(str(self.file))
        self.assertRegexpMatches(ref, expression1)
        self.assertRegexpMatches(ref, expression2)
        last = ref.split(expression1)[-1]
        return self.assertRegexpMatches(last, expression2)


