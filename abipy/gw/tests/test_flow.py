from __future__ import division, print_function, unicode_literals
import os
from pymatgen.util.testing import PymatgenTest
from pymatgen.io.abinit.flows import Flow
from pymatgen.io.abinit.tasks import TaskManager, Task
from pymatgen.io.abinit.works import Work

__author__ = 'setten'

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", 'test_files')


class FlowTest(PymatgenTest):
    def test_flow(self):
        """
        Testing flow creation and task registering
        """
        flow = Flow(workdir=test_dir, manager=TaskManager.from_file(os.path.join(test_dir, "taskmanager.yml")))
        inp = {}
        flow.register_task(input=inp)
        flow.allocate()
        self.assertTrue(flow.allocated)
        self.assertIsInstance(flow[0], Work)
        self.assertIsInstance(flow[0][0], Task)
        self.assertEqual(flow.check_status(), None)
