"""Tests for dfpt_works module"""
from __future__ import print_function, division, unicode_literals, absolute_import

import abipy.data as abidata
import abipy.flowtk as flowtk

from abipy.core.testing import AbipyTest
#from abipy.flowtk import dfpt_works


class TestDfptWorks(AbipyTest):

    def test_nscfddkswork(self):
        """Testing NscfDdksWork."""
        scf_task = self.get_gsinput_si(as_task=True)
        work = flowtk.NscfDdksWork.from_scf_task(scf_task, ddk_ngkpt=[8, 8, 8],
            ddk_shiftk=[0, 0, 0], ddk_nband=10)
        assert len(work) == 4
        self.abivalidate_work(work)

    def test_elastic_work(self):
        """Testing ElasticWork."""
        scf_task = self.get_gsinput_si(as_task=True)
        scf_input = scf_task.input
        den_deps = {scf_task: "DEN"}
        tolerances = dict(nscf={"tolwfr": 1.0e-10}, ddk={"tolwfr": 1.0e-12}, strain={"tolvrs": 1.0e-10})
        work = flowtk.ElasticWork.from_scf_input(scf_input,
	    with_relaxed_ion=True, with_piezo=True, with_dde=True, tolerances=tolerances,
            den_deps=den_deps, manager=None)
        self.abivalidate_work(work)

        #assert len(work) == 4
        assert work[0].input["iscf"] == -2
        assert work[0].input["tolwfr"] == tolerances["nscf"]["tolwfr"]
        assert isinstance(work[0], flowtk.NscfTask)

        for task in work[1:]:
            assert task.input["kptopt"] == 2
        assert all(isinstance(task, flowtk.ElasticTask) for task in work[-6:])
        for task in work[-6:]:
            assert task.input["tolvrs"] == tolerances["strain"]["tolvrs"]
