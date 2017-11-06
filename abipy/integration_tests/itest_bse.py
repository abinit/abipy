"""Integration tests for BSE flows."""
from __future__ import print_function, division, unicode_literals, absolute_import

import pytest
import abipy.data as abidata
import abipy.abilab as abilab
import abipy.flowtk as flowtk

from abipy.core.testing import has_abinit, has_matplotlib

# TODO
#def itest_bse_with_mdf(fwp, tvars):
#    pseudos = abidata.pseudos("14si.pspnc")
#    structure = abilab.Structure.from_file(abidata.cif_file("si.cif"))
#
#    kppa = scf_kppa = 1
#    nscf_nband = 6
#    nscf_ngkpt = [4,4,4]
#    nscf_shiftk = [0.1, 0.2, 0.3]
#    bs_loband = 2
#    bs_nband = nscf_nband
#    soenergy = 0.7
#    mdf_epsinf = 12
#    max_ncpus = 1
#    ecuteps = 2
#
#    extra_abivars = dict(
#        ecut=12,
#        istwfk="*1",
#    )
#
#    flow = flowtk.Flow(workdir=fwp.workdir, manager=fwp.manager)
#
#    # BSE calculation with model dielectric function.
#    from pymatgen.io.abinit.calculations import bse_with_mdf
#    work = bse_with_mdf(structure, pseudos, scf_kppa, nscf_nband, nscf_ngkpt, nscf_shiftk,
#                       ecuteps, bs_loband, bs_nband, soenergy, mdf_epsinf,
#                       accuracy="normal", spin_mode="unpolarized", smearing=None,
#                       charge=0.0, scf_solver=None, **extra_abivars)
#
#    flow.register_work(work)
#    flow.allocate()
#    flow.build_and_pickle_dump(abivalidate=True)
#
#    for task in flow:
#        task.start_and_wait()
#        assert task.status == task.S_DONE
#
#    flow.check_status()
#    flow.show_status()
#    assert all(work.finalized for work in flow)
#    assert flow.all_ok
#    assert flow.validate_json_schema()
