"""Configuration file for pytest."""
from __future__ import print_function, division, unicode_literals

import os
import pytest
import abipy.abilab as abilab
import abipy.data as abidata
from abipy.htc.factories import ebands_input
from abipy.data.benchmark_structures import simple_semiconductors, simple_metals

from fireworks import LaunchPad, FWorker

from pymatgen.matproj.rest import MPRester


TESTDB_NAME = 'fireworks_unittest'
MODULE_DIR = os.path.dirname(os.path.abspath(__file__))


@pytest.fixture(scope="module")
def lp(request):
    lp = LaunchPad(name=TESTDB_NAME, strm_lvl='ERROR')
    lp.reset(password=None, require_password=False)

    def fin():
        lp.connection.drop_database(TESTDB_NAME)

    request.addfinalizer(fin)
    return lp


@pytest.fixture(scope="module")
def fworker():
    return FWorker()


@pytest.fixture(scope="function")
def cleandb(request, lp):
    def fin():
        lp.reset(password=None, require_password=False)
    request.addfinalizer(fin)


@pytest.fixture(scope="function")
def input_scf_si_low():
    pseudos = abidata.pseudos("14si.pspnc")
    cif_file = abidata.cif_file("si.cif")
    structure = abilab.Structure.from_file(cif_file)

    return ebands_input(structure, pseudos, kppa=100, ecut=6).split_datasets()[0]


@pytest.fixture(scope="function", params=simple_semiconductors)
def benchmark_input_scf(request):
    pseudos = abidata.pseudos("14si.pspnc", "6c.pspnc", "3li.pspnc", "9f.pspnc",
                              "12mg.pspnc", "8o.pspnc", "31ga.pspnc", "7n.pspnc")
    rest = MPRester()
    structure = rest.get_structure_by_material_id(request.param)
    try:
        return ebands_input(structure, pseudos, kppa=100, ecut=6).split_datasets()[0]
    except:
        #to deal with missing pseudos
        pytest.skip('Cannot create input for material {}.'.format(request.param))


@pytest.fixture()
def fwp(tmpdir):
    """
    Parameters used to initialize Flows.
    """
    # Temporary working directory
    fwp.workdir = str(tmpdir)

    # Create the TaskManager.
    fwp.manager = abilab.TaskManager.from_file(os.path.join(os.path.dirname(__file__), "manager.yml"))

    fwp.scheduler = abilab.PyFlowScheduler.from_file(os.path.join(os.path.dirname(__file__), "scheduler.yml"))

    return fwp



