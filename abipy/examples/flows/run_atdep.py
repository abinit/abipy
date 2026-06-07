from pathlib import Path

import abipy.data as abidata
from abipy.abio.inputs import AtdepInput, AnaddbInput
from abipy import flowtk

def main():

    workdir = Path(__file__).name.replace(".py", "").replace("run_", "flow_")

    structure = abidata.structure_from_mpid('mp-1265')  # MgO
    ddb_file = abidata.ref_file('MgO_dfpt_aimd/MgO_DDB.nc')
    hist_file = abidata.ref_file('MgO_dfpt_aimd/MgO_HIST.nc')

    flow = get_atdep_flow(workdir, structure, hist_file, ddb_file)
    flow.build_and_pickle_dump()


def get_atdep_flow(workdir, structure, hist_file, ddb_file=None):

    atdep_inp = get_atdep_input(structure)
    anaddb_inp = get_anaddb_input(structure)

    anaddb_inp.set_vars(
        ngqpt = atdep_inp['ngqpt1'],
        rifcsph = atdep_inp['rcut'],
        )

    atdep_task = flowtk.AtdepTask(atdep_inp, hist_node=hist_file, ddb_node=ddb_file)
    anaddb_task = flowtk.AnaddbTask(anaddb_inp, ddb_node=atdep_task)

    flow = flowtk.Flow(workdir=workdir)

    flow.register_task(atdep_task)
    flow.register_task(anaddb_task, deps={atdep_task: "DDB"})

    return flow


def get_atdep_input(structure):

    inp = AtdepInput(structure)
    inp.set_vars(
        multiplicity=[[-2,2,2],[2,-2,2],[2,2,-2]],
        nstep_max=20,
        nstep_min=1,
        temperature=300,
        rcut=7.0,
        ngqpt1=[2,2,2],
        ngqpt2=[1,1,1],
        )

    return inp

def get_anaddb_input(structure):

    inp = AnaddbInput.phbands_and_dos(
        structure, ngqpt=[2,2,2], ndivsm=40, line_density=None,
        nqsmall=10, qppa=None, q1shft=(0, 0, 0), qptbounds=None,
        asr=2, chneut=1, dipdip=1, dipquad=0, quadquad=0,
        dos_method='tetra', lo_to_splitting='automatic',
        with_ifc=True, anaddb_kwargs={}, spell_check=False)

    return inp


if __name__ == "__main__":
    main()
