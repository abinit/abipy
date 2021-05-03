import sys
import os
import abipy.data as abidata
import abipy.abilab as abilab
import abipy.flowtk as flowtk
from abipy.abilab import abiopen

def build_flow(options):
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    if not options.workdir:
        options.workdir = os.path.basename(sys.argv[0]).replace(".py", "").replace("run_", "flow_")

    # Initialize the flow
    flow = flowtk.Flow(workdir=options.workdir, manager=options.manager)

    # Initialize the tmesh, sigma_kerange, sigma_erange and dense meshes for the eph integrations
    tmesh = [300,300,1]
    sigma_kerange = [0, "0.5 eV"]
    sigma_erange = [0, "0.25 eV"]
    dense_meshes = [[30,30,30],
                    [40,40,40]]

    # We first compute the DDB and DVDB files
    structure = abidata.structure_from_ucell("AlAs")
    pseudos = abidata.pseudos("13al.981214.fhi", "33as.pspnc")

    # Ground-state computation for 1) the phonons and 2) the WFK generation
    scf_input = abilab.AbinitInput(structure, pseudos=pseudos)

    scf_input.set_vars(
        nband=8,
        ecut=2.0,
        ngkpt=[4, 4, 4],
        nshiftk=1,
        shiftk=[0, 0, 0],
        tolvrs=1.0e-10,
        diemac=9.0,
        prtden=1,
        #iomode=3,
    )

    ddb_ngqpt = [4, 4, 4]

    # Add the ground-state work to the flow
    work_scf = flow.register_scf_task(scf_input)

    # Band structure calculation to make sure everything is ok
    bs_input = scf_input.make_ebands_input(tolwfr=1e-12, ndivsm=10, nb_extra=4)
    work_scf.register_nscf_task(bs_input, deps={work_scf[0]: "DEN"})

    # Add the phonon work to the flow
    ph_work = flowtk.PhononWork.from_scf_task(work_scf[0], qpoints=ddb_ngqpt, is_ngqpt=True, with_becs=True)
    flow.register_work(ph_work)

    # NSCF input for the WFK needed to interpolate with kerange
    nscf_input = abilab.AbinitInput(structure, pseudos)
    nscf_input.set_vars(
            ecut=2,
            nband=8,
            iscf=-2,
            tolwfr=1e-20,
            prtwf=1,
            ngkpt=[16, 16, 16], # Should be dense enough so that the kerange interpolation works
            shiftk=[0.0, 0.0, 0.0],
    )

    work_nscf = flowtk.Work()
    work_nscf.register_nscf_task(nscf_input, deps={work_scf[0]: "DEN"})
    flow.register_work(work_nscf)

    # We loop over the dense meshes
    for i, sigma_ngkpt in enumerate(dense_meshes):
        # Use the kerange trick to generate a WFK file
        kerange_input, wfk_input = nscf_input.make_wfk_kerange_input(sigma_erange=sigma_kerange, sigma_ngkpt=sigma_ngkpt).split_datasets()

        work_eph = flowtk.Work()
        work_eph.register_kerange_task(kerange_input, deps={work_nscf[0]: "WFK"})
        work_eph.register_nscf_task(wfk_input, deps={work_scf[0]: "DEN", work_eph[0]: "KERANGE.nc"})

        # Generate the input file for the transport calculation
        eph_input = wfk_input.make_eph_transport_input(ddb_ngqpt=ddb_ngqpt, sigma_erange=sigma_erange,
                                                       tmesh=tmesh, eph_ngqpt_fine=sigma_ngkpt)

        # We compute the phonon dispersion to be able to check they are ok
        if i==0:
            eph_input.set_qpath(20)

        work_eph.register_eph_task(eph_input, deps={work_eph[1]: "WFK", ph_work: ["DDB", "DVDB"]})

        flow.register_work(work_eph)
        

    flow.allocate(use_smartio=True)

    return flow


# This block generates the thumbnails in the AbiPy gallery.
# You can safely REMOVE this part if you are using this script for production runs.
if os.getenv("READTHEDOCS", False):
    __name__ = None
    import tempfile
    options = flowtk.build_flow_main_parser().parse_args(["-w", tempfile.mkdtemp()])
    build_flow(options).graphviz_imshow()


@flowtk.flow_main
def main(options):
    """
    This is our main function that will be invoked by the script.
    flow_main is a decorator implementing the command line interface.
    Command line args are stored in `options`.
    """
    
    return build_flow(options)


if __name__ == "__main__":
    sys.exit(main())
