#!/usr/bin/env python
r"""
Flow for ZSISA calculations
===========================

This examples shows how to compute temperature-dependent lattice parameters,
thermal expansion coefficients and relaxed-ions temperature-dependent elastic constants
using the ZSISA-QHA method. External pressure can be included as well.

By incorporating second-order derivatives of the vibrational free energy with respect
to lattice degrees of freedom, one can significantly reduce the number of
required phonon band structure calculations. See <https://arxiv.org/pdf/2503.09531>.

The AbiPy flow performs the following operations:

1) The initial structure is relaxed and a minimum set of deformed configurations is generated.
2) The atomic positions in the deformed structures are relaxed at fixed cell,
   and the relaxed configuration is then used to compute phonons, Born-effective charges,
   electronic dielectric tensort (eps_inf) and dynamical, quadrupoles with DFPT.
3) Phonon DOSes are computed for all the deformed structures on a much denser q-mesh
   and the results are used to compute thermal stresses.
4) An iterative process for determining lattice parameters at temperature T and external pressure P.
   The process begins with an initial guess for the lattice configuration [R],
   followed by the computation of thermal and BO stresses.
   A target stress is defined based on thermal stress and P, and the lattice and atomic positions
   are relaxed iteratively until the target stress matches the BO stress and the BO forces are zeroed,
   ensuring convergence to the optimal lattice configuration.
5) Finally, relaxed-ions elastic constants for the thermal-relaxed configurations are computed with DFPT.

The results of the calculations are stored in the outdata directory of the flow.
Two files are produced:

zsisa_data.csv:

    lattice parameters, thermal expansion, and elastic tensor components for each T and P in CSV format.

ZsisaResults.json:

    JSON file with the location of the different files (GSR.nc, DDB) on the filesystem.
    To recostruct a python object from file, use:

        from abipy.zsisa import ZsisaResults
        data = ZsisaResults.json_load("ABSPATH_TO_FILE")

        # Post-processing methods.
        data.get_dataframe()
        data.plot_lattice_vs_temp()
        data.plot_thermal_expansion()
"""
import sys
import os
import abipy.abilab as abilab
import abipy.data as abidata

from abipy import flowtk
from abipy.flowtk.zsisa import ZsisaFlow


def build_flow(options):
    """
    Create a `ZsisaFlow` for QHA calculations within the ZSISA method
    """
    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    if not options.workdir:
        __file__ = os.path.join(os.getcwd(), "run_qha_zsisa.py")
        options.workdir = os.path.basename(__file__).replace(".py", "").replace("run_", "flow_")

    # Initialize structure and pseudos for ZnO.
    structure = abilab.Structure.from_abistring("""
natom 4
ntypat 2
typat
1 1 2
2
znucl 30 8
xred
   0.0000000000    0.0000000000   -0.0025137620
   0.3333333333    0.6666666667    0.4974862380
   0.0000000000    0.0000000000    0.3835201241
   0.3333333333    0.6666666667    0.8835201241
acell    1.0    1.0    1.0
rprim
   6.3016720624    0.0000000000    0.0000000000
  -3.1508360312    5.4574080923    0.0000000000
   0.0000000000    0.0000000000    9.7234377918
""")

    # Initialize structure and pseudos
    # FIXME: This is just to make the computation faster.
    structure = abilab.Structure.from_file(abidata.cif_file("si.cif"))

    # IMPORTANT:
    # The Zsisa code assumes structure in primitive standard settings.
    # This call to abi_sanitize enforces the correct settings.
    structure = structure.abi_sanitize(primitive_standard=True)

    # Use NC PBEsol pseudos from pseudodojo v0.4.
    from abipy.flowtk.psrepos import get_oncvpsp_pseudos
    pseudos = get_oncvpsp_pseudos(xc_name="PBEsol", version="0.4")

    scf_input = abilab.AbinitInput(structure, pseudos)

    # Set other important variables for the SCF run.
    # Assuming non-magnetic semiconductor with scalar wavefunctions.
    scf_input.set_vars(
        nband=scf_input.num_valence_electrons // 2,
        nstep=100,
        ecutsm=1.0,
        tolvrs=1.0e-8,      # SCF stopping criterion.
        paral_kgb=0,
    )

    # Select k-mesh for electrons and q-mesh for phonons.
    # NB: The q-mesh should divide ngkpt.
    ngkpt = [2, 2, 2]; ngqpt = [1, 1, 1]
    #ngkpt = [6, 6, 4]; ngqpt = [1, 1, 1]
    #ngkpt = [4, 4, 4]; ngqpt = [2, 2, 2]

    scf_input.set_kmesh(ngkpt=ngkpt, shiftk=[0, 0, 0])

    # Flags to activate the computation of Born effective charges,
    # and dynamical quadrupoles with DFPT.
    with_becs = False
    #with_becs = True
    with_quad = False
    #with_quad = not structure.has_zero_dynamical_quadrupoles

    # List of temperatures in Kelvin and pressures in Gpa.
    #temperatures = [10, 50, 100, 200, 300]
    temperatures = [10, 100, 200]
    pressures_gpa = [0]
    #pressures_gpa = [0, 5]

    eps = 0.005  # Strain magnitude to be applied to the reference lattice.
    mode = "TEC" # "TEC" for thermal expansion only,
    mode = "ECs" # "ECs" to include T-dependent elastic constants
                 #  (more calculations are usually needed)

    # Q-mesh for the computation of the phonon DOS. See docstring.
    nqsmall_or_qppa = np.max(ngqpt) * 10

    # Relaxation with thermal stress may require several iterations.
    # Here we set programmatically the maximum number of restarts to 20.
    from abipy.flowtk.tasks import set_user_config_taskmanager_attrs
    set_user_config_taskmanager_attrs(max_num_launches=20)

    flow = ZsisaFlow.from_scf_input(options.workdir, scf_input, eps, mode, ngqpt,
                                    with_becs, with_quad, temperatures, pressures_gpa,
                                    nqsmall_or_qppa=nqsmall_or_qppa,
                                    )
    return flow


# This block generates the thumbnails in the Abipy gallery.
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
