.. _post-processing-howto:

***********************
Post-processing How-To
***********************

This page collects frequently asked questions about the AbiPy scripts.
Feel free to suggest new entries!

.. contents::
   :backlinks: top


Preliminary considerations
---------------------------

The AbiPy scripts detect the file type from the file extension,
so **do not change the file extension**.
Remember also that you can get the documentation of any script simply by adding
the ``--help`` option to the command line.
For example:

.. code-block:: shell

    abistruct.py --help

shows the documentation and usage examples for the abistruct.py_ script, while:

.. code-block:: shell

    abistruct.py COMMAND --help

prints the documentation and the options supported by ``COMMAND``.


Get information about a generic ``FILE``
----------------------------------------

Use::

    abiopen.py FILE --print

to print information about a file in the terminal, or::

    abiopen.py FILE --expose

to generate a set of matplotlib figures that depend on the type of FILE.

Use ``--verbose`` or ``-v`` to increase the verbosity level.
The option can be given multiple times, e.g. ``-vv``.

Get all file extensions supported by ``abiopen.py``
---------------------------------------------------

Use::

    abiopen.py --help

.. command-output:: abiopen.py --help


Convert the structure stored in ``FILE`` to a different format
--------------------------------------------------------------

Use::

    abistruct.py convert FILE

to read the structure from ``FILE`` and generate a CIF_ file (the default behaviour).

Most of the netcdf_ files produced by Abinit contain structural information,
so this command works with netcdf output files as well as with Abinit input/output
files and all the other formats supported by pymatgen, e.g. POSCAR files.
Other output formats can be selected with the ``-f`` option.
For example::

    abistruct.py convert FILE -f abivars

prints the Abinit input variables, while::

    abistruct.py convert FILE -f xsf > out.xsf

exports the structure in the ``xsf`` format (xcrysden_) and saves it to a file.

Use::

    abistruct.py convert --help

to list the supported formats.

Check if my Abinit input file is OK
-----------------------------------

First of all, use::

    abiopen.py ../abipy/data/refs/si_ebands/run.abi -p

to print the crystalline structure and determine the space group with the spglib_ library.

If the structure looks good, validate the input file with Abinit using the ``validate`` command of the abinp.py_ script::

    abinp.py validate run.abi

This requires a ``manager.yml`` file and, of course, Abinit.

The script provides other commands that invoke Abinit
to get space group information, the list of k-points in the IBZ,
the list of atomic perturbations for phonons or the list of autoparal configurations.
See ``abinp.py --help`` for more information.

Print the warnings in the log file
----------------------------------

Use::

    abiopen.py run.log -p

to get::

    Events found in /Users/gmatteo/git_repos/abipy/abipy/examples/flows/develop/flow_from_files/w0/t0/run.log

    [1] <AbinitWarning at m_nctk.F90:568>
        netcdf lib does not support MPI-IO and: NetCDF: Parallel operation on file opened for non-parallel access

    [2] <AbinitWarning at m_nctk.F90:588>
        The netcdf library does not support parallel IO, see message above
        Abinit won't be able to produce files in parallel e.g. when paral_kgb==1 is used.
        Action: install a netcdf4+HDF5 library with MPI-IO support.

    [3] <AbinitWarning at m_hdr.F90:4258>
        input kptrlatt= 0 0 0 0 0 0 0 0 0  /= disk file kptrlatt=8 0 0 0 8 0 0 0 8

    [4] <AbinitWarning at m_hdr.F90:4261>
        input kptopt= -2  /= disk file kptopt= 1

    num_errors: 0, num_warnings: 4, num_comments: 0, completed: True

A similar interface is also available with::

    abiview.py log run.log


Get a quick look at a file
--------------------------

The abiview.py_ script is designed specifically for this task.
The syntax is ``abiview.py COMMAND FILE``, where ``COMMAND`` is either
the Abinit file extension (without ``.nc``, if any) or the AbiPy object to visualize.

To get a quick look at the DDB file, use::

    abiview.py ddb out_DDB

This command invokes anaddb to compute the phonon bands and DOS from the DDB and produces matplotlib_ plots.

If ``FILE`` contains electronic band energies, use e.g.::

    abiview.py ebands out_GSR.nc

to plot the KS eigenvalues (the same command works for other files such as ``WFK.nc``, ``DEN.nc``, etc.).

Note that abiview.py_ visualizes the data according to a predefined logic.
Some options let you tune parameters or export the data in different formats,
but exposing the full AbiPy API from the command line is not practical.

For a more flexible interface, use::

    abiopen.py FILE

to start an ipython_ shell in which you can interact with the Python object directly.

If jupyter_ is installed on your machine/cluster and you have a web browser, use::

    abiopen.py FILE -nb

to automatically generate a predefined jupyter notebook for that file type.

Visualize a structure
---------------------

Structure visualization is delegated to external graphical applications
that must be installed on your machine.
AbiPy extracts the structure from ``FILE``, converts it to one of the formats
supported by the graphical application and then invokes the executable.
If vesta_ is installed in one of the standard
locations on your machine, simply run::

    abistruct.py visualize FILE

in the terminal.
Other applications can be selected with the ``--application`` option.
At present, AbiPy supports vesta_, ovito_, xcrysden_, avogadro_, and v_sim_.

To visualize the crystalline structure inside a jupyter_ notebook, you may want to
try the nbjsmol_ jupyter extension.

Get a high-symmetry kpath for a given structure
-----------------------------------------------

Use the `kpath` command with a FILE that provides structural information::

    abistruct.py kpath FILE

to generate a template with the input variables that define the k-path:

.. code-block:: shell

     # Abinit Structure
     natom 2
     ntypat 1
     typat 1 1
     znucl 14
     xred
        0.0000000000    0.0000000000    0.0000000000
        0.2500000000    0.2500000000    0.2500000000
     acell    1.0    1.0    1.0
     rprim
        0.0000000000    5.1085000000    5.1085000000
        5.1085000000    0.0000000000    5.1085000000
        5.1085000000    5.1085000000    0.0000000000

     # K-path in reduced coordinates:
     # tolwfr 1e-20 iscf -2 getden ??
     ndivsm 10
     kptopt -11
     kptbounds
        +0.00000  +0.00000  +0.00000 # $\Gamma$
        +0.50000  +0.00000  +0.50000 # X
        +0.50000  +0.25000  +0.75000 # W
        +0.37500  +0.37500  +0.75000 # K
        +0.00000  +0.00000  +0.00000 # $\Gamma$
        +0.50000  +0.50000  +0.50000 # L
        +0.62500  +0.25000  +0.62500 # U
        +0.50000  +0.25000  +0.75000 # W
        +0.50000  +0.50000  +0.50000 # L
        +0.37500  +0.37500  +0.75000 # K
        +0.62500  +0.25000  +0.62500 # U
        +0.50000  +0.00000  +0.50000 # X


Re-symmetrize a structure when Abinit reports fewer symmetries than expected
----------------------------------------------------------------------------

Crystalline structures saved in text format (e.g. CIF files downloaded from
the Materials Project website) may not have enough significant digits.
Since the default tolerance for symmetry detection in Abinit is tight (tolsym = 1e-8),
Abinit may then find a different space group from the one reported by the source.

In this case, use the `abispg` command of abistruct.py to compute the space group
with Abinit using a tolerance larger than the default value::

    abistruct.py abispg problematic.cif --tolsym=1e-6

Hopefully, the code will detect the correct space group, re-symmetrize
the initial lattice vectors and atomic positions, and print the new symmetrized structure to the terminal.


Get neighbors for each atom in the unit cell out to a distance radius
---------------------------------------------------------------------

To analyze the environment (nearest neighbours) of the atoms in the unit cell
and their coordination, use::

    abistruct.py neighbors sic_relax_HIST.nc

.. code-block:: shell

    Finding neighbors for each atom in the unit cell, out to a distance 2 [Angstrom]

    [0] site PeriodicSite: C (0.0000, -0.0000, 0.0000) [-0.0000, 0.0000, -0.0000] has 4 neighbors:
             PeriodicSite: Si (1.0836, -1.0836, -1.0836) [-0.7500, 0.2500, 0.2500]  at distance 1.8767766107
             PeriodicSite: Si (-1.0836, 1.0836, -1.0836) [0.2500, -0.7500, 0.2500]  at distance 1.8767766107
             PeriodicSite: Si (-1.0836, -1.0836, 1.0836) [0.2500, 0.2500, -0.7500]  at distance 1.8767766107
             PeriodicSite: Si (1.0836, 1.0836, 1.0836) [0.2500, 0.2500, 0.2500]  at distance 1.8767766107

    [1] site PeriodicSite: Si (1.0836, 1.0836, 1.0836) [0.2500, 0.2500, 0.2500] has 4 neighbors:
             PeriodicSite: C (0.0000, 0.0000, 0.0000) [0.0000, 0.0000, 0.0000]  at distance 1.8767766107
             PeriodicSite: C (2.1671, 2.1671, 0.0000) [0.0000, 0.0000, 1.0000]  at distance 1.8767766107
             PeriodicSite: C (2.1671, 0.0000, 2.1671) [0.0000, 1.0000, 0.0000]  at distance 1.8767766107
             PeriodicSite: C (0.0000, 2.1671, 2.1671) [1.0000, 0.0000, 0.0000]  at distance 1.8767766107


Search on the Materials Project database for structures
-------------------------------------------------------

Use::

    abistruct.py mp_search LiF

to search the `materials project`_ database for structures matching a
chemical system or formula, e.g. ``Fe2O3``, ``Li-Fe-O``, or
``Ir-O-*`` for wildcard pattern matching.

The script prints the results to the terminal in tabular form:

.. code-block:: bash

    # Found 2 structures in materials project database (use `verbose` to get full info)
               pretty_formula  e_above_hull  energy_per_atom  \
    mp-1138               LiF      0.000000        -4.845142
    mp-1009009            LiF      0.273111        -4.572031

                formation_energy_per_atom  nsites     volume spacegroup.symbol  \
    mp-1138                     -3.180439       2  17.022154             Fm-3m
    mp-1009009                  -2.907328       2  16.768040             Pm-3m

                spacegroup.number  band_gap  total_magnetization material_id
    mp-1138                   225    8.7161                  0.0     mp-1138
    mp-1009009                221    7.5046                 -0.0  mp-1009009


.. important::

    The script connects to the Materials Project server,
    so you need a ``~/.pmgrc.yaml`` configuration file in your home directory
    containing the authentication token **PMG_MAPI_KEY**.
    For more information, see the
    `pymatgen documentation <http://pymatgen.org/usage.html#pymatgen-matproj-rest-integration-with-the-materials-project-rest-api>`_

The script provides other commands to get (experimental) structures from the COD_ database,
find matching structures on the `materials project`_ website, and generate phase diagrams.
See ``abistruct.py --help`` for more examples.

Compare my structure with the Materials Project database
--------------------------------------------------------

Suppose you have performed a structural relaxation and want
to compare your results with the Materials Project data.
Use the abicomp.py_ script to extract the structure from the HIST.nc_
file and compare it with the database::

    abicomp.py mp_structure ../abipy/data/refs/sic_relax_HIST.nc

To select only the structures with the same space group number as the input structure, use::

    abicomp.py mp_structure ../abipy/data/refs/sic_relax_HIST.nc --same-spgnum

which produces:

.. code-block:: ipython

    Lattice parameters:
            formula  natom  angle0  angle1  angle2         a         b         c  \
    mp-8062  Si1 C1      2    60.0    60.0    60.0  3.096817  3.096817  3.096817
    this     Si1 C1      2    60.0    60.0    60.0  3.064763  3.064763  3.064763

                volume abispg_num spglib_symb  spglib_num
    mp-8062  21.000596       None       F-43m         216
    this     20.355222       None       F-43m         216

    Use --verbose to print atomic positions.

The HIST.nc_ file can be replaced by any other file that provides a structure object.

.. important::

    The Materials Project structures have been obtained with the GGA-PBE functional
    and may include a U term in the Hamiltonian.
    Take these settings into account when comparing structural relaxations.


Visualize the iterations of the SCF cycle
-----------------------------------------

Use::

    abiview.py abo run.abo

to plot the SCF iterations, the steps of a structural relaxation or the DFPT SCF cycles,
depending on the content of run.abo.

You can also use::

    abiview.py log run.log

to print the warnings/comments/errors reported in the Abinit log file ``run.log``.

Export bands to xmgrace format
------------------------------

Both |ElectronBands| and |PhononBands| provide a ``to_xmgrace`` method to produce xmgrace_ files.
To export the data to xmgrace, use abiview.py_ with the ``--xmgrace`` option.
For electrons, use::

    abiview.py ebands out_GSR.nc --xmgrace

and::

    abiview.py phbands out_PHBST.nc --xmgrace

for phonons.

Visualize the Fermi surface
---------------------------

Use::

    abiview.py ebands out_GSR.nc --bxsf

to export the band energies in BXSF format,
which is suitable for visualizing the Fermi surface with xcrysden_.
Then use::

    xcrysden --bxsf BXSF_FILE

to visualize the Fermi surface.
The same file can be produced from Python with:

.. code-block:: ipython

    abifile.ebands.to_bxsf("mgb2.bxsf")

.. important::

    This option requires k-points in the irreducible wedge and a Gamma-centered k-mesh.

Visualize phonon displacements
------------------------------

AbiPy is interfaced with the phononwebsite_ project.
If you have installed the Python package from `GitHub <https://github.com/henriquemiranda/phononwebsite>`_,
you can export the ``PHBST.nc`` file to JSON and then load it in the web interface.
Alternatively, the entire procedure can be automated with the abiview.py_ script.

Use::

    abiview.py phbands out_PHBST.nc -web

to start a local web server and open the HTML page in the default browser
(use the ``--browser`` option to select a different browser).

You can also visualize the phonon modes directly from a DDB_ file with::

    abiview.py ddb out_DDB -web

In this case, AbiPy invokes anaddb to produce the ``PHBST.nc`` file along an automatically
generated q-path and then starts the web server.

Visualize the results of a structural relaxation
------------------------------------------------

The quickest way is to use::

    abiview.py hist out_HIST.nc

to plot the results with matplotlib, or::

    abiopen.py out_HIST.nc -p

to print the most important results to the terminal.

You can also generate an ``XDATCAR`` file with::

    abiview.py hist out_HIST.nc --xdatcar

and visualize the evolution of the crystalline structure with ovito_::

    abiview.py hist out_HIST.nc --appname=ovito

.. important::

    The XDATCAR format assumes a fixed unit cell, so changes in the
    lattice vectors will not be visible in ovito.


Plot results stored in a text file in tabular format
----------------------------------------------------

Use::

    abiview.py data FILE_WITH_COLUMNS

to plot all the columns in the file with matplotlib_.
By default, the first column is used for the x-axis;
use the `--use-index` option to use the row index instead.
Multiple datasets, i.e. blocks of data separated by blank lines, are supported.

To compare multiple files, use::

    abicomp.py data FILE1 FILE2

Standard tools such as gnuplot_ and xmgrace_ work as well, of course, but
the AbiPy scripts are handy for a quick analysis of the results.

Compare multiple files
----------------------

The abicomp.py_ script is designed specifically for this task.
It operates on multiple files (usually with the same extension) and
either produces matplotlib_ plots or creates AbiPy robots with methods
to analyze the results, perform convergence studies and build pandas DataFrames_.

The ``COMMAND`` argument specifies the quantity to compare and is followed by a list of filenames.

For example, to compare the structure in an Abinit input file with the structure
stored in a GSR.nc_ file, use::

    abicomp.py structure run.abi out_GSR.nc

.. note::

    In this example, we can mix files of different types because
    both provide a Structure object. The same philosophy applies to other commands as well:
    everything works as long as AbiPy can extract the quantity of interest from the file.

To plot multiple electronic structures on a grid, use::

    abicomp.py ebands *_GSR.nc out2_WFK.nc -p

Remember that you can use the shell syntax ``*_GSR.nc`` to select all files with a given extension.
For nested directories, use the Unix ``find`` command to scan the directory tree for files matching a pattern.
For example::

    abicomp.py ebands `find . -name *_GSR.nc`

finds all the ``GSR.nc`` files in the current working directory and its subdirectories,
and passes them to the abicomp.py_ script.

.. note::

    Note the `backtick syntax <https://unix.stackexchange.com/questions/27428/what-does-backquote-backtick-mean-in-commands>`_
    used in the command.

Profile the scripts
-------------------

All AbiPy scripts can be run in profile mode by prepending the ``prof`` keyword
to the command-line arguments.
This is useful if a script seems slow and you want to understand why.

Use::

    abiopen.py prof FILE

or::

    abistruct.py prof COMMAND FILE

if the script requires a ``COMMAND`` argument.

Get the description of a variable
---------------------------------

The abidoc.py_ script provides a simplified interface to the Abinit documentation.

Use::

    abidoc.py man ecut

to print the official documentation for the ``ecut`` variable to the terminal.

To list all the variables that depend on the ``natom`` dimension, use::

    abidoc.py withdim natom

More options are available; see ``abidoc.py --help``.

Avoid transferring files from the cluster to localhost just to use matplotlib
-----------------------------------------------------------------------------

Use `SSHFS <https://www.digitalocean.com/community/tutorials/how-to-use-sshfs-to-mount-remote-file-systems-over-ssh>`_
to mount the remote file system over SSH.
You can then run the AbiPy scripts in a terminal on your local machine
to open and visualize the files stored on the cluster.

