.. _post-processing-howto:

***********************
Post-processing How-To
***********************

This is a list of Frequently Asked Questions about the usage of the AbiPy scripts. 
Feel free to suggest new entries!

.. contents::
   :backlinks: top

.. important::

    The AbiPy scripts detect the file type by looking at the file extension to so don't change it.

    Also, remember that it is possible to get the documentation 
    by just adding ``--help`` after the scripts name.

    For example::

        abistruct.py --help
        
    gives the documentation and usage examples for the ``abistruct.py`` script while::

        abistruct.py COMMAND --help
    
    prints the documentation for ``COMMAND`` and the options supported by ``COMMAND``.


Get information about a generic ``FILE``
----------------------------------------

Use::

    abiopen.py FILE --print 

to print information about a file inside the terminal or

    abiopen.py FILE --expose

to generate multiple matplotlib figures depending on FILE.

Use ``--verbose`` or ``-v`` to increase verbosity level. 
The option can be can be supplied multiple times e.g. ``-vv``.

Get all file extensions supported by ``abiopen.py``
---------------------------------------------------

Use::

    abiopen.py --help

.. command-output:: abiopen.py --help


Convert the structure stored in ``FILE`` to a different format 
--------------------------------------------------------------

Use::

    abistruct.py convert FILE

to read the structure from ``FILE`` and generate a CIF_ file (default behaviour).

The majority of the netcdf_ files produced by Abinit contain structural information
so this command can be used with netcdf ouput files as well as Abinit input/output 
files and all the other formats supported by pymatgen e.g. POSCAR files.
Other formats can be specified with the ``-f`` option.
For example::

    abistruct.py convert FILE -f abivars

gives the Abinit input variables while::

    abistruct.py convert FILE -f xsf > out.xsf

exports the structure to the ``xsf`` format (xcrysden_) and save it to file.

Use::

    abistruct.py convert --help

to list the formats supported.

Check if my Abinit input file is OK
-----------------------------------

First of all, one can use::

    abiopen.py ../abipy/data/refs/si_ebands/run.abi -p

to print the crystalline structure and find the space group with the spglib_ library.

If the structure looks good, use the :ref:`abinp.py` script with the ``validate`` command as in::

    abinp.py validate run.abi       
    
to validate the input file with Abinit (requires ``manager.yml`` and, obviously, Abinit).

The script provides other options to invoke Abinit 
to get space group information, the list of k-points in the IBZ.
the list of atomic perturbations for phonons or the list of autoparal configurations.
See ``abinp.py --help`` for futher info.

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

A similar interface is also available via::

    abiview.py log run.log


Get a quick look to a file
--------------------------

The :ref:`abiview.py` script is especially designed for this task.
The syntax is ``abistruct.py COMMAND FILE`` where ``COMMAND`` is either 
the Abinit file extension (without ``.nc``, if any) or the AbiPy object we want to visualize.

To get a quick look at the DDB file, use::

    abiview.py ddb out_DDB

This command invokes anaddb to compute phonon bands and DOS from the DDB and produces matplotlib_ plots.

If ``FILE`` contains electronic band energies, use e.g.::

    abiview.py ebands out_GSR.nc

to plot the KS eigenvalues (the same command works for other files such as ``WFK.nc``, ``DEN.nc`` etcetera. 

Note that :ref:`abiview.py` uses a predefined logic to visualize the data.
There are options to tune some parameters and/or export data in different formats
but exposing the AbiPy API from the command line is not easy.

For a more flexible interface, we suggest to use::

    abiopen.py FILE

to start an ipython_ shell that will allow you to interact with the python object directly.

If you have installed jupyter_ on your machine/cluster and you have a web browser, use::

    abiopen.py FILE -nb

to generate automatically a predefined jupyter notebook associated to the file type.

Visualize a structure
---------------------

The visualization of the structure is delegated to external graphical applications
that must be istalled on your machine. 
AbiPy will extract the structure from ``FILE``, convert it to one of the formats 
supported by the graphical application and finally invoke the executable.
If you have vesta_ installed in one of the standard 
locations of your machine, you can simply execute::

    abistruct.py visualize FILE

inside the terminal. 
Other applications can be specified with the ``--application`` option.
At present, AbiPy supports vesta_, ovito_, xcrysden_, avogadro_, and v_sim_.

To visualize the crystalline structure inside the jupyter_ notebook, you may want to
try the nbjsmol_ jupyter extension.

Get a high-symmetry kpath for a given structure
-----------------------------------------------

Use the `kpath` command and pass a FILE providing structural info::

    abistruct.py kpath FILE

to generate a template with the input variables defining the k-path

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


Get neighbors for each atom in the unit cell out to a distance radius
---------------------------------------------------------------------

If we are interested in the environment/nearest neighbours of the atoms in the unit cell,
we can analyze the different coordinations with::

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

to search on the `materials project`_ database for structures corresponding to a 
chemical system or formula e.g. ``Fe2O3`` or ``Li-Fe-O`` or
``Ir-O-*`` for wildcard pattern matching.

The script prints the results to terminal in tabular form:

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

    The script will try to connect to the materials project server.
    You need a ``~/.pmgrc.yaml`` configuration file inside your home directory
    with the authentication token **PMG_MAPI_KEY**.
    For further info please refer to the 
    `pymatgen documentation <http://pymatgen.org/usage.html#pymatgen-matproj-rest-integration-with-the-materials-project-rest-api>`_

The script provides other commands to get (experimental) structures from the COD_ database,
find matching structures on the `materials project`_ website and generate phase diagrams.
See ``abistruct.py --help`` for more examples.

Compare my structure with the Materials Project database
--------------------------------------------------------

Let's assume we have performed a structural relaxation and we want
to compare our results with the Materials Project data.
One can use the :ref:`abicomp.py` structure to extract the structure from the HIST.nc_
file and compare the data with the database::

    abicomp.py mp_structure ../abipy/data/refs/sic_relax_HIST.nc

It's also possible to select only the structures with the same space group number as the input structure with::

    abicomp.py mp_structure ../abipy/data/refs/sic_relax_HIST.nc --same-spgnum

that produces

.. code-block:: ipython

    Lattice parameters:
            formula  natom  angle0  angle1  angle2         a         b         c  \
    mp-8062  Si1 C1      2    60.0    60.0    60.0  3.096817  3.096817  3.096817
    this     Si1 C1      2    60.0    60.0    60.0  3.064763  3.064763  3.064763

                volume abispg_num spglib_symb  spglib_num
    mp-8062  21.000596       None       F-43m         216
    this     20.355222       None       F-43m         216

    Use --verbose to print atomic positions.

Note that one can replace the HIST.nc_ file with any other file providing a structure object.

.. important::

    The structures of the materials project have been obtained with the GGA-PBE functional
    and they might include the U term in the Hamiltonian.
    One should take into account these different settings when comparing structural relaxations.


Visualize the iterations of the SCF cycle
-----------------------------------------

Use::

    abiview.py abo run.abo

to plot the SCF iterations or the steps of the structural relaxations or the DFPT SCF cycles
(depending on the content of run.abo).

Note that one can also use::
    
    abiview.py log run.log

to print the warnings/comments/errors reported in the Abinit log file ``run.log``.

Export bands to xmgrace format 
------------------------------

But |ElectronBands| and |PhononBands| provide a ``to_xmgrace`` method to produce xmgrace_ files.
To export the data to xmgrace, use :ref:`abiview.py` with the ``--xmgrace`` option.
For electrons, use::

    abiview.py ebands out_GSR.nc --xmgrace

and::

    abiview.py phbands out_PHBST.nc -xmgrace 

for phonons.

Visualize the Fermi surface
---------------------------

Use::

    abiview.py ebands out_GSR.nc --bxsf

to export a set of band energies in BXSF format
suitable for the visualization of the Fermi surface with xcrysden_.
Then use::

    xcrysden --bxsf BXSF_FILE

to visualize the Fermi surface with xcrysden_

.. code-block:: ipython

    abifile.ebands.to_bxsf("mgb2.bxsf")    

.. important::

    This option requires k-points in the irreducible wedge and a gamma-centered k-mesh.

Visualize phonon displacements
------------------------------

AbiPy is interfaced with the phononwebsite_ project 
If you have already installed the python package from `github <https://github.com/henriquemiranda/phononwebsite>`_
it's possbile to export the ``PHBST.nc`` to JSON and then load the file via the web-interface.
Alternatively, it's possible to automate the entire procedure with the :ref:`abiview.py` script.

Use::

    abiview.py phbands out_PHBST.nc -web

to start a local webserver and open the HTML page inside the default browser 
(the browser can be changed with the ``--browser`` option).

It is also possible to visualize the phonon modes starting directly from a DDB_ file with::

    abiview.py ddb -web

In this case, AbiPy will invoke anaddb to produce the ``PHBST.nc`` file on an automatically 
generated q-path and then start the webserver.

Visualize the results of a structural relaxation
------------------------------------------------

The quickest way is to use::

    abiview hist out_HIST.nc

to plot the results with matplotlib or:: 

    abiopen.py out_HIST.nc -p
    
to print the most important results to terminal.

Note that it's possible to generate a ``XDATCAR`` file with::

    abiview hist out_HIST.nc --xdatcar

and visualize the crystalline structure with ovito_::

    abiview hist out_HIST.nc -appname=ovito

.. important::

    The XDATCAR format assumes a fixed unit cell so you won't be able
    to visualize the modifications of the unit cell lattice vectors in ovito.

Compare multiple files
----------------------

The :ref:`abicomp.py` script is explicitely designed for this kind of task.
It operates on multiple files (usually files with the same extension) and 
either produces matplotlib_ plots or creates AbiPy robots providing methods
to analyze the results, perform convergence studies and build pandas DataFrames_.

The ``COMMAND`` defines the quantity to be compared, followed by a list of filenames.

To compare e.g. the structure given in one Abinit input file with the structure
coming from a GSR.nc_ file, use::

    abicomp.py structure run.abi out_GSR.nc

.. note::

    In this example, we can use files of different type because they
    both have a Structure object. This philosophy can be applied to other commands as well:
    everything works as long as AbiPy is able to extract the quantity of interest from the file.

To plot multiple electronic structures on a grid, use::

    abicomp.py ebands *_GSR.nc out2_WFK.nc -p

Remember that it is possible to use the shell syntax ``*_GSR.nc`` to select all files with a given extension.
If you have nested directories, use unix ``find`` to scan the directory tree for files matching a given pattern
For example::

    abicomp.py ebands `find . -name *_GSR.nc` 

finds all ``GSR.nc`` files contained withing the current working directory.
The output of ``find`` is then passed to the :ref:`abicomp.py` script.

.. note::

    Note the `backticks syntax <https://unix.stackexchange.com/questions/27428/what-does-backquote-backtick-mean-in-commands>`_
    used in the command.

Profile the scripts
-------------------

All AbiPy script can be executed in profile mode by just prepending the ``prof`` keyword  
to the command line arguments. 
This option could be useful if the script seems to be slow and you need to understand what's happening.

Use::

    abiopen.py prof FILE

or::

    abistruct.py prof COMMAND FILE

if the script requires a ``COMMAND`` argument.

Get the description of a variable 
---------------------------------

The :ref:`abidoc.py` script provides a simplified interface to the Abinit documentation.

Use::

    abidoc.py man ecut

to print the official documentation for the ``ecut`` variable to the terminal.

To list all the variables depending on the ``natom`` dimension, use::

    abidoc.py withdim natom

More options are available. See ``abidoc.py --help``.
