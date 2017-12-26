.. Links to websites
.. _Sphinx: http://sphinx.pocoo.org
.. _Abinit: https://www.abinit.org
.. _abiconfig: https://github.com/abinit/abiconfig
.. _abiflows: https://github.com/abinit/abiflows
.. _abitutorials: https://github.com/abinit/abitutorials
.. _abiconda: https://github.com/abinit/abiconda
.. _pseudo-dojo: http://www.pseudo-dojo.org/
.. _pseudo-dojo repository: https://github.com/abinit/pseudo_dojo
.. _pymatgen: http://www.pymatgen.org
.. _fireworks: https://materialsproject.github.io/fireworks/
.. _mongodb: https://www.mongodb.com/
.. _`materials project`: https://materialsproject.org/
.. _conda: https://conda.io/docs/
.. _Anaconda: https://continuum.io/downloads
.. _abinit-channel: https://anaconda.org/abinit
.. _matsci: http://materials.sh/
.. _spack: https://github.com/LLNL/spack
.. _matplotlib: http://matplotlib.org
.. _pandas: http://pandas.pydata.org
.. _scipy: https://www.scipy.org/
.. _numpy: http://www.numpy.org/
.. _seaborn: https://seaborn.pydata.org/
.. _simpy: https://simpy.readthedocs.io/en/latest/
.. _networkx: https://networkx.github.io/
.. _pytest: https://docs.pytest.org/en/latest/contents.html
.. _netcdf4-python: http://unidata.github.io/netcdf4-python/
.. _nbformat: https://github.com/jupyter/nbformat
.. _pip: https://pypi.python.org/pypi/pip
.. _ipython: https://ipython.org/index.html
.. _jupyter: http://jupyter.org/
.. _Python: http://www.python.org/
.. _spglib: https://atztogo.github.io/spglib/
.. _ase: https://wiki.fysik.dtu.dk/ase/
.. _COD: http://www.crystallography.net/cod/
.. _CIF: http://www.iucr.org/resources/cif
.. _vesta: http://jp-minerals.org/vesta/en/
.. _xcrysden: http://www.xcrysden.org/
.. _xmgrace: http://plasma-gate.weizmann.ac.il/Grace/
.. _gnuplot: http://www.gnuplot.info/
.. _ovito: https://ovito.org/
.. _v_sim: http://inac.cea.fr/L_Sim/V_Sim/
.. _mayavi: http://docs.enthought.com/mayavi/mayavi/
.. _avogadro: https://avogadro.cc/
.. _nbjsmol: https://github.com/gmatteo/nbjsmol
.. _phononwebsite: http://henriquemiranda.github.io/phononwebsite/
.. _netcdf: https://www.unidata.ucar.edu/software/netcdf/docs/faq.html#whatisit
.. _YAML: https://en.wikipedia.org/wiki/YAML
.. _JSON: https://en.wikipedia.org/wiki/JSON
.. _slurm: https://slurm.schedmd.com/
.. _pbspro: http://pbspro.org/
.. _sge: http://gridscheduler.sourceforge.net/howto/GridEngineHowto.html
.. _torque: http://www.adaptivecomputing.com/products/open-source/torque/
.. _moab: http://www.adaptivecomputing.com/products/hpc-products/moab-hpc-basic-edition/
.. _loadleveler: https://www.ibm.com/support/knowledgecenter/en/SSFJTW

.. Links to important python objects.
.. _POSCAR: http://cms.mpi.univie.ac.at/vasp/guide/node59.html
.. _DataFrame: https://pandas.pydata.org/pandas-docs/stable/generated/pandas.DataFrame.html
.. _DataFrames: https://pandas.pydata.org/pandas-docs/stable/generated/pandas.DataFrame.html

.. Links to jupyter notebooks associated to AbiPy files available at
   https://nbviewer.jupyter.org/github/abinit/abitutorials/blob/master/abitutorials/index.ipynb?flush_cache=true
.. _AbipyStructure: https://nbviewer.jupyter.org/github/abinit/abitutorials/blob/master/abitutorials/structure.ipynb
.. _AbinitInput: https://nbviewer.jupyter.org/github/abinit/abitutorials/blob/master/abitutorials/abinit_input.ipynb
.. _GSR.nc: https://nbviewer.jupyter.org/github/abinit/abitutorials/blob/master/abitutorials/gsr.ipynb
.. _HIST.nc: https://nbviewer.jupyter.org/github/abinit/abitutorials/blob/master/abitutorials/hist.ipynb
.. _FATBANDS.nc: https://nbviewer.jupyter.org/github/abinit/abitutorials/blob/master/abitutorials/efatbands.ipynb
.. _DDB: https://nbviewer.jupyter.org/github/abinit/abitutorials/blob/master/abitutorials/ddb.ipynb
.. _SIGRES.nc: https://nbviewer.jupyter.org/github/abinit/abitutorials/blob/master/abitutorials/sigres.ipynb

.. Important Abipy objects.
.. |AttrDict| replace:: :class:`monty.collections.AttrDict`
.. |Function1D| replace:: :class:`abipy.core.func1d.Function1D`
.. |Mesh3d| replace:: :class:`abipy.core.mesh3d.Mesh3d`
.. |GSphere| replace:: :class:`abipy.core.gsphere.GSphere`
.. |Kpoint| replace:: :class:`abipy.core.kpoints.Kpoint`
.. |KpointList| replace:: :class:`abipy.core.kpoints.KpointList`
.. |Kpath| replace:: :class:`abipy.core.kpoints.Kpath`
.. |IrredZone| replace:: :class:`abipy.core.kpoints.IrredZone`
.. |KpointStar| replace:: :class:`abipy.core.kpoints.KpointStar`
.. |Structure| replace:: :class:`abipy.core.structure.Structure`
.. |pymatgen-Structure| replace:: :class:`pymatgen.core.structure.Structure`
.. |Lattice| replace:: :class:`pymatgen.core.lattice.Lattice`
.. |AbinitInput| replace:: :class:`abipy.abio.inputs.AbinitInput`
.. |MultiDataset| replace:: :class:`abipy.abio.inputs.MultiDataset`
.. |ElectronBands| replace:: :class:`abipy.electrons.ebands.ElectronBands`
.. |SkwInterpolator| replace:: :class:`abipy.core.skw.SkwInterpolator`
.. |ElectronDos| replace:: :class:`abipy.electrons.ebands.ElectronDos`
.. |PhononBands| replace:: :class:`abipy.dfpt.phonons.PhononBands`
.. |ScfTask| replace:: :class:`pymatgen.io.abinit.tasks.ScfTask`
.. |NscfTask| replace:: :class:`pymatgen.io.abinit.tasks.NscfTask`
.. |TaskManager| replace:: :class:`pymatgen.io.abinit.tasks.TaskManager`
.. |GsrFile| replace:: :class:`abipy.electrons.gsr.GsrFile`
.. |GsrRobot| replace:: :class:`abipy.electrons.gsr.GsrRobot`
.. |DdbFile| replace:: :class:`abipy.dfpt.ddb.DdbFile`
.. |DdbRobot| replace:: :class:`abipy.dfpt.ddb.DdbRobot`
.. |PhbstFile| replace:: :class:`abipy.dfpt.phonons.PhbstFile`
.. |PhdosFile| replace:: :class:`abipy.dfpt.phonons.PhdosFile`
.. |PhononDos| replace:: :class:`abipy.dfpt.phonons.PhononDos`
.. |PhononBandsPlotter| replace:: :class:`abipy.dfpt.phonons.PhononBandsPlotter`
.. |PhononDosPlotter| replace:: :class:`abipy.dfpt.phonons.PhononDosPlotter`
.. |Pseudo| replace:: :class:`pymatgen.io.abinit.pseudos.Pseudo`
.. |PseudoTable| replace:: :class:`pymatgen.io.abinit.pseudos.PseudoTable`
.. |Visualizer| replace:: :class:`abipy.iotools.visualizer.Visualizer`

.. Important objects provided by libraries.
.. |matplotlib-Figure| replace:: :class:`matplotlib.figure.Figure`
.. |matplotlib-Axes| replace:: :class:`matplotlib.axes.Axes`
.. |pandas-DataFrame| replace:: :class:`pandas.DataFrame`
.. |pandas-DataFrames| replace:: :class:`pandas.DataFrame`
.. |numpy-array| replace:: :class:`numpy.ndarray`
