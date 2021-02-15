========
Overview
========

AbiPy is a Python_ package to analyze the results produced by Abinit_
an open-source program for the ab-initio calculations of the physical properties of materials
within Density Functional Theory and Many-Body perturbation theory.
AbiPy also provides tools to generate input files and workflows to automate
ab-initio calculations and typical convergence studies.
AbiPy is interfaced with pymatgen_ thus allowing users to
benefit from the different tools and python objects available in the pymatgen ecosystem.

AbiPy can be used in conjunction with  matplotlib_, pandas_, seaborn_,
ipython_ and jupyter_ thus providing a powerful and user-friendly environment for data analysis and visualization.
Check out the list of plotting scripts available in our :ref:`plot-gallery`.
To learn more about the integration between jupyter_ and AbiPy, visit `our collection of notebooks
<https://nbviewer.jupyter.org/github/abinit/abitutorials/blob/master/abitutorials/index.ipynb>`_
or click the **Launch Binder** badge to start a Docker image with Abinit, AbiPy and all the other python dependencies
required to run the code inside the jupyter notebooks.
The notebook will be opened in your browser after building.

.. AbiPy supports both Python 2.7 as well as Python >= 3.4.
.. Note however that Python 2.7 is more intensively tested than py3k, especially at the level of workflows.
.. We hence still recommend py2.7 if you plan to run automatic calculations with AbiPy.

Note also that the majority of the post-processing tools available in AbiPy require Abinit output files in
netcdf_ format so we strongly suggest to compile Abinit with netcdf support
(use ``--with_trio_flavor="netcdf-fallback"`` at configure time to activate the internal netcdf library,
to link Abinit against an external netcdf library please consult the configuration examples
provided by the abiconfig_ package.

AbiPy is free to use. However, we also welcome your help to improve this library by making your own contributions.
Please report any bugs and issues at AbiPy's `Github page <https://github.com/abinit/abipy>`_.
