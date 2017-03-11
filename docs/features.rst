================
Feature Overview
================

AbiPy is a Python library to analyze the results produced by `ABINIT <http://www.abinit.org>`_,
an open-source program for the ab-initio calculations of the physical properties of materials
within Density Functional Theory and Many-Body perturbation theory.
It also provides tools to generate input files and workflows to automate
ab-initio calculations and typical convergence studies.
AbiPy is interfaced with `Pymatgen <http://www.pymatgen.org>`_ and this allows users to
benefit from the different tools and python objects available in the pymatgen ecosystem.

AbiPy can be used in conjunction with  `matplotlib <http://matplotlib.org>`_, `pandas <http://pandas.pydata.org>`,
`ipython <https://ipython.org/index.html>`_ and `jupyter <http://jupyter.org/>`_
thus providing a powerful and user-friendly environment for data analysis and visualization.
Check out the list of plotting scripts available in our :doc:`gallery </examples/index>`.
To learn more about the integration between jupyter and AbiPy, visit our collection of `notebooks
<http://nbviewer.ipython.org/github/abinit/abipy/blob/master/abipy/examples/notebooks/index.ipynb>`_
and the
`AbiPy lessons <http://nbviewer.ipython.org/github/abinit/abipy/blob/master/abipy/examples/notebooks/lessons/index.ipynb>`_.

AbiPy supports both Python 2.7 as well as Python >= 3.4.
Note however that Python 2.7 is more intensively tested than py3k especially at the level of workflows
so we still recommend py2.7 if you plan to run automatic calculations with AbiPy.

Note also that the majority of the post-processing tools available in AbiPy require output files in
``netcdf`` format so we strongly suggest to compile Abinit with netcdf support
(use ``--with_trio_flavor="netcdf-fallback"`` at configure time to activate the internal netcdf library,
to link Abinit against an external netcdf library please consult the configuration examples
provided by the `abiconfig package <https://github.com/abinit/abiconfig>`_.


AbiPy is free to use. However, we also welcome your help to improve this library by making your own contributions.
Please report any bugs and issues at AbiPy's `Github page <https://github.com/abinit/abipy>`_.
