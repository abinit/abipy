========
Overview
========

AbiPy is a Python_ package for analyzing the results produced by Abinit_,
an open-source program for the ab-initio calculation of the physical properties of materials
within Density Functional Theory and Many-Body Perturbation Theory.
AbiPy also provides tools to generate input files and to build workflows that automate
ab-initio calculations and typical convergence studies.
Since AbiPy is interfaced with pymatgen_, users can also take advantage of
the many tools and Python objects available in the pymatgen ecosystem.

AbiPy works well together with matplotlib_, pandas_, seaborn_,
ipython_ and jupyter_, providing a powerful and user-friendly environment for data analysis and visualization.
Check out the plotting scripts available in our :ref:`plot-gallery`.
To learn more about how AbiPy integrates with jupyter_, browse `our collection of notebooks
<https://nbviewer.jupyter.org/github/abinit/abitutorials/blob/master/abitutorials/index.ipynb>`_
or click the **Launch Binder** badge to start a Docker image with Abinit, AbiPy and all the Python dependencies
needed to run the notebooks.
Once the image is built, the notebook opens in your browser.

Note that most of the post-processing tools in AbiPy require Abinit output files in
netcdf_ format, so we strongly recommend compiling Abinit with netcdf support.
Use ``--with_trio_flavor="netcdf-fallback"`` at configure time to activate the internal netcdf library.
To link Abinit against an external netcdf library, consult the configuration examples
provided by the abiconfig_ package.

AbiPy is free to use, and we welcome contributions that help improve the library.
Please report bugs and issues on AbiPy's `GitHub page <https://github.com/abinit/abipy>`_.
