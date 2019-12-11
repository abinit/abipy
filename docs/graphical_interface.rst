.. _graphical-interface:

Graphical interface
===================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

AbiPy provides interactive dashboards that can be used either as a standalone web apps 
with the `bokeh server <http://docs.bokeh.org/>`_ or inside jupyter notebooks.
This document explains how to install the required dependencies and how to
generate dashboards either with the command line interface or inside jupyter notebooks.

.. important::

    Note that you will need a running python kernel to execute the callbacks triggerered 
    by the GUI hence the examples given in this page are only meant to show how to build the the GUI.

Installation
------------

Install the `panel <https://panel.pyviz.org/>`_  package either from pip with:

.. code-block:: bash

  pip install panel

or with conda:

.. code-block:: bash

  conda install panel -c conda-forge

If you want to work with JupyterLab, you will also need to install 
the optional PyViz JupyterLab extension:

.. code-block:: bash

    conda install -c conda-forge jupyterlab
    jupyter labextension install @pyviz/jupyterlab_pyviz


Basic Usage
-----------

The AbiPy structure and many AbiPY files provide a ``get_panel`` method that returns 
a dashboard that can be used inside the jupyter notebook.
To enable the integration with ``panel`` inside a jupyter notebook, execute the below code:

.. jupyter-execute::

    # Import panel and activate extensions
    import panel as pn
    pn.extension()

Now one can start to construct AbiPy objects and use the ``get_panel`` method to generate graphical interfaces. 
In our first example, we use the ``abiopen`` function to open a ``GSR`` file 
and then we call ``get_panel`` to build a set of widgets that allows us to interact with the object:


.. jupyter-execute::

    # Import AbiPy modules.
    from abipy import abilab
    import abipy.data as abidata

    filename = abidata.ref_file("si_nscf_GSR.nc")
    gsr = abilab.abiopen(filename)

    gsr.get_panel()


The same approach can be used with a ``DDB`` file.
In this case, we get more tabs and options because one can use the GUI 
to set the input parameters, invoke ``anaddb`` and visualize the results:


.. jupyter-execute::

    # Open DDB file with abiopen and invoke get_panel method.
    #ddb_path = abidata.ref_file("mp-1009129-9x9x10q_ebecs_DDB")

    #abilab.abiopen(ddb_path).get_panel()


Calling ``get_structure`` with an AbiPy structure, creates a set of widgets 
to facilitate common operations such as exporting to a different format or
generating a Abinit input file for GS calculations:

.. jupyter-execute::

    gsr.structure.get_panel()


There are, however, cases in which you don't need the interactive environment provided 
by jupyter notebooks as you are mainly interested in the visualization of the results.
In this case, it is possible to use the command line interface to automatically generate 
a dashboard with widgets without having to start a jupyter-lab application.

To build a dashboard for a ``Structure`` object extract from ``FILE``, use:

.. code-block:: shell

    abistruct.py panel FILE


where ``FILE`` is any file providing a ``Structure`` object e.g. netcdf file, cif files, abi, abo etc.

To build a dashboard associated to one of the AbiPy file, use the syntax:

.. code-block:: shell

    abiopen.py FILE --panel


where ``FILE`` is one of the Abinit files supported by ``abiopen.py``.
For instance, one can create a dashboard to interact with a ``DDB`` file with:

.. code-block:: shell

    abiopen.py out_DDB --panel

.. important::

    To build a dashboard for an AbiPy Flow use:

        abirun.py FLOWDIR panel

    or, alternatively:

        abiopen.py FLOWDIR/__AbinitFlow__.pickle --panel

.. jupyter-execute::

    import numpy as np
    from matplotlib import pyplot
    %matplotlib inline

    x = np.linspace(1E-3, 2 * np.pi)

    pyplot.plot(x, np.sin(x) / x)
    pyplot.plot(x, np.cos(x))
    pyplot.grid()
