.. _graphical-interface:

Graphical interface
===================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

AbiPy provides interactive dashboards that can be used either as a standalone web applications
(**dashboards**) with the `bokeh server <http://docs.bokeh.org/>`_ or inside jupyter notebooks.
This document explains how to install the required dependencies and how to
generate dashboards/GUIs either with the command line interface (CLI) or inside jupyter notebooks.

.. important::

    Please note that one needs a **running python backend**
    to execute the callbacks triggerered by the GUI widgets.
    This part, indeed, is implemented in HTML/CSS/JS code executed
    by the frontend (i.e. **your browser**) that sends the signal
    to the python server (the **backend**).
    The python server is supposed to process the data
    and send the results back to the frontend for visualization purposes

    Don't be surprised if you start to click buttons and **nothing happens**!
    The examples provided in this page are only meant to show how to build GUI
    or dashboards with AbiPy.


Installation
------------

Install the `panel <https://panel.pyviz.org/>`_  package either from pip with:

.. code-block:: bash

    pip install panel

or with conda (**recommended**) using:

.. code-block:: bash

    conda install panel -c conda-forge

If you plan to use panel within JupyterLab, you will also need to install
the PyViz JupyterLab extension and activate it with:

.. code-block:: bash

    conda install -c conda-forge jupyterlab
    jupyter labextension install @pyviz/jupyterlab_pyviz


Basic Usage
-----------

Several AbiPy objects provide a ``get_panel`` method returning
an object that can be displayed inside the jupyter notebook or inside a Bokeh server.
When running inside a jupyter notebook, remember enable the integration
with the ``panel`` infrastructure by executing:

.. jupyter-execute::

    from abipy import abilab
    abilab.abipanel();

**before calling** any AbiPy ``get_panel`` method.

.. note::

    The ``abipanel`` function is needed to load extensions and javascript packages
    required by AbiPy.
    This function is just a small wrapper around the official panel API:

    .. code-block:: bash

        import panel as pn
        pn.extension()


At this point, we can start to construct AbiPy objects.
For our first example, we use the ``abiopen`` function to open a ``GSR`` file,
then we call ``get_panel`` to build a set of widgets that allows us to interact
with the `GsrFile`:

.. jupyter-execute::

    from abipy import abilab
    import abipy.data as abidata

    filename = abidata.ref_file("si_nscf_GSR.nc")
    gsr = abilab.abiopen(filename)

    gsr.get_panel()

The **summary** tab provides a string representation of the file
but there is not widget to interact with it.
If you select the **e-Bands** tab, you will see several widgets and a button
that activates the visualization of the KS band energies.
Again, in this HTML page there is no python server running in background so
if you click the **Plot e-bands** button nothing happens (this is not a bug!).

The advantage of this notebook-based approach is that it is possible to mix
the panel GUIs with python code that can be used to perform
more advanced tasks not supported by the GUI.

Obviously it is possible to have multiple panels running in the same notebook.
Calling ``get_structure`` with an AbiPy structure, for instance, creates a set of widgets
to facilitate common operations such as exporting the structure to a different format or
generating a basic Abinit input file for e.g. GS calculations:

.. jupyter-execute::

    gsr.structure.get_panel()

.. note::

    At present, not all the AbiPy objects support the ``get_panel`` protocol
    but we plan to gradually support more objects, especially the most important
    netcdf files produced by Abinit

To generate a notebook from the command line, use the abiopen.py_ script:

.. code-block:: bash

    abiopen.py si_nscf_GSR.nc -nb  # short for --notebook

that will automatically open the notebook inside jupyterlab.
If you prefer classic jupyter notebooks, use the ``-nb --classic-notebook`` options

If you do not need to execute python code, you may want to generate a panel dashboard with:

.. code-block:: bash

    abiopen.py si_nscf_GSR.nc -pn  # short for --panel

The same approach can be used with a ``DDB`` file.
In this case, we get more tabs and options because one can use the GUI
to set the input parameters, invoke ``anaddb`` and visualize the results:

.. jupyter-execute::

    # Open DDB file with abiopen and invoke get_panel method.
    ddb_path = abidata.ref_file("mp-1009129-9x9x10q_ebecs_DDB")
    abilab.abiopen(ddb_path).get_panel()

The same result can be obtained from the CLI with

.. code-block:: bash

    abiopen.py mp-1009129-9x9x10q_ebecs_DDB -nb

There are, however, cases in which you don't need the interactive environment provided
by jupyter notebooks as you are mainly interested in the visualization of the results.
In this case, it is possible to use the command line interface to automatically generate
a dashboard with widgets without having to start a jupyter-lab application.

To build a dashboard for a ``Structure`` object extracted from ``FILE``, use::

    abistruct.py panel FILE

where ``FILE`` is **any** file providing a ``Structure`` object
e.g. netcdf files, cif files, abi, abo files etc.

To build a dashboard associated to one of the AbiPy file, use the syntax::


    abiopen.py FILE --panel

where ``FILE`` is one of the Abinit files supported by ``abiopen.py``.
For instance, one can create a dashboard to interact with a ``DDB`` file with::

    abiopen.py out_DDB --panel


.. important::

    To build a dashboard for an AbiPy Flow use::

        abirun.py FLOWDIR panel

    or alternatively::

        abiopen.py FLOWDIR/__AbinitFlow__.pickle --panel
