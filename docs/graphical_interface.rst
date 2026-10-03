.. _graphical-interface:

Graphical interface
===================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

AbiPy provides interactive dashboards that can be used either as standalone web applications
served by the `bokeh server <http://docs.bokeh.org/>`_ or inside jupyter notebooks.
This page explains how to install the required dependencies and how to
generate dashboards from the command line interface (CLI) or inside jupyter notebooks.

.. important::

    The callbacks triggered by the widgets in the HTML page require a **running Python backend**.
    The widgets are implemented in HTML/CSS/JS code executed
    by the frontend (i.e. **your browser**), which sends signals
    to the Python server (the **backend**).
    The server processes the data
    and sends the results back to the frontend for visualization.

    This static page has no backend, so don't be surprised if **nothing happens** when you click the buttons!
    The examples on this page are only meant to show how to build GUIs
    and dashboards with AbiPy.


Installation
------------

Install the `panel <https://panel.pyviz.org/>`_ package with pip:

.. code-block:: bash

    pip install panel

or with conda (**recommended**):

.. code-block:: bash

    conda install panel -c conda-forge

If you plan to use panel within JupyterLab, you also need to install
and activate the PyViz JupyterLab extension:

.. code-block:: bash

    conda install -c conda-forge jupyterlab
    jupyter labextension install @pyviz/jupyterlab_pyviz


Basic Usage
-----------

Several AbiPy objects provide a ``get_panel`` method that returns
an object that can be served to a web browser or displayed inside a jupyter notebook.
When working in a jupyter notebook, remember to enable the integration
with ``panel`` by executing:

.. jupyter-execute::

    from abipy import abilab
    abilab.abipanel();

**before calling** any AbiPy ``get_panel`` method.

.. note::

    The ``abipanel`` function loads the extensions and javascript packages
    required by AbiPy.
    It is just a small wrapper around the panel API:

    .. code-block:: bash

        import panel as pn
        pn.extension()


We can now start building AbiPy objects.
In our first example, we use the abiopen function to open a ``GSR`` file
and then call ``get_panel`` to build a set of widgets for interacting
with the |GsrFile|:

.. jupyter-execute::

    from abipy import abilab
    import abipy.data as abidata

    filename = abidata.ref_file("si_nscf_GSR.nc")
    gsr = abilab.abiopen(filename)

    gsr.get_panel()

The **Summary** tab shows a string representation of the file
and has no interactive widgets.
The **e-Bands** tab, on the other hand, contains several widgets and a button
that plots the KS band energies.
Since no Python server is running behind this HTML page,
clicking the **Plot e-bands** button does nothing (this is not a bug!).

The advantage of the notebook-based approach is that you can mix
the panel GUIs with Python code to perform
more advanced tasks not supported by the GUI.

You can also have multiple panels in the same notebook.
For instance, calling ``get_panel`` on an AbiPy structure creates a set of widgets
for common operations such as exporting the structure to a different format or
generating a basic Abinit input file, e.g. for GS calculations:

.. jupyter-execute::

    gsr.structure.get_panel()

.. note::

    Not all AbiPy objects support the ``get_panel`` protocol yet,
    but we plan to gradually extend it to more objects, starting from the most important
    netcdf files produced by Abinit.

To generate a notebook from the command line, use the abiopen.py_ script:

.. code-block:: bash

    abiopen.py si_nscf_GSR.nc -nb  # short for --notebook

which automatically opens the notebook in JupyterLab.
If you prefer the classic jupyter notebook, use the ``-nb --classic-notebook`` options.

If you do not need to execute Python code, you can generate a panel dashboard instead with:

.. code-block:: bash

    abiopen.py si_nscf_GSR.nc -pn  # short for --panel

The same approach works with a ``DDB`` file.
In this case, there are more tabs and options, since the GUI can be used
to set the input parameters, invoke ``anaddb`` and visualize the results:

.. jupyter-execute::

    # Open DDB file with abiopen and invoke get_panel method.
    ddb_path = abidata.ref_file("mp-1009129-9x9x10q_ebecs_DDB")
    abilab.abiopen(ddb_path).get_panel()

The same result can be obtained from the CLI with:

.. code-block:: bash

    abiopen.py mp-1009129-9x9x10q_ebecs_DDB -nb

Sometimes, however, you do not need the interactive environment provided
by jupyter notebooks because you are mainly interested in visualizing the results.
In this case, you can use the command line interface to generate
a dashboard with widgets without starting a notebook.

To build a dashboard for a |Structure| object extracted from ``FILE``, use:

.. code-block:: bash

    abistruct.py panel FILE

where ``FILE`` is **any** file that provides a ``Structure`` object,
e.g. netcdf files, CIF files, abi and abo files, etc.

To build a dashboard for one of the files supported by AbiPy, use:

.. code-block:: bash

    abiopen.py FILE --panel

where ``FILE`` is any Abinit file supported by ``abiopen.py``.
For instance, to create a dashboard for a ``DDB`` file, use:

.. code-block:: bash

    abiopen.py out_DDB --panel

To build a dashboard for an AbiPy Flow, use:

.. code-block:: bash

        abirun.py FLOWDIR panel

or, equivalently:

.. code-block:: bash

        abiopen.py FLOWDIR/__AbinitFlow__.pickle --panel

Serving dashboards from a remote server
---------------------------------------

In all the examples so far, we assumed that AbiPy and the web browser
run on the same machine.
Calculations, however, are often performed on clusters where a web browser
is not available or where the connection is too slow to use one.
You could copy the files from the cluster to your local machine with scp
or mount the remote file system with sshfs, but neither approach is ideal.
Ideally, we would like to run AbiPy and Abinit on the remote cluster
and visualize the results directly on our local machine.

This section explains how to start a web server on the remote cluster and
connect to it from your local machine.
The procedure is inspired by `this blog post
<https://ljvmiranda921.github.io/notebook/2018/01/31/running-a-jupyter-notebook/>`_.

In what follows, ``localuser`` and ``localhost`` denote the local user and host,
while ``remoteuser`` and ``remotehost`` denote the remote user and host.
Make sure that AbiPy and all its dependencies are installed on ``remotehost``,
including the ``manager.yml`` configuration file.

**Step 1: Start the server on the remote machine**

Log in to the remote machine via ssh as usual with ``ssh remoteuser@remotehost``, then run:

.. code-block:: bash

    abiopen.py FILE --panel --no-browser --port 49412

The ``--no-browser`` option starts the server without opening a browser.
The server listens on port 49412 of ``remotehost``.
If this port is already in use, choose another one, but remember that ports below 1024 are
reserved.

**Step 2: Forward the remote port to your local machine**

The server is now running on port ``REMOTE_PORT`` (49412 in our example) of the remote host.
Next, forward this port to a port ``LOCAL_PORT`` of your local machine so that you can
connect to the server from your browser.
On your local machine, run:

.. code-block:: bash

    localuser@localhost: ssh -N -f -L localhost:LOCAL_PORT:localhost:REMOTE_PORT remoteuser@remotehost

where the options have the following meaning:

- ``-N``: do not execute a remote command (typically used for port forwarding).
- ``-f``: send ssh to the background before executing the command.
- ``-L``: bind ``LOCAL_PORT`` on the local machine to ``REMOTE_PORT`` on the remote machine.
  The argument has the form ``local_socket:remote_socket``.

**Step 3: Open the dashboard in your local browser**

Start the web browser on your local machine and type the following in the address bar::

    localhost:LOCAL_PORT

If everything worked, you should see the dashboard.
At the same time, the terminal on the remote machine should show log messages
as you interact with the page.

**Closing the connections**

To close the connections, stop the server on the remote machine with ``CTRL + C``,
then find the ssh process listening on ``LOCAL_PORT`` on your local machine:

.. code-block:: bash

    localuser@localhost: sudo netstat -lpn | grep :LOCAL_PORT

This command shows the process ID (PID) of the process bound to ``LOCAL_PORT``.
Kill it with:

.. code-block:: bash

    localuser@localhost: kill PID
