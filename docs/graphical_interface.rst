.. _graphical-interface:

Graphical interface
===================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

AbiPy provides interactive dashboards that can be used either as a standalone web applications
with the `bokeh server <http://docs.bokeh.org/>`_ or inside jupyter notebooks.
This document explains how to install the required dependencies and how to
generate dashboards either from the command line interface (CLI) or inside jupyter notebooks.

.. important::

    Please note that one needs a **running python backend**
    to execute the callbacks triggerered by the widgets in the HTML page.
    This part, indeed, is implemented by HTML/CSS/JS code executed
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
an object that can be served by a web browser or displayed inside the jupyter notebook.
When running inside a jupyter notebook, remember to enable the integration 
with the ``panel`` infrastructure by executing:

.. jupyter-execute::

    from abipy import abilab
    abilab.abipanel();

**before calling** any AbiPy ``get_panel`` method.

.. note::

    The ``abipanel`` function is needed to load extensions and javascript packages
    required by AbiPy.
    This function is just a small wrapper around the panel API:

    .. code-block:: bash

        import panel as pn
        pn.extension()


At this point, we can start to construct AbiPy objects.
For our first example, we use the `abiopen` function to open a `GSR` file,
then we call `get_panel` to build a set of widgets that allows us to interact
with the `GsrFile`:

.. jupyter-execute::

    from abipy import abilab
    import abipy.data as abidata

    filename = abidata.ref_file("si_nscf_GSR.nc")
    gsr = abilab.abiopen(filename)

    gsr.get_panel()

The **summary** tab provides a string representation of the file
but there is no widget to interact with it.
If you select the **e-Bands** tab, you will see several widgets and a button
that activates the visualization of the KS band energies.
Again, in this HTML page there is no python server running in background so
if you click the **Plot e-bands** button nothing happens (this is not a bug!).

The advantage of this notebook-based approach is that it is possible to mix
the panel GUIs with python code that can be used to perform
more advanced tasks not supported by the GUI.

Obviously it is possible to have multiple panels running in the same notebook.
Calling ``get_panel`` with an AbiPy structure, for instance, creates a set of widgets
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
a dashboard with widgets without having to start a notebook.

To build a dashboard for a ``Structure`` object extracted from ``FILE``, use:

.. code-block:: bash

    abistruct.py panel FILE

where ``FILE`` is **any** file providing a ``Structure`` object
e.g. netcdf files, cif files, abi, abo files etc.

To build a dashboard associated to one of the AbiPy file, use the syntax:

.. code-block:: bash

    abiopen.py FILE --panel

where ``FILE`` is one of the Abinit files supported by ``abiopen.py``.
For instance, one can create a dashboard to interact with a ``DDB`` file with:

.. code-block:: bash

    abiopen.py out_DDB --panel

To build a dashboard for an AbiPy Flow use:

.. code-block:: bash

        abirun.py FLOWDIR panel

or alternatively:

.. code-block:: bash

        abiopen.py FLOWDIR/__AbinitFlow__.pickle --panel

Serving dashboards from a remote server
---------------------------------------

inspired to https://ljvmiranda921.github.io/notebook/2018/01/31/running-a-jupyter-notebook/

In all the examples presented so far we assumed that both AbiPy and the web browser 
are running on the same machine.
In many cases, however, calculations are performed on clusters in which web browsers are
not always available or, even if the browser is installed, the connection may be too slow..
Obviously, one can always copy files from the remote cluster to the local machine with scp
or mount the remote file system with sshfs but both approaches is far from optimal.
In principle, we would like to be able to execute AbiPy and Abinit on the remote cluster 
and visualize the results directly in our local machine.

In this section, we discuss how to start a web server on the remote cluster and how 
to connect to it from our local machine.

Before starting, let us introduce some notation. 
Let us define the local user and host as localuser and localhost, respectively. 
Similarly, let us define the remote user and remote host as remoteuser and remotehost. 
Needless to say, make sure that Abipy and all its dependencies are installed on the remotehost,
including the `taskmanager.yml` configuration file.

Step 1: 

Run Jupyter Notebook from remote machine
Log-in to your remote machine as usual via ssh.
Once the console shows, type the following:

.. code-block:: bash

    abiopen.py FILE --panel --no-browser --port 49412

jupyter notebook: simply fires up your notebook

--no-browser starts the notebook without opening the browser on the remote host
--port=XXXX sets the port for starting the web server.
When it’s occupied, it finds the next available port.

Step 2: 

Forward port XXXX to YYYY and listen to it
In your remote, the notebook is now running at the port XXXX that you specified. 
What you’ll do next is forward this to port YYYY of your machine so that you can listen 
and run it from your browser. 
To achieve this, we use the following command:

localuser@localhost: ssh -N -f -L localhost:YYYY:localhost:XXXX remoteuser@remotehost

ssh: your handy ssh command. See man page for more info

-N: suppresses the execution of a remote command. Pretty much used in port forwarding.
-f: this requests the ssh command to go to background before execution.
-L: this argument requires an input in the form of local_socket:remote_socket. 
Here, we’re specifying our port as YYYY which will be binded to the port XXXX from your remote connection.

Step 3: Fire-up Jupyter Notebook
To open up the Jupyter notebook from your remote machine, simply start your web browser on your local machine 
and type the following in your address bar:

localhost:YYYY

Again, the reason why we’re opening it at YYYY and not at XXXX is because the latter 
is already being forwarded to the former. 
XXXX and YYYY can be the “same” number (not the same port, technically) because they are from different machines.

If you’re successful, you should see the typical Jupyter Notebook home screen in the directory 
where you ran the command in Step 1. At the same time, if you look in your remote terminal, 
you should see some log actions happening as you perform some tasks.

In your first connection, you may be prompted to enter an Access Token as typical to most Jupyter notebooks. 
Normally, I’d just copy-paste it from my terminal, but to make things easier for you, 
you can set-up your own notebook password.

Closing all connections
To close connections, I usually stop my notebook from remote via CTRL + C then Y, and kill the process on YYYY via:

localuser@localhost: sudo netstat -lpn |grep :YYYY

# This will show the process ID (PID), e.g. ABCDEF of the one running in YYYY,
# you can kill it by simply typing

localuser@localhost: kill ABCDEF
