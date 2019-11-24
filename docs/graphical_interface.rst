.. _graphical-interface:

Graphical interface
===================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

AbiPy provides interactive dashboards that can be used
either as a standalone web apps or inside jupyter notebooks.
The implementation is based on `panel <https://panel.pyviz.org/>`_

This document explains how to install the dependencies and 
use the AbiPy dashboard either with the command line interface or inside jupyter notebooks.


Installation
------------

Get panel from pip:

.. code-block:: bash

  pip install panel

or conda:

.. code-block:: bash

  conda install panel -c conda-forge

If you want to work with JupyterLab, you will also need to install the optional PyViz JupyterLab extension:

.. code-block:: bash

    conda install -c conda-forge jupyterlab
    jupyter labextension install @pyviz/jupyterlab_pyviz

Basic Usage
-----------

You can use the ``jupyter-execute`` directive to embed code into the document::

  .. jupyter-execute::

    name = 'world'
    print('hello ' + name + '!')

The above is rendered as follows:

.. jupyter-execute::

  name = 'world'
  print('hello ' + name + '!')

Note that the code produces *output* (printing the string 'hello world!'), and the output
is rendered directly after the code snippet.

Because all code cells in a document are run in the same kernel, cells later in the document
can use variables and functions defined in cells earlier in the document:

.. jupyter-execute::

    a = 1
    print('first cell: a = {}'.format(a))

.. jupyter-execute::

    a += 1
    print('second cell: a = {}'.format(a))

Because jupyter-sphinx uses the machinery of ``nbconvert``, it is capable of rendering
any rich output, for example plots:

.. jupyter-execute::

    import numpy as np
    from matplotlib import pyplot
    %matplotlib inline

    x = np.linspace(1E-3, 2 * np.pi)

    pyplot.plot(x, np.sin(x) / x)
    pyplot.plot(x, np.cos(x))
    pyplot.grid()

LaTeX output:

.. jupyter-execute::

  from IPython.display import Latex
  Latex('∫_{-∞}^∞ e^{-x²}dx = \sqrt{π}')

or even full-blown javascript widgets:

.. jupyter-execute::

    import ipywidgets as w
    from IPython.display import display

    a = w.IntSlider()
    b = w.IntText()
    w.jslink((a, 'value'), (b, 'value'))
    display(a, b)

It is also possible to include code from a regular file by passing the filename as argument
to ``jupyter-execute``::

  .. jupyter-execute:: some_code.py

``jupyter-execute`` may also be used in docstrings within your Python code, and will be executed
when they are included with Sphinx autodoc.

.. jupyter-execute::
    
    from abipy import abilab
    import abipy.data as abidata
    import panel as pn
    pn.extension()

.. jupyter-execute::

    filename = abidata.ref_file("si_nscf_GSR.nc")
    abifile = abilab.abiopen(filename)
    abifile.structure.get_panel()

.. jupyter-execute::

    abifile.get_panel()


.. jupyter-execute::

    import nglview as nv
    view = nv.show_pymatgen(abifile.structure)
    view.add_unitcell()
    view

.. jupyter-execute::

    view = nv.demo(gui=True)
