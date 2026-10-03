.. _documenting-abipy:

Documenting AbiPy
=================

.. contents::
   :backlinks: top

Organization of documentation
-----------------------------

The AbiPy documentation is written in reStructuredText and built with the Sphinx_ documentation generator.
The sources are in the :file:`~/docs/` directory of the repository.
The main items are:

* index.rst - the top-level document of the AbiPy docs
* api - placeholders used to generate the API documentation automatically
* scripts - documentation for the scripts
* flow_gallery/gallery - generated automatically by sphinx-gallery
* workflows - documentation on Abinit flows/works/tasks and the TaskManager
* README.rst - this file
* conf.py - the Sphinx configuration file
* _static - used by the Sphinx build system

The main entry point is :file:`docs/index.rst`, which pulls in the files
for the user guide, the developer guide and the API reference.
The reStructuredText files for the API of the subpackages, the scripts and the workflows are kept
in :file:`docs/api`, :file:`docs/scripts` and :file:`docs/workflows`, respectively.

To add a file to one of the guides, include its base
name (the ``.rst`` extension is not needed) in the table of contents.
You can also include other documents with an include
directive, such as::

  .. include:: ../../TODO

The Sphinx output can be configured by editing the :file:`conf.py` file in :file:`docs/`.
Before building the documentation, install the Sphinx extensions listed
in the ``docs`` optional-dependency group of :file:`abipy/pyproject.toml` with::

    cd abipy
    python -m pip install --editable ".[docs]"

To build the HTML documentation, type ``make``, which executes::

    sphinx-build -b html -d _build/doctrees . _build/html

Remember to set::

    export READTHEDOCS=1

before running ``make`` to generate the thumbnails for :file:`abipy/examples/flows`.

The documentation is written to :file:`_build/html`.
Use::

	open _build/html/index.html

to open the homepage in your browser.

Run ``make help`` to list all the available make targets.

.. _formatting-abipy-docs:

Formatting
----------

The Sphinx website has plenty of documentation_ on reST markup and
on working with Sphinx in general.
Here are a few additional points to keep in mind:

* Familiarize yourself with the Sphinx directives for `inline markup`_.
  AbiPy's documentation makes heavy use of cross-referencing and other semantic markup.
  Several aliases are defined in :file:`abipy/docs/links.rst` and are automatically
  included in every ``rst`` file via `rst_epilog <https://www.sphinx-doc.org/en/stable/config.html#confval-rst_epilog>`_.

* Mathematical expressions are rendered in HTML with `mathjax <https://www.mathjax.org/>`_.
  For example:

  ``:math:`\sin(x_n^2)``` yields :math:`\sin(x_n^2)`, and::

    .. math::

      \int_{-\infty}^{\infty}\frac{e^{i\phi}}{1+x^2\frac{e^{i\phi}}{1+x^2}}

  yields:

  .. math::

    \int_{-\infty}^{\infty}\frac{e^{i\phi}}{1+x^2\frac{e^{i\phi}}{1+x^2}}

* BibTeX citations are supported via the
  `sphinxcontrib-bibtex extension <https://sphinxcontrib-bibtex.readthedocs.io/en/latest/>`_.
  The BibTeX entries are declared in the :file:`abipy/docs/abiref.bib` file.
  For example::

    See :cite:`Gonze2016` for a brief description of recent developments in ABINIT.

  yields: See :cite:`Gonze2016` for a brief description of recent developments in ABINIT.

  To add a new BibTeX entry to the database, use the :program:`doi2bibtex` tool
  provided by the `betterbib package <https://github.com/nschloe/betterbib>`_::

    doi2bibtex https://doi.org/10.1103/PhysRevB.33.7017 >> abiref.bib

  then change the BibTeX identifier to the name of the first author followed by the publication year.

* Interactive ipython_ sessions can be shown in the documentation with the following directive::

    .. sourcecode:: ipython

      In [69]: lines = plot([1, 2, 3])

  which yields:

  .. sourcecode:: ipython

    In [69]: lines = plot([1, 2, 3])

* Use the *note* and *warning* directives, sparingly, to draw attention to important comments::

    .. note::
       Here is a note

  yields:

  .. note::
     Here is a note

  Similarly:

  .. warning::
     Here is a warning

* Use the *deprecated* directive when appropriate::

    .. deprecated:: 0.98
       This feature is obsolete, use something else.

  yields:

  .. deprecated:: 0.98
     This feature is obsolete, use something else.

* The *versionadded* and *versionchanged* directives have a syntax similar
  to *deprecated*::

    .. versionadded:: 0.2
       The transforms have been completely revamped.

  yields:

  .. versionadded:: 0.2
     The transforms have been completely revamped.

* The autodoc extension handles index entries for the API, but any additional
  index entries must be added explicitly.

.. _documentation: http://www.sphinx-doc.org/en/master/
.. _`inline markup`: http://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html?highlight=inline#inline-markup

Docstrings
----------

In addition to the formatting suggestions above:

* Docstrings follow the
  `Google Python Style Guide <http://google.github.io/styleguide/pyguide.html>`_.
  The `napoleon <https://sphinxcontrib-napoleon.readthedocs.io/en/latest/>`_ extension
  converts Google-style docstrings to reStructuredText before Sphinx parses them.

Dynamically generated figures
-----------------------------

Figures can be generated automatically from scripts and included in the docs
with `sphinx-gallery <https://github.com/sphinx-gallery/sphinx-gallery>`_.
There is no need to save the figure explicitly in the script: this is done
automatically at build time, which also ensures that the included code runs and produces the advertised figure.

Plots specific to the documentation should be added to the :file:`examples/plot/` directory and committed to git.
