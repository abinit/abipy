.. _documenting-abipy:

Documenting AbiPy
==================

.. contents::
   :backlinks: top

Organization of documentation
-----------------------------

The documentation for AbiPy is generated from ReStructured Text using the Sphinx_ documentation generation tool. 
The documentation sources are found in the :file:`~/docs/` directory in the repository. Major items : 

* index.rst - the top level include document for AbiPy docs
* api - placeholders to automatically generate the api documentation
* scripts - documentation for scripts
* workflows - documentation for workflows 
* README.rst - the present file
* conf.py - the sphinx configuration
* _static - used by the sphinx build system
* _templates - used by the sphinx build system

The main entry point is :file:`docs/index.rst`, which pulls in 
the files for the users guide, developers guide, api reference. 
The actual ReStructured Text files for the APIs of the subpackages, for the scripts and for the workflows (resp.) are kept 
in :file:`docs/api`, :file:`docs/scripts` and :file:`docs/workflows` (resp.). 

Additional files can be added to the various guides by including their base
file name (the ``.rst`` extension is not necessary) in the table of contents.
It is also possible to include other documents through the use of an include
statement, such as::

  .. include:: ../../TODO

The output produced by Sphinx can be configured by editing the :file:`conf.py` file located in the :file:`docs/`.
Before building the documentation, you need to install the sphinx extensions listed 
in :file:`abipy/docs/requirements.txt` with::

    pip install -r abipy/docs/requirements.txt

To build the HTML documentation, install sphinx then type ``make html`` that will execute::

    sphinx-build -b html -d _build/doctrees . _build/html

Remember to issue::

    export READTHEDOCS=1

before running ``make`` to activate the generation of the thumbnails in :file:`abipy/examples/flows`.

The documentation is produced in :file:`_build/html`.

You can run ``make help`` to see information on all possible make targets.

Use::

   ./ghp_import.py _build/html/ -n -p

to deploy to gh-pages.


.. _formatting-abipy-docs:

Formatting
----------

The Sphinx website contains plenty of documentation_ concerning ReST markup and
working with Sphinx in general. 
Here are a few additional things to keep in mind:

* Please familiarize yourself with the Sphinx directives for `inline markup`_. 
  Abipy's documentation makes heavy use of cross-referencing and other semantic markup. 
  Several aliases are defined in :file:`abipy/docs/links.rst` and are automatically
  included in each ``rst`` file via `rst_epilog <http://www.sphinx-doc.org/en/stable/config.html#confval-rst_epilog>`_

* Mathematical expressions can be rendered with `mathjax <https://www.mathjax.org/>`_ in html.
  For example:

  ``:math:`\sin(x_n^2)``` yields: :math:`\sin(x_n^2)`, and::

    .. math::

      \int_{-\infty}^{\infty}\frac{e^{i\phi}}{1+x^2\frac{e^{i\phi}}{1+x^2}}

  yields:

  .. math::

    \int_{-\infty}^{\infty}\frac{e^{i\phi}}{1+x^2\frac{e^{i\phi}}{1+x^2}}

* Bibtex citations are supported via the 
  `sphinxcontrib-bibtex extension <https://sphinxcontrib-bibtex.readthedocs.io/en/latest/>`_
  The bibtext entries are declared in the :file:`abipy/docs/refs.bib` file.
  For example::

    See :cite:`Gonze2016` for a brief description of recent developments in ABINIT.

  yelds: See :cite:`Gonze2016` for a brief description of recent developments in ABINIT.

  To add a new bibtex entry to the database, please use the :program:`doi2bibtex` tool
  provided by the `betterbib package <https://github.com/nschloe/betterbib>`_::

    doi2bibtex https://doi.org/10.1103/PhysRevB.33.7017 >> refs.bib

  then change the bibtex identifier (use the name of the first author and the publication year).

* Interactive ipython_ sessions can be illustrated in the documentation using the following directive::

    .. sourcecode:: ipython

      In [69]: lines = plot([1, 2, 3])

  which would yield:

  .. sourcecode:: ipython

    In [69]: lines = plot([1, 2, 3])

* Use the *note* and *warning* directives, sparingly, to draw attention to important comments::

    .. note::
       Here is a note

  yields:

  .. note::
     here is a note

  also:

  .. warning::
     here is a warning

* Use the *deprecated* directive when appropriate::

    .. deprecated:: 0.98
       This feature is obsolete, use something else.

  yields:

  .. deprecated:: 0.98
     This feature is obsolete, use something else.

* Use the *versionadded* and *versionchanged* directives, which have similar
  syntax to the *deprecated* role::

    .. versionadded:: 0.2
       The transforms have been completely revamped.

  .. versionadded:: 0.2
     The transforms have been completely revamped.

* The autodoc extension will handle index entries for the API, but additional
  entries in the index_ need to be explicitly added.

.. _documentation: http://sphinx.pocoo.org/contents.html
.. _`inline markup`: http://sphinx.pocoo.org/markup/inline.html
.. _index: http://sphinx.pocoo.org/markup/para.html#index-generating-markup

Docstrings
----------

In addition to the aforementioned formatting suggestions:

* Docstrings are written following the 
  `Google Python Style Guide <http://google.github.io/styleguide/pyguide.html>`_.
  We use the `napoleon <https://sphinxcontrib-napoleon.readthedocs.io/en/latest/>`_ extension
  to convert Google style docstrings to reStructuredText before Sphinx attempts to parse them.

* Please limit the text width of docstrings to (around) 90 characters.

* Keyword arguments should be described using a definition list.

Dynamically generated figures
-----------------------------

Figures can be automatically generated from scripts and included in the docs.  
It is not necessary to explicitly save the figure in the script, this will be done 
automatically at build time to ensure that the code that is included runs and produces the advertised figure.

Any plots specific to the documentation should be added to the :file:`examples/plot/` directory and committed to git.  

`sphinx-gallery <https://github.com/sphinx-gallery/sphinx-gallery>`_
