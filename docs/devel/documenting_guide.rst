.. _documenting-abipy:

*****************
Documenting abipy
*****************

Getting started
===============

The documentation for abipy is generated from ReStructured Text using the Sphinx_ documentation generation tool. 
Sphinx-1.0 or later is required.

The documentation sources are found in the :file:`docs/` directory in the trunk.  
To build the users guide in html format, cd into :file:`docs/` and do::

  make html 

The output produced by Sphinx can be configured by editing the :file:`conf.py` file located in the :file:`docs/`.


Organization of AbiPy's documentation
=====================================

The actual ReStructured Text files are kept in :file:`docs/users`, :file:`docs/devel`, :file:`docs/api`. 
The main entry point is :file:`docs/index.rst`, which pulls in the :file:`index.rst` 
file for the users guide, developers guide, api reference. 

Additional files can be added to the various guides by including their base
file name (the .rst extension is not necessary) in the table of contents.
It is also possible to include other documents through the use of an include
statement, such as::

  .. include:: ../../TODO


.. _formatting-abipy-docs:

Formatting
==========

The Sphinx website contains plenty of documentation_ concerning ReST markup and
working with Sphinx in general. Here are a few additional things to keep in mind:

* Please familiarize yourself with the Sphinx directives for `inline
  markup`_. Abipy's documentation makes heavy use of cross-referencing and
  other semantic markup. For example, when referring to external files, use the
  ``:file:`` directive.

* Function arguments and keywords should be referred to using the *emphasis*
  role. This will keep abipy's documentation consistant with Python's
  documentation::

    Here is a description of *argument*

  Please do not use the `default role`::

    Please do not describe `argument` like this.

  nor the ``literal`` role::

    Please do not describe ``argument`` like this.

* Sphinx does not support tables with column- or row-spanning cells for
  latex output. Such tables can not be used when documenting matplotlib.

* Mathematical expressions can be rendered as png images in html, and in the
  usual way by latex. For example:

  ``:math:`\sin(x_n^2)``` yields: :math:`\sin(x_n^2)`, and::

    .. math::

      \int_{-\infty}^{\infty}\frac{e^{i\phi}}{1+x^2\frac{e^{i\phi}}{1+x^2}}

  yields:

  .. math::

    \int_{-\infty}^{\infty}\frac{e^{i\phi}}{1+x^2\frac{e^{i\phi}}{1+x^2}}

* Interactive IPython sessions can be illustrated in the documentation using the following directive::

    .. sourcecode:: ipython

      In [69]: lines = plot([1,2,3])

  which would yield:

  .. sourcecode:: ipython

    In [69]: lines = plot([1,2,3])

* Footnotes [#]_ can be added using ``[#]_``, followed later by::

    .. rubric:: Footnotes

    .. [#]

  .. rubric:: Footnotes

  .. [#] For example.

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

.. _Sphinx: http://sphinx.pocoo.org
.. _documentation: http://sphinx.pocoo.org/contents.html
.. _`inline markup`: http://sphinx.pocoo.org/markup/inline.html
.. _index: http://sphinx.pocoo.org/markup/para.html#index-generating-markup

Docstrings
----------

In addition to the aforementioned formatting suggestions:

* Please limit the text width of docstrings to 70 characters.

* Keyword arguments should be described using a definition list.


Figures
=======

Dynamically generated figures
-----------------------------

Figures can be automatically generated from scripts and included in the docs.  
It is not necessary to explicitly save the figure in the script, this will be done 
automatically at build time to ensure that the code that is included runs and produces the advertised figure.

Any plots specific to the documentation should be added to the ``examples/plot/`` directory and committed to git.  
