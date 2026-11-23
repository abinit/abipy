.. _coding-guide:

Coding guide
============

.. contents::
   :backlinks: top

Committing changes
------------------

When committing changes to AbiPy, keep the following points in mind:

* If your changes are non-trivial, add an entry to :file:`CHANGELOG.rst`.
  The changelog follows the format used by the
  `releases <https://github.com/bitprophet/releases>`_ Sphinx extension.

* If you change the API, document the modifications in the docstring with the ``versionadded`` role::

    .. versionadded:: 0.2
       Add new argument ``foobar``

* Do your changes pass the automatic tests?

* Can you add a test for your changes?

* If you have added new files or directories, or reorganized existing
  ones, are the new files matched by the patterns in :file:`MANIFEST.in`?
  This file determines what goes into the source distribution.

Importing and name spaces
-------------------------

For numpy_, use::

  import numpy as np
  a = np.array([1,2,3])

For matplotlib_, **avoid** the high-level interface, as in::

  import matplotlib.pyplot as plt
  plt.plot(x, y)

and use the object-oriented API provided by |matplotlib-Axes| instead::

    from abipy.tools.plotting import get_ax_fig_plt
    ax, fig, plt = get_ax_fig_plt(ax=None)
    ax.plot(x, y)

Plotting methods should accept an Axes ``ax`` argument and
use the ``add_fig_kwargs`` decorator::

    @add_fig_kwargs
    def plot(self, ax=None, **kwargs):
        """
        Plot the object ...

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.

        Returns: |matplotlib-Figure|
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        ax.plot(self.xvals, self.yvalsm, **kwargs)
        return fig

Naming, spacing, and formatting conventions
-------------------------------------------

In general, we try to follow as closely as possible the standard
Python coding guidelines written by Guido van Rossum in `PEP0008 <https://www.python.org/dev/peps/pep-0008>`_.

* functions and class methods: ``lower`` or ``lower_underscore_separated``
* attributes and variables: ``lower``
* classes: ``Upper`` or ``MixedCase``

Prefer the shortest names that are still readable.

Configure your editor to use spaces, not hard tabs.
The standard indentation unit is always four spaces;
a file with tabs or a different number of spaces is a bug, so please fix it.

Keep docstrings uniformly indented as in the example above, with nothing to the left of the triple quotes.

Limit line length to around 90 characters.
If you wonder why we deviate from PEP8, which specifies a maximum line length of 79 characters,
check out this video by Raymond Hettinger:

.. youtube:: wf-BqAjZb8M

If a logical line needs to be longer, break it using parentheses rather than an escaped newline.
Sometimes it is better to introduce a temporary variable and replace a single
long line with two shorter, more readable ones.

Please do not commit lines with trailing whitespace, as they add noise to diffs.

Writing examples
----------------

The examples live in subdirectories of :file:`abipy/examples` and are automatically
run when the website is built, so that they appear in both the :file:`examples`
and :file:`gallery` sections of the website.

Many people find these examples on the website and do not have easy access to the
:file:`examples` directory itself.
Any data required by an example should therefore be added to the :file:`abipy/data` directory.

Testing
-------

AbiPy has a testing infrastructure based on :mod:`unittest` and pytest_.

Common test support is provided by :mod:`abipy.core.testing`.
Data files are stored in :file:`abipy/data`; in particular, :file:`abipy/data/refs`
contains several output files that can be used to write unit tests and examples.

To install pytest together with useful plugins, run::

    python -m pip install --editable ".[tests]"

in the top-level directory of the package.
To run the tests of the abio.inputs module, use::

    pytest -v abio/tests/test_inputs.py

while::

    pytest -v

runs the entire test suite.
