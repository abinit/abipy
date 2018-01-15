"""
Global variables used to initialize AbiPy environment in notebooks.
"""
from __future__ import print_function, division, unicode_literals

from monty.termcolor import cprint

import os
import tempfile


__IN_NOTEBOOK = False

def in_notebook():
    """True if we are running inside a jupyter notebook (and enable_notebook has been called)."""
    return __IN_NOTEBOOK


def disable_notebook():
    """Set ``in_notebook`` flag to False."""
    global __IN_NOTEBOOK
    __IN_NOTEBOOK = False


def enable_notebook(with_seaborn=True):
    """
    Set ``in_notebook`` flag to True and activate seaborn settings for notebooks if ``with_seaborn``.
    """
    global __IN_NOTEBOOK
    __IN_NOTEBOOK = True

    # Use seaborn settings for plots (optional)
    if with_seaborn:
        import seaborn as sns
        sns.set(context='notebook', style='darkgrid', palette='deep',
                font='sans-serif', font_scale=1, color_codes=False, rc=None)


def get_abinb_workdir():
    """
    Return the absolute path of the scratch directory used to produce
    and save temporary files when we are runnning inside a jupyter_ notebook.

    .. note:

        Due to web-browser policy, files used in the notebook must be within the current working directory.
    """
    wdir = os.path.join(os.getcwd(), "__abinb_workdir__")
    if not os.path.exists(wdir): os.mkdir(wdir)
    return wdir


def abinb_mkstemp(force_abinb_workdir=False, use_relpath=False, **kwargs):
    """
    Invoke mkstep with kwargs, return the (fd, name) of the temporary file.
    kwargs are passed to ``mkstemp`` except for ``dir`` if we are inside a jupyter notebook.

    Args:
        use_abipy_nbworkdir:
        use_relpath: Return relative path (os.path.relpath) if True else absolute (default)
            Relative paths are required if we are gonna use the temporary file in
            notebooks or in web browers.
            In this case, the caller is responsbile for calling the function with the correct flag.

    .. example:

        _, filename = abinb_mkstep(suffix="." + ext, text=True)
    """
    if in_notebook() or force_abinb_workdir:
        d = kwargs.pop("dir", None)
        if d is not None:
            cprint("Files should be created inside abipy_nbworkdir if we are inside jupyter or force_abinb_workdir",
                    "yellow")
        fd, path = tempfile.mkstemp(dir=get_abinb_workdir(), **kwargs)
    else:
        fd, path = tempfile.mkstemp(**kwargs)

    if use_relpath:
        path = os.path.relpath(path)

    return fd, path
