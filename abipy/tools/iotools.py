# coding: utf-8
"""IO related utilities."""
from __future__ import annotations

import os
import codecs
import errno
import tempfile
import datetime
import ruamel.yaml as yaml
import pandas as pd

from pathlib import Path
from contextlib import ExitStack
from subprocess import call
from typing import Any
from monty.termcolor import cprint
from monty.string import list_strings
from abipy.tools.typing import PathLike


def make_executable(filepath: PathLike) -> None:
    """Make file executable"""
    mode = os.stat(filepath).st_mode
    mode |= (mode & 0o444) >> 2    # copy R bits to X
    os.chmod(filepath, mode)


def try_files(filepaths: list[PathLike]) -> Path:
    """
    Return the first existent file in filepaths
    or raise RuntimeError.
    """
    for path in filepaths:
        path = Path(str(path))
        if path.exists(): return path

    raise RuntimeError("Cannot find {filepaths=}")


def file_with_ext_indir(ext: str, directory: PathLike) -> Path:
    """
    Find file with extension `ext` inside directory.
    Raise RuntimeError if no file can be found.
    """
    directory = Path(str(directory))
    for path in directory.listdir():
        if path.is_dir(): continue
        if path.suffix == ext:
            return path.absolute()

    raise RuntimeError(f"Cannot find file with extension {ext} in {directory=})")


def yaml_safe_load(string: str) -> Any:
    """Load Yaml string"""
    return yaml.YAML(typ='safe', pure=True).load(string)


def yaml_safe_load_path(filepath: str) -> Any:
    """Load Yaml document from filepath"""
    with open(os.path.expanduser(filepath), "rt") as fh:
        return yaml.YAML(typ='safe', pure=True).load(fh.read())


def dataframe_from_filepath(filepath: str, **kwargs) -> pd.DataFrame:
    """
    Try to read a dataframe from an external file according to the file extension.
    """
    _, ext = os.path.splitext(filepath)
    if ext == 'csv': return pd.read_csv(filepath, **kwargs)
    if ext == "json": return pd.read_json(filepath, **kwargs)
    if ext in ("xls", "xlsx"): return pd.read_excel(filepath, **kwargs)

    raise ValueError(f"Don't know how to construct DataFrame from file {filepath} with extension: {ext}")


class ExitStackWithFiles(ExitStack):
    """
    Context manager for dynamic management of a stack of file-like objects.
    Mainly used in a callee that needs to return files to the caller

    Usage example:

    .. code-block:: python

        exit_stack = ExitStackWithFiles()
        exit_stack.enter_context(phbst_file)
        return exit_stack
    """
    def __init__(self):
        self.files = []
        super().__init__()

    def enter_context(self, myfile):
        # If my file is None, we add it to files but without registering the callback.
        self.files.append(myfile)
        if myfile is not None:
            return super().enter_context(myfile)

    def __iter__(self):
        return self.files.__iter__()

    def __next__(self):
        return self.files.__next__()

    def __getitem__(self, slice):
        return self.files.__getitem__(slice)


def get_input(prompt: str):
    """
    Wraps python builtin input so that we can easily mock it in unit tests using:

    Example:

        from unittest.mock import patch
        with patch('abipy.tools.iotools.get_input', return_value='no'):
            do_something_that_uses_get_input
    """
    return input(prompt)


def ask_yes_no(prompt: str, default=None):  # pragma: no cover
    """
    Ask a question and return a boolean (y/n) answer.

    If default is given (one of 'y','n'), it is used if the user input is
    empty. Otherwise the question is repeated until an answer is given.

    An EOF is treated as the default answer.  If there is no default, an
    exception is raised to prevent infinite loops.

    Valid answers are: y/yes/n/no (match is not case sensitive).
    """
    # Fixes py2.x
    answers = {'y': True, 'n': False, 'yes': True, 'no': False}
    ans = None
    while ans not in answers.keys():
        try:
            ans = get_input(prompt + ' ').lower()
            if not ans:
                # response was an empty string
                ans = default
        except KeyboardInterrupt:
            pass
        except EOFError:
            if default in answers.keys():
                ans = default
                print("")
            else:
                raise

    return answers[ans]


def _user_wants_to_exit(): # pragma: no cover
    try:
        answer = get_input("Do you want to continue [Y/n]")
    except EOFError:
        return True

    if answer.lower().strip() in ["n", "no"]: return True
    return False


class EditorError(Exception):
    """Base class for exceptions raised by `Editor`"""


class Editor:  # pragma: no cover
    DEFAULT_EDITOR = "vi"

    Error = EditorError

    def __init__(self, editor=None):
        if editor is None:
            self.editor = os.getenv("EDITOR", self.DEFAULT_EDITOR)
        else:
            self.editor = str(editor)

    def edit_file(self, filepath):
        retcode = call([self.editor, filepath])
        if retcode != 0:
            cprint("Retcode %s while editing file: %s" % (retcode, filepath), "red")
        return retcode

    def edit_files(self, filepaths, ask_for_exit=True):
        for idx, fname in enumerate(list_strings(filepaths)):
            exit_status = self.edit_file(fname)

            if exit_status != 0:
                return exit_status

            if ask_for_exit and idx != len(filepaths) - 1 and _user_wants_to_exit():
                break

        return 0


def input_from_editor(message=None):  # pragma: no cover
    if message is not None:
        print(message, end="")

    from tempfile import mkstemp
    fd, fname = mkstemp(text=True)

    Editor().edit_file(fname)
    with open(fname, "rt") as fileobj:
        return fileobj.read()


def ask_yesno(question: str, default=True):
    """
    Args:
        question ():
        default ():
    Returns:
    """
    try:
        answer = input(question)
        return answer.lower().strip() in ["y", "yes"]
    except EOFError:
        return default


umask = os.umask(0)
os.umask(umask)

def _maketemp(name: str, createmode=None) -> str:
    """
    Create a temporary file with the filename similar the given ``name``.
    The permission bits are copied from the original file or ``createmode``.
    Returns: the name of the temporary file.
    """
    d, fn = os.path.split(name)
    fd, tempname = tempfile.mkstemp(prefix=f".{fn}-", dir=d)
    os.close(fd)

    # Temporary files are created with mode 0600, which is usually not
    # what we want. If the original file already exists, just copy its mode.
    # Otherwise, manually obey umask.
    try:
        st_mode = os.lstat(name).st_mode & 0o777
    except OSError as err:
        if err.errno != errno.ENOENT:
            raise
        st_mode = createmode
        if st_mode is None:
            st_mode = ~umask
        st_mode &= 0o666
    os.chmod(tempname, st_mode)

    return tempname


class AtomicFile:
    """
    This is a straight port of Alexander Saltanov's atomicfile package.
    Writeable file object that atomically writes a file.
    All writes will go to a temporary file.
    Call ``close()`` when you are done writing, and AtomicFile will rename
    the temporary copy to the original name, making the changes visible.
    If the object is destroyed without being closed, all your writes are
    discarded.
    If an ``encoding`` argument is specified, codecs.open will be called to open
    the file in the wanted encoding.
    """

    def __init__(self, name, mode="w+b", createmode=None, encoding=None):
        """
        Args:
            name ():
            mode ():
            createmode ():
            encoding ():
        """
        self.__name = name  # permanent name
        self._tempname = _maketemp(name, createmode=createmode)
        if encoding:
            self._fp = codecs.open(self._tempname, mode, encoding)  # pylint: disable=R1732
        else:
            self._fp = open(self._tempname, mode)  # pylint: disable=R1732

        # delegated methods
        self.write = self._fp.write
        self.fileno = self._fp.fileno

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, exc_tb):
        if exc_type:
            return
        self.close()

    def close(self) -> None:
        """
        Close the file.
        """
        if not self._fp.closed:
            self._fp.close()
            # This to avoid:
            #   FileExistsError: [WinError 183] Cannot create a file when that file already exists:
            # On Windows, if dst already exists, OSError will be raised even if it is a file;
            # there may be no way to implement an atomic rename when dst names an existing file.
            if os.name == "nt" and os.path.exists(self.__name):
                os.remove(self.__name)
            os.rename(self._tempname, self.__name)

    def discard(self) -> None:
        """
        Discard the file.
        """
        if not self._fp.closed:
            try:
                os.unlink(self._tempname)
            except OSError:
                pass
            self._fp.close()

    def __del__(self):
        if getattr(self, "_fp", None):  # constructor actually did something
            self.discard()


def workdir_with_prefix(workdir, prefix, exist_ok=False) -> Path:
    """
    if workdir is None, create temporary directory with prefix else check that workdir does not exist.
    If exist_ok is False (the default), a FileExistsError is raised if the target directory already exists.
    """
    if workdir is None:
        if prefix is not None:
            prefix += datetime.datetime.now().strftime("dT%H-%M-%S-")
        workdir = tempfile.mkdtemp(prefix=prefix, dir=os.getcwd())
    else:
        workdir = str(workdir)
        if os.path.exists(workdir) and not exist_ok:
            raise RuntimeError(f"{workdir=} already exists!")
        os.makedirs(workdir, exist_ok=exist_ok)

    return Path(workdir).absolute()


def change_ext_from_top(top: PathLike, old_ext: str, new_ext: str) -> int:
    """
    Change the extension of all the files with extension old_ext with new_ext.

    Args:
        top: Walk the file system starting from top.
        old_ext: Old file extension.
        new_ext: New file extension.

    Return: Number of files whose extension has been changed.
    """
    count = 0
    for root, dirs, files in os.walk(str(top)):
        root = Path(root)
        for filepath in files:
            filepath = root / Path(filepath)
            if not filepath.name.endswith(old_ext): continue
            new_name = filepath.name[:-len(old_ext)] + new_ext
            filepath.rename(root / new_name)
            count += 1

    return count


class _Script:
    """
    Base class for Script objects.
    """

    def __init__(self, filepath: str):
        self.filepath = filepath

        self.text = f"""\
#!/usr/bin/env python
# This script has been automatically generated by AbiPy.
from __future__ import annotations

import numpy as np
import pandas as pd

if False:
    import seaborn as sns
    sns.set(context="paper", style='darkgrid', palette='deep',
            font='sans-serif', font_scale=0.8, color_codes=False, rc=None)
"""

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Activated at the end of the with statement."""
        self.write()

    def add_text(self, text: str):
        """Add `text` to the script."""
        self.text += "\n" + text
        return self

    def write(self):
        """
        Write python script and json file with the list of files in the Robot.
        """
        with open(self.filepath, "wt") as fh:
            fh.write(self.text)
        make_executable(self.filepath)
        return self


class PythonScript(_Script):
    """
    Small object used to generate python scripts programmatically.
    Client code can then add additional logic to the script and write it to disk.

    Example:

        with PythonScript("script.py") as script:
            script.add_text("a = 1")
    """

    def __init__(self, filepath: str):
        super().__init__(filepath)

        self.text = f"""\
#!/usr/bin/env python
# This script has been automatically generated by AbiPy.
from __future__ import annotations

import numpy as np
import pandas as pd

if False:
    import seaborn as sns
    sns.set(context="paper", style='darkgrid', palette='deep',
            font='sans-serif', font_scale=0.8, color_codes=False, rc=None)
"""

    def add_main(self):
        """Add main section"""
        self.text += \
"""

if __name__ == "__main__":
    main()
"""
        return self



class ShellScript(_Script):
    """
    Small object used to generate a shell scripts programmatically.
    Client code can then add additional logic to the script and write it to disk.

    Example:

        with ShellScript("script.sh") as script:
            script.add_text("a = 1")
    """

    def __init__(self, filepath: str):
        super().__init__(filepath)

        self.text = f"""\
#!/bin/bash
# This script has been automatically generated by AbiPy.

"""
