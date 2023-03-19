# coding: utf-8
"""IO related utilities."""
from __future__ import annotations

import os
import codecs
import errno
import tempfile
import ruamel.yaml as yaml

from contextlib import ExitStack
from subprocess import call
from typing import Any
from monty.termcolor import cprint
from monty.string import list_strings


def yaml_safe_load(string: str) -> Any:
    """Load Yaml string"""
    return yaml.YAML(typ='safe', pure=True).load(string)


def yaml_safe_load_path(filepath: str) -> Any:
    """Load Yaml document from filepath"""
    with open(filepath, "rt") as fh:
        return yaml.YAML(typ='safe', pure=True).load(fh.read())


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


def get_input(prompt):
    """
    Wraps python builtin input so that we can easily mock it in unit tests using:

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


class Editor(object):  # pragma: no cover
    DEFAULT_EDITOR = "vi"

    Error = EditorError

    def __init__(self, editor=None):
        if editor is None:
            self.editor = os.getenv("EDITOR", self.DEFAULT_EDITOR)
        else:
            self.editor = str(editor)

    def edit_file(self, fname):
        retcode = call([self.editor, fname])
        if retcode != 0:
            cprint("Retcode %s while editing file: %s" % (retcode, fname), "red")
        return retcode

    def edit_files(self, fnames, ask_for_exit=True):
        for idx, fname in enumerate(list_strings(fnames)):
            exit_status = self.edit_file(fname)

            if exit_status != 0:
                return exit_status

            if ask_for_exit and idx != len(fnames) - 1 and _user_wants_to_exit():
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


def ask_yesno(question, default=True):
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

def _maketemp(name, createmode=None):
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

    def close(self):
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

    def discard(self):
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
