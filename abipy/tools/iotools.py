# coding: utf-8
"""IO related utilities."""
import os

from contextlib import ExitStack
from subprocess import call
from monty.termcolor import cprint
from monty.string import list_strings


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


def ask_yes_no(prompt, default=None):  # pragma: no cover
    """
    Ask a question and return a boolean (y/n) answer.

    If default is given (one of 'y','n'), it is used if the user input is
    empty. Otherwise the question is repeated until an answer is given.

    An EOF is treated as the default answer.  If there is no default, an
    exception is raised to prevent infinite loops.

    Valid answers are: y/yes/n/no (match is not case sensitive).
    """
    # Fixes py2.x
    try:
        my_input = raw_input
    except NameError:
        my_input = input

    answers = {'y': True, 'n': False, 'yes': True, 'no': False}
    ans = None
    while ans not in answers.keys():
        try:
            ans = my_input(prompt + ' ').lower()
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
    # Fixes py2.x
    try:
        my_input = raw_input
    except NameError:
        my_input = input

    try:
        answer = my_input("Do you want to continue [Y/n]")
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
    with open(fname) as fileobj:
        return fileobj.read()
