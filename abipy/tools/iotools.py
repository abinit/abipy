# coding: utf-8
"""IO related utilities."""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import tempfile

from subprocess import call
from monty.termcolor import cprint
from monty.string import list_strings


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
            cprint("Retcode %s while editing file: %s" % (retcode, fname) ,"red")
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
