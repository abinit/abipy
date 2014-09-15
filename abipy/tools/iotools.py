"""IO related utilities."""
from __future__ import print_function, division, unicode_literals

import os
import tempfile

from subprocess import call
from .text import list_strings


def file_read(filename):
    """Read a file and close it.  Return the file source using read()."""
    fobj = open(filename, 'r')
    source = fobj.read()
    fobj.close()
    return source


def file_readlines(filename):
    """Read a file and close it.  Return the file source using readlines()."""
    fobj = open(filename, 'r')
    lines = fobj.readlines()
    fobj.close()
    return lines


def diff_files(fname1, fname2, no_dupspace=True):
    """
    Return the `Differ`-style delta of two files.

    If no_duspace == True, multiple spaces are removed.
    """
    from difflib import ndiff

    lines1 = file_readlines(fname1)
    lines2 = file_readlines(fname2)

    if no_dupspace:
        for idx, l1 in enumerate(lines1): lines1[idx] = " ".join(l1.split()) + "\n"
        for idx, l2 in enumerate(lines2): lines2[idx] = " ".join(l2.split()) + "\n"
        # Could use re.sub("\s+", " ", s) maybe faster but leading space is not removed

    return ndiff(lines1, lines2) # diff generator.


def write_temp_file(src, suffix='', prefix='tmp', dir=None, text=True):
    """
    Make a temporary file, return filename and filehandle.

    :arg src: string or list of strings (no need for ending newlines if list)
     Source code to be written to the file.

    If 'suffix' is specified, the file name will end with that suffix,
    otherwise there will be no suffix.

    If 'prefix' is specified, the file name will begin with that prefix,
    otherwise a default prefix is used.

    If 'dir' is specified, the file will be created in that directory,
    otherwise a default directory is used.

    If 'text' is specified and False, the file is opened in binary mode.
    Else (the default) the file is opened in text mode.
    On some operating systems, this makes no difference.

    The file is readable and writable only by the creating user ID.
    If the operating system uses permission bits to indicate whether a
    file is executable, the file is executable by no one.
    The file descriptor is not inherited by children of this process

    Returns:
        (filename, open filehandle).
        It is the caller's responsibility to close the open file and unlink it.
    """
    fname = tempfile.mkstemp(suffix, prefix, dir, text)[1]
    f = open(fname, 'w')
    f.write(src)
    f.flush()
    return fname, f


def raw_input_ext(prompt='', ps2='... '):
    """Similar to raw_input(), but accepts extended lines if input ends with \\."""

    line = raw_input(prompt)
    while line.endswith('\\'):
        line = line[:-1] + raw_input(ps2)
    return line


def ask_yes_no(prompt, default=None):
    """
    Ask a question and return a boolean (y/n) answer.

    If default is given (one of 'y','n'), it is used if the user input is
    empty. Otherwise the question is repeated until an answer is given.

    An EOF is treated as the default answer.  If there is no default, an
    exception is raised to prevent infinite loops.

    Valid answers are: y/yes/n/no (match is not case sensitive).
    """

    answers = {'y': True, 'n': False, 'yes': True, 'no': False}
    ans = None
    while ans not in answers.keys():
        try:
            ans = raw_input(prompt + ' ').lower()
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


def stream_has_colours(stream):
    """Python cookbook, #475186"""
    if not hasattr(stream, "isatty"):
        return False

    if not stream.isatty():
        # auto color only on TTYs
        return False

    try:
        import curses
        curses.setupterm()
        return curses.tigetnum("colors") > 2

    except:
        # guess false in case of error.
        return False


class StringColorizer(object):
    COLORS = {"default": "",
               "blue": "\x1b[01;34m",
               "cyan": "\x1b[01;36m",
               "green": "\x1b[01;32m",
               "red": "\x1b[01;31m",
               # lighting colors.
               #"lred":    "\x1b[01;05;37;41m"
    }

    def __init__(self, stream):
        self.has_colours = stream_has_colours(stream)

    def __call__(self, string, colour):
        if self.has_colours:
            code = self.COLORS.get(colour, "")
            if code:
                return code + string + "\x1b[00m"
            else:
                return string
        else:
            return string


def _user_wants_to_exit():
    try:
        answer = raw_input("Do you want to continue [Y/n]")
    except EOFError:
        return True

    if answer.lower().strip() in ["n", "no"]: return True
    return False


class EditorError(Exception):
    """Base class for exceptions raised by `Editor`"""


class Editor(object):
    DEFAULT_EDITOR = "vi"

    Error = EditorError

    def __init__(self, editor=None):
        if editor is None:
            self.editor = os.getenv("EDITOR", self.DEFAULT_EDITOR)
        else:
            self.editor = str(editor)

    def edit_file(self, fname):
        retcode = call([self.editor, fname])
        #if retcode != 0:
        #    raise self.Error("Received retcode %s while editing file: %s" % (retcode, fname))
        return retcode

    def edit_files(self, fnames, ask_for_exit=True):
        for (idx, fname) in enumerate(list_strings(fnames)):
            exit_status = self.edit_file(fname)

            if exit_status != 0:
                return exit_status

            if ask_for_exit and idx != len(fnames) - 1 and _user_wants_to_exit():
                break

        return 0


def input_from_editor(message=None):
    if message is not None:
        print(message, end="")

    from tempfile import mkstemp
    fd, fname = mkstemp(text=True)

    Editor().edit_file(fname)

    with open(fname) as fileobj:
        return fileobj.read()

