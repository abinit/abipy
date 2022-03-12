# coding: utf-8
import os
import tempfile
import warnings

import logging
logger = logging.getLogger(__file__)


def profile(statement, global_vars, local_vars):
    """
    Run statement under profiler, supplying your own globals and locals

    Example::

        stats = profile("main()", global_vars=globals(), local_vars=locals())
    """
    import pstats
    import cProfile
    _, filename = tempfile.mkstemp()
    cProfile.runctx(statement, global_vars, local_vars, filename=filename)

    s = pstats.Stats(filename)
    s.strip_dirs().sort_stats("time")
    s.print_stats()
    os.remove(filename)
    return s


class HtmlDiff(object):
    """
    This object produces diff files in HTML format and displays them in the browser.

    Usage example:

    .. code-block:: python

        HtmlDiff(filepaths).open_browser()
    """
    def __init__(self, filepaths):
        if len(filepaths) < 2:
            raise ValueError("You need more than one file to compare!")
        self.filepaths = filepaths

    def open_browser(self, diffmode="difflib", **kwargs):
        """
        Generate diff with ``diffmode``, open browser, return exit code.
        """
        try:
            func = getattr(self, diffmode)
        except AttributeError:
            raise AttributeError("Unsupported value for diffmode: `%s`" % str(diffmode))

        return func(**kwargs)

    def _launch_browser(self, tmpname):  # pragma: no cover
        """Open ``tmpname`` file in the default browser."""
        # warning: This code is not portable since we should pass a url.
        if not tmpname.startswith("file://"): tmpname = "file://" + tmpname
        import webbrowser
        try:
            return int(webbrowser.open(tmpname))
        except webbrowser.Error as exc:
            # Warn the user and ignore the exception.
            warnings.warn(str(exc))
            return 1

    def difflib(self, **kwargs):
        """
        Use difflib to generate a HTML file with the diff.
        Open file in the browser.
        """
        with open(self.filepaths[0], 'rt') as fh:
            fromlines = fh.readlines()

        diffs = []
        for path in self.filepaths[1:]:
            with open(path, 'rt') as fh:
                tolines = fh.readlines()

            _, tmpname = tempfile.mkstemp(suffix=".html", text=True)
            import difflib
            #diff = difflib.HtmlDiff().make_table(fromlines, tolines,
            diff = difflib.HtmlDiff().make_file(fromlines, tolines,
                                                self.filepaths[0], path) #context=options.c, numlines=n)
            with open(tmpname, "wt") as fh:
                fh.writelines(diff)

            return self._launch_browser(tmpname)

    def pygmentize(self):
        """
        Execute ``diff`` and ``pygmentize`` in a subprocess to generate a HTML file with the diff.
        Open file in the browser.
        """
        for file2 in self.filepaths[1:]:
            _, tmpname = tempfile.mkstemp(suffix=".html", text=True)

            # https://stackoverflow.com/questions/641055/diff-to-html-diff2html-program
            #cmd = "/usr/bin/diff -y %s %s | pygmentize -l diff -f html -O full -o %s" % (
            #    self.filepaths[0], file2, tmpname)
            cmd = "diff -U9999999 %s %s | pygmentize -l diff -f html -O full -o %s" % (
                self.filepaths[0], file2, tmpname)

            retcode = os.system(cmd)
            if retcode != 0:
                print("Non-zero exit status returned by: `%s", cmd)
                print("Please install pygmentize with `pip install pygmentize` or `conda install pygmentize`")
                return retcode

            return self._launch_browser(tmpname)


def display_top(snapshot, key_type='lineno', limit=3):
    """
    Profile memory usage in Python.
    Taken from https://stackoverflow.com/questions/552744/how-do-i-profile-memory-usage-in-python

    Example::

        tracemalloc.start()
        main()
        snapshot = tracemalloc.take_snapshot()
        display_top(snapshot)
    """
    import tracemalloc
    import linecache
    snapshot = snapshot.filter_traces((
        tracemalloc.Filter(False, "<frozen importlib._bootstrap>"),
        tracemalloc.Filter(False, "<unknown>"),
    ))
    top_stats = snapshot.statistics(key_type)

    print("Top %s lines" % limit)
    for index, stat in enumerate(top_stats[:limit], 1):
        frame = stat.traceback[0]
        # replace "/path/to/module/file.py" with "module/file.py"
        filename = os.sep.join(frame.filename.split(os.sep)[-2:])
        print("#%s: %s:%s: %.1f KiB"
              % (index, filename, frame.lineno, stat.size / 1024))
        line = linecache.getline(frame.filename, frame.lineno).strip()
        if line:
            print('    %s' % line)

    other = top_stats[limit:]
    if other:
        size = sum(stat.size for stat in other)
        print("%s other: %.1f KiB" % (len(other), size / 1024))
    total = sum(stat.size for stat in top_stats)
    print("Total allocated size: %.1f KiB" % (total / 1024))
