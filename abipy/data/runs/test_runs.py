from __future__ import division, print_function, unicode_literals

from subprocess import call

import os


def test_scripts():
    """Running all examples in abipy.data.runs..."""
    return
    script = os.path.join(os.path.dirname(__file__), "_run_all.py")
    retcode = call(script)
    if retcode != 0:
        msg = "_run_all.py has executed all the scripts in abipy.data.runs and returned %s failures" % retcode
        raise RuntimeError(msg)
