from __future__ import print_function, division

import os

from subprocess import call, Popen

#def apathin(apath, basename):
#    if os.path.isdir(apath)
#        return os.path.join(apath, basename)
#    else:
#        # file
#        return os.path.join(os.path.dirname(apath), basename)

def test_scripts():
    script = os.path.join(os.path.dirname(__file__), "_run_all.py")
    retcode = call(script)
    assert retcode == 0
