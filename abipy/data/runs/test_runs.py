from __future__ import print_function, division

from subprocess import call, Popen

def test_scripts():
   retcode = call("_run_all.py")
   assert retcode == 0
