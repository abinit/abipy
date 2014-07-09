from __future__ import print_function, division

from pymatgen.util.io_utils import FileLock

import os
import re
import subprocess


def number_of_cpus():
    """
    Number of virtual or physical CPUs on this system, i.e.
    user/real as output by time(1) when called with an optimally scaling userspace-only program
    Return -1 if ncpus cannot be detected
    taken from:
    http://stackoverflow.com/questions/1006289/how-to-find-out-the-number-of-cpus-in-python
    """

    # Python 2.6+
    try:
        import multiprocessing
        return multiprocessing.cpu_count()
    except (ImportError, NotImplementedError):
        pass

    # POSIX
    try:
        res = int(os.sysconf('SC_NPROCESSORS_ONLN'))
        if res > 0: return res
    except (AttributeError, ValueError):
        pass

    # Windows
    try:
        res = int(os.environ['NUMBER_OF_PROCESSORS'])
        if res > 0: return res
    except (KeyError, ValueError):
        pass

    # jython
    try:
        from java.lang import Runtime
        runtime = Runtime.getRuntime()
        res = runtime.availableProcessors()
        if res > 0: return res
    except ImportError:
        pass

    # BSD
    try:
        sysctl = subprocess.Popen(['sysctl', '-n', 'hw.ncpu'], stdout=subprocess.PIPE)
        scStdout = sysctl.communicate()[0]
        res = int(scStdout)
        if res > 0: return res
    except (OSError, ValueError):
        pass

    # Linux
    try:
        res = open('/proc/cpuinfo').read().count('processor\t:')
        if res > 0: return res
    except IOError:
        pass

    # Solaris
    try:
        pseudoDevices = os.listdir('/devices/pseudo/')
        expr = re.compile('^cpuid@[0-9]+$')
        res = 0
        for pd in pseudoDevices:
            if expr.match(pd) is not None:
                res += 1
        if res > 0: return res
    except OSError:
        pass

    # Other UNIXes (heuristic)
    try:
        try:
            dmesg = open('/var/run/dmesg.boot').read()
        except IOError:
            dmesgProcess = subprocess.Popen(['dmesg'], stdout=subprocess.PIPE)
            dmesg = dmesgProcess.communicate()[0]

        res = 0
        while '\ncpu' + str(res) + ':' in dmesg:
            res += 1

        if res > 0: return res
    except OSError:
        pass

    return -1
    #raise Exception('Cannot determine number of CPUs on this system')

import json


class NumPyArangeEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist() # or map(int, obj)

        return json.JSONEncoder.default(self, obj)


def profile(statement, global_vars=None, local_vars=None):
    """
    Run statement under profiler, supplying your own globals and locals,
    """
    import pstats
    import cProfile
    import tempfile
    global_vars = globals() if global_vars is None else global_vars
    local_vars = locals() if local_vars is None else local_vars

    _, filename = tempfile.mkstemp()
    cProfile.runctx(statement, global_vars, local_vars, filename=filename)

    s = pstats.Stats(filename)
    s.strip_dirs().sort_stats("time").print_stats()

    os.remove(filename)
