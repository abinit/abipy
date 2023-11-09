# coding: utf-8
from __future__ import annotations

import os
import tempfile
import warnings


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


class HtmlDiff:
    """
    This object produces diff files in HTML format and displays them in the browser.

    Usage example:

    .. code-block:: python

        HtmlDiff(filepaths).open_browser()
    """
    def __init__(self, filepaths: list[str]):
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



def get_size(bytes, suffix="B"):
    """
    Scale bytes to its proper format
    e.g:
        1253656 => '1.20MB'
        1253656678 => '1.17GB'
    """
    factor = 1024
    for unit in ["", "K", "M", "G", "T", "P"]:
        if bytes < factor:
            return f"{bytes:.2f}{unit}{suffix}"
        bytes /= factor


def print_hardware_system_info() -> None:
    """
    Taken from <https://thepythoncode.com/article/get-hardware-system-information-python>
    """
    import platform
    uname = platform.uname()
    print("="*40, "System Information", "="*40)
    print(f"System: {uname.system}")
    print(f"Node Name: {uname.node}")
    print(f"Release: {uname.release}")
    print(f"Version: {uname.version}")
    print(f"Machine: {uname.machine}")
    print(f"Processor: {uname.processor}")

    import psutil
    # let's print CPU information
    print("="*40, "CPU Info", "="*40)
    # number of cores
    print("Physical cores:", psutil.cpu_count(logical=False))
    print("Total cores:", psutil.cpu_count(logical=True))
    # CPU frequencies
    cpufreq = psutil.cpu_freq()
    print(f"Max Frequency: {cpufreq.max:.2f}Mhz")
    print(f"Min Frequency: {cpufreq.min:.2f}Mhz")
    print(f"Current Frequency: {cpufreq.current:.2f}Mhz")
    # CPU usage
    print("CPU Usage Per Core:")
    for i, percentage in enumerate(psutil.cpu_percent(percpu=True, interval=1)):
        print(f"Core {i}: {percentage}%")
    print(f"Total CPU Usage: {psutil.cpu_percent()}%")

    # Memory Information
    print("="*40, "Memory Information", "="*40)
    # get the memory details
    svmem = psutil.virtual_memory()
    print(f"Total: {get_size(svmem.total)}")
    print(f"Available: {get_size(svmem.available)}")
    print(f"Used: {get_size(svmem.used)}")
    print(f"Percentage: {svmem.percent}%")
    print("="*20, "SWAP", "="*20)
    # get the swap memory details (if exists)
    swap = psutil.swap_memory()
    print(f"Total: {get_size(swap.total)}")
    print(f"Free: {get_size(swap.free)}")
    print(f"Used: {get_size(swap.used)}")
    print(f"Percentage: {swap.percent}%")

    # Disk Information
    print("="*40, "Disk Information", "="*40)
    print("Partitions and Usage:")
    # get all disk partitions
    partitions = psutil.disk_partitions()
    for partition in partitions:
        print(f"=== Device: {partition.device} ===")
        print(f"  Mountpoint: {partition.mountpoint}")
        print(f"  File system type: {partition.fstype}")
        try:
            partition_usage = psutil.disk_usage(partition.mountpoint)
        except PermissionError:
            # this can be catched due to the disk that
            # isn't ready
            continue
        print(f"  Total Size: {get_size(partition_usage.total)}")
        print(f"  Used: {get_size(partition_usage.used)}")
        print(f"  Free: {get_size(partition_usage.free)}")
        print(f"  Percentage: {partition_usage.percent}%")

    # get IO statistics since boot
    disk_io = psutil.disk_io_counters()
    print(f"Total read: {get_size(disk_io.read_bytes)}")
    print(f"Total write: {get_size(disk_io.write_bytes)}")

    # Network information
    print("="*40, "Network Information", "="*40)
    # get all network interfaces (virtual and physical)
    if_addrs = psutil.net_if_addrs()
    for interface_name, interface_addresses in if_addrs.items():
        for address in interface_addresses:
            print(f"=== Interface: {interface_name} ===")
            if str(address.family) == 'AddressFamily.AF_INET':
                print(f"  IP Address: {address.address}")
                print(f"  Netmask: {address.netmask}")
                print(f"  Broadcast IP: {address.broadcast}")
            elif str(address.family) == 'AddressFamily.AF_PACKET':
                print(f"  MAC Address: {address.address}")
                print(f"  Netmask: {address.netmask}")
                print(f"  Broadcast MAC: {address.broadcast}")
    # get IO statistics since boot
    net_io = psutil.net_io_counters()
    print(f"Total Bytes Sent: {get_size(net_io.bytes_sent)}")
    print(f"Total Bytes Received: {get_size(net_io.bytes_recv)}")
