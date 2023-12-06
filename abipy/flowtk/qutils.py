# coding: utf-8
"""
Collection of low-level tools to facilitate the interface with resource managers.

The preferred way of importing this module is:

    import qutils as qu
"""
from __future__ import annotations

import os
import json

from monty.string import is_string
from pymatgen.core.units import Time, Memory
from abipy.tools.typing import PathLike
from abipy.tools import duck
from abipy.tools.text import rm_multiple_spaces


def slurm_parse_timestr(s: str) -> Time:
    """
    A slurm time parser. Accepts a string in one the following forms:

        # "days-hours",
        # "days-hours:minutes",
        # "days-hours:minutes:seconds".
        # "minutes",
        # "minutes:seconds",
        # "hours:minutes:seconds",

    Returns: Time in seconds.

    Raises:
        `ValueError` if string is not valid.
    """
    days, hours, minutes, seconds = 0, 0, 0, 0

    if duck.is_number_like(s):
        return Time(s, "s")

    if '-' in s:
        # "days-hours",
        # "days-hours:minutes",
        # "days-hours:minutes:seconds".
        days, s = s.split("-")
        days = int(days)

        if ':' not in s:
            hours = int(float(s))
        elif s.count(':') == 1:
            hours, minutes = map(int, s.split(':'))
        elif s.count(':') == 2:
            hours, minutes, seconds = map(int, s.split(':'))
        else:
            raise ValueError("More that 2 ':' in string!")

    else:
        # "minutes",
        # "minutes:seconds",
        # "hours:minutes:seconds",
        if ':' not in s:
            minutes = int(float(s))
        elif s.count(':') == 1:
            minutes, seconds = map(int, s.split(':'))
        elif s.count(':') == 2:
            hours, minutes, seconds = map(int, s.split(':'))
        else:
            raise ValueError("More than 2 ':' in string!")

    return Time((days*24 + hours)*3600 + minutes*60 + seconds, "s")


def time2slurm(timeval: float, unit="s") -> str:
    """
    Convert a number representing a time value in the given unit (Default: seconds)
    to a string following the slurm convention: "days-hours:minutes:seconds".

    >>> assert time2slurm(61) == '0-0:1:1' and time2slurm(60*60+1) == '0-1:0:1'
    >>> assert time2slurm(0.5, unit="h") == '0-0:30:0'
    """
    d, h, m, s = 24*3600, 3600, 60, 1

    timeval = Time(timeval, unit).to("s")
    days, hours = divmod(timeval, d)
    hours, minutes = divmod(hours, h)
    minutes, secs = divmod(minutes, m)

    return "%d-%d:%d:%d" % (days, hours, minutes, secs)


def time2pbspro(timeval: float, unit="s") -> str:
    """
    Convert a number representing a time value in the given unit (Default: seconds)
    to a string following the PbsPro convention: "hours:minutes:seconds".

    >>> assert time2pbspro(2, unit="d") == '48:0:0'
    """
    h, m, s = 3600, 60, 1

    timeval = Time(timeval, unit).to("s")
    hours, minutes = divmod(timeval, h)
    minutes, secs = divmod(minutes, m)

    return "%d:%d:%d" % (hours, minutes, secs)


def time2loadlever(timeval: float, unit="s") -> str:
    """
    Convert a number representing a time value in the given unit (Default: seconds)
    to a string following the LoadLever convention. format hh:mm:ss (hours:minutes:seconds)

    >>> assert time2loadlever(2, unit="d") == '48:00:00'
    """
    h, m, s = 3600, 60, 1

    timeval = Time(timeval, unit).to("s")
    hours, minutes = divmod(timeval, h)
    minutes, secs = divmod(minutes, m)

    return "%d:%02d:%02d" % (hours, minutes, secs)


def timelimit_parser(s):
    """Convert a float or a string into time in seconds."""
    try:
        return Time(float(s), "s")
    except ValueError:
        return slurm_parse_timestr(s)


def any2mb(s):
    """Convert string or number to memory in megabytes."""
    if is_string(s):
        return int(Memory.from_str(s).to("Mb"))
    else:
        return int(s)


def slurm_get_jobs(username=None) -> dict[int, dict]:
    """
    Invoke squeue, parse output and return list of dictionaries with job info indexed by job id.
    """
    # Based on https://gist.github.com/stevekm/7831fac98473ea17d781330baa0dd7aa
    username = os.getlogin() if username is None else username
    import subprocess as sp
    process = sp.Popen(['squeue', '-u',  username, "-o", '%all'],
                       stdout=sp.PIPE, stderr=sp.PIPE, shell=False, universal_newlines=True)
    proc_stdout, proc_stderr = process.communicate()

    lines = proc_stdout.split('\n')
    header_line = lines.pop(0)
    header_cols = header_line.split('|')
    entries = []
    error_lines = [] # do something with this later
    for line in lines:
        parts = line.split('|')
        if len(parts) != len(header_cols):
            error_lines.append((len(parts), line, parts))
        else:
            d = {}
            for i, key in enumerate(header_cols):
                d[key] = parts[i]
                if key == "JOBID":
                    d[key] = int(d[key])
            entries.append(d)

    return {e["JOBID"]: e for e in entries}



class SlurmJobArray:
    """

    Example:

    header = '''\
#!/bin/bash

#SBATCH --account=battab
#SBATCH --job-name=abiml_md
#SBATCH --time=0-16:0:0
#SBATCH --partition=batch
#SBATCH --nodes=1             # 1 node has 128 cores
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1

conda activate env3.10
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
ulimit -s unlimited
'''

command = "abiml.py md"
arr_options = ["--help", "--version"]
job_array = SlurmJobArray(header, command, arr_options)
print(job_array)
queue_id = job_array.sbatch("job_array.sh")
    """

    def __init__(self, header: str, command: str, arr_options: list[str]):
        self.command = command
        if not self.command.endswith(" "): self.command += " "
        self.header = header
        self.arr_options = arr_options
        self.arr_options_str = rm_multiple_spaces("\n".join(arr_options))

    def __str__(self):
        # Add slurm array section.
        lines = self.header.splitlines()
        for il, line in enumerate(lines):
            if line.startswith("#SBATCH"):
                break
        else:
            raise ValueError("Cannot find line starting with #SBATCH")

        lines.insert(il, f"#SBATCH --array=0-{len(self.arr_options)-1}")
        header = "\n".join(lines)

        select_opts = r"""
# multiline string
OPTS_STRING="
%s
"

# Convert the multiline string to an array
IFS=$'\n' read -rd '' -a OPTS_LIST <<< "$OPTS_STRING"

# Index of the entry you want (0-based)
index=${SLURM_ARRAY_TASK_ID}

# Check if the index is within the range of the array
OPTS="${OPTS_LIST[index]}"
echo "Running entry at index $index:\nwith OPTS=$OPTS"

env
""" % (self.arr_options_str)

        end = f"""
{self.command} ${{OPTS}} > job_${{index}}.log 2> job_${{index}}.err

# Remove the file with the Slurm job id
me=$(basename "$0")
rm ${{me}}.qid
"""

        return header + select_opts + end

    def sbatch(self, slurm_filepath: PathLike) -> int:
        """
        Write slurm submission script to slurm_filepath and submits it.
        Return Slurm JOB id.
        """
        # Make sure no slurm job is already running by checking for a .qid file.
        path_qid = slurm_filepath + ".qid"
        if os.path.exists(path_qid):
            with open(path_qid, "rt") as fh:
                queue_id = int(fh.read().split("#"))
                err_msg = f"Found slurm job ID {queue_id} in {path_qid}" + \
                          "This usually indicates that a similar array job is already running\n" + \
                          f"If this not the case, please remove {path_qid} and rerun the script."
                raise RuntimeError(err_msg)

        with open(slurm_filepath, "wt") as fh:
            fh.write(str(self))

        queue_id = slurm_sbatch(slurm_filepath)

        print("Saving slurm job ID in:", path_qid)
        with open(path_qid, "wt") as fh:
            fh.write(str(queue_id) + " # Slurm job id")

        return queue_id


def slurm_sbatch(script_file) -> int:
    """
    Submit a job script to the queue with sbatch. Return Slurm JOB ID.
    """
    from subprocess import Popen, PIPE
    # need string not bytes so must use universal_newlines
    process = Popen(['sbatch', script_file], stdout=PIPE, stderr=PIPE, universal_newlines=True)
    out, err = process.communicate()

    # grab the returncode. SLURM returns 0 if the job was successful
    if process.returncode == 0:
        try:
            # output should of the form '2561553.sdb' or '352353.jessup' - just grab the first part for job id
            queue_id = int(out.split()[3])
            print(f"Job submission was successful and {queue_id=}")
            return queue_id
        except Exception as exc:
            # probably error parsing job code
            print('Could not parse job id following slurm...')
            raise exc
    else:
        raise RuntimeError(f"Error while submitting {script_file=}")
