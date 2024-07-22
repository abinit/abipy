# coding: utf-8
"""
Collection of low-level tools to facilitate the interface with resource managers.

The preferred way of importing this module is:

    import qutils as qu
"""
from __future__ import annotations

import os

from subprocess import Popen, PIPE, run
from monty.string import is_string
from pymatgen.core.units import Time, Memory, UnitError
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
        try:
            # latest pymatgen version (as of july 2024)
            mem = int(Memory.from_str(s.upper()).to("MB"))
        except (KeyError, UnitError):  # For backward compatibility with older pymatgen versions
            try:
                mem = int(Memory.from_str(s.replace("B", "b")).to("Mb"))
            except AttributeError:  # For even older pymatgen versions
                mem = int(Memory.from_string(s.replace("B", "b")).to("Mb"))
        return mem
    else:
        return int(s)


def slurm_get_jobs() -> dict[int, dict]:
    """
    Invoke squeue, parse output and return list of dictionaries with job info indexed by job id.
    """
    # Based on https://gist.github.com/stevekm/7831fac98473ea17d781330baa0dd7aa
    process = Popen(["squeue", "--me", "-o", '%all'],
                       stdout=PIPE, stderr=PIPE, shell=False, universal_newlines=True)
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
        return queue_id


def slurm_write_and_sbatch(script_filepath: str, slurm_script_str: str) -> int:
    """
    Write job script and submit it to the queue with Slurm sbatch. Return Slurm JOB ID.
    """
    with open(script_filepath, "wt") as fh:
        fh.write(slurm_script_str)
        return slurm_sbatch(script_filepath)


def slurm_sbatch(slurm_filepath: PathLike) -> int:
    """
    Submit a job script to the queue with Slurm sbatch. Return Slurm JOB ID.
    """
    from monty.os import cd
    dirpath = os.path.dirname(slurm_filepath)
    #print("dirpath", dirpath)
    with cd(dirpath):

        # need string not bytes so must use universal_newlines
        slurm_filepath = str(slurm_filepath)
        process = Popen(['sbatch', slurm_filepath], stdout=PIPE, stderr=PIPE, universal_newlines=True)
        out, err = process.communicate()

        # grab the returncode. SLURM returns 0 if the job was successful
        if process.returncode == 0:
            try:
                # output should of the form '2561553.sdb' or '352353.jessup' - just grab the first part for job id
                queue_id = int(out.split()[3])
                path_qid = slurm_filepath + ".qid"
                print(f"Job submission was successful and queue_id: {queue_id}")
                print("Saving slurm job ID in:", path_qid)
                with open(path_qid, "wt") as fh:
                    fh.write(str(queue_id) + " # Slurm job id")
                return queue_id

            except Exception as exc:
                # probably error parsing job code
                print('Could not parse job id following slurm...')
                raise exc
        else:
            raise RuntimeError(f"Error while submitting {slurm_filepath=} with {process.returncode=},\n{out=}\n{err=}")


def get_sacct_info():
    """
    Run the sacct command to get the job information
    """
    try:

        result = run(['sacct', '--format=JobID,JobName,Partition,Account,AllocCPUS,State,ExitCode', '--noheader'],
                      stdout=PIPE, stderr=PIPE, text=True)

        # Check if the command was successful
        if result.returncode != 0:
            print(f"Error running sacct: {result.stderr}")
            return None

        # Process the output
        jobs_info = result.stdout.strip().split('\n')
        jobs = [dict(zip(['JobID', 'JobName', 'Partition', 'Account', 'AllocCPUS', 'State', 'ExitCode'], job.split())) 
                for job in jobs_info]
        return jobs

    except Exception as e:
        print(f"An error occurred: {e}")
        return None


def get_completed_job_info(job_id: int | str):
    try:
        # Define the fields we want to retrieve
        fields = "JobID,JobName,Partition,Account,AllocCPUS,State,ExitCode,Start,End,Elapsed,TotalCPU,MaxRSS"

        # Run the sacct command with the specified fields for the given job ID
        result = run(
            ['sacct', '--jobs', str(job_id), '--format', fields, '--noheader', '--parsable2'],
            stdout=PIPE, stderr=PIPE, text=True
        )

        # Check if the command was successful
        if result.returncode != 0:
            print(f"Error running sacct: {result.stderr}")
            return None

        # Process the output
        lines = result.stdout.strip().split('\n')
        jobs = [dict(zip(fields.split(','), line.split('|'))) for line in lines]
        return jobs

    except Exception as e:
        print(f"An error occurred: {e}")
        return None


def get_slurm_template(body: str) -> str:
    """
    Return template for slurm submission that is supposed to be customized by the user.
    """
    return f"""\
#!/bin/bash
#
# Please CUSTOMIZE this section according to your cluster and the type of calculation
#
#SBATCH --job-name=my_job
#SBATCH --output=%j_%x.slurm.out
#SBATCH --error=%j_%x.slurm.err
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=2:00:00
#SBATCH --account=htforft

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# ------------------------------------------------------------------------------
# Printing some information
# ------------------------------------------------------------------------------

echo "------------------- Job info -------------------"
echo "job_id             : $SLURM_JOB_ID"
echo "jobname            : $SLURM_JOB_NAME"
echo "queue              : $SLURM_JOB_PARTITION"
echo "qos                : $SLURM_JOB_QOS"
echo "account            : $SLURM_JOB_ACCOUNT"
echo "submit dir         : $SLURM_SUBMIT_DIR"
echo "number of mpi tasks: $SLURM_NTASKS tasks"
echo "OMP_NUM_THREADS    : $OMP_NUM_THREADS"
echo "number of gpus     : $SLURM_GPUS_ON_NODE"
echo "------------------- Node list ------------------"
echo $SLURM_JOB_NODELIST

echo "---------------- Printing limits ---------------"
ulimit -s unlimited
ulimit -a

# ------------------------------------------------------------------------------
# Setting up the environment
# ------------------------------------------------------------------------------

echo "----------------- Environment ------------------"

source $HOME/vasp.6.2.1/modules.sh
module list

echo "--------------- Running the code ---------------"
echo -n "This run started on: "
date

{body}

echo -n "This run completed on: "
date
"""


def get_custodian_template() -> str:
    return """\
#!/usr/bin/env python

from custodian.custodian import Custodian
from custodian.vasp.jobs import VaspJob
from custodian.vasp.handlers import VaspErrorHandler, UnconvergedErrorHandler

# List of error handlers
handlers = [
    VaspErrorHandler(),        # Handles common VASP errors
    UnconvergedErrorHandler()  # Handles unconverged calculations
]

# Define the VASP job with appropriate command and parameters
jobs = [
    VaspJob(
        #["mpirun", "vasp_std"],  # Replace NCORES with the number of cores
        #["mpirun", "-np", "1", "vasp_std"],  # Replace NCORES with the number of cores
        ["mpirun", "vasp_std"],  # Replace NCORES with the number of cores
        auto_npar=False,
        final=True
    )
]

# Create the Custodian instance with handlers and jobs
c = Custodian(handlers, jobs, max_errors=5)

# Run the Custodian job
c.run()
"""