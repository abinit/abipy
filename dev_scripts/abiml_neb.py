cat run_neb.py
#!/usr/bin/env python

from pathlib import Path
import sys
import os
import json
import numpy as np
import pandas as pd

from monty.os.path import find_exts
from abipy.tools.iotools import make_executable

from abipy.flowtk.qutils import SlurmJobArray

NN_NAMES = ["chgnet", "m3gnet", "matgl"]

TOP = Path(os.path.dirname(os.path.realpath(__file__)))

VERBOSE = 0


def run():
    # IMPORTANT: Please customize the slurm options according to your machine.
    conda_env = os.environ['CONDA_DEFAULT_ENV']
    print(f"Slurm script will be executed in {conda_env=}")

    header = f"""\
#!/bin/bash

#SBATCH --account=battab
#SBATCH --job-name=abiml_neb
#SBATCH --time=0-1:0:00
#SBATCH --partition=batch
#SBATCH --nodes=1             # 1 node has 128 cores
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
ulimit -s unlimited

# Important: This is the conda environment providing the NN potential(s)
conda activate {conda_env}
"""
    npaths = 3

    print(f"""\
Performing NEB calculations with the following parameters:

{npaths=}
{NN_NAMES=}
""")

    arr_options = []
    for nn_name in NN_NAMES:
        for i in range(npaths):
            dirpath = Path(f"path{i+1}")
            workdir = dirpath / nn_name
            ini_vasp = str(dirpath / "ini.vasp")
            fin_vasp = str(dirpath / "fin.vasp")
            opts = f"{ini_vasp} {fin_vasp} --nn-name {nn_name} -w {str(workdir)}"
            arr_options.append(opts)

    command = "abiml.py neb"
    job_array = SlurmJobArray(header, command, arr_options)
    job_array.sbatch("job_array.sh")


def process():
    """Post-process the results."""
    json_paths = find_exts(TOP, "neb_data.json")
    neb_data_list = []
    for path in json_paths:
        if VERBOSE: print("About to read json data from", path)
        parts = path.split(os.sep)
        path_index, nn_name = parts[-3], parts[-2]
        path_index = int(path_index.replace("path", ""))
        if VERBOSE: print(parts, "\n", path_index, nn_name)
        with open(path, "rt") as fh:
            d = json.load(fh)
            d["path"] = path_index
            d["nn_name"] = nn_name
            neb_data_list.append(d)

    # Sort dict by path.
    neb_data_list = sorted(neb_data_list, key=lambda d: d["path"])

    keys = [
        "nn_name",
        "path",
        "barrier_with_fit",
        "barrier_without_fit",
        "energy_change_with_fit",
        "energy_change_without_fit",
    ]

    d_list = []
    for d in neb_data_list:
        d_list.append({k: d[k] for k in keys})

    df = pd.DataFrame(d_list)
    #df.to_csv("XXX_ML_barriers.csv")
    print(df)

    from abipy.tools.plotting import get_axarray_fig_plt
    ax_list, fig, plt = get_axarray_fig_plt(
        None, nrows=1, ncols=3, sharex=True, sharey=True, squeeze=False)
    ax_list = ax_list.ravel()
    cmap = plt.get_cmap("jet")
    fontsize = 8

    for ix, (ax, nn_name) in enumerate(zip(ax_list, NN_NAMES)):
        my_data_list = [d for d in neb_data_list if d["nn_name"] == nn_name]
        my_data_list = sorted(my_data_list, key=lambda d: d["path"])

        for i, data in enumerate(my_data_list):
            enes = np.array(data["energies_images"])
            ax.plot(enes - enes[0], label=f"path{i+1}", color=cmap(i/len(my_data_list)))

        ax.set_title(nn_name)
        ax.set_xlabel('Image index', fontsize=fontsize)
        ax.set_ylabel('$\Delta$ energy [eV]', fontsize=fontsize)
        ax.legend(loc="best", shadow=True, fontsize=fontsize)

    fig.suptitle("K2Cu3Ge5O14")
    plt.show()


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "md_gen":
        md_generate(sys.argv[2])
    elif cmd == "run":
        run()
    elif cmd == "process":
        process()
    else:
        raise ValueError(f"Invalid {cmd=}")
