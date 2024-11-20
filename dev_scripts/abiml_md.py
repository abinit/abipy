#!/usr/bin/env python
import sys
import os

from abipy.flowtk.qutils import SlurmJobArray


def get_slurm_header(conda_env):
    # IMPORTANT: You need customize the slurm options below according to your machine.
    #conda_env = os.environ['CONDA_DEFAULT_ENV']
    #print(f"Slurm script will be executed in {conda_env=}")

    slurm_header = f"""\
#!/bin/bash

#SBATCH --account=battab
#SBATCH --job-name=abiml_md
#SBATCH --time=0-1:0:00
#SBATCH --partition=batch
#SBATCH --nodes=1             # 1 node has 128 cores
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
ulimit -s unlimited

# Important: This is the conda environment providing the NN potential(s) nn_names
source $HOME/.bashrc
conda activate {conda_env}
"""
    return slurm_header


def main():
    filepath = "LLZO_cubic.vasp"
    nn_names_env = {"mace_mp": "mace_mp"}
    temperature_list = [600, 800, 1000, 1200]
    steps = 4000
    steps = 100

    print(f"""\
Performing MD calculations with the following parameters:

{filepath=}
{nn_names_env=}
{temperature_list=}
{steps=}
""")

    arr_options = []
    for nn_name, conda_env in nn_names_env.items():
        slurm_header = get_slurm_header(conda_env)
        for temperature in temperature_list:
            workdir = f"nn-{nn_name}_T-{temperature}"
            opts = f"{filepath} --nn-name {nn_name} --temperature {temperature} --timestep 1 \
            --loginterval 100 --steps {steps} --ensemble nvt -w {workdir}"
            arr_options.append(opts)

        job_array = SlurmJobArray(slurm_header, "abiml.py md", arr_options)
        job_array.sbatch("job_array.sh")


if __name__ ==  "__main__":
    main()
