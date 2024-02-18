#!/usr/bin/env python

import os
from abipy.flowtk.qutils import SlurmJobArray

def main():
    # IMPORTANT: You need customize the slurm options according to your machine.
    conda_env = os.environ['CONDA_DEFAULT_ENV']
    print(f"Slurm script will be executed in {conda_env=}")

    header = f"""\
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
conda activate {conda_env}
"""
    filepath = "LLZO_cubic.vasp"
    nn_names = ["m3gnet", "chgnet"]
    temp_list = [600, 800, 1000, 1200]
    steps = 4000

    print(f"""\
Performing MD calculations with the following parameters:

{filepath=}
{nn_names=}
{temp_list=}
{steps=}
""")

    arr_options = []
    for nn_name in nn_names:
        for temperature in temp_list:
            workdir = f"nn-{nn_name}_T-{temperature}"
            opts = f"{filepath} --nn-name {nn_name} --temperature {temperature} --timestep 1 \
            --loginterval 100 --steps {steps} --ensemble nvt -w {workdir}"
            arr_options.append(opts)

    job_array = SlurmJobArray(header, "abiml.py md", arr_options)
    job_array.sbatch("job_array.sh")


if __name__ ==  "__main__":
    main()
