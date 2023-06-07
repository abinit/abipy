#!/bin/bash
#set -e  # exit on first error
#set -ev  # exit on first error, print each command
#set -x # or set -o xtrace expands variables and prints a little + sign before the line

#echo "Running preliminary tests for Abipy with Abinit..."
#echo "PMG_MAPI_KEY: foobar" > ${HOME}/.pmgrc.yaml

env_name="abigui"
echo "Creating ${env_name} conda environment from scratch..."

conda create -n ${env_name} -y python=3.9
#source activate ${env_name} 
conda activate ${env_name} 

pip install -r ./requirements.txt
pip install -r ./requirements-optional.txt
python setup.py install

#conda install abinit -c conda-forge

list_of_repos=( 
  ONCVPSP-PBE-SR-PDv0.4
  ATOMPAW-LDA-JTHv1.1
)

for repo in "${list_of_repos[@]}"; do
    echo "Installing pseudopotential repo: " ${repo}
    abips.py install ${repo}
done

echo "Installation completed."
serve_cmd="abigui.py --num_procs 4 --has-remote-server --max_size_mb=150"
echo "Now you can start the sever with e.g.:\n\n\t${serve_cmd}"
