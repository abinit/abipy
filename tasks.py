"""
Deployment file to facilitate AbiPy releases.
Use invoke --list to get list of tasks
"""
from __future__ import annotations

import shutil
import json
import os
import re
import subprocess
import webbrowser
from datetime import datetime, timezone
from typing import TYPE_CHECKING

import requests
from invoke import task
from monty.os import cd


if TYPE_CHECKING:
    from invoke import Context

#from abipy.core.release import __version__ as CURRENT_VER

ABIPY_ROOTDIR = os.path.dirname(__file__)
DOCS_DIR = os.path.join(ABIPY_ROOTDIR, "docs")


@task
def pull(ctx):
    """"Execute `git stash && git pull --recurse-submodules && git stash apply && makemake`"""
    ctx.run("git stash")
    ctx.run("git pull --recurse-submodules")
    ctx.run("git stash apply")


@task
def submodules(ctx):
    """Update submodules."""
    with cd(ABIPY_ROOTDIR):
        # https://stackoverflow.com/questions/1030169/easy-way-to-pull-latest-of-all-git-submodules
        ctx.run("git submodule update --remote --init", pty=True)
        ctx.run("git submodule update --recursive --remote", pty=True)


@task
def make_doc(ctx):
    """Build the website"""
    with cd(DOCS_DIR):
        ctx.run("make clean")
        ctx.run("make", env=dict(READTHEDOCS="1"), pty=True)
        open_doc(ctx)


@task
def open_doc(ctx):
    """Open the index.html in docs/_build/html/."""
    import webbrowser
    webbrowser.open_new_tab("file://" + os.path.join(ABIPY_ROOTDIR, "docs/_build/html/index.html"))


@task
def twine(ctx):
    """Upload new release with twine."""
    with cd(ABIPY_ROOTDIR):
        ctx.run("rm dist/*.*", warn=True)
        ctx.run("python setup.py register sdist bdist_wheel")
        ctx.run("twine upload dist/*")


@task
def pytest(ctx):
    """Execute pytest."""
    pytest_cmd = r"""\
pytest -n 2 --cov-config=.coveragerc --cov=abipy -v --doctest-modules abipy \
    --ignore=abipy/integration_tests --ignore=abipy/data/refs --ignore=abipy/scripts/ \
    --ignore=abipy/examples/plot --ignore=abipy/examples/flows --ignore=abipy/gui
"""
    with cd(ABIPY_ROOTDIR):
        ctx.run(pytest_cmd, pty=True)


@task
def style(ctx):
    """Execute pycodestyle."""
    with cd(ABIPY_ROOTDIR):
        ctx.run("pycodestyle 2>&1 | tee style.log", pty=True)
        #ctx.run("pydocstyle abipy | tee -a style.log", pty=True)


@task
def flake(ctx):
    """Execute flake8."""
    with cd(ABIPY_ROOTDIR):
        ctx.run("flake8 --count --show-source --statistics | tee -a style.log", pty=True)


@task
def plots(ctx):
    """Run the scripts in the examples/plots directory."""
    with cd(os.path.join(ABIPY_ROOTDIR, "abipy", "examples")):
        ctx.run("_runplots.py", pty=True)


@task
def flows(ctx):
    """Run the scripts in the examples/flows directory."""
    with cd(os.path.join(ABIPY_ROOTDIR, "abipy", "examples")):
        ctx.run("_runflows.py", pty=True)


@task
def pygrep(ctx, pattern):
    """
    Grep for `pattern` in all py files contained in the project.
    """
    # grep -r -i --include \*.h
    # Syntax notes:
    #    -r - search recursively
    #    -i - case-insensitive search
    #    --include=\*.${file_extension} - search files that match the extension(s) or file pattern only
    with cd(os.path.join(ABIPY_ROOTDIR, "abipy",)):
        cmd = 'grep -r -i --color --include "*.py" "%s" .' % pattern
        print("Executing:", cmd)
        ctx.run(cmd, pty=True)


@task
def update_vars(ctx, abinit_repo_path):
    """Update the database of Abinit variables."""
    abinit_repo_path = os.path.abspath(abinit_repo_path)
    dir_with_pyfiles = os.path.join(ABIPY_ROOTDIR, "abipy", "abio", "abivar_database")

    with cd(ABIPY_ROOTDIR):
        local_files = [f for f in os.listdir(dir_with_pyfiles) if f.startswith("variables_")]
        for local_file in local_files:
            local_file = os.path.join(dir_with_pyfiles, local_file)
            # "vimdiff $ABINIT_REPOPATH/abimkdocs/variables_abinit.py variables_abinit.py
            source = os.path.join(abinit_repo_path, "abimkdocs", local_file)
            cmd = f"vimdiff {source} {local_file}"
            print(f"Executing: {cmd}")
            os.system(cmd)


@task
def pyclean(ctx):
    """remove all pyc files and all __pycache__ directory."""
    def rm_pycaches(top: str) -> None:
        """remove __pycache__ directory."""
        count = 0
        for dirpath, dirnames, filenames in os.walk(top):
            for d in dirnames:
                if not d ==  "__pycache__": continue
                path = os.path.join(dirpath, d)
                #print("Will remove %s" % path)
                shutil.rmtree(path)
                count += 1

        print("Removed %d __pycache__ directories" % count)

    def rm_pycfiles(top: str) -> int:
        """remove all pyc files."""
        count = 0
        for dirpath, dirnames, filenames in os.walk(top):
            for f in filenames:
                if not f.endswith(".pyc"): continue
                path = os.path.join(dirpath, f)
                #print("Will remove %s" % path)
                os.remove(path)
                count += 1

        print("Removed %d .pyc files" % count)
        return count

    top = ABIPY_ROOTDIR
    rm_pycaches(top)
    rm_pycfiles(top)


@task
def tuna(ctx: Context) -> None:
    """Execute tuna import profiler."""
    cmd = 'python -X importtime -c "import abipy" 2> __abipy_import.log'
    print("Executing:", cmd)
    ctx.run(cmd, pty=True)
    cmd = "tuna __abipy_import.log"
    ctx.run(cmd, pty=True)


def generate_rst(package_path, output_dir):
    """
    Generate .rst files for all Python modules in a package.

    Parameters:
        package_path: Path to the root of the Python package.
        output_dir: Directory where .rst files will be generated.
    """
    package_path = os.path.join(ABIPY_ROOTDIR, "abipy")
    output_dir = os.path.join(DOCS_DIR, "api")
    #if not os.path.exists(output_dir):
    #    os.makedirs(output_dir)

    for root, dirs, files in os.walk(package_path):
        rel_dir = os.path.relpath(root, package_path)
        if rel_dir == ".":
            rel_dir = ""
        rst_dir = os.path.join(output_dir, rel_dir)
        #if not os.path.exists(rst_dir):
        #    os.makedirs(rst_dir)

        # Process Python modules and packages
        for file in files:
            if file.endswith(".py") and file != "__init__.py":
                module_name = file[:-3]  # Remove .py extension
                rst_file = os.path.join(rst_dir, f"{module_name}.rst")
                #with open(rst_file, "w") as f:
                f = sys.stdout
                full_module = (
                    f"{os.path.basename(package_path)}.{rel_dir.replace('/', '.')}.{module_name}"
                    if rel_dir
                    else f"{os.path.basename(package_path)}.{module_name}"
                )
                f.write(f"{module_name}\n{'=' * len(module_name)}\n\n")
                f.write(f".. automodule:: {full_module}\n")
                f.write("    :members:\n")
                f.write("    :undoc-members:\n")
                f.write("    :show-inheritance:\n")

        # Process subpackages
        for dir_ in dirs:
            if os.path.isfile(os.path.join(root, dir_, "__init__.py")):
                package_name = dir_
                rst_file = os.path.join(rst_dir, f"{package_name}.rst")
                #with open(rst_file, "w") as f:
                f = sys.stdout
                full_package = (
                    f"{os.path.basename(package_path)}.{rel_dir.replace('/', '.')}.{package_name}"
                    if rel_dir
                    else f"{os.path.basename(package_path)}.{package_name}"
                )
                f.write(f"{package_name}\n{'=' * len(package_name)}\n\n")
                f.write(f".. automodule:: {full_package}\n")
                f.write("    :members:\n")
                f.write("    :undoc-members:\n")
                f.write("    :show-inheritance:\n")

        print(f"Documentation files have been generated in '{output_dir}'.")
