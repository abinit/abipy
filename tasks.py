"""
Deployment file to facilitate AbiPy releases.
Use invoke --list to get list of tasks
"""
from __future__ import annotations

import shutil
import json
import sys
import os
import re
import subprocess
import webbrowser

from datetime import datetime, timezone
from typing import TYPE_CHECKING
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
        ctx.run("touch api/index.rst", warn=True)
        ctx.run("sphinx-apidoc --implicit-namespaces -M -d 1 -o api -f ../abipy ../**/tests/* ../abipy/benchmarks ../abipy/data ../abipy/integration_tests ../abipy/test_files ../abipy/examples")

        rst_files = sorted([f for f in os.listdir(os.path.join(DOCS_DIR, "api")) if f.endswith(".rst") and f != "modules.rst"])
        rst_files = [3 * " " + f for f in rst_files]
        #print(rst_files)

        header = """\
.. _api-index:

=================
API documentation
=================

.. toctree::
   :maxdepth: 1

""" + "\n".join(rst_files)

        with open(os.path.join(DOCS_DIR, "api", "index.rst"), "wt") as fh:
            fh.write(header)

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
    print(f"{abinit_repo_path=}")

    with cd(ABIPY_ROOTDIR):
        local_files = [f for f in os.listdir(dir_with_pyfiles) if f.startswith("variables_")]
        for local_file in local_files:
            target_file = os.path.join(dir_with_pyfiles, local_file)
            # "vimdiff $ABINIT_REPOPATH/abimkdocs/variables_abinit.py variables_abinit.py
            source = os.path.join(abinit_repo_path, "abimkdocs", local_file)
            cmd = f"vimdiff {source} {target_file}"
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

        print(f"Removed {count} __pycache__ directories")

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

        print(f"Removed {count} .pyc files")
        return count

    rm_pycaches(ABIPY_ROOTDIR)
    rm_pycfiles(ABIPY_ROOTDIR)


@task
def tuna(ctx: Context) -> None:
    """Execute tuna import profiler."""
    cmd = 'python -X importtime -c "import abipy" 2> __abipy_import.log'
    print("Executing:", cmd)
    ctx.run(cmd, pty=True)
    cmd = "tuna __abipy_import.log"
    ctx.run(cmd, pty=True)


@task
def git_info(ctx: Context, top_n=20) -> None:
    """Scan Git history for largest files"""

    def get_git_objects():
        """Return list of all Git objects (hash, path)."""
        result = subprocess.run(
            ['git', 'rev-list', '--objects', '--all'],
            stdout=subprocess.PIPE,
            text=True,
            check=True
        )
        objects = []
        for line in result.stdout.splitlines():
            parts = line.split(' ', 1)
            if len(parts) == 2:
                objects.append((parts[0], parts[1]))
        return objects

    def get_blob_sizes(hashes):
        """Return a dict of {hash: (size_in_bytes, path)} for blobs."""
        input_text = '\n'.join(hashes)
        result = subprocess.run(
            ['git', 'cat-file', '--batch-check=%(objectname) %(objecttype) %(objectsize)'],
            input=input_text,
            stdout=subprocess.PIPE,
            text=True,
            check=True
        )

        sizes = {}
        for line in result.stdout.splitlines():
            obj_hash, obj_type, obj_size = line.split()
            if obj_type == "blob":
                sizes[obj_hash] = int(obj_size)
        return sizes

    print("Scanning Git history for largest files...")

    objects = get_git_objects()
    hashes = [obj[0] for obj in objects]
    paths = {obj[0]: obj[1] for obj in objects}

    sizes = get_blob_sizes(hashes)

    sorted_blobs = sorted(
        ((size, paths[_hash], _hash) for _hash, size in sizes.items() if _hash in paths),
        reverse=True
    )

    print(f"\nTop {top_n} largest files ever committed:")
    for size, path, obj_hash in sorted_blobs[:top_n]:
        print(f"{size / (1024*1024):7.2f} MB\t{path}")

    ctx.run("git count-objects -vH", pty=True)
