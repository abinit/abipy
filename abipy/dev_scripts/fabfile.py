#!/usr/bin/env python
"""Fabric file providing useful tools for the synchronization of the code on the CECI clusters."""
from __future__ import division, print_function

import os

from cStringIO import StringIO
from fabric.api import local, settings, abort, run, cd, env, prefix, put
from contextlib import contextmanager 

env.user = "gmatteo"

# The directory of the virtual environment.
VENV = "~/VENV-2.7"

#env.password = ""
def all_hosts():
    env.hosts = [
        #"green.cism.ucl.ac.be",
        #"manneback.cism.ucl.ac.be",
        "hmem.cism.ucl.ac.be",
        "lemaitre2.cism.ucl.ac.be",
        "vega.ulb.ac.be",
        "dragon1.umons.ac.be",
        "hercules.ptci.unamur.be",
    ]

git_reposdir = '~gmatteo/git_repos/'

git_urls = {
    "abipy":     "https://github.com/gmatteo/abipy.git",
    "pymatgen":   "https://github.com/gmatteo/pymatgen.git",
    #"abipy": "git@github.com:gmatteo/abipy.git",
    #"pymatgen": "git@github.com:gmatteo/pymatgen.git",

}

git_repospaths = [os.path.join(git_reposdir, dirpath) for dirpath in git_urls]

@contextmanager
def virtualenv(venv_dir):
    with prefix("source %s" % os.path.join(venv_dir, "bin", "activate")):
        yield

def _exists(path):
    with settings(warn_only=True):
        return not run('test -e %s' % path).failed

def git_deploy():
    """Synchronize the git branches with the master branch located on github."""
    # Create ~/git_repos and clone the repositories if this is the first time.
    if not _exists(git_reposdir):
        run("mkdir %s" % git_reposdir)
        with cd(git_reposdir):
            for name, url in git_urls.items():
                run("git clone %s" % url)

    # Pull from github.
    for apath in git_repospaths:
        with cd(apath):
            run("git pull")
            with virtualenv(VENV):
                run("python setup.py clean")
                run("python setup.py install")

def pytest():
    """Run the test suite with py.test"""
    for apath in git_repospaths:
        if "pymatgen" in apath: continue
        with cd(apath):
            with virtualenv(VENV):
                #run("nosetests -v")
                run("py.test -v")

def pip_install(*options):
    """
    Execute `pip install options` on the remote hosts. 
    (use a virtual environment).
    """
    with virtualenv(VENV):
        run("pip install %s" % " ".join(o for o in options))

def py_version():
    out = run("python --version", quiet=True)
    py_version = out.split()[1]

    out = StringIO()
    run("module available", stdout=out)

    py_modules = []
    for o in out.getvalue().splitlines():
        py_modules.extend([s for s in o.split() if "python" in s])
            
    print("default python:", py_version, "python modules:", py_modules)


def bzr_deploy():
    url = "bzr+ssh://forge.abinit.org/abinit/gmatteo/7.5.4-private"
    to_location = os.path.join(*url.split("/")[-2:]).replace("/", "_")

    # Create ~/bzr_repos and clone the repositories if this is the first time.
    bzr_reposdir = '~gmatteo/bzr_repos/'
    repo_path = os.path.join(bzr_reposdir, to_location)
    print("repo_path", repo_path, "to_location", to_location)

    with virtualenv(VENV):
        if not _exists(bzr_reposdir):
            run("mkdir %s" % bzr_reposdir)

        if not _exists(repo_path):
            with cd(bzr_reposdir):
                run("bzr get %s %s" % (url, to_location))

        # Pull the branch
        with cd(repo_path):
            run("bzr pull")


def upload_file(local_path, remote_path, mode=None):
    """
    Copy local_path to host:remote_path. 
    If mode is not None `chmod mode remote_path` is executed.
    """
    put(local_path, remote_path)

    if mode is not None:
        run("chmod %s %s" % (mode, remote_path))


def upload_rsa():
    if not _exists("~/.ssh/"): run("mkdir ~/.ssh")
    put("~/.ssh/id_rsa.ceci", "~/.ssh/id_rsa.ceci")
    run("chmod go-rwx ~/.ssh/id_rsa.ceci")

    put("~/.ssh/gmatteo_dsa", "~/.ssh/gmatteo_dsa")
    run("chmod go-rwx ~/.ssh/gmatteo_dsa")


def use_ssh_agent():
    run("eval $(ssh-agent)")
    run("ssh-add ~/.ssh/id_rsa.ceci")
