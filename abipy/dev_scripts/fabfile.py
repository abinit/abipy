#!/usr/bin/env python
"""
Fabric file providing useful tools for the synchronization of the code on the CECI clusters.
"""
from __future__ import division, print_function

import os
import cStringIO as StringIO

from fabric.api import local, settings, abort, run, cd, env, prefix, put
from contextlib import contextmanager as _contextmanager 

env.user = "gmatteo"

USER_HOME = "~" + env.user

# The directory of the virtual environment.
VENV = "~/VENV-2.7"

#env.password = ""

# We store our git branches in this directory 
GIT_REPOSDIR = USER_HOME + "/git_repos"

git_urls = {
    "abipy":    "https://github.com/gmatteo/abipy.git",
    "pymatgen": "https://github.com/gmatteo/pymatgen.git",
    "matplotlib": "https://github.com/matplotlib/matplotlib.git",
    #"abipy": "git@github.com:gmatteo/abipy.git",
    #"pymatgen": "git@github.com:gmatteo/pymatgen.git",
    #"matplotlib": "git@github.com:matplotlib/matplotlib.git",
    #git://git.gnome.org/pygtk
}

git_repospaths = [os.path.join(GIT_REPOSDIR, dirpath) for dirpath in git_urls]

# We store our bzr branches in this directory 
bzr_reposdir = USER_HOME + "/bzr_repos"

bzr_branchurl = "bzr+ssh://forge.abinit.org/abinit/gmatteo/7.5.4-private"

to_location = os.path.join(*bzr_branchurl.split("/")[-2:]).replace("/", "_")

repo_path = os.path.join(bzr_reposdir, to_location)
build_path = os.path.join(repo_path, "build")


@_contextmanager
def _virtualenv(venv_dir):
    with prefix("source %s" % os.path.join(venv_dir, "bin", "activate")):
        yield

def _exists(path):
    """True if path exists on the remote host."""
    return run('test -e %s' % path, warn_only=True, quiet=True).succeeded

def _cpu_count():
    """Returns the number of CPUs in the remote host."""
    with _virtualenv(VENV):
        cmd = run("python -c 'import multiprocessing, sys; sys.exit(multiprocessing.cpu_count())'") 
        return cmd.return_code


def all_hosts():
    """
    Used to run the command on all the CECI clusters.
    example: fab all_hosts pip_install abipy
    """
    env.hosts = [
        #"green.cism.ucl.ac.be",
        #"manneback.cism.ucl.ac.be",
        "hmem.cism.ucl.ac.be",
        "lemaitre2.cism.ucl.ac.be",
        "vega.ulb.ac.be",
        "dragon1.umons.ac.be",
        "hercules.ptci.unamur.be",
    ]

def git_pull():
    """
    Synchronize the git branches with the master branches located on github.
    """
    # Create ~/git_repos and clone the repositories if this is the first time.
    if not _exists(GIT_REPOSDIR):
        run("mkdir %s" % GIT_REPOSDIR)

    with cd(GIT_REPOSDIR):
        for name, url in git_urls.items():
            if not _exists(name):
                run("git clone %s" % url)

    # Pull from github.
    for apath in git_repospaths:
        with cd(apath):
            run("git pull")

def git_install():
    """Install the code of the git branches."""
    git_pull()
    for apath in git_repospaths:
        with cd(apath), _virtualenv(VENV):
            run("python setup.py clean")
            #run("rm -rf build")
            #run("rm -rf sdist")
            run("python setup.py install")


def pytest(opts=""):
    """
    Run the test suite with py.test. Usage: pytests:"-v"
    """
    for apath in git_repospaths:
        if "pymatgen" in apath: continue
        with cd(apath), _virtualenv(VENV):
            run("py.test %s" % opts)

def pip_install(*options):
    """
    Execute `pip install options` on the remote hosts. 
    Examples: 

        pip_install foo bar

    .. note:
        use the virtual environment VENV
    """
    with _virtualenv(VENV):
        run("pip install %s" % " ".join(o for o in options))

def py_version():
    """
    Show the default python version and the python modules 
    available on the remote machines.
    """
    out = run("python --version", quiet=True)
    py_version = out.split()[1]

    out = StringIO.StringIO()
    run("module available", stdout=out)

    py_modules = []
    for o in out.getvalue().splitlines():
        py_modules.extend([s for s in o.split() if "python" in s])
            
    print("default python:", py_version, "python modules:", py_modules)

def bzr_pull():
    """
    Upload the bzr repository on the remote clusters.
    """
    # Create ~/bzr_repos and clone the repositories if this is the first time.
    print("repo_path", repo_path, "to_location", to_location)

    with _virtualenv(VENV):
        if not _exists(bzr_reposdir):
            run("mkdir %s" % bzr_reposdir)

        if not _exists(repo_path):
            with cd(bzr_reposdir):
                run("bzr get %s %s" % (bzr_branchurl, to_location))

        # Pull the branch
        with cd(repo_path):
            run("bzr pull")
            run("bzr status")


def upload_file(local_path, remote_path, mode=None):
    """
    Copy local_path to host:remote_path. 
    If mode is not None `chmod mode remote_path` is executed.
    """
    put(local_path, remote_path)

    if mode is not None:
        run("chmod %s %s" % (mode, remote_path))


def wget(urls):
    """
    Download files on the remove host from URLS to ~/Downloads

    fab wget:"foo.tar.gz bar.tar.gz"
    """
    downloads_dir = "Downloads"

    with cd (USER_HOME):
        if not _exists(downloads_dir):
            run("mkdir %s" % downloads_dir)

        with cd(downloads_dir):
            for u in urls.split():
                run("wget %s " % u)


def upload_key(key):
    """
    key: path of the SSH key to upload on the cluster.
    """
    if not _exists("~/.ssh/"): run("mkdir ~/.ssh")
    filename = os.path.basename(key)
    remotepath = "~/.ssh/" + filename
    put(key, remotepath)
    run("chmod go-rwx %s" % remotepath)

#
#def use_ssh_agent():
#    run("eval $(ssh-agent)")
#    run("ssh-add ~/.ssh/id_rsa.ceci")

def abinit_version():
    """Show the version of Abinit available on the remote host."""
    if run("abinit --version", warn_only=True, quiet=True).return_code:
        run("mpirun abinit --version") # Try to prepend mpirun

def abinit_build(self):
    """Show the Abinit build parameters."""
    if run("abinit --build", warn_only=True, quiet=True).return_code:
        run("mpirun abinit --build") # Try to prepend mpirun

def abinit_makemake():
    """Run abinit makemake to generate configure script."""
    with cd(repo_path):
        run("./config/scripts/makemake")

def abinit_build():
    with cd(build_path):
        conf_file = os.path.join(bzr_reposdir, "build.ac")
        run("../configure --with-config-file=%s" % conf_file)
        run("make clean")
        run("make -j8")

def abinit_makeall(self):
    """Runs abinit makemake and build."""
    abinit_makemake()
    abinit_build()

def abinit_runtests(opts=""):
    """
    Run the abinit automatic tests. Usage: abinit_runtests:"-j8 -k GW"
    """
    with cd(build_path):
        run("../tests/runtests.py %s" % opts)
