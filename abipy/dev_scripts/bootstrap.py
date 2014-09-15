#!/usr/bin/env python
"""
Taken from http://eli.thegreenplace.net/2013/04/20/bootstrapping-virtualenv/
See also http://stackoverflow.com/questions/4324558/whats-the-proper-way-to-install-pip-virtualenv-and-distribute-for-python
"""
from __future__ import print_function, division, unicode_literals
import sys
import subprocess

VENV_VERSION = '1.9.1'
#VENV_VERSION = '1.10.1'
PYPI_VENV_BASE = 'https://pypi.python.org/packages/source/v/virtualenv'

PYTHON = 'python'
#PYTHON = 'python2.7'
INITIAL_ENV = 'py-env0'

ENV_OPTS = []
#ENV_OPTS = ["--no-site-packages",]
ENV_OPTS = " ".join(ENV_OPTS)


def shellcmd(cmd, echo=True):
    """
    Run 'cmd' in the shell and return its standard out.
    """
    if echo: print('[cmd] {0}'.format(cmd))
    out = subprocess.check_output(cmd, stderr=sys.stderr, shell=True)
    if echo: print(out)
    return out

def main():

    dirname = 'virtualenv-' + VENV_VERSION
    tgz_file = dirname + '.tar.gz'

    # Fetch virtualenv from PyPI
    venv_url = PYPI_VENV_BASE + '/' + tgz_file
    shellcmd('curl -O {0}'.format(venv_url))

    # Untar
    shellcmd('tar xzf {0}'.format(tgz_file))

    # Create the initial env
    shellcmd('{0} {1}/virtualenv.py {2} {3}'.format(PYTHON, dirname, ENV_OPTS, INITIAL_ENV))

    # Install the virtualenv package itself into the initial env
    shellcmd('{0}/bin/pip install {1}'.format(INITIAL_ENV, tgz_file))

    # Cleanup
    shellcmd('rm -rf {0} {1}'.format(dirname, tgz_file))

    return 0

if __name__ == "__main__":
    sys.exit(main())
