#!/usr/bin/env python
from __future__ import print_function

import sys
import subprocess

def shellcmd(cmd, echo=True):
    """
    Run 'cmd' in the shell and return its standard out.
    """
    if echo: print('[cmd] {0}'.format(cmd))
    out = subprocess.check_output(cmd, stderr=sys.stderr, shell=True)
    if echo: print(out)
    return out


# Netcdf4+HDF5 tarballs.
nc4_tarballs = dict(
    hdf5="ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4/hdf5-1.8.11.tar.gz",
    netcdf4="ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4.3.0.tar.gz",
    netcdf4_f90="ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-fortran-4.2.tar.gz",
)

# Useful packages for development.
dev_tarballs = dict(
    ctags="http://prdownloads.sourceforge.net/ctags/ctags-5.8.tar.gz",
    python2.7="http://www.python.org/ftp/python/2.7.6/Python-2.7.6.tar.xz",
)


#matplotlib = dict(
#"http://download.savannah.gnu.org/releases/freetype/freetype-2.5.1.tar.gz",
#)

# BLAS + LAPACK (if the host does not provide them).
#linalg_tarballs = dict(
#)

# PyPy software stack. 
pypy_tarballs = dict(
    pypy="https://bitbucket.org/pypy/pypy/downloads/pypy-2.2-src.tar.bz2",
)


def download_tarballs():
    hdf5_url = "ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4/hdf5-1.8.11.tar.gz"
    nc4_url = "ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4.3.0.tar.gz"
    ncf4_url = "ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-fortran-4.2.tar.gz"

    shellcmd('curl -O %s' % hdf5_url)
    shellcmd('curl -O %s' % nc4_url)
    #shellcmd('curl -O %s' % ncf4_url)

    return 0

def main():
    test = False
    download = True

    if test:
        try:
            import netCDF4
            return 0
        except ImportError:
            return 1


    if download:
        download_tarballs()

    return 0

if __name__ == '__main__':
    sys.exit(main())
