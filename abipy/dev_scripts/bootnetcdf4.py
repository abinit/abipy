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
