# coding: utf-8
"""DDB File."""
from __future__ import print_function, division, unicode_literals



import os
import numpy as np

from monty.collections import AttrDict
from monty.functools import lazy_property
from abipy.core.structure import Structure


class DdbFile(object):

    def __init__(self, filepath):
        self.filepath = os.path.abspath(filepath)

    def __enter__(self):
        return self.open()

    def open(self, mode="r"):
        self._file = open(self.filepath, mode=mode)
        return self

    def __iter__(self):
        return iter(self._file)

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Activated at the end of the with statement. It automatically closes the file."""
        self.close()

    def close(self):
        try:
            self._file.close()
        except:
            pass

    def seek(self, offset, whence=0):
        """
        seek(offset[, whence])
        Set the file's current position, like stdio's fseek(). 
        The whence argument is optional and defaults to 0 (absolute file positioning); other values
        are 1 (seek relative to the current position) and 2 (seek relative to the file's end). 
        There is no return value. Note that if the file is opened
        for appending (mode 'a' or 'a+'), any seek() operations will be undone at the next write. 
        If the file is only opened for writing in append mode (mode 'a'), this method is essentially a no-op, 
        but it remains useful for files opened in append mode with reading enabled (mode 'a+'). If the
        file is opened in text mode (without 'b'), only offsets returned by tell() are legal. 
        Use of other offsets causes undefined behavior.
        """
        self._file.seek(offset, whence)

    @lazy_property
    def structure(self):
        return Structure.from_abivars(**self.header)

    @lazy_property
    def header(self):
        """
        Dictionary with the values reported in the header section. 
        Use ddb.header.ecut to access its values
        """
        return self._parse_header()

    def _parse_header(self):
        """Parse the header sections. Returns :class:`AttrDict` dictionary."""
        #ixc         7
        #kpt  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00
        #     0.25000000000000D+00  0.00000000000000D+00  0.00000000000000D+00
        self.seek(0)
        keyvals, in_header = [], False
        for i, line in enumerate(self):
            line = line.strip()
            if not line: continue
            if "Version" in line:
                # +DDB, Version number    100401
                version = int(line.split()[-1])

            if line == "Description of the potentials (KB energies)":
                # Skip section with psps info.
                break
            if i == 6: in_header = True

            if in_header:
                # Python does not support exp format with D 
                line = line.replace("D+", "E+").replace("D-", "E-")
                tokens = line.split()
                key = None
                try:
                    float(tokens[0])
                    parse = float if "." in tokens[0] else int
                    keyvals[-1][1].extend(map(parse, tokens))
                except ValueError:
                    # We have a new key
                    key = tokens.pop(0)
                    parse = float if "." in tokens[0] else int
                    keyvals.append((key, map(parse, tokens)))

        h = AttrDict(version=version)
        for key, value in keyvals:
            if len(value) == 1: value = value[0]
            h[key] = value

        # Convert to array. Note that znucl is converted into nucl 
        # to avoid problems with pymatgen routines that expect integral Z
        # This of course will break any calculation that uses alchemical mixing
        arrays = {
            "kpt": dict(shape=(h.nkpt, 3), dtype=np.double),
            "rprim": dict(shape=(3, 3), dtype=np.double),
            "symrel": dict(shape=(h.nsym, 3, 3), dtype=np.int),
            "tnons": dict(shape=(h.nsym, 3), dtype=np.double),
            "xred":  dict(shape=(h.natom, 3), dtype=np.double),
            "znucl": dict(shape=(-1,), dtype=np.int),
        }

        for k, ainfo in arrays.items():
            h[k] = np.reshape(np.array(h[k], dtype=ainfo["dtype"]), ainfo["shape"])

        return h

    @lazy_property
    def qpoints(self):
        """`ndarray` with the list of q-points in reduced coordinates."""
        return self._read_qpoints()

    def _read_qpoints(self):
        """Read the list q-points from the DDB file. Returns `ndarray`"""
        # 2nd derivatives (non-stat.)  - # elements :      36
        # qpt  2.50000000E-01  0.00000000E+00  0.00000000E+00   1.0

        # Since there are multiple occurrences of qpt in the DDB file
        # we use seen to remove duplicates.
        self.seek(0)
        tokens, seen = [], set()

        for line in self:
            line = line.strip()
            if line.startswith("qpt") and line not in seen:
                seen.add(line)
                tokens.append(line.replace("qpt", ""))

        qpoints, weights = [], []
        for tok in tokens:
            nums = map(float, tok.split())
            qpoints.append(nums[:3])
            weights.append(nums[3])

        return np.reshape(qpoints, (-1,3))

    #def calc_phbands_and_dos(self):
    #def calc_thermo(self):

#from pymatgen.io.abinitio.tasks import AnaddbTask, TaskManager
#task = AnaddbTask(Mock(), "out_DDB", workdir="test", manager=TaskManager.from_user_config())
#task.start_and_wait(autoparal=False)

class Mock(object):
    anaddb_input = """
!Flags
ifcflag   1     ! Interatomic force constant flag

!Wavevector grid number 1 (coarse grid, from DDB)
ngqpt   1  1  1   ! Monkhorst-Pack indices
nqshft  1         ! number of q-points in repeated basic q-cell
q1shft  0.25 0 0 

!Effective charges
asr   1     ! Acoustic Sum Rule. 1 => imposed asymetrically
chneut   1  ! Charge neutrality requirement for effective charges.

!Interatomic force constant info
dipdip  1   ! Dipole-dipole interaction treatment

nqpath 2
qpath 0.25 0 0    
       0.0 0 0 
ndivsmall 1
"""
    def add_extra_abivars(self, extra):
        pass

    def make_input(self):
        return self.anaddb_input
