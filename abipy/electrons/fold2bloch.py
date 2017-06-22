#!/usr/bin/env python
# coding: utf-8
"""Python interface to Fold2Bloch."""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import numpy as np

from collections import OrderedDict
from pymatgen.core.lattice import Lattice
from abipy.tools.plotting import set_axlims, add_fig_kwargs, get_ax_fig_plt
from abipy.electrons.ebands import ElectronsReader


class Fold2Bloch(object):

    @classmethod
    def from_wfkpath(cls, wfkpath, mult, workdir=None, manager=None):
        # Run task in workdir
        # Usage: $fold2Bloch file_WFK x:y:z (folds)
        #manager = TaskManager.as_manager(manager).to_shell_manager(mpi_procs=1)

        # Build a simple manager to run the job in a shell subprocess
        #import tempfile
        #workdir = tempfile.mkdtemp() if workdir is None else workdir

        #mrgddb = wrappers.Mrgddb(manager=manager, verbose=0)
        #mrgddb.merge(self.outdir.path, ddb_files, out_ddb=out_ddb, description=desc)

        return cls.from_directory(cls, workdir=workdir)

    @classmethod
    def from_directory(cls, workdir="."):
        """Build object from a directory containing fold2bloch output files."""
        return cls(self, workdir)

    def __init__(self, workdir):
        self.workdir = os.path.abspath(workdir)

        tail = "_foldbloch.nc"
        files = [f for f os.listdir(self.workdir) if f.endswith(tail)]
        if not files:
            raise ValueError("Cannot find netcdf file ending with `%s` in workdir `%s`" % (tail, self.workdir))
        if len(files) > 1:
            raise ValueError("Found multiple netcdf files ending with `%s` in workdir `%s`" % (tail, self.workdir)
        seedname = files[:-len(tail)]

        # Get info from netcdf file.
        with ElectronsReader(os.path.join(self.workdir, files[0])) as r:
	    self.supercell = r.read_structure()
	    self.nsppol, self.nspinor = r.read_value("nsppol"), r.read_value("nspinor")
	    #self.nband = r.read_value("mband")
            #self.fermie
	    self.mult = self.read_value("mult")

        # Direct lattice of the primitive cell
        self.pc_lattice = Lattice((self.supercell.lattice.matrix.T * self.mult).T.copy())

        # Names of output files (depend on nsppol, nspinor)
        datafiles = []
	if self.nsppol == 1
            if self.nspinor == 1:
                datafiles.append(seedname + ".f2b")
            else:
                # Weights for spin up/down spinor component.
                datafiles.append(seedname + "_SPOR_1.f2b")
                datafiles.append(seedname + "_SPOR_2.f2b")
        elif self.nsppol == 2:
            # Weights for spin up/down
            datafiles.append(seedname + "_UP.f2b")
            datafiles.append(seedname + "_DOWN.f2b")

        datafiles = [os.path.join(self.workdir, p) for p in datafiles]

        # Now Read output files.
	# Columns 1-3 correspond to kx, ky and kz of the unfolded bands;
	# the 4th column is the energy eigenvalue in [Ha] and the 5th column corresponds
	# to a spectral weight of the k-point after unfolding.
        data = []
        for isp in range(len(datafiles)):
            od = OrderedDict()
            with open(datafiles[isp], "rt") as fh:
                for line in fh:
                    tokens = line.split()
                    kstr, eig, w = " ".join(tokens[:3]), float(tokens[3]), float(tokens[4])
                    if kstr not in od: od[kstr] = []
                    od[kstr].append((eig, w))
            data.append(od)

	# Build numpy arrays.
        for spin, od in enumerate(data):
            kpoints = np.reshape([tuple(map(float, t.split())) for k in od]), (-1, 3))
            if spin == 1:
                self.kpoints = kpoints
                self.nkpt = len(self.kpoints)
                self.eigens = np.empty((self.nsppol, self.nkpt, self.nband))
                self.weights = np.empty((self.nsppol * self.nspinor, self.nkpt, self.nband))
            else:
                if not np.all(self.kpoints == kpoints):
                    print("self.kpoints with shape:", self.kpoints.shape, "\nfrac_coords:\n", self.kpoints)
                    print("kpoints with shape:", kpoints.shape, "\nfrac_coords:\n", kpoints)
                    raise RuntimeError("Something wrong in the k-point parser")

            for ik, (kstr, tuple_list) in enumerate(od.items()):
                self.weights[spin, ik] = [t[1] for t in tuple_list]
                if spin == 2 and self.nspinor == 2: continue
                self.eigens[spin, ik] = [t[0] for t in tuple_list]

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose=0):
        """String representation"""
        lines = []
        app = lines.append
        app("nk_unfolded: %s, nsppol: %d, nspinor: %s, nband: %s" % (self.nkpt, self.nsppol, self.nspinor, self.nband))
        app("mult: %s" % (str(self.mult)))
        app("Supercell:")
        app(str(self.supercell))
        app("Direct lattice of the primitive cell:")
        app(str(self.pc_lattice))

        return "\n".join(lines)

    @add_fig_kwargs
    def plot(self, klabels=None, ylims=None, ax=None, **kwargs):
	"""
        Args:
            klabels: dictionary whose keys are tuple with the reduced coordinates of the k-points.
                The values are the labels. e.g. `klabels = {(0.0,0.0,0.0): "$\Gamma$", (0.5,0,0): "L"}`.
            ylims: Set the data limits for the y-axis. Accept tuple e.g. `(left, right)`
                   or scalar e.g. `left`. If left (right) is None, default values are used
            ax: matplotlib :class:`Axes` or None if a new figure should be created.

        Returns:
            `matplotlib` figure
	"""
        # Find k-points close to the input path.
        dist_tol = 0.001
        klines = []
        for i, k0 in enumerate(kbounds[:-1]):
            k1 = kbounds[i + 1]
            linds = self._find_points_close_to_line(k0, k1, dist_tol)
            if linds:
                klines.append(linds)
            else:
                print("Cannot find k-points close to line %s --> %s with dist_tol: %s" % (str(k0), str(k1), dist_tol))
        if not klines: return None

        # Compute ascissas so that proportions between segments is preserved.
        lmin = ls.min()
        xs = []
        for linds, l in zip(klines, ls):
            x0 = 0.0 if not xs else xs[-1]
            nn = len(linds)
            lps = x0 + (j / nn) * np.array([(l / lmin) for i in range(nn)])
            xs.extend(lps.flat)

        klines = klines.flatten()

        fact = 1.0
        e0 = self.fermie
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        for spin in range(self.nsppol):
            for band in range(self.nband):
                ys = self.eigens[spin, klines, band] - e0
                ws = self.weights[spin, klines, band] * fact
                ax.scatter(xs, ys, s=ws, marker="^", label="spin %s" % spin)

        #self.decorate_ax(ax, klabels=klabels)
        ax.grid(True)
        ax.set_ylabel('Energy [eV]')
        #set_axlims(ax, ylims, "y")
        ax.legend(loc="best")

        return fig


if __name__ == "__main__":
    import sys
    f = Fold2Bloch(sys.argv[1])
    print(f)
