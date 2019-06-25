# coding: utf-8
"""
Object to plot DFPT potentials in the phonon mode representation.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import numpy as np

from collections import OrderedDict
from monty.string import marquee
from monty.functools import lazy_property
from abipy.tools.plotting import add_fig_kwargs, get_axarray_fig_plt
from abipy.core.mixins import AbinitNcFile, Has_Structure, NotebookWriter
from abipy.core.kpoints import KpointList, Kpoint
from abipy.tools import duck
from abipy.iotools import Visualizer, xsf, ETSF_Reader #, cube


class V1qnuFile(AbinitNcFile, Has_Structure, NotebookWriter):

    def __init__(self, filepath):
        super(V1qnuFile, self).__init__(filepath)
        self.reader = r = ETSF_Reader(filepath)
        # Read dimensions.
        #self.nfft = r.read_dimvalue("nfft")
        #self.nspden = r.read_dimvalue("nspden")
        #self.natom3 = len(self.structure) * 3
        #self.symv1scf = r.read_value("symv1scf")
        #self.has_dielt_zeff = bool(r.read_value("has_dielt_zeff"))
        #self.add_lr_part = bool(r.read_value("add_lr_part"))

    @lazy_property
    def structure(self):
        """|Structure| object."""
        return self.reader.read_structure()

    def close(self):
        self.reader.close()

    @lazy_property
    def params(self):
        """:class:`OrderedDict` with parameters that might be subject to convergence studies."""
        return {}

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose=0):
        """String representation."""
        lines = []; app = lines.append
        app(marquee("File Info", mark="="))
        app(self.filestat(as_string=True))
        app("")
        app(self.structure.to_string(verbose=verbose, title="Structure"))
        app("")
        #app("has_dielt_zeff: %s" % self.has_dielt_zeff)
        #app("add_lr_part: %s" % self.add_lr_part)
        #app("symv1scf: %s" % self.symv1scf)
        for iq, qpt in enumerate(self.qpoints):
            app("[%d] %s" % (iq, repr(qpt)))

        return "\n".join(lines)

    @lazy_property
    def qpoints(self):
        return KpointList(self.structure.reciprocal_lattice, frac_coords=self.reader.read_value("qlist"))

    def _find_iqpt_qpoint(self, qpoint):
        if duck.is_intlike(qpoint):
            iq = qpoint
            qpoint = self.qpoints[iq]
        else:
            qpoint = Kpoint.as_kpoint(qpoint, self.structure.reciprocal_lattice)
            iq = self.qpoints.index(qpoint)

        return iq, qpoint

    def visualize_qpoint_nu(self, qpoint, nu, appname="vesta"):
        iq, qpoint = self._find_iqpt_qpoint(qpoint)

        # Fortran array nctkarr_t("v1qnu", "dp", "two, nfft, nspden, natom3, nqlist")])
        v1qnu = self.reader.read_variable("v1qnu")[iq]
        v1qnu = v1qnu[..., 0] + 1j * v1qnu[..., 1]  
        #fact = np.sqrt(2 * self.reader.read_variable["phfreqs"][nu]

        #visu = Visualizer.from_name(appname)
        #from abipy.core.globals import abinb_mkstemp
        #_, filename = abinb_mkstemp(suffix="." + ext, text=True)

        #with open(filename, mode="wt") as fh:
        #    if ext == "xsf":
        #        #xsf.xsf_write_structure(fh, self.structure)
        #        #xsf.xsf_write_data(fh, self.structure, datar, add_replicas=True)
        #    else:
        #        raise NotImplementedError("extension %s is not supported." % ext)

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        """
        maxnq = 3
        for iq, qpoint in enumerate(self.qpoints):
            if iq > maxnq:
                print("Only the first %d q-points are show..." % maxnq)
                break
            yield self.plot_diff_at_qpoint(qpoint=iq, **kwargs, show=False)

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("ncfile = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(ncfile)"),
        ])

        for iq, qpoint in enumerate(self.qpoints):
            nb.cells.append(nbv.new_code_cell("ncfile.plot_diff_at_qpoint(qpoint=%d);" % iq))

        return self._write_nb_nbpath(nb, nbpath)
