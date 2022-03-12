# coding: utf-8
"""
Object to analyze the results stored in the WR.nc file
"""
import numpy as np

from monty.string import marquee
from monty.functools import lazy_property
from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt #, get_axarray_fig_plt
from abipy.core.mixins import AbinitNcFile, Has_Structure, NotebookWriter
from abipy.iotools import ETSF_Reader


class WrNcFile(AbinitNcFile, Has_Structure, NotebookWriter):

    def __init__(self, filepath):
        super().__init__(filepath)
        self.reader = r = ETSF_Reader(filepath)

        # Read dimensions.
        self.nfft = r.read_dimvalue("nfft")
        self.nspden = r.read_dimvalue("nspden")
        self.natom3 = len(self.structure) * 3
        self.method = r.read_value("method")
        assert self.method == 0
        self.ngqpt = r.read_value("ngqpt")
        self.rpt = r.read_value("rpt")
        self.nrpt = len(self.rpt)
        # FFT mesh.
        self.ngfft = r.read_value("ngfft")

    def create_xsf(self, iatom=0, red_dir=(1, 0, 0), u=1.0, ispden=0):

        nfft, nrpt = self.nfft, self.nrpt

        nx, ny, nz = self.ngfft
        nqx, nqy, nqz = self.ngqpt
        box_shape = self.ngqpt * self.ngfft
        box_size = np.product(box_shape)
        print("ngqpt:", self.ngqpt)
        print("nrpt:", self.nrpt)
        print("Unit cell FFT shape:", self.ngfft)
        print("Big box shape:", box_shape)
        print("rpt:\n", self.rpt)

        # Get FFT points in reduced coordinates of the microcell.
        # ix is the fastest index here because we are gonna access
        # FFT values produced by Fortran via ifft
        fft_inds = np.empty((nfft, 3), dtype=np.int64)
        ifft = -1
        for iz in range(nz):
            for iy in range(ny):
                for ix in range(nx):
                    ifft += 1
                    fft_inds[ifft, :] = [ix, iy, iz]

        def ig2gfft(ig, ng):
            # Use the following indexing (N means ngfft of the adequate direction)
            # 0 1 2 3 ... N/2    -(N-1)/2 ... -1    <= gc
            # 1 2 3 4 ....N/2+1  N/2+2    ...  N    <= index ig

            #if ( ig <= 0 or ig > ng):
            #  # Wrong ig, returns huge. Parent code will likely crash with SIGSEV.
            #  gc = huge(1)
            #  return

            #if (ig  > ng/2 + 1):
            #  gc = ig - ng -1
            #else
            #  gc = ig -1
            #return gc
            raise NotImplementedError("Foo")

        # Find index of (r, R) in the bigbox
        # Use C notation (z runs faster)
        #box2fr = np.full((box_size, 2), np.inf, dtype=np.int64)

        #iffr2box = np.empty((nfft, nrpt), dtype=np.int64)
        print("Building r-R dictionary")
        d = {}
        for ir, rpt in enumerate(self.rpt):
            for ifft, fft_ijk in enumerate(fft_inds):
                ijk = (fft_ijk - rpt * self.ngfft) % box_shape
                key = tuple(map(int, ijk))
                d[key] = (ifft, ir)
                #i, j, k = ijk
                #box_iloc = k + j * box_shape[2] + i * (box_shape[2] * box_shape[1])
                #print(box_iloc, "(i, j, k): ", fft_ijk, "rpt:", rpt)
                #box2fr[int(box_iloc)] = [ifft, ir]
                #iffr2box[ifft, ir] = box_iloc
        print("Done")

        #ip = 0
        #idir = ip % 3
        #iatom = (ip - idir) // 3 # + 1

        # nctkarr_t("v1scf_rpt_sr", "dp", "two, nrpt, nfft, nspden, natom3")
        # use iorder = "f" to transpose the last 3 dimensions since ETSF
        # stores data in Fortran order while AbiPy uses C-ordering.
        # (z,y,x) --> (x,y,z)
        # datar = transpose_last3dims(datar)
        wsr_var = self.reader.read_variable("v1scf_rpt_sr")
        wlr_var = self.reader.read_variable("v1scf_rpt_lr")

        wsr = np.zeros(nfft, nrpt)
        wlr = np.zeros(nfft, nrpt)
        for idir, red_comp in enumerate(red_dir):
            ip = idir + 3 * iatom
            wsr += u * red_comp * wsr_var[ip, ispden, :, :, 0]
            wlr += u * red_comp * wlr_var[ip, ispden, :, :, 0]

        print("wsr.shape:", wsr.shape)
        print("Max |Re Wsr|:", np.max(np.abs(wsr.real)), "Max |Im Wsr|:", np.max(np.abs(wsr.imag)))
        #print("Max |Re Wlr|:", np.max(np.abs(wlr.real)), "Max |Im Wlr|:", np.max(np.abs(wlr.imag)))

        r0 = np.array([0, 0, 0], dtype=np.int)
        qgrid = np.where(self.ngqpt > 2, self.ngqpt, 0)
        r0 = - (self.ngqpt - 1) // 2
        print("Origin of datagrid set at R0:", r0)

        # Build datagrid in the supercell using C indexing
        # This is what xsf_write_data expects.
        miss = []
        data_lr = np.empty(box_shape)
        data_sr = np.empty(box_shape)

        print("Filling data array")
        for ix in range(box_shape[0]):
            for iy in range(box_shape[1]):
                for iz in range(box_shape[2]):
                    y = np.array((ix, iy, iz), dtype=np.int)
                    x = (y + r0) % box_shape
                    key = tuple(map(int, x))
                    try:
                        ifft, ir = d[key]
                    except KeyError:
                        #print("Cannot find r - R with key:", key)
                        miss.append(key)
                        continue

                    data_lr[ix, iy, iz] = wlr[ifft, ir]
                    #data_sr[ix, iy, iz] = wsr[ifft, ir]

        if miss:
            #print(d.keys())
            raise RuntimeError("Cannot find r-R points! nmiss:", len(miss))

        super_structure = self.structure * self.ngqpt

        def dump_xsf(filename, data):
            from abipy.iotools import xsf
            xsf.xsf_write_structure_and_data_to_path(filename, super_structure, data, cplx_mode="abs")

        dump_xsf("foo_lr.xsf", data_lr)
        dump_xsf("foo_sr.xsf", data_sr)

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

        return "\n".join(lines)

    @add_fig_kwargs
    def plot_maxw(self, scale="semilogy", ax=None, fontsize=8, **kwargs):
        """
        Plot the decay of max_{r,idir,ipert} `|W(R,r,idir,ipert)|`
        for the long-range and the short-range part.

        Args:
            scale: "semilogy", "loglog" or "plot".
            ax: |matplotlib-Axes| or None if a new figure should be created.
            fontsize: fontsize for legends and titles

        Return: |matplotlib-Figure|
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        f = {"plot": ax.plot, "semilogy": ax.semilogy, "loglog": ax.loglog}[scale]

        rmod = self.reader.read_value("rmod")

        # Plot short-range part.
        # normalize wrt the R=0 value
        # Fortran array: nctkarr_t("maxw_sr", "dp", "nrpt, natom3")
        maxw_sr = self.reader.read_value("maxw_sr")
        data = np.max(maxw_sr, axis=0)
        #data = data / data[0]
        f(rmod, data, marker="o", ls=":", lw=0, label="SR", **kwargs)

        # Plot long-range part.
        maxw_lr = self.reader.read_value("maxw_lr")
        data = np.max(maxw_lr, axis=0)
        #data = data / data[0]
        f(rmod, data, marker="x", ls="-", lw=0, label="LR", **kwargs)

        # Plot the ratio
        data = np.max(maxw_lr, axis=0) / np.max(maxw_sr, axis=0)
        #f(rmod, data, marker="x", ls="-", lw=0, label="LR/SR", **kwargs)

        #rmod = self.reader.read_value("rmod_lrmodel")
        #maxw = self.reader.read_value("maxw_lrmodel")
        #data = np.max(maxw, axis=0)
        #data = data / data[0]
        #f(rmod, data, marker="x", ls="-", lw=0, label="W_LR_only", **kwargs)

        ax.grid(True)
        ax.set_ylabel(r"$Max_{({\bf{r}}, idir, ipert)} \| W({\bf{r}}, {\bf{R}}, idir, ipert) \|$")
        ax.set_xlabel(r"$\|{\bf{R}}\|$ (Bohr)")
        ax.legend(loc="best", fontsize=fontsize, shadow=True)

        #if kwargs.pop("with_title", True):
        #    ax.set_title("dvdb_add_lr %d, qdamp: %s, symv1scf: %d" % (self.dvdb_add_lr, self.qdamp, self.symv1scf),
        #                 fontsize=fontsize)
        return fig

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        """
        yield self.plot_maxw(scale="semilogy")

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

        #nb.cells.append(nbv.new_code_cell("ncfile.plot_diff_at_qpoint(qpoint=%d);" % iq))

        return self._write_nb_nbpath(nb, nbpath)


if __name__ == "__main__":
    import sys
    ncfile = WrNcFile.from_file(sys.argv[1])

    #print(ncfile)
    ncfile.plot_maxw(scale="semilogy", ax=None, fontsize=8)
    #ncfile.create_xsf(iatom=0, red_dict=(-1, +1, +1), u=0.1)
