"""Post-processing tools for GWAN.nc files containing the e-ph vertex in real space."""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
from typing import TYPE_CHECKING

import numpy as np
from monty.string import marquee

import abipy.core.abinit_units as abu
from abipy.core.mixins import AbinitNcFile, Has_ElectronBands, Has_Header, Has_Structure
from abipy.eph.common import BaseEphReader
from abipy.tools.plotting import add_fig_kwargs, get_axarray_fig_plt, set_grid_legend

if TYPE_CHECKING:
    from abipy.core.structure import Structure
    from abipy.electrons.ebands import ElectronBands
    from abipy.tools.typing import Figure, PathLike


@dataclass(frozen=True)
class GwanSpin:
    """Real-space Hamiltonian and e-ph vertex for one spin channel."""

    spin: int
    r_h: np.ndarray
    r_e: np.ndarray
    r_p: np.ndarray
    ndegen_h: np.ndarray
    ndegen_e: np.ndarray
    ndegen_p: np.ndarray
    hwan_r: np.ndarray
    grpe_wwp: np.ndarray

    @property
    def nwan(self) -> int:
        """Number of Wannier functions."""
        return self.grpe_wwp.shape[1]

    @property
    def nr_e(self) -> int:
        """Number of electronic lattice vectors."""
        return len(self.r_e)

    @property
    def nr_p(self) -> int:
        """Number of phonon lattice vectors."""
        return len(self.r_p)

    def get_decay(self) -> tuple[np.ndarray, np.ndarray]:
        """
        Return the maximum vertex magnitude associated with each real-space vector.

        ``decay_e[ire]`` is ``max_{Rp,m,n,nu} |g(m,n,nu; Re,Rp)|`` and
        ``decay_p[irp]`` is ``max_{Re,m,n,nu} |g(m,n,nu; Re,Rp)|``.
        The values have the same units as ``grpe_wwp`` (Ha/Bohr in ABINIT).
        """
        # In Python: grpe_wwp(natom3, nwan, nwan, nr_e, nr_p).
        abs_g = np.abs(self.grpe_wwp)
        decay_e = np.max(abs_g, axis=(0, 1, 2, 4))
        decay_p = np.max(abs_g, axis=(0, 1, 2, 3))
        return decay_e, decay_p

    def get_hwan_decay(self) -> np.ndarray:
        """Return ``max_{m,n}|H_mn(Rh)|`` in Hartree for each Hamiltonian lattice vector."""
        # In Python: hwan_r(nwan, nwan, nr_h).
        return np.max(np.abs(self.hwan_r), axis=(0, 1))


class GwanFile(AbinitNcFile, Has_Header, Has_Structure, Has_ElectronBands):
    """File containing the e-ph vertex in the real-space Wannier representation."""

    @classmethod
    def from_file(cls, filepath: PathLike) -> GwanFile:
        """Initialize the object from a NetCDF file."""
        return cls(filepath)

    def __init__(self, filepath: PathLike):
        """Open ``filepath`` and initialize the GWAN reader."""
        super().__init__(filepath)
        self.r = GwanReader(filepath)

    @cached_property
    def ebands(self) -> ElectronBands:
        """Electronic bands stored in the root group."""
        return self.r.read_ebands()

    @property
    def structure(self) -> Structure:
        """Crystalline structure."""
        return self.ebands.structure

    @cached_property
    def params(self) -> dict:
        """Parameters used by generic convergence-analysis interfaces."""
        return {}

    @cached_property
    def gwan_spin(self) -> tuple[GwanSpin, ...]:
        """Real-space Wannier data for each spin channel."""
        return tuple(self.r.read_gwan_spin(spin) for spin in range(self.r.nsppol))

    def close(self) -> None:
        """Close the NetCDF reader."""
        self.r.close()

    def __str__(self) -> str:
        return self.to_string()

    def to_string(self, verbose: int = 0) -> str:
        """Return a string summary with verbosity level ``verbose``."""
        lines = [
            marquee("File Info", mark="="),
            self.filestat(as_string=True),
            "",
            self.structure.to_string(verbose=verbose, title="Structure"),
            "",
            self.ebands.to_string(with_structure=False, verbose=verbose, title="Electronic Bands"),
            "",
            f"Number of spin channels: {self.r.nsppol}",
        ]
        for data in self.gwan_spin:
            lines.append(
                f"Spin {data.spin}: nwan={data.nwan}, nr_h={len(data.r_h)}, "
                f"nr_e={data.nr_e}, nr_p={data.nr_p}"
            )
        if verbose > 1:
            lines.extend(["", self.hdr.to_string(verbose=verbose, title="Abinit Header")])
        return "\n".join(lines)

    @add_fig_kwargs
    def plot_decay(
        self,
        spin: int = 0,
        ax_mat=None,
        yscale: str = "log",
        marker: str = "o",
        markersize: float = 3,
        fontsize: int = 8,
        **kwargs,
    ) -> Figure:
        """
        Plot the real-space decay of the e-ph vertex for one spin channel.

        The left panel shows ``max_{Rp,m,n,nu}|g(Re,Rp)|`` versus ``|Re|``;
        the right panel shows ``max_{Re,m,n,nu}|g(Re,Rp)|`` versus ``|Rp|``.
        This is the same electronic decay measure written by ABINIT to
        ``*_spinN_GWAN.txt``, complemented by the phonon-lattice decay.

        Args:
            spin: Zero-based spin index.
            ax_mat: Two matplotlib axes or None to create them.
            yscale: Scale for the vertical axes, typically ``"log"`` or ``"linear"``.
            marker: Matplotlib marker for individual real-space vectors.
            markersize: Marker size.
            fontsize: Font size used for labels and legends.
            **kwargs: Additional options accepted by the plotting decorator.
        """
        if not 0 <= spin < self.r.nsppol:
            raise ValueError(f"Invalid {spin=}; expected 0 <= spin < {self.r.nsppol}")

        data = self.gwan_spin[spin]
        decay_e, decay_p = data.get_decay()
        lattice = self.structure.lattice.matrix
        rmod_e = np.linalg.norm(np.matmul(data.r_e, lattice), axis=1)
        rmod_p = np.linalg.norm(np.matmul(data.r_p, lattice), axis=1)

        ax_mat, fig, _plt = get_axarray_fig_plt(
            ax_mat, nrows=1, ncols=2, sharey=True, squeeze=False, rescale_fig=True
        )
        axes = ax_mat.ravel()
        for ax, radii, values, symbol, title in zip(
            axes,
            (rmod_e, rmod_p),
            (decay_e, decay_p),
            (r"R_e", r"R_p"),
            ("Electronic lattice vectors", "Phonon lattice vectors"),
            strict=True,
        ):
            order = np.argsort(radii)
            ax.plot(radii[order], values[order], linestyle="none", marker=marker, markersize=markersize)
            ax.set_yscale(yscale)
            ax.set_title(title)
            set_grid_legend(ax, fontsize, xlabel=rf"$|{symbol}|$ (Å)", ylabel=r"max $|g|$ (Ha/Bohr)")

        fig.suptitle(f"Real-space decay of the e-ph vertex, spin {spin}")
        return fig

    @add_fig_kwargs
    def plot_hwan_decay(
        self,
        spin: int = 0,
        ax=None,
        units: str = "eV",
        yscale: str = "log",
        marker: str = "o",
        markersize: float = 3,
        fontsize: int = 8,
        **kwargs,
    ) -> Figure:
        """
        Plot the real-space decay of the Hamiltonian in the Wannier representation.

        The plotted quantity is ``max_{m,n}|H_mn(Rh)|`` for each Hamiltonian
        lattice vector. Distances are reported in Å.

        Args:
            spin: Zero-based spin index.
            ax: Matplotlib axis or None to create one.
            units: Energy units, either ``"eV"`` or ``"Ha"``.
            yscale: Scale for the vertical axis, typically ``"log"`` or ``"linear"``.
            marker: Matplotlib marker for individual real-space vectors.
            markersize: Marker size.
            fontsize: Font size used for labels and legends.
            **kwargs: Additional options accepted by the plotting decorator.
        """
        if not 0 <= spin < self.r.nsppol:
            raise ValueError(f"Invalid {spin=}; expected 0 <= spin < {self.r.nsppol}")

        data = self.gwan_spin[spin]
        values = data.get_hwan_decay()
        if units == "eV":
            values = values * abu.Ha_eV
        elif units != "Ha":
            raise ValueError(f"Invalid {units=}; expected 'eV' or 'Ha'")

        rmod_h = np.linalg.norm(np.matmul(data.r_h, self.structure.lattice.matrix), axis=1)
        ax_mat, fig, _plt = get_axarray_fig_plt(ax, nrows=1, ncols=1)
        ax = np.asarray(ax_mat).ravel()[0]
        order = np.argsort(rmod_h)
        ax.plot(rmod_h[order], values[order], linestyle="none", marker=marker, markersize=markersize)
        ax.set_yscale(yscale)
        ax.set_title(f"Wannier Hamiltonian decay, spin {spin}")
        set_grid_legend(ax, fontsize, xlabel=r"$|R_h|$ (Å)", ylabel=rf"max $|H_{{mn}}|$ ({units})")
        return fig

    def yield_figs(self, **kwargs):  # pragma: no cover
        """Generate the standard GWAN figures."""
        for spin in range(self.r.nsppol):
            yield self.plot_decay(spin=spin, show=False)
            yield self.plot_hwan_decay(spin=spin, show=False)


class GwanReader(BaseEphReader):
    """NetCDF reader for the root group and per-spin GWAN groups."""

    def __init__(self, filepath: PathLike):
        """Open ``filepath`` and read the root-level dimensions."""
        super().__init__(filepath)
        self.nsppol = self.read_dimvalue("number_of_spins")

    @cached_property
    def path2group(self) -> dict:
        """Map NetCDF group names to group objects."""
        return self.rootgrp.groups

    def read_gwan_spin(self, spin: int) -> GwanSpin:
        """Read the real-space arrays for zero-based spin index ``spin``."""
        if not 0 <= spin < self.nsppol:
            raise ValueError(f"Invalid {spin=}; expected 0 <= spin < {self.nsppol}")

        path = f"gwan_spin{spin + 1}"
        return GwanSpin(
            spin=spin,
            r_h=self.read_value("r_h", path=path),
            r_e=self.read_value("r_e", path=path),
            r_p=self.read_value("r_p", path=path),
            ndegen_h=self.read_value("ndegen_h", path=path),
            ndegen_e=self.read_value("ndegen_e", path=path),
            ndegen_p=self.read_value("ndegen_p", path=path),
            hwan_r=self.read_value("hwan_r", path=path, cmode="c"),
            grpe_wwp=self.read_value("grpe_wwp", path=path, cmode="c"),
        )
