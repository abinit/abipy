# coding: utf-8
"""Tools to compute speed of sound."""
from __future__ import print_function, division, absolute_import  # unicode_literals,

import os
import math
import numpy as np
import pandas as pd
import abipy.core.abinit_units as abu

from abipy.core.mixins import Has_Structure, NotebookWriter
from abipy.dfpt.ddb import DdbFile
from abipy.dfpt.phonons import PhononBands, get_dyn_mat_eigenvec, match_eigenvectors
from abipy.abio.inputs import AnaddbInput
from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt
from pymatgen.core.units import bohr_to_angstrom, eV_to_Ha


class SoundVelocity(Has_Structure, NotebookWriter):
    """
    Compute the speed of sound by fitting phonon frequencies along selected directions
    by linear least-squares fit.
    """
    def __init__(self, directions, sound_velocities, mode_types, structure,
                 labels=None, phfreqs=None, qpts=None):
        """
        Args:
            directions: list of qpoints identifying the directions for the calculation
                of the speed of sound. In fractional coordinates.
            sound_velocities: array of shape (len(directions), 3) with the values of the
                sound velocities in SI units.
            mode_types: array of shape (len(directions), 3) with the type of modes (transverse,
                longitudinal). None if not known.
            structure: |Structure| object.
            labels: list with the same length as direction with labels for the directions.
            phfreqs: array with shape (len(directions), 3, num_points) with the acoustic phonon
                frequencies along the directions.
            qpts: array with shape (len(directions), num_points, 3) with the coordinates of
                the qpoints in fractional used to fit the phonon frequencies.
        """
        self.directions = directions
        self.sound_velocities = np.array(sound_velocities)
        self.mode_types = mode_types
        self.labels = labels
        self.phfreqs = np.array(phfreqs) if phfreqs else None
        self.qpts = np.array(qpts) if qpts else None
        self._structure = structure

    @property
    def structure(self):
        """|Structure| object"""
        return self._structure

    @property
    def n_directions(self):
        """Number of directions."""
        return len(self.directions)

    @classmethod
    def from_ddb(cls, ddb_path, directions=None, labels=None, num_points=20, qpt_norm=0.1,
                 ignore_neg_freqs=True, asr=2, chneut=1, dipdip=1, ngqpt=None, spell_check=True,
                 anaddb_kwargs=None, verbose=0, mpi_procs=1, workdir=None, manager=None):
        """
        Creates and instance of the object. Runs anaddb along the specified
        directions or the standard directions in the standard paths given
        in :cite:`Setyawan2010`. The values of the speed of sound
        will be calculated as the slope of the linear fits along each direction.

        Args:
            ddb_path (str): path to the ddb file.
            directions (list): list of qpoints identifying the directions for the calculation
                of the speed of sound. In fractional coordinates.
            labels (list): list of string with the name of the directions.
            num_points (int): number of points calculated along each direction and used to
                fit the speed of sound.
            qpt_norm (float): Norm of the largest point in fractional coordinates for
                each of the directions considered.
            ignore_neg_freqs (bool): if True points with negative frequencies will not be
                considered in the fit, in order to ignore inaccuracies in the long range
                behavior.
            asr, chneut, dipdip: Anaddb input variable. See official documentation.
            ngqpt: Number of divisions for the q-mesh in the DDB file. Auto-detected if None (default).
            anaddb_kwargs: additional kwargs for anaddb.
            verbose: verbosity level. Set it to a value > 0 to get more information.
            mpi_procs: Number of MPI processes to use.
            workdir: Working directory. If None, a temporary directory is created.
            manager: |TaskManager| object. If None, the object is initialized from the configuration file.
        long.

        Returns:
            an instance of SoundVelocity
        """
        with DdbFile(ddb_path) as ddb:
            if ngqpt is None: ngqpt = ddb.guessed_ngqpt

            inp = AnaddbInput(ddb.structure, comment="ANADDB input for speed of sound",
                              anaddb_kwargs=anaddb_kwargs, spell_check=spell_check)

            q1shft = [[0, 0, 0]]
            inp.set_vars(
                ifcflag=1,
                ngqpt=np.array(ngqpt),
                q1shft=q1shft,
                nqshft=len(q1shft),
                asr=asr,
                chneut=chneut,
                dipdip=dipdip,
            )

            if not directions:
                hs = ddb.structure.hsym_kpath
                kpath = hs.kpath

                directions = []
                labels = []

                for chunk in kpath["path"]:
                    for i, q in enumerate(chunk):
                        if "Gamma" in q:
                            if i > 0 and q not in labels:
                                new_q = kpath["kpoints"][chunk[i - 1]]
                                directions.append(new_q)
                                labels.append(chunk[i - 1])
                            if i < len(chunk) - 1 and q not in labels:
                                new_q = kpath["kpoints"][chunk[i + 1]]
                                directions.append(new_q)
                                labels.append(chunk[i + 1])

            qpts = []
            for q in directions:
                q = qpt_norm * q / np.linalg.norm(q)
                steps = q / num_points
                qpts.extend((steps[:, None] * np.arange(num_points)).T)

            n_qpoints = len(qpts)
            qph1l = np.zeros((n_qpoints, 4))

            qph1l[:, :-1] = qpts
            qph1l[:, -1] = 1

            inp['qph1l'] = qph1l.tolist()
            inp['nph1l'] = n_qpoints

            task = ddb._run_anaddb_task(inp, mpi_procs=mpi_procs, workdir=workdir, manager=manager,
                                        verbose=verbose)

            phbst_path = os.path.join(task.workdir, "run.abo_PHBST.nc")

            return cls.from_phbst(phbst_path, ignore_neg_freqs=ignore_neg_freqs, labels=labels)

    @classmethod
    def from_phbst(cls, phbst_path, ignore_neg_freqs=True, labels=None):
        """
        Creates an instance of the object starting interpolating the acoustic frequencies
        from a PHBST netcdf file.
        The file should contain a series of directions starting from gamma and with the
        same number of points for each direction, as the one produced in the from_ddb method.

        Args:
            phbst_path: path to the PHBST netcdf file.
            ignore_neg_freqs (bool): if True points with negative frequencies will not be
                considered in the fit, in order to ignore inaccuracies in the long range
                behavior.
            labels (list): list of string with the name of the directions.

        Returns:
            an instance of SoundVelocity
        """
        phb = PhononBands.from_file(phbst_path)
        structure = phb.structure

        rlatt = structure.lattice.reciprocal_lattice
        # q points in cartesian coordinate in 1/bohr, the original units are 1/A
        qpt_cart_coords = [rlatt.get_cartesian_coords(c) * bohr_to_angstrom for c in
                           phb.qpoints.frac_coords]
        qpt_cart_norms = np.linalg.norm(qpt_cart_coords, axis=1)

        # find the indices of the gamma points
        gamma_ind = []
        for i, q in enumerate(phb.qpoints.frac_coords):
            if np.array_equal(q, [0,0,0]):
                gamma_ind.append(i)

        n_directions = len(gamma_ind)

        n_points = len(phb.qpoints) / n_directions
        if not n_points.is_integer():
            raise ValueError('Error extracting information from {}'.format(phbst_path))
        n_points = int(n_points)

        phfreqs = phb.phfreqs
        eigdisp = phb.phdispl_cart

        sound_velocities = []
        mode_types = []
        directions = []

        all_acoustic_freqs = []
        all_qpts = []

        for i in range(n_directions):
            start = n_points * i
            # index of the end point used for the slice
            # (the position of the last point is actually end-1)
            end = n_points * (i + 1)
            dir_freqs = phfreqs[start: end]
            dir_displ = eigdisp[start: end]

            # matching bands
            dir_eigv = get_dyn_mat_eigenvec(dir_displ, structure, amu=phb.amu)
            n_freqs = 3 * len(structure)
            ind_match = np.zeros((n_points, n_freqs), dtype=np.int)
            ind_match[0] = range(n_freqs)

            for j in range(1, n_points):
                k = j - 1
                match = match_eigenvectors(dir_eigv[k], dir_eigv[j])
                ind_match[j] = [match[m] for m in ind_match[k]]

            acoustic_freqs = (dir_freqs[np.arange(n_points)[:, None], ind_match])[:, 0:3]
            acoustic_displ = (dir_displ[np.arange(n_points)[:, None], ind_match])[:, 0:3]

            direction = phb.qpoints[end - 1].frac_coords
            direction = np.array(direction) / np.linalg.norm(direction)
            directions.append(direction)

            # identify the first (not gamma) qpoint with all positive frequencies
            first_positive_freq_ind = None
            for j in range(1, n_points):
                if min(acoustic_freqs[j]) > 0:
                    first_positive_freq_ind = j
                    break

            if first_positive_freq_ind is None or first_positive_freq_ind - n_points / 2 > 0:
                raise ValueError("too many negative frequencies along direction {}".format(direction))

            sv = []
            mt = []

            cart_versor = qpt_cart_coords[end -1] / np.linalg.norm(qpt_cart_coords[end -1])
            for k in range(3):
                slope, se, _, _ = np.linalg.lstsq(qpt_cart_norms[start:end][:, np.newaxis],
                                                  acoustic_freqs[:, k] * eV_to_Ha, rcond=None)
                sv.append(slope[0] * abu.velocity_at_to_si)

                # identify the type of the mode (longitudinal/transversal) based on the
                # scalar product between the eigendisplacement and the direction.
                disp_0 = acoustic_displ[first_positive_freq_ind + 1, k, 0:3]
                disp_0 = disp_0 / np.linalg.norm(disp_0)

                scalar_prod = np.abs(np.dot(disp_0, cart_versor))
                if scalar_prod > 0.9:
                    mt.append("longitudinal")
                elif scalar_prod < 0.1:
                    mt.append("transversal")
                else:
                    mt.append(None)

            # sort the lists based on the sound velocites
            sv, mt, freqs = zip(*sorted(zip(sv, mt, acoustic_freqs.T)))

            sound_velocities.append(sv)
            mode_types.append(mt)
            all_acoustic_freqs.append(freqs)
            all_qpts.append(phb.qpoints.frac_coords[start:end])

        return cls(directions=directions, sound_velocities=sound_velocities, mode_types=mode_types,
                   structure=structure, labels=labels, phfreqs=all_acoustic_freqs, qpts=all_qpts)

    def get_dataframe(self):
        """
        Return a |pandas-DataFrame| with the data of the speed of sound.
        """
        columns = ["direction", "label", "sound velocity (m/s)", "mode type"]
        rows = []
        for i in range(self.n_directions):
            for m in range(3):
                rows.append([
                    tuple(np.round(self.directions[i], decimals=5)),
                    self.labels[i] if self.labels else "",
                    self.sound_velocities[i][m],
                    self.mode_types[i][m]
                ])

        return pd.DataFrame(rows, columns=columns).set_index(["direction", "label"])

    @add_fig_kwargs
    def plot_fit_freqs_dir(self, idir, ax=None, units="eV", **kwargs):
        """
        Plots the phonon frequencies, if available, along the specified direction.
        The line representing the fitted value will be shown as well.

        Args:
            idir: index of the direction.
            ax: |matplotlib-Axes| or None if a new figure should be created.
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz"). Case-insensitive.

        Returns:
            |matplotlib-Figure|
        """
        if self.phfreqs is None or self.qpts is None:
            raise ValueError("The plot requires the phonon frequencies and the qpoints.")

        ax, fig, plt = get_ax_fig_plt(ax=ax)

        ax.margins(x=0, y=0)
        rlatt = self.structure.lattice.reciprocal_lattice
        freqs = self.phfreqs[idir]
        qpt_cart_coords = np.array([np.linalg.norm(rlatt.get_cartesian_coords(c)) for c in self.qpts[idir]])
        slope = self.sound_velocities[idir] / abu.velocity_at_to_si * bohr_to_angstrom / eV_to_Ha

        units_factor = abu.phfactor_ev2units(units)

        title = "{:.4f}, {:.4f}, {:.4f}".format(*self.directions[idir])
        if self.labels:
            title += " - {}".format(self.labels[idir])

        for i, c in enumerate(["r", "b", "g"]):
            ax.scatter(qpt_cart_coords, freqs[i] * units_factor, color=c, marker="x")
            ax.plot(qpt_cart_coords, slope[i] * qpt_cart_coords * units_factor, color=c, ls="-")

        ax.set_title(title)
        ax.set_xlabel("Wave Vector")
        ax.set_ylabel(abu.wlabel_from_units(units))

        return fig

    @add_fig_kwargs
    def plot(self, units="eV", **kwargs):
        """
        Plots the phonon frequencies, if available, along all the directions.
        The lines representing the fitted values will be shown as well.

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.
            units: Units for phonon plots. Possible values in ("eV", "meV", "Ha", "cm-1", "Thz"). Case-insensitive.

        Returns:
            |matplotlib-Figure|
        """
        ax, fig, plt = get_ax_fig_plt(ax=None)
        from matplotlib.gridspec import GridSpec

        nrows, ncols = math.ceil(self.n_directions / 2),  2
        gspec = GridSpec(nrows=nrows, ncols=ncols, wspace=0.15, hspace=0.25)

        for i in range(self.n_directions):
            axi = plt.subplot(gspec[i])
            self.plot_fit_freqs_dir(i, ax=axi, units=units, show=False)

        return fig

    def yield_figs(self, **kwargs):   # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        """
        for i in range(self.n_directions):
            yield self.plot_fit_freqs_dir(i)

    def write_notebook(self, nbpath=None):
        """
        Write an jupyter_ notebook to nbpath. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        tmpfile = self.pickle_dump()

        nb.cells.extend([
            nbv.new_code_cell("sv = abilab.SoundVelocity.pickle_load('{}')".format(tmpfile)),
            nbv.new_code_cell("frame = sv.get_dataframe()\ndisplay(frame)")
        ])
        if self.phfreqs is not None and self.qpts is not None:
            for i in range(self.n_directions):
                nb.cells.append(nbv.new_code_cell("sv.plot_fit_freqs_dir({});".format(i)))

        return self._write_nb_nbpath(nb, nbpath)
