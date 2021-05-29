import numpy as np
import abipy.core.abinit_units as abu
from abipy.iotools import ETSF_Reader
from abipy.core.func1d import Function1D
from abipy.tools.plotting import add_fig_kwargs, get_ax_fig_plt
from collections import namedtuple


PowderIntensity = namedtuple("PowderIntensity", ("paral", "perp", "tot"))


class Raman:
    """
    An object that allows to obtain the Raman intensities based on the calculated
    susceptibilities. The expressions used to extract the different values of the
    Raman intensities can be found for example in PRB 71, 214307, Physics of the
    Earth and Planetary Interiors 174, 113 or PRB 71, 125107.
    """

    def __init__(self, susceptibility, phfreqs, non_anal_susceptibility=None,
                 non_anal_phfreqs=None, non_anal_directions=None):
        """
        Args:
            susceptibility: a numpy array with shape (n modes, 3 3) containing the
                Raman susceptibilities of transverse zone-center phonon modes.
            phfreqs: a numpy array with the phonon frequencies at Gamma, without non
                analytical contributions.
            non_anal_susceptibility: a numpy array with shape (n directions, n modes, 3 3)
                containing the Raman susceptibilities of zone-center phonons, with non-analyticity
                in  different directions.
            non_anal_phfreqs: a numpy array with the phonon frequencies with shape
                (n directions, n modes), with non analytical contributions along the different
                directions.
            non_anal_directions: a numpy array with shape (n directions, 3) with the directions
                along which the non analytical contribution has been calculated.
        """
        self.susceptibility = susceptibility
        self.phfreqs = phfreqs
        self.non_anal_susceptibility = non_anal_susceptibility
        self.non_anal_phfreqs = non_anal_phfreqs
        self.non_anal_directions = non_anal_directions

    @classmethod
    def from_file(cls, filepath):
        """
        Create the object from an anaddb.nc netcdf file.

        Args:
            filepath: path to the netcdf file.

        Returns:
            An instance of Raman.
        """

        with ETSF_Reader(filepath) as r:
            try:
                susceptibility = r.read_value("raman_sus").T
                phfreqs = r.read_value("gamma_phonon_modes")
                non_anal_susceptibility = r.read_value("non_analytical_raman_sus", default=None)
                if non_anal_susceptibility is not None:
                    non_anal_susceptibility = non_anal_susceptibility.T
                non_anal_phfreqs = r.read_value("non_analytical_phonon_modes", default=None)
                non_anal_directions = r.read_value("non_analytical_directions", default=None)
            except Exception:
                import traceback
                msg = traceback.format_exc()
                msg += ("Error while trying to read Raman from file.\n"
                        "Verify that the required variables are used in anaddb: nlflag\n")
                raise ValueError(msg)

            return cls(susceptibility=susceptibility, phfreqs=phfreqs, non_anal_susceptibility=non_anal_susceptibility,
                       non_anal_phfreqs=non_anal_phfreqs, non_anal_directions=non_anal_directions)

    def get_modes_intensities(self, temp, laser_freq, non_anal_dir=None, relative=False, units="eV",
                              pol_in=None, pol_out=None):
        """
        Calculates the Raman intensities for each mode in arbitrary units. It is possible to use the
        susceptibilities from the transverse modes only or to specify one of the directions with non
        analytical contributions. By default it returns an array with shape (n modes, 3, 3) where each
        component of the 3x3 matrix are those with different incoming (the first index) and outgoing
        (the second index) polarization. These are along the components along the standard axes (e.g.
        "xx" or "yz" components). If pol_in and pol_out an array with shape (n modes) is returned
        corresponding to the selected polarizations.

        Args:
            temp: temperature in K.
            laser_freq: frequency of the incident laser. The units are determined the "units"
                argument.
            non_anal_dir: index of the direction along which the non analytical contribution
                has been calculated. Corresponds to the indices in the non_anal_directions attribute.
            relative: if True the intensities will be rescaled so that the largest value is 1.
            units: the units in which the input and the output frequencies will be given.
                Possible values in ("eV", "meV", "Ha", "cm-1", "Thz")
            pol_in: the polarization of the incoming photon. If not None can be either either a string
                with one of the cartesian components i.e. "x", "y", "z" or an array with 3 elements
                representing the polarization direction. If not None pol_out can not be None.
            pol_out: the polarization of the outgoing photon. If not None can be either either a string
                with one of the cartesian components i.e. "x", "y", "z" or an array with 3 elements
                representing the polarization direction. If not None pol_in can not be None.

        Returns:
            An array with the Raman intensities. If pol_in==pol_out==None has shape (n modes, 3, 3)
            with all the components. Otherwise an array with size (n modes) with the intensities of
            the selected polarizations.
        """

        if non_anal_dir is None:
            w = self.phfreqs
            sus = self.susceptibility
        else:
            w = self.non_anal_phfreqs[non_anal_dir]
            sus = self.non_anal_susceptibility[non_anal_dir]

        laser_freq = laser_freq / abu.phfactor_ev2units(units)

        c = self._get_prefactor(w=w, temp=temp, laser_freq=laser_freq)
        if pol_in is None and pol_out is None:
            i = c[:, np.newaxis, np.newaxis] * sus**2
            # this will make the indices of the i,j component such that the first
            # will refer to the polarization of the incoming photon and the second
            # to the polarization of the created one.
            np.transpose(i, axes=(0, 2, 1))
        else:
            if pol_in is None or pol_out is None:
                raise ValueError("pol_in and pol_out should be either both None or both defined")
            dxyz = {"x": [1, 0, 0], "y": [0, 1, 0], "z": [0, 0, 1]}
            if isinstance(pol_in, str):
                pol_in = dxyz[pol_in.lower()]
            if isinstance(pol_out, str):
                pol_out = dxyz[pol_out.lower()]

            pol_in = np.array(pol_in) / np.linalg.norm(pol_in)
            pol_out = np.array(pol_out) / np.linalg.norm(pol_out)

            i = c * np.einsum("ijk, j, k -> i", sus, pol_out, pol_in) ** 2

        if relative:
            i /= i.max()

        return i

    @staticmethod
    def _get_prefactor(w, temp, laser_freq):
        """
        Helper method to calculate the coefficient for the Raman intensities.

        Args:
            w: the selected frequencies in eV.
            temp: the temperature in K.
            laser_freq: the frequency of the laser in eV

        Returns:
            An array with shape (n modes) with the coefficient for the Raman intensities.
        """

        c = np.zeros_like(w)
        ind = np.where(w > 1e-5)

        bose_factor = 1 / (1 - np.exp(-w[ind] / (abu.kb_eVK * temp)))

        c[ind] = (w[ind] - laser_freq) ** 4 / (2 * w[ind]) * bose_factor

        return c

    def _get_lorentz_freqs_and_factor(self, intensity, non_anal_dir, min_freq, max_freq, num, width, units):
        """
        Helper method to get the list of frequencies and the main spread factors to
        calculate the broadened Raman intensities with a Lorentz distribution.

        Args:
            intensity: the Raman intensities at specified frequencies.
            non_anal_dir: ndex of the direction along which the non analytical contribution
                has been calculated. Corresponds to the indices in the non_anal_directions attribute.
            min_freq: minimum frequency considered. If None it will be given by the minimum
                frequency with non zero intensities minus 10 times the width of the distribution.
            max_freq: maximum frequency considered. If None it will be given by the maximum
                frequency with non zero intensities plus 10 times the width of the distribution.
            width: the width of the Lorentz distribution.
            units: the units in which the input and the output frequencies will be given.
                Possible values in ("eV", "meV", "Ha", "cm-1", "Thz")

        Returns:
            Tuple with list of "num" frequencies in eV and factors for the Lorentz broadening
            with shape (n modes, num).
        """

        if non_anal_dir is None:
            w = self.phfreqs
        else:
            w = self.non_anal_phfreqs[non_anal_dir]

        units_factor = abu.phfactor_ev2units(units)

        width = width / units_factor

        if min_freq is None:
            min_ind = np.where(intensity/intensity.max() > 1e-10)[0].min()
            min_freq = w[min_ind] - 10 * width
        else:
            min_freq = min_freq / units_factor

        if max_freq is None:
            max_ind = np.where(intensity/intensity.max() > 1e-10)[0].max()
            max_freq = w[max_ind] + 10 * width
        else:
            max_freq = max_freq / units_factor

        freqs = np.linspace(min_freq, max_freq, num)

        lorentz = width / ((freqs - w.reshape((-1, 1)))**2 + width**2) / np.pi

        return freqs, lorentz

    def get_lorentz_intensity(self, temp, laser_freq, width, non_anal_dir=None, min_freq=None, max_freq=None,
                              num=1000, relative=False, units="eV", pol_in=None, pol_out=None):
        """
        Calculates the broadened Raman intensities in arbitrary units for frequencies in an interval. It is
        possible to use the susceptibilities from the transverse modes only or to specify one of the directions
        with non analytical contributions. By default it returns a 3x3 matrix where each component is a
        Function1D object with the Raman intensities with different incoming (the first index) and outgoing
        (the second index) polarization. These are along the components along the standard axes (e.g. "xx" or
        "yz" components). If pol_in and pol_out a single Function1D is returned corresponding to the selected
        polarizations.

        Args:
            temp: temperature in K.
            laser_freq: frequency of the incident laser. The units are determined the "units"
                argument.
            width: the width of the Lorentz distribution. The units are determined the "units"
                argument.
            non_anal_dir: index of the direction along which the non analytical contribution
                has been calculated. Corresponds to the indices in the non_anal_directions attribute.
            min_freq: minimum frequency considered. If None it will be given by the minimum
                frequency with non zero intensities minus 10 times the width of the distribution.
                If given the units are determined by the "units" argument.
            max_freq: maximum frequency considered. If None it will be given by the maximum
                frequency with non zero intensities plus 10 times the width of the distribution.
                If given the units are determined by the "units" argument.
            num: number of frequencies in the interval (min_freq, max_freq).
            relative: if True the intensities will be rescaled so that the largest value is 1.
            units: the units in which the input and the output frequencies will be given.
                Possible values in ("eV", "meV", "Ha", "cm-1", "Thz")
            pol_in: the polarization of the incoming photon. If not None can be either either a string
                with one of the cartesian components i.e. "x", "y", "z" or an array with 3 elements
                representing the polarization direction. If not None pol_out can not be None.
            pol_out: the polarization of the outgoing photon. If not None can be either either a string
                with one of the cartesian components i.e. "x", "y", "z" or an array with 3 elements
                representing the polarization direction. If not None pol_in can not be None.

        Returns:
            If pol_in==pol_out==None a 3x3 list with a Function1D corresponding to the different
            components of the intensities. Otherwise a single Function1D with the  intensities of
            the selected polarizations. Each Function1D has "num" points.
        """

        i = self.get_modes_intensities(temp=temp, laser_freq=laser_freq, non_anal_dir=non_anal_dir,
                                       units=units, pol_in=pol_in, pol_out=pol_out)

        freqs, lorentz = self._get_lorentz_freqs_and_factor(intensity=i, non_anal_dir=non_anal_dir, min_freq=min_freq,
                                                            max_freq=max_freq, num=num, width=width, units=units)

        # convert the frequencies to the desired units for the output
        x = freqs * abu.phfactor_ev2units(units)

        if pol_in is not None and pol_out is not None:
            li = np.dot(i, lorentz)
            if relative:
                li /= li.max()

            return Function1D(x, li)

        else:
            li = np.einsum("ij, ikl -> jkl", lorentz, i)

            li_func = [[None]*3]*3

            for i in range(3):
                for j in range(3):
                    y = li[:, i, j]
                    if relative:
                        y /= y.max()
                    li_func[i][j] = Function1D(x, y)

            return li_func

    def get_powder_intensity(self, temp, laser_freq, non_anal_dir=None, relative=False, units="eV"):
        """
        Calculates the Raman intensities in arbitrary units for each mode integrated over all possible
        orientation to reproduce the powder measurements. It is possible to use the susceptibilities from
        transverse modes only or to specify one of the directions with non analytical contributions.

        Args:
            temp: temperature in K.
            laser_freq: frequency of the incident laser. The units are determined the "units"
                argument.
            non_anal_dir: index of the direction along which the non analytical contribution
                has been calculated. Corresponds to the indices in the non_anal_directions attribute.
            relative: if True the intensities will be rescaled so that the largest value of the
                total intensity is 1.
            units: the units in which the input and the output frequencies will be given.
                Possible values in ("eV", "meV", "Ha", "cm-1", "Thz")

        Returns:
            A PowderIntensity with the parallel, perpendicular and total components of the powder
            intensities. Each one is an array with length n modes.
        """

        if non_anal_dir is None:
            w = self.phfreqs
            sus = self.susceptibility
        else:
            w = self.non_anal_phfreqs[non_anal_dir]
            sus = self.non_anal_susceptibility[non_anal_dir]

        g0 = np.trace(sus, axis1=1, axis2=2)**2 / 3
        g1 = ((sus[:, 0, 1] - sus[:, 1, 0])**2 + (sus[:, 0, 2] - sus[:, 2, 0])**2 + (sus[:, 2, 1] - sus[:, 1, 2])**2) / 2
        g2 = ((sus[:, 0, 1] + sus[:, 1, 0])**2 + (sus[:, 0, 2] + sus[:, 2, 0])**2 + (sus[:, 2, 1] + sus[:, 1, 2])**2) / 2 + \
             ((sus[:, 0, 0] - sus[:, 1, 1])**2 + (sus[:, 0, 0] - sus[:, 2, 2])**2 + (sus[:, 1, 1] - sus[:, 2, 2])**2) / 3

        laser_freq = laser_freq / abu.phfactor_ev2units(units)

        c = self._get_prefactor(w=w, temp=temp, laser_freq=laser_freq)

        paral = c * (10 * g0 + 4 * g2)
        perp = c * (5 * g1 + 3 * g2)
        tot = paral + perp
        if relative:
            m = tot.max()
            paral /= m
            perp /= m
            tot /= m

        return PowderIntensity(paral, perp, tot)

    def get_powder_lorentz_intensity(self, temp, laser_freq, width, non_anal_dir=None, min_freq=None,
                                     max_freq=None, num=1000, relative=False, units="eV"):
        """
        Calculates the broadened Raman intensities in arbitrary units integrated over all possible
        orientation to reproduce the powder measurements for frequencies in an interval. It is possible to
        use the susceptibilities from the transverse modes only or to specify one of the directions with non
        analytical contributions.

        Args:
            temp: temperature in K.
            laser_freq: frequency of the incident laser. The units are determined the "units"
                argument.
            width: the width of the Lorentz distribution. The units are determined the "units"
                argument.
            non_anal_dir: index of the direction along which the non analytical contribution
                has been calculated. Corresponds to the indices in the non_anal_directions attribute.
            min_freq: minimum frequency considered. If None it will be given by the minimum
                frequency with non zero intensities minus 10 times the width of the distribution.
                If given the units are determined by the "units" argument.
            max_freq: maximum frequency considered. If None it will be given by the maximum
                frequency with non zero intensities plus 10 times the width of the distribution.
                If given the units are determined by the "units" argument.
            num: number of frequencies in the interval (min_freq, max_freq).
            relative: if True the intensities will be rescaled so that the largest value of the
                total intensity is 1.
            units: the units in which the input and the output frequencies will be given.
                Possible values in ("eV", "meV", "Ha", "cm-1", "Thz")

        Returns:
            A PowderIntensity with the parallel, perpendicular and total components of the powder
            intensities. Each one is a Function1D with "num" points.
        """

        pi = self.get_powder_intensity(temp=temp, laser_freq=laser_freq, non_anal_dir=non_anal_dir, units=units)

        freqs, lorentz = self._get_lorentz_freqs_and_factor(intensity=pi.tot, non_anal_dir=non_anal_dir, min_freq=min_freq,
                                                            max_freq=max_freq, num=num, width=width, units=units)

        lpi = np.array([i.dot(lorentz) for i in pi])
        if relative:
            lpi /= lpi[2].max()

        # now convert the frequencies to the desired units for the output
        x = freqs * abu.phfactor_ev2units(units)

        return PowderIntensity(*(Function1D(x, y) for y in lpi))

    @add_fig_kwargs
    def plot_intensity(self, temp, laser_freq, width, value, non_anal_dir=None, min_freq=None, max_freq=None,
                       num=1000, relative=False, units="eV", ax=None, plot_phfreqs=False, **kwargs):
        """
        Plot one representation of the broadened Raman intensities.

        Args:
            temp: temperature in K.
            laser_freq: frequency of the incident laser. The units are determined the "units"
                argument.
            width: the width of the Lorentz distribution. The units are determined the "units"
                argument. If None or 0 a plot of only the frequencies for each mode will be given.
            value: a string describing the value that should be plotted. Can be "powder" or
                a string of the type "xz" with the polarization of the incoming and outgoing
                phonon. All the combinations of "x", "y" and "z" are accepted.
            non_anal_dir: index of the direction along which the non analytical contribution
                has been calculated. Corresponds to the indices in the non_anal_directions attribute.
            min_freq: minimum frequency considered. If None it will be given by the minimum
                frequency with non zero intensities minus 10 times the width of the distribution.
                If given the units are determined by the "units" argument.
            max_freq: maximum frequency considered. If None it will be given by the maximum
                frequency with non zero intensities plus 10 times the width of the distribution.
                If given the units are determined by the "units" argument.
            num: number of frequencies in the interval (min_freq, max_freq).
            relative: if True the intensities will be rescaled so that the largest value of the
                total intensity is 1.
            units: the units in which the input and the output frequencies will be given.
                Possible values in ("eV", "meV", "Ha", "cm-1", "Thz")
            ax: |matplotlib-Axes| or None if a new figure should be created.
            plot_phfreqs: if True vertical dashed lines are added to the figure for all the
                phonon modes.
            **kwargs: arguments passed to the plot function.

        Returns:
            |matplotlib-Figure|
        """

        ax, fig, plt = get_ax_fig_plt(ax=ax)

        if width:
            if value == "powder":
                f = self.get_powder_lorentz_intensity(temp=temp, laser_freq=laser_freq, width=width,
                                                      non_anal_dir=non_anal_dir, min_freq=min_freq, max_freq=max_freq,
                                                      num=num, relative=relative, units=units).tot
            else:
                pol_in = value[0]
                pol_out = value[1]
                f = self.get_lorentz_intensity(temp=temp, laser_freq=laser_freq, width=width, non_anal_dir=non_anal_dir,
                                               min_freq=min_freq, max_freq=max_freq, num=num, relative=relative,
                                               units=units, pol_in=pol_in, pol_out=pol_out)

            f.plot(ax=ax, **kwargs)

            if plot_phfreqs:
                if non_anal_dir is None:
                    w = self.phfreqs
                else:
                    w = self.non_anal_phfreqs[non_anal_dir]

                w = w * abu.phfactor_ev2units(units)

                min_freq = f.mesh[0]
                max_freq = f.mesh[-1]

                for wi in w:
                    if min_freq < wi < max_freq:
                        ax.axvline(x=wi, ls="--", color="k", lw=0.5)

        else:
            if value == "powder":
                ri = self.get_powder_intensity(temp=temp, laser_freq=laser_freq, non_anal_dir=non_anal_dir,
                                               relative=relative, units=units)

                i = ri.tot
            else:
                if len(value) != 2:
                    raise ValueError("The value should contain the ingoing and outgoing polarizations.")
                pol_in = value[0]
                pol_out = value[1]
                i = self.get_modes_intensities(temp=temp, laser_freq=laser_freq, non_anal_dir=non_anal_dir,
                                               relative=relative, units=units, pol_in=pol_in, pol_out=pol_out)

            if non_anal_dir is None:
                w = self.phfreqs * abu.phfactor_ev2units(units)
            else:
                w = self.non_anal_phfreqs[non_anal_dir] * abu.phfactor_ev2units(units)

            ax.stem(w, i, **kwargs)

        return fig
