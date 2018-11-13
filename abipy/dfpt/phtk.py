# coding: utf-8
"""
Phonon Toolkit: This module gathers low-level tools to operate on phonons.
"""
from __future__ import print_function, division, absolute_import  # unicode_literals,

import warnings
import numpy as np
import abipy.core.abinit_units as abu

from monty.functools import lazy_property
from pymatgen.core.periodic_table import Element
from abipy.core.mixins import Has_Structure
from abipy.iotools import ETSF_Reader


# TODO: amu should become mandatory.
def get_dyn_mat_eigenvec(phdispl, structure, amu=None, amu_symbol=None):
    """
    Converts the phonon displacements to the orthonormal eigenvectors of the dynamical matrix.
    Small discrepancies with the original values may be expected due to the different values of the atomic masses in
    abinit and pymatgen.

    .. note::

        These eigenvectors are orthonormalized and should be very close to the ones computed by Abinit in a.u.
        Note, however, that the output vectors are given in atomic units so dividing then by the sqrt(Mass)
        won't give the dipl_cart used in PhononBands that are in Angstrom.

    Args:
        phdispl: a numpy array containing the displacements in cartesian coordinates. The last index should have
            size 3*(num atoms), but the rest of the shape is arbitrary. If qpts is not None the first dimension
            should match the q points.
        structure: |Structure| object.
        amu: dictionary that associates the atomic numbers present in the structure to the values of the atomic
            mass units used for the calculation. Incompatible with amu_sumbol. If None and amu_symbol is None, values
            from pymatgen will be used.  Note that this will almost always lead to inaccuracies in the conversion.
        amu_symbol: dictionary that associates the symbol present in the structure to the values of the atomic
            mass units used for the calculation. Incompatible with amu. If None and amu_symbol is None, values from
            pymatgen will be used. that this will almost always lead to inaccuracies in the conversion.

    Returns:
        A |numpy-array| of the same shape as phdispl containing the eigenvectors of the dynamical matrix
    """
    eigvec = np.zeros(np.shape(phdispl), dtype=np.complex)

    if amu is not None and amu_symbol is not None:
        raise ValueError("Only one between amu and amu_symbol should be provided!")

    if amu is not None:
        amu_symbol = {Element.from_Z(n).symbol: v for n, v in amu.items()}

    if amu_symbol is None:
        warnings.warn("get_dyn_mat_eigenvec has been called with amu=None. Eigenvectors may not be orthonormal.")
        amu_symbol = {e.symbol: e.atomic_mass for e in structure.composition.elements}

    for j, a in enumerate(structure):
        eigvec[...,3*j:3*(j+1)] = phdispl[...,3*j:3*(j+1)] * np.sqrt(amu_symbol[a.specie.symbol]*abu.amu_emass) / abu.Bohr_Ang

    return eigvec


def match_eigenvectors(v1, v2):
    """
    Given two list of vectors, returns the pair matching based on the complex scalar product.
    Returns the indices of the second list that match the vectors of the first list in ascending order.
    """
    prod = np.absolute(np.dot(v1, v2.transpose().conjugate()))

    indices = np.zeros(len(v1), dtype=np.int)
    missing_v1 = [True] * len(v1)
    missing_v2 = [True] * len(v1)
    for m in reversed(np.argsort(prod, axis=None)):
        i, j = np.unravel_index(m, prod.shape)
        if missing_v1[i] and missing_v2[j]:
            indices[i] = j
            missing_v1[i] = missing_v2[j] = False
            if not any(missing_v1):
                if any(missing_v2):
                    raise RuntimeError('Something went wrong in matching vectors: {} {}'.format(v1, v2))
                break

    return indices


class NonAnalyticalPh(Has_Structure):
    """
    Phonon data at gamma including non analytical contributions
    Read from anaddb.nc
    """

    def __init__(self, structure, directions, phfreqs, phdispl_cart, amu=None):
        """
        Args:
            structure: |Structure| object.
            directions: Cartesian directions along which the non analytical frequencies have been calculated
            phfreqs: Phonon frequencies with non analytical contribution in eV along directions
            phdispl_cart: Displacement in Angstrom in Cartesian coordinates with non analytical contribution
                along directions
            amu: dictionary that associates the atomic species present in the structure to the values of the atomic
                mass units used for the calculation
        """
        self._structure = structure
        self.directions = directions
        self.phfreqs = phfreqs
        self.phdispl_cart = phdispl_cart
        self.amu = amu
        self.amu_symbol = None
        if amu is not None:
            self.amu_symbol = {}
            for z, m in amu.items():
                el = Element.from_Z(int(z))
                self.amu_symbol[el.symbol] = m

    @classmethod
    def from_file(cls, filepath):
        """
        Reads the non analytical directions, frequencies and displacements from the anaddb.nc file specified.
        Non existence of displacements is accepted for compatibility with abinit 8.0.6
        Raises an error if the other values are not present in anaddb.nc.
        """
        with ETSF_Reader(filepath) as r:
            directions = r.read_value("non_analytical_directions")
            phfreq = r.read_value("non_analytical_phonon_modes")

            # need a default as the first abinit version including IFCs in the netcdf doesn't have this attribute
            phdispl_cart = r.read_value("non_analytical_phdispl_cart", cmode="c", default=None)

            structure = r.read_structure()

            amu_list = r.read_value("atomic_mass_units", default=None)
            if amu_list is not None:
                # ntypat arrays
                atomic_numbers = r.read_value("atomic_numbers")
                amu = {at: a for at, a in zip(atomic_numbers, amu_list)}
            else:
                amu = None

            return cls(structure=structure, directions=directions, phfreqs=phfreq, phdispl_cart=phdispl_cart, amu=amu)

    @lazy_property
    def dyn_mat_eigenvect(self):
        """
        [ndirection, 3*natom, 3*natom] array with the orthonormal eigenvectors of the dynamical matrix.
        in Cartesian coordinates.
        """
        return get_dyn_mat_eigenvec(self.phdispl_cart, self.structure, amu=self.amu)

    @property
    def structure(self):
        """|Structure| object."""
        return self._structure

    def index_direction(self, direction, cartesian=False):
        """
        Returns: the index of direction. Raises: `ValueError` if not found.

        Args:
            direction: a 3 element list indicating the direction. Can be a generic vector
            cartesian: if True the direction are already in cartesian coordinates, if False it
                will be converted to match the internal description of the directions.
        """
        if not cartesian:
            direction = self.structure.lattice.reciprocal_lattice_crystallographic.get_cartesian_coords(direction)
        else:
            direction = np.array(direction)
        direction = direction / np.linalg.norm(direction)

        for i, d in enumerate(self.directions):
            d = d / np.linalg.norm(d)
            if np.allclose(d, direction):
                return i

        raise ValueError("Cannot find direction: `%s` with cartesian: `%s` in non_analytical cartesian directions:\n%s" %
                (str(direction), cartesian, str(self.directions)))

    def has_direction(self, direction, cartesian=False):
        """
        Checks if the input direction is among those available.

        Args:
            direction: a 3 element list indicating the direction. Can be a generic vector
            cartesian: if True the direction are already in cartesian coordinates, if False it
                will be converted to match the internal description of the directions.
        """
        try:
            self.index_direction(direction, cartesian=cartesian)
            return True
        except ValueError:
            return False


def open_file_phononwebsite(filename, port=8000,
                            website="http://henriquemiranda.github.io/phononwebsite",
                            host="localhost", browser=None): # pragma: no cover
    """
    Take a file, detect the type and open it on the phonon website
    Based on a similar function in <https://github.com/henriquemiranda/phononwebsite/phononweb.py>

    Args:
        filename: file with phonon data in phononwebsite format.
        port: Initial port.
        website: Website URL
        host: localhost name.
        browser: Open webpage in ``browser``. Use default if $BROWSER if None.
    """
    if filename.endswith(".json"):
        filetype = "json"
    elif filename.endswith(".yaml"):
        filetype = "yaml"
    else:
        filetype = "rest"

    try:
        from http.server import HTTPServer, SimpleHTTPRequestHandler
    except ImportError:
        from BaseHTTPServer import HTTPServer
        # python 2 requires internal implementation
        from abipy.tools.SimpleHTTPServer import SimpleHTTPRequestHandler

    # Add CORS header to the website
    class CORSRequestHandler (SimpleHTTPRequestHandler):
        def end_headers (self):
            #self.send_header('Access-Control-Allow-Origin', website)
            self.send_header('Access-Control-Allow-Origin', "http://henriquemiranda.github.io")
            SimpleHTTPRequestHandler.end_headers(self)
        def log_message(self, format, *args):
            return

    # Initialize http server thread
    print('Starting HTTP server at port %d ...' % port, end=" ")
    trial, max_ntrial = 0, 50
    while trial < max_ntrial:
        try:
            server = HTTPServer(('', port), CORSRequestHandler)
            #print("got port:", port)
            break
        except OSError:
            trial += 1
            port += 1
            print(port, end=", ")
    else:
        raise RuntimeError("Cannot find available port after %s attempts" % max_ntrial)

    # Create threads python
    server.url = 'http://{}:{}'.format(host, server.server_port)
    from threading import Thread
    t = Thread(target=server.serve_forever)
    t.daemon = True
    t.start()

    # Open website with the file
    try:
        from urllib.parse import quote
    except ImportError:
        from urllib import quote

    url_filename = 'http://{}:{}/{}'.format(host, server.server_port, quote(filename))
    url = '%s/phonon.html?%s=%s' % (website, filetype, url_filename)
    print("\nOpening URL:", url)
    print("Using default browser, if the webpage is not displayed correctly",
          "\ntry to change browser either via command line options or directly in the shell with e.g:\n\n"
          "     export BROWSER=firefox\n")
    print('Press Ctrl+C to terminate the HTTP server')
    import webbrowser
    webbrowser.get(browser).open_new_tab(url)

    # Quit application when SIGINT is received
    def signal_handler(signal, frame):
        sys.exit(0)

    import signal
    signal.signal(signal.SIGINT, signal_handler)
    signal.pause()

