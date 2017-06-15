"""
Objects used to extract and plot results from output files in text format.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import numpy as np

from collections import OrderedDict
from monty.string import is_string
from monty.functools import lazy_property
from monty.termcolor import cprint
from pymatgen.core.units import bohr_to_ang
from abipy.core.structure import Structure, frames_from_structures, AbinitSpaceGroup
from abipy.core.kpoints import has_timrev_from_kptopt
from abipy.core.mixins import TextFile, AbinitNcFile, NotebookWriter
from abipy.abio.inputs import GEOVARS
from abipy.abio.timer import AbinitTimerParser
from abipy.flowtk import EventsParser, NetcdfReader, GroundStateScfCycle, D2DEScfCycle


class AbinitTextFile(TextFile):
    """
    Class for the ABINIT main output file and the log file.
    """
    @property
    def events(self):
        """
        List of ABINIT events reported in the file.
        """
        # Parse the file the first time the property is accessed or when mtime is changed.
        stat = os.stat(self.filepath)
        if stat.st_mtime != self._last_mtime or not hasattr(self, "_events"):
            self._events = EventsParser().parse(self.filepath)
        return self._events

    def get_timer(self):
        """
        Timer data.
        """
        timer = AbinitTimerParser()
        timer.parse(self.filepath)
        return timer


class AbinitLogFile(AbinitTextFile, NotebookWriter):
    """Class representing the Abinit log file."""

    def write_notebook(self, nbpath=None):
        """
        Write an ipython notebook to nbpath. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("abilog = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(abilog.events)"),
        ])

        return self._write_nb_nbpath(nb, nbpath)


class AbinitOutputFile(AbinitTextFile, NotebookWriter):
    """
    Class representing the main Abinit output file.
    """

    def __init__(self, filepath):
        super(AbinitOutputFile, self).__init__(filepath)
        self._parse()

    def _parse(self):
        """
        header: String with the input variables
        footer: String with the output variables
        datasets: Dictionary mapping dataset index to list of strings.
        """
        # Get code version and find magic line signaling that the output file is completed.
        self.version, self.run_completed = None, False
        with open(self.filepath) as fh:
            for line in fh:
                if self.version is None and line.startswith(".Version"):
                    self.version = line.split()[1]
                if " Calculation completed." in line:
                    self.run_completed = True
                    break

        # Parse header to get important dimensions and variables
        self.header, self.footer, self.datasets = [], [], OrderedDict()
        where = "in_header"
        verbose = 0

        with open(self.filepath, "rt") as fh:
            for line in fh:
                if "== DATASET" in line:
                    # Save dataset number
                    # == DATASET  1 ==================================================================
                    where = int(line.replace("=", "").split()[-1])
                    assert where not in self.datasets
                    self.datasets[where] = []
                elif "== END DATASET(S) " in line:
                    where = "in_footer"

                if where == "in_header":
                    self.header.append(line)
                elif where == "in_footer":
                    self.footer.append(line)
                else:
                    self.datasets[where].append(line)

        if not self.datasets:
            raise NotImplementedError("Empty dataset sections.")

        self.header = "".join(self.header)
        if verbose:
            print("header")
            print(self.header)

        #if " jdtset " in self.header:
        #    raise NotImplementedError("jdtset is not supported")
        #if " udtset " in self.header:
        #    raise NotImplementedError("udtset is not supported")

        for key, data in self.datasets.items():
            if verbose: print("data")
            self.datasets[key] = "".join(data)
            if verbose: print(self.datasets[key])

        self.footer = "".join(self.footer)
        if verbose:
            print("footer")
            print(self.footer)

        self.ndtset = len(self.datasets)
        self.initial_vars_global, self.initial_vars_dataset = self._parse_variables("header")
        self.final_vars_global, self.final_vars_dataset = None, None
        if self.run_completed:
            self.final_vars_global, self.final_vars_dataset = self._parse_variables("footer")

    def _parse_variables(self, what):
        vars_global = OrderedDict()
        vars_dataset = OrderedDict([(k, OrderedDict()) for k in self.datasets.keys()])
        #print("keys", vars_dataset.keys())

        lines = getattr(self, what).splitlines()
        if what == "header":
            magic_start = " -outvars: echo values of preprocessed input variables --------"
        else:
            magic_start = " -outvars: echo values of variables after computation  --------"
        magic_stop = "================================================================================"

        # Select relevant portion with variables.
        for i, line in enumerate(lines):
            if magic_start in line:
                break
        else:
            raise ValueError("Cannot find magic_start line: %s" % magic_start)
        lines = lines[i+1:]

        for i, line in enumerate(lines):
            if magic_stop in line:
                break
        else:
            raise ValueError("Cannot find magic_stop line: %s" % magic_stop)
        lines = lines[:i]

        # Parse data. Assume format:
        #   timopt          -1
        #    tnons      0.0000000  0.0000000  0.0000000     0.2500000  0.2500000  0.2500000
        #               0.0000000  0.0000000  0.0000000     0.2500000  0.2500000  0.2500000
        def get_dtindex_key_value(line):
            tokens = line.split()
            s, value = tokens[0], " ".join(tokens[1:])
            l = []
            for i, c in enumerate(s[::-1]):
                if c.isalpha():
                    key = s[:len(s)-i]
                    break
                l.append(c)
            else:
                raise ValueError("Cannot find dataset index in token: %s" % s)

            #print(line, "\n", l)
            dtindex = None
            if l:
                l.reverse()
                dtindex = int("".join(l))
            return dtindex, key, value

        # (varname, dtindex), [line1, line2 ...]
        stack_var, stack_lines = None, []
        def pop_stack():
            if stack_lines:
                key, dtidx = stack_var
                value = " ".join(stack_lines)
                if dtidx is None:
                    vars_global[key] = value
                else:
                    vars_dataset[dtidx][key] = value

        for line in lines:
            if not line: continue
            # Ignore first char
            line = line[1:].lstrip().rstrip()
            #print(line)
            if line[0].isalpha():
                pop_stack()
                stack_lines = []
                dtidx, key, value = get_dtindex_key_value(line)
                stack_var = (key, dtidx)
                stack_lines.append(value)
            else:
                stack_lines.append(line)

        pop_stack()

        return vars_global, vars_dataset

    def _get_structures(self, what):
        if what == "header":
            vars_global, vars_dataset = self.initial_vars_global, self.initial_vars_dataset
        else:
            vars_global, vars_dataset = self.final_vars_global, self.final_vars_dataset

        #print("global", vars_global["acell"])
        from abipy.abio.abivars import is_abiunit
        inigeo = {k: vars_global[k] for k in GEOVARS if k in vars_global}

        spgvars = ("spgroup", "symrel", "tnons", "symafm")
        spgd_global = {k: vars_global[k] for k in spgvars if k in vars_global}
        global_kptopt = vars_global.get("kptopt", 1)

        structures = []
        for i in self.datasets.keys():
            # This code breaks down if there are conflicting GEOVARS in globals and dataset.
            d = inigeo.copy()
            d.update({k: vars_dataset[i][k] for k in GEOVARS if k in vars_dataset[i]})

            for key, value in d.items():
                # Must handle possible unit.
                fact = 1.0
                tokens = [t.lower() for t in value.split()]
                if is_abiunit(tokens[-1]):
                    tokens, unit = tokens[:-1], tokens[-1]
                    if unit in ("angstr", "angstrom", "angstroms"):
                        fact = 1.0 / bohr_to_ang
                    elif unit in ("bohr", "bohrs", "au"):
                        fact = 1.0
                    else:
                        raise ValueError("Don't know how to handle unit: %s" % unit)

                s = " ".join(tokens)
                dtype = np.float if key not in ("ntypat", "typat", "natom") else np.int
                try:
                    #print(key, s)
                    value = np.fromstring(s, sep=" ", dtype=dtype)
                    #print(key, value)
                    if fact != 1.0: value *= fact # Do not change integer arrays e.g typat!
                    d[key] = value
                except ValueError as exc:
                    print(key, s)
                    raise exc

            if "rprim" not in d and "angdeg" not in d: d["rprim"] = np.eye(3)
            if "natom" in d and d["natom"] == 1 and all(k not in d for k in ("xred", "xcart", "xangst")):
                d["xred"] = np.zeros(3)
            #print(d)
            abistr = Structure.from_abivars(d)

            # Extract Abinit spacegroup.
            spgd = spgd_global.copy()
            spgd.update({k: vars_dataset[i][k] for k in spgvars if k in vars_dataset[i]})

            spgid = int(spgd.get("spgroup", 0))
            if "symrel" not in spgd:
                symrel = np.reshape(np.eye(3, 3, dtype=np.int), (1, 3, 3))
                spgd["symrel"] = " ".join((str(i) for i in symrel.flatten()))
            else:
                symrel = np.reshape(np.array([int(n) for n in spgd["symrel"].split()], dtype=np.int), (-1, 3, 3))
            nsym = len(symrel)
            assert nsym == spgd.get("nsym", nsym) #; print(symrel.shape)

            if "tnons" in spgd:
                tnons = np.reshape(np.array([float(t) for t in spgd["tnons"].split()], dtype=np.float), (nsym, 3))
            else:
                tnons = np.zeros((nsym, 3))

            if "symafm" in spgd:
                symafm = np.array([int(n) for n in spgd["symafm"].split()], dtype=np.int)
                symafm.shape = (nsym,)
            else:
                symafm = np.ones(nsym, dtype=np.int)

            try:
                has_timerev = has_timrev_from_kptopt(vars_dataset[i].get("kptopt", global_kptopt))
                abi_spacegroup = AbinitSpaceGroup(spgid, symrel, tnons, symafm, has_timerev, inord="C")
                abistr.set_abi_spacegroup(abi_spacegroup)
            except Exception as exc:
                print("Cannot build AbinitSpaceGroup from the variables reported in file!\n", str(exc))

            structures.append(abistr)

        return structures

    @lazy_property
    def initial_structures(self):
        """List of initial structures."""
        return self._get_structures("header")

    @property
    def has_same_initial_structures(self):
        """True if all initial structures are equal."""
        return all(self.initial_structures[0] == s for s in self.initial_structures)

    @lazy_property
    def final_structures(self):
        """List of final structures."""
        if self.run_completed:
            return self._get_structures("footer")
        else:
            print("Cannot extract final structures from file.\n %s" % str(exc))
            return []

    @lazy_property
    def initial_structure(self):
        """
        The structure defined in the output file.

        If the input file contains multiple datasets **AND** the datasets
        have different structures, this property returns None.
        In this case, one has to access the structure of the individual datasets.
        For example:

            self.initial_structures[0]

        gives the structure of the first dataset.
        """
        if not self.has_same_initial_structures:
            print("Datasets have different structures. Returning None. Use initial_structures[0]")
            return None
        return self.initial_structures[0]

    @property
    def has_same_final_structures(self):
        """True if all initial structures are equal."""
        return all(self.final_structures[0] == s for s in self.final_structures)

    @lazy_property
    def final_structure(self):
        """
        The structure defined in the output file.

        If the input file contains multiple datasets **AND** the datasets
        have different structures, this property returns None.
        In this case, one has to access the structure of the individual datasets.
        For example:

            self.final_structures[0]

        gives the structure of the first dataset.
        """
        if not self.has_same_final_structures:
            print("Datasets have different structures. Returning None. Use final_structures[0]")
            return None
        return self.final_structures[0]

    def diff_datasets(self, dt_list1, dt_list2, with_params=True, differ="html", dryrun=False):
        """
        Compare datasets
        """
        if not isinstance(dt_list1, (list, tuple)): dt_list1 = [dt_list1]
        if not isinstance(dt_list2, (list, tuple)): dt_list2 = [dt_list2]

        dt_lists = [dt_list1, dt_list2]
        import tempfile
        tmp_names = []
        for i in range(2):
            _, tmpname = tempfile.mkstemp(text=True)
            tmp_names.append(tmpname)
            with open(tmpname, "wt") as fh:
                if with_params: fh.write(self.header)
                for idt in dt_lists[i]:
                    fh.write(self.datasets[idt])
                if with_params: fh.write(self.footer)

        if differ == "html":
            from abipy.tools.devtools import HtmlDiff
            diff = HtmlDiff(tmp_names)
            if dryrun:
                return diff
            else:
                return diff.open_browser()
        else:
            cmd = "%s %s %s" % (differ, tmp_names[0], tmp_names[1])
            if dryrun:
                return cmd
            else:
                return os.system(cmd)

    def to_string(self, verbose=0):
        lines = ["ndtset: %d, completed: %s" % (self.ndtset, self.run_completed)]
        app = lines.append

        # Different cases depending whether final structures are available
        # and whether structures are equivalent.
        if self.run_completed:
            if self.has_same_final_structures:
                if self.initial_structure != self.final_structure:
                    # Structural relaxation.
                    df = frames_from_structures([self.initial_structure, self.final_structure],
                                                index=["initial", "final"])
                    app("Lattice parameters:")
                    app(str(df.lattice))
                    app("Atomic coordinates:")
                    app(str(df.coords))
                else:
                    # initial == final. Print final structure.
                    app(self.final_structure.to_string(verbose=verbose))
        else:
            # Final structures are not available.
            if self.has_same_initial_structures:
                app(self.initial_structure.to_string(verbose=verbose))
            else:
                df = frames_from_structures(self.initial_structures,
                                            index=[i+1 for i in range(self.ndtset)])
                app("Lattice parameters:")
                app(str(df.lattice))
                app("Atomic coordinates:")
                app(str(df.coords))

        return "\n".join(lines)

    def next_gs_scf_cycle(self):
        """
        Return the next :class:`GroundStateScfCycle` in the file. None if not found.
        """
        return GroundStateScfCycle.from_stream(self)

    def next_d2de_scf_cycle(self):
        """
        Return :class:`GroundStateScfCycle` with information on the GS iterations. None if not found.
        """
        return D2DEScfCycle.from_stream(self)

    def compare_gs_scf_cycles(self, others, show=True):
        """
        Produce and returns a list of `matplotlib` figure comparing the GS self-consistent
        cycle in self with the ones in others.

        Args:
            others: list of `AbinitOutputFile` objects or strings with paths to output files.
            show: True to diplay plots.
        """
        # Open file here if we receive a string. Files will be closed before returning
        close_files = []
        for i, other in enumerate(others):
            if is_string(other):
                others[i] = self.__class__.from_file(other)
                close_files.append(i)

        fig, figures = None, []
        while True:
            cycle = self.next_gs_scf_cycle()
            if cycle is None: break

            fig = cycle.plot(show=False)
            for i, other in enumerate(others):
                other_cycle = other.next_gs_scf_cycle()
                if other_cycle is None: break
                last = (i == len(others) - 1)
                fig = other_cycle.plot(axlist=fig.axes, show=show and last)
                if last:
                    fig.tight_layout()
                    figures.append(fig)

        self.seek(0)
        for other in others: other.seek(0)

        if close_files:
            for i in close_files: others[i].close()

        return figures

    def compare_d2de_scf_cycles(self, others, show=True):
        """
        Produce and returns a `matplotlib` figure comparing the DFPT self-consistent
        cycle in self with the ones in others.

        Args:
            others: list of `AbinitOutputFile` objects or strings with paths to output files.
            show: True to diplay plots.
        """
        # Open file here if we receive a string. Files will be closed before returning
        close_files = []
        for i, other in enumerate(others):
            if is_string(other):
                others[i] = self.__class__.from_file(other)
                close_files.append(i)

        fig, figures = None, []
        while True:
            cycle = self.next_d2de_scf_cycle()
            if cycle is None: break

            fig = cycle.plot(show=False)
            for i, other in enumerate(others):
                other_cycle = other.next_d2de_scf_cycle()
                if other_cycle is None: break
                last = (i == len(others) - 1)
                fig = other_cycle.plot(axlist=fig.axes, show=show and last)
                if last:
                    fig.tight_layout()
                    figures.append(fig)

        self.seek(0)
        for other in others: other.seek(0)

        if close_files:
            for i in close_files: others[i].close()

        return figures

    def write_notebook(self, nbpath=None):
        """
        Write an ipython notebook to nbpath. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("abo = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(abo.events)"),
            nbv.new_code_cell("""
abo.seek(0); icycle = -1
while True:
    gs_cycle = abo.next_gs_scf_cycle()
    if gs_cycle is None: break
    icycle += 1
    gs_cycle.plot(title="SCF cycle no %d" % icycle, tight_layout=True)"""),

        nbv.new_code_cell("""
abo.seek(0); icycle = -1
while True:
    d2de_cycle = abo.next_d2de_scf_cycle()
    if d2de_cycle is None: break
    icycle += 1
    d2de_cycle.plot(title="DFPT cycle no %d" % icycle, tight_layout=True)"""),

       nbv.new_code_cell("""
abo.seek(0)
timer = abo.get_timer()
if timer:
    timer.plot_all()
"""),
        ])

        return self._write_nb_nbpath(nb, nbpath)


class OutNcFile(AbinitNcFile):
    """
    Class representing the _OUT.nc file containing the dataset results
    produced at the end of the run. The netcdf variables can be accessed
    via instance attribute e.g. `outfile.ecut`. Provides integration with ipython.
    """
    def __init__(self, filepath):
        super(OutNcFile, self).__init__(filepath)
        self.reader = NetcdfReader(filepath)
        self._varscache= {k: None for k in self.reader.rootgrp.variables}

    def __dir__(self):
        """Ipython integration."""
        return sorted(list(self._varscache.keys()))

    def __getattribute__(self, name):
        try:
            return super(OutNcFile, self).__getattribute__(name)
        except AttributeError:
            # Look in self._varscache
            varscache = super(OutNcFile, self).__getattribute__("_varscache")
            if name not in varscache:
                raise AttributeError("Cannot find attribute %s" % name)
            reader = super(OutNcFile, self).__getattribute__("reader")
            if varscache[name] is None:
                varscache[name] = reader.read_value(name)
            return varscache[name]

    def close(self):
        self.reader.close()

    def get_allvars(self):
        """
        Read all netcdf variables present in the file.
        Return dictionary varname --> value
        """
        for k, v in self._varscache.items():
            if v is not None: continue
            self._varscache[k] = self.reader.read_value(k)
        return self._varscache


def validate_output_parser(abitests_dir=None, output_files=None):
    """
    Validate/test Abinit output parser.

    Args:
        dirpath: Abinit tests directory.
        output_files: List of Abinit output files.

    Return: Exit code.
    """
    def is_abinit_output(path):
        """
        True if path is one of the output files used in the Abinit Test suite.
        """
        if path.endswith(".abo"): return True
        if not path.endswith(".out"): return False

        with open(path, "rt") as fh:
            for i, line in enumerate(fh):
                if i == 1:
                    line = line.rstrip().lower()
                    if line.endswith("abinit"): return True
                    return False
            return False

    # Files are collected in paths.
    paths = []

    if abitests_dir is not None:
        print("Analyzing directory %s for input files" % abitests_dir)

        for dirpath, dirnames, filenames in os.walk(abitests_dir):
            for fname in filenames:
                path = os.path.join(dirpath, fname)
                if is_abinit_output(path): paths.append(path)

    if output_files is not None:
        print("Analyzing files ", str(output_files))
        for arg in output_files:
            if is_abinit_output(arg): paths.append(arg)

    nfiles = len(paths)
    if nfiles == 0:
        cprint("Empty list of input files.", "red")
        return 0

    print("Found %d Abinit output files" % len(paths))
    errpaths = []
    for path in paths:
        print(path + ": ", end="")
        try:
            out = AbinitOutputFile.from_file(path)
            s = out.to_string(verbose=2)
            assert out.run_completed
            cprint("OK", "green")
        except Exception as exc:
            if not isinstance(exc, NotImplementedError):
                cprint("FAILED", "red")
                errpaths.append(path)
                import traceback
                print(traceback.format_exc())
                #print("[%s]: Exception:\n%s" % (path, str(exc)))
                #with open(path, "rt") as fh:
                #    print(10*"=" + "Input File" + 10*"=")
                #    print(fh.read())
                #    print()
            else:
                cprint("NOTIMPLEMENTED", "magenta")

    if errpaths:
        cprint("failed: %d/%d [%.1f%%]" % (len(errpaths), nfiles, 100 * len(errpaths)/nfiles), "red")
        for i, epath in enumerate(errpaths):
            cprint("[%d] %s" % (i, epath), "red")
    else:
        cprint("All input files successfully parsed!", "green")

    return len(errpaths)
