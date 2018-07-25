"""
Objects used to extract and plot results from output files in text format.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import numpy as np
import pandas as pd

from collections import OrderedDict
from six.moves import cStringIO
from monty.string import is_string, marquee
from monty.functools import lazy_property
from monty.termcolor import cprint
from pymatgen.core.units import bohr_to_ang
from abipy.core.symmetries import AbinitSpaceGroup
from abipy.core.structure import Structure, dataframes_from_structures
from abipy.core.kpoints import has_timrev_from_kptopt
from abipy.core.mixins import TextFile, AbinitNcFile, NotebookWriter
from abipy.abio.inputs import GEOVARS
from abipy.abio.timer import AbinitTimerParser
from abipy.abio.robots import Robot
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
    """
    Class representing the Abinit log file.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: AbinitLogFile
    """

    def to_string(self, verbose=0):
        return str(self.events)

    def plot(self, **kwargs):
        """Empty placeholder."""
        return None

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        """
        yield None

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter_ notebook to ``nbpath``. If nbpath is None, a temporay file in the current
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

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: AbinitOutputFile
    """
    # TODO: Extract number of errors and warnings.

    def __init__(self, filepath):
        super(AbinitOutputFile, self).__init__(filepath)
        self.debug_level = 0
        self._parse()

    def _parse(self):
        """
        header: String with the input variables
        footer: String with the output variables
        datasets: Dictionary mapping dataset index to list of strings.
        """
        # Get code version and find magic line signaling that the output file is completed.
        self.version, self.run_completed = None, False
        self.overall_cputime, self.overall_walltime = 0.0, 0.0
        self.proc0_cputime, self.proc0_walltime = 0.0, 0.0

        with open(self.filepath) as fh:
            for line in fh:
                if self.version is None and line.startswith(".Version"):
                    self.version = line.split()[1]

                if line.startswith("- Proc."):
                    #- Proc.   0 individual time (sec): cpu=         25.5  wall=         26.1
                    tokens = line.split()
                    self.proc0_walltime = float(tokens[-1])
                    self.proc0_cputime = float(tokens[-3])

                if line.startswith("+Overall time"):
                    #+Overall time at end (sec) : cpu=         25.5  wall=         26.1
                    tokens = line.split()
                    self.overall_cputime = float(tokens[-3])
                    self.overall_walltime = float(tokens[-1])

                if " Calculation completed." in line:
                    self.run_completed = True

        # Parse header to get important dimensions and variables
        self.header, self.footer, self.datasets = [], [], OrderedDict()
        where = "in_header"

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
                    # dataset number --> lines
                    self.datasets[where].append(line)

        self.header = "".join(self.header)
        if self.debug_level: print("header:\n", self.header)
        # Output files produced in dryrun_mode contain the following line:
        # abinit : before driver, prtvol=0, debugging mode => will skip driver
        self.dryrun_mode = "debugging mode => will skip driver" in self.header
        #print("dryrun_mode:", self.dryrun_mode)

        #if " jdtset " in self.header: raise NotImplementedError("jdtset is not supported")
        #if " udtset " in self.header: raise NotImplementedError("udtset is not supported")

        self.ndtset = len(self.datasets)
        if not self.datasets:
            #raise NotImplementedError("Empty dataset sections.")
            self.ndtset = 1
            self.datasets[1] = "Empty dataset"

        for key, data in self.datasets.items():
            if self.debug_level: print("data")
            self.datasets[key] = "".join(data)
            if self.debug_level: print(self.datasets[key])

        self.footer = "".join(self.footer)
        if self.debug_level: print("footer:\n", self.footer)

        self.initial_vars_global, self.initial_vars_dataset = self._parse_variables("header")
        self.final_vars_global, self.final_vars_dataset = None, None
        if self.run_completed:
            if self.dryrun_mode:
                # footer is not present. Copy values from header.
                self.final_vars_global, self.final_vars_dataset = self.initial_vars_global, self.initial_vars_dataset
            else:
                self.final_vars_global, self.final_vars_dataset = self._parse_variables("footer")

    def _parse_variables(self, what):
        vars_global = OrderedDict()
        vars_dataset = OrderedDict([(k, OrderedDict()) for k in self.datasets.keys()])
        #print("keys", vars_dataset.keys())

        lines = getattr(self, what).splitlines()
        if what == "header":
            magic_start = " -outvars: echo values of preprocessed input variables --------"
        elif what == "footer":
            magic_start = " -outvars: echo values of variables after computation  --------"
        else:
            raise ValueError("Invalid value for what: `%s`" % str(what))

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
            if not line: continue
            #print("line", line)
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
        elif what == "footer":
            vars_global, vars_dataset = self.final_vars_global, self.final_vars_dataset
        else:
            raise ValueError("Invalid value for what: `%s`" % str(what))

        #print("global", vars_global["acell"])
        from abipy.abio.abivars import is_abiunit
        inigeo = {k: vars_global[k] for k in GEOVARS if k in vars_global}

        spgvars = ("spgroup", "symrel", "tnons", "symafm")
        spgd_global = {k: vars_global[k] for k in spgvars if k in vars_global}
        global_kptopt = vars_global.get("kptopt", 1)

        structures = []
        for i in self.datasets:
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
        """List of initial |Structure|."""
        return self._get_structures("header")

    @property
    def has_same_initial_structures(self):
        """True if all initial structures are equal."""
        return all(self.initial_structures[0] == s for s in self.initial_structures)

    @lazy_property
    def final_structures(self):
        """List of final |Structure|."""
        if self.run_completed:
            return self._get_structures("footer")
        else:
            print("Cannot extract final structures from file.\n %s" % str(exc))
            return []

    @lazy_property
    def initial_structure(self):
        """
        The |Structure| defined in the output file.

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
        The |Structure| defined in the output file.

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

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose=0):
        """String representation."""
        lines = ["ndtset: %d, completed: %s" % (self.ndtset, self.run_completed)]
        app = lines.append

        # Different cases depending whether final structures are available
        # and whether structures are equivalent.
        if self.run_completed:
            if self.has_same_final_structures:
                if self.initial_structure != self.final_structure:
                    # Structural relaxation.
                    df = dataframes_from_structures([self.initial_structure, self.final_structure],
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
                df = dataframes_from_structures(self.initial_structures,
                                                index=[i+1 for i in range(self.ndtset)])
                app("Lattice parameters:")
                app(str(df.lattice))
                app("Atomic coordinates:")
                app(str(df.coords))

        # Print dataframe with dimensions.
        dims_dataset, spginfo_dataset = self.get_dims_spginfo_dataset(verbose=verbose)
        rows = []
        for dtind, dims in dims_dataset.items():
            d = OrderedDict()
            d["dataset"] = dtind
            d.update(dims)
            d.update(spginfo_dataset[dtind])
            rows.append(d)

        from abipy.tools.printing import print_dataframe
        df = pd.DataFrame(rows, columns=list(rows[0].keys()) if rows else None)
        df = df.set_index('dataset')
        strio = cStringIO()
        print_dataframe(df, file=strio)
        strio.seek(0)
        app("")
        app(marquee("Dimensions of calculation", mark="="))
        app("".join(strio))

        return "\n".join(lines)

    def get_dims_spginfo_dataset(self, verbose=0):
        """
        Parse the section with the dimensions of the calculation.

        Args:
            verbose: Verbosity level.

        Return: (dims_dataset, spginfo_dataset)
            where dims_dataset[i] is an OrderedDict with the dimensions of dataset `i`
            spginfo_dataset[i] is a dictionary with space group information.
        """
        # If single dataset, we have to parse
        #
        #  Symmetries : space group Fd -3 m (#227); Bravais cF (face-center cubic)
        # ================================================================================
        #  Values of the parameters that define the memory need of the present run
        #      intxc =       0    ionmov =       0      iscf =       7    lmnmax =       6
        #      lnmax =       6     mgfft =      18  mpssoang =       3    mqgrid =    3001
        #      natom =       2  nloc_mem =       1    nspden =       1   nspinor =       1
        #     nsppol =       1      nsym =      48    n1xccc =    2501    ntypat =       1
        #     occopt =       1   xclevel =       2
        # -    mband =           8        mffmem =           1         mkmem =          29
        #        mpw =         202          nfft =        5832          nkpt =          29
        # ================================================================================
        # P This job should need less than                       3.389 Mbytes of memory.
        #   Rough estimation (10% accuracy) of disk space for files :
        # _ WF disk file :      0.717 Mbytes ; DEN or POT disk file :      0.046 Mbytes.
        # ================================================================================

        # If multi datasets we have to parse:

        #  DATASET    2 : space group F-4 3 m (#216); Bravais cF (face-center cubic)
        # ================================================================================
        #  Values of the parameters that define the memory need for DATASET  2.
        #      intxc =       0    ionmov =       0      iscf =       7    lmnmax =       2
        #      lnmax =       2     mgfft =      12  mpssoang =       3    mqgrid =    3001
        #      natom =       2  nloc_mem =       1    nspden =       1   nspinor =       1
        #     nsppol =       1      nsym =      24    n1xccc =    2501    ntypat =       2
        #     occopt =       1   xclevel =       1
        # -    mband =          10        mffmem =           1         mkmem =           2
        #        mpw =          69          nfft =        1728          nkpt =           2
        # ================================================================================
        # P This job should need less than                       1.331 Mbytes of memory.
        #   Rough estimation (10% accuracy) of disk space for files :
        # _ WF disk file :      0.023 Mbytes ; DEN or POT disk file :      0.015 Mbytes.
        # ================================================================================

        magic = "Values of the parameters that define the memory need"
        memory_pre = "P This job should need less than"
        magic_exit = "------------- Echo of variables that govern the present computation"
        filesizes_pre = "_ WF disk file :"
        #verbose = 1

        def parse_spgline(line):
            """Parse the line with space group info, return dict."""
            # Could use regular expressions ...
            i = line.find("space group")
            spg_str, brav_str = line[i:].replace("space group", "").split(";")
            toks = spg_str.split()
            return {
                "spg_symbol": "".join(toks[:-1]),
                "spg_number": int(toks[-1].replace("(", "").replace(")", "").replace("#", "")),
                "bravais": brav_str.strip(),
            }

        from abipy.tools.numtools import grouper
        dims_dataset, spginfo_dataset = OrderedDict(), OrderedDict()
        inblock = 0
        with open(self.filepath, "rt") as fh:
            for line in fh:
                line = line.strip()
                if verbose > 1: print("inblock:", inblock, " at line:", line)

                if line.startswith(magic_exit):
                    break
                if (not line or line.startswith("===") or line.startswith("---")
                    or line.startswith("Rough estimation") or line.startswith("PAW method is used")):
                    continue

                if line.startswith("DATASET") or line.startswith("Symmetries :"):
                    # Get dataset index, parse space group and lattice info, init new dims dict.
                    inblock = 1
                    if line.startswith("Symmetries :"):
                        # No multidataset
                        dtindex = 1
                    else:
                        tokens = line.split()
                        dtindex = int(tokens[1])

                    dims_dataset[dtindex] = dims = OrderedDict()
                    spginfo_dataset[dtindex] = parse_spgline(line)
                    continue

                if inblock == 1 and line.startswith(magic):
                    inblock = 2
                    continue

                if inblock == 2:
                    # Lines with data.
                    if line.startswith(memory_pre):
                        dims["mem_per_proc_mb"] = float(line.replace(memory_pre, "").split()[0])
                    elif line.startswith(filesizes_pre):
                        tokens = line.split()
                        mbpos = [i - 1 for i, t in enumerate(tokens) if t.startswith("Mbytes")]
                        assert len(mbpos) == 2
                        dims["wfk_size_mb"] = float(tokens[mbpos[0]])
                        dims["denpot_size_mb"] = float(tokens[mbpos[1]])
                    elif line.startswith("Pmy_natom="):
                        dims.update(my_natom=int(line.replace("Pmy_natom=", "").strip()))
                        #print("my_natom", dims["my_natom"])
                    else:
                        if line and line[0] == "-": line = line[1:]
                        tokens = grouper(2, line.replace("=", "").split())
                        if verbose > 1: print("tokens:", tokens)
                        dims.update([(t[0], int(t[1])) for t in tokens])

            return dims_dataset, spginfo_dataset

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

    def plot(self, tight_layout=True, with_timer=False, show=True):
        """
        Plot GS/DFPT SCF cycles and timer data found in the output file.

        Args:
            with_timer: True if timer section should be plotted
        """
        from abipy.tools.plotting import MplExpose
        with MplExpose(slide_mode=False, slide_timeout=5.0) as e:
            e(self.yield_figs(tight_layout=tight_layout, with_timer=with_timer))

    # TODO: Use header and vars to understand if we have SCF/DFPT/Relaxation
    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        """
        tight_layout = kwargs.pop("tight_layout", True)
        with_timer = kwargs.pop("with_timer", True)

        self.seek(0)
        icycle = -1
        while True:
            gs_cycle = self.next_gs_scf_cycle()
            if gs_cycle is None: break
            icycle += 1
            yield gs_cycle.plot(title="SCF cycle #%d" % icycle, tight_layout=tight_layout, show=False)

        self.seek(0)
        icycle = -1
        while True:
            d2de_cycle = self.next_d2de_scf_cycle()
            if d2de_cycle is None: break
            icycle += 1
            yield d2de_cycle.plot(title="DFPT cycle #%d" % icycle, tight_layout=tight_layout, show=False)

        if with_timer:
            self.seek(0)
            try:
                yield self.get_timer().plot_all(tight_layout=tight_layout, show=False)
            except Exception:
                print("Abinit output files does not contain timopt data")


    def compare_gs_scf_cycles(self, others, show=True):
        """
        Produce and returns a list of matplotlib_ figure comparing the GS self-consistent
        cycle in self with the ones in others.

        Args:
            others: list of :class:`AbinitOutputFile` objects or strings with paths to output files.
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
                fig = other_cycle.plot(ax_list=fig.axes, show=show and last)
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
        Produce and returns a matplotlib_ figure comparing the DFPT self-consistent
        cycle in self with the ones in others.

        Args:
            others: list of :class:`AbinitOutputFile` objects or strings with paths to output files.
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
                fig = other_cycle.plot(ax_list=fig.axes, show=show and last)
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
        Write a jupyter_ notebook to nbpath. If ``nbpath`` is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        nb.cells.extend([
            nbv.new_code_cell("abo = abilab.abiopen('%s')" % self.filepath),
            nbv.new_code_cell("print(abo.events)"),
            nbv.new_code_cell("abo.plot()"),
        ])

        return self._write_nb_nbpath(nb, nbpath)


def validate_output_parser(abitests_dir=None, output_files=None):  # pragma: no cover
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
        if not path.endswith(".abo"): return False
        if not path.endswith(".out"): return False

        with open(path, "rt") as fh:
            for i, line in enumerate(fh):
                if i == 1:
                    return line.rstrip().lower().endswith("abinit")
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
        print("Analyzing files:", str(output_files))
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


class AboRobot(Robot):
    """
    This robot analyzes the results contained in multiple Abinit output files.
    Can compare dimensions, SCF cycles, analyze timers.

    .. rubric:: Inheritance Diagram
    .. inheritance-diagram:: AboRobot
    """
    EXT = "abo"

    def get_dims_dataframe(self, with_time=True, index=None):
        """
        Build and return |pandas-DataFrame| with the dimensions of the calculation.

        Args:
            with_time: True if walltime and cputime should be added
            index: Index of the dataframe. Use relative paths of files if None.
        """
        rows, my_index = [], []
        for i, abo in enumerate(self.abifiles):
            try:
                dims_dataset, spg_dataset = abo.get_dims_spginfo_dataset()
            except Exception as exc:
                cprint("Exception while trying to get dimensions from %s\n%s" % (abo.relpath, str(exc)), "yellow")
                continue

            for dtindex, dims in dims_dataset.items():
                dims = dims.copy()
                dims.update({"dtset": dtindex})
                # Add walltime and cputime in seconds
                if with_time:
                    dims.update(OrderedDict([(k, getattr(abo, k)) for k in
                        ("overall_cputime", "proc0_cputime", "overall_walltime", "proc0_walltime")]))
                rows.append(dims)
                my_index.append(abo.relpath if index is None else index[i])

        return pd.DataFrame(rows, index=my_index, columns=list(rows[0].keys()))

    def get_dataframe(self, with_geo=True, with_dims=True, abspath=False, funcs=None):
        """
        Return a |pandas-DataFrame| with the most important results and the filenames as index.

        Args:
            with_geo: True if structure info should be added to the dataframe
            with_dims: True if dimensions should be added
            abspath: True if paths in index should be absolute. Default: Relative to getcwd().
            funcs: Function or list of functions to execute to add more data to the DataFrame.
                Each function receives a |GsrFile| object and returns a tuple (key, value)
                where key is a string with the name of column and value is the value to be inserted.
        """
        rows, row_names = [], []
        for label, abo in self.items():
            row_names.append(label)
            d = OrderedDict()

            if with_dims:
                dims_dataset, spg_dataset = abo.get_dims_spginfo_dataset()
                if len(dims_dataset) > 1:
                    cprint("Multiple datasets are not supported. ARGH!", "yellow")
                d.update(dims_dataset[1])

            # Add info on structure.
            if with_geo and abo.run_completed:
                d.update(abo.final_structure.get_dict4pandas(with_spglib=True))

            # Execute functions
            if funcs is not None: d.update(self._exec_funcs(funcs, abo))
            rows.append(d)

        row_names = row_names if not abspath else self._to_relpaths(row_names)
        return pd.DataFrame(rows, index=row_names, columns=list(rows[0].keys()))

    def get_time_dataframe(self):
        """
        Return a |pandas-DataFrame| with the wall-time, cpu time in seconds and the filenames as index.
        """
        rows, row_names = [], []

        for label, abo in self.items():
            row_names.append(label)
            d = OrderedDict([(k, getattr(abo, k)) for k in
                ("overall_cputime", "proc0_cputime", "overall_walltime", "proc0_walltime")])
            rows.append(d)

        return pd.DataFrame(rows, index=row_names, columns=list(rows[0].keys()))

    # TODO
    #def gridplot_timer(self)

    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        """
        yield None

    def write_notebook(self, nbpath=None):
        """
        Write a jupyter_ notebook to nbpath. If nbpath is None, a temporay file in the current
        working directory is created. Return path to the notebook.
        """
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

        args = [(l, f.filepath) for l, f in self.items()]
        nb.cells.extend([
            #nbv.new_markdown_cell("# This is a markdown cell"),
            nbv.new_code_cell("robot = abilab.AboRobot(*%s)\nrobot.trim_paths()\nrobot" % str(args)),
            nbv.new_code_cell("# robot.get_dims_dataframe()"),
            nbv.new_code_cell("robot.get_dataframe()"),
        ])

        # Mixins
        nb.cells.extend(self.get_baserobot_code_cells())

        return self._write_nb_nbpath(nb, nbpath)


class OutNcFile(AbinitNcFile):
    """
    Class representing the _OUT.nc file containing the dataset results
    produced at the end of the run. The netcdf variables can be accessed
    via instance attribute e.g. ``outfile.ecut``. Provides integration with ipython_.
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

    @lazy_property
    def params(self):
        """:class:`OrderedDict` with parameters that might be subject to convergence studies."""
        return {}

    def close(self):
        """Close the file."""
        self.reader.close()

    def get_allvars(self):
        """
        Read all netcdf_ variables present in the file.
        Return dictionary varname --> value
        """
        for k, v in self._varscache.items():
            if v is not None: continue
            self._varscache[k] = self.reader.read_value(k)
        return self._varscache
