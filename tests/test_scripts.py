# coding: utf-8
"""Test abipy command line scripts."""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import abipy.data as abidata
import abipy.flowtk as flowtk

from scripttest import TestFileEnvironment
from monty.inspect import all_subclasses
from pymatgen.io.abinit.qadapters import QueueAdapter
from abipy.core.testing import AbipyTest
from abipy import abilab


script_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "abipy", "scripts"))


def test_if_all_scripts_are_tested():
    """Testing if all scripts are tested"""
    tested_scripts = set(os.path.basename(c.script) for c in all_subclasses(ScriptTest))
    all_scripts = set(f for f in os.listdir(script_dir) if f.endswith(".py"))
    not_tested = all_scripts.difference(tested_scripts)

    if not_tested:
        print("The following scripts are not tested")
        for i, s in enumerate(not_tested):
            print("[%d] %s" % (i, s))

    assert len(not_tested) == 0


class ScriptTest(AbipyTest):
    loglevel = "--loglevel=ERROR"
    verbose = "--verbose"

    expect_stderr = True   # else tests fail due to warnings and deprecation messages

    def get_env(self, check_help_version=True):
        #import tempfile
        #env = TestFileEnvironment(tempfile.mkdtemp(suffix='', prefix='test_' + script))
        env = TestFileEnvironment()

        # Use Agg backend for plots.
        #env.writefile("matplotlibrc", "backend: Agg")
        with open(os.path.join(env.base_path, "matplotlibrc"), "wt") as fh:
            fh.write("backend: Agg\n")

        if check_help_version:
            # Start with --help. If this does not work...
            r = env.run(self.script, "--help")
            assert r.returncode == 0

            # Script must provide a version option
            r = env.run(self.script, "--version", expect_stderr=self.expect_stderr)
            assert r.returncode == 0
            print("stderr", r.stderr)
            print("stdout", r.stdout)
            verstr = r.stderr.strip()
            if not verstr: verstr = r.stdout.strip()  # py3k
            assert verstr == abilab.__version__

        return env


class TestAbidoc(ScriptTest):
    script = os.path.join(script_dir, "abidoc.py")

    def test_abidoc(self):
        """Testing abidoc.py script"""
        env = self.get_env()
        r = env.run(self.script, "man", "ecut", self.loglevel, self.verbose, expect_stderr=self.expect_stderr)
        r = env.run(self.script, "apropos", "test", self.loglevel, self.verbose, expect_stderr=self.expect_stderr)
        r = env.run(self.script, "find", "paw", self.loglevel, self.verbose, expect_stderr=self.expect_stderr)
        r = env.run(self.script, "list", self.loglevel, self.verbose, expect_stderr=self.expect_stderr)
        r = env.run(self.script, "withdim", "natom", self.loglevel, self.verbose, expect_stderr=self.expect_stderr)


class TestAbinp(ScriptTest):
    script = os.path.join(script_dir, "abinp.py")

    def test_abinp(self):
        """Testing abinp.py script"""
        env = self.get_env()
        runabi = abidata.ref_file("refs/si_ebands/run.abi")
        # Commands operating on input files.
        r = env.run(self.script, "validate", runabi, self.loglevel, self.verbose, expect_stderr=self.expect_stderr)
        r = env.run(self.script, "autoparal", runabi, self.loglevel, self.verbose, expect_stderr=self.expect_stderr)
        r = env.run(self.script, "ibz", runabi, self.loglevel, self.verbose, expect_stderr=self.expect_stderr)
        r = env.run(self.script, "phperts", runabi, self.loglevel, self.verbose, expect_stderr=self.expect_stderr)

        # Commands generating input files.
        gan2_cif = abidata.cif_file("gan2.cif")
        r = env.run(self.script, "gs", gan2_cif, self.loglevel, self.verbose, expect_stderr=self.expect_stderr)
        r = env.run(self.script, "ebands", gan2_cif, self.loglevel, self.verbose, expect_stderr=self.expect_stderr)
        r = env.run(self.script, "phonons", gan2_cif, self.loglevel, self.verbose, expect_stderr=self.expect_stderr)
        r = env.run(self.script, "g0w0", gan2_cif, self.loglevel, self.verbose, expect_stderr=self.expect_stderr)
        r = env.run(self.script, "anaph", gan2_cif, self.loglevel, self.verbose, expect_stderr=self.expect_stderr)


class TestAbiopen(ScriptTest):
    script = os.path.join(script_dir, "abiopen.py")

    def test_abiopen(self):
        """Testing abiopen.py script"""
        env = self.get_env()
        gan2_cif = abidata.cif_file("gan2.cif")
        r = env.run(self.script, gan2_cif, "-p", self.loglevel, self.verbose, expect_stderr=self.expect_stderr)
        for f in ["tgw1_9o_DS4_SIGRES.nc", "si_scf_WFK.nc"]:
            path = abidata.ref_file(f)
            r = env.run(self.script, path, "-p", self.loglevel, self.verbose, expect_stderr=self.expect_stderr)


class TestAbistruct(ScriptTest):
    script = os.path.join(script_dir, "abistruct.py")

    def test_spglib(self):
        """Testing abistruct spglib"""
        ncfile = abidata.ref_file("tgw1_9o_DS4_SIGRES.nc")
        env = self.get_env()
        r = env.run(self.script, "spglib", ncfile, self.loglevel, self.verbose,
                    expect_stderr=self.expect_stderr)

    def test_abispg(self):
        """Testing abistruct abispg"""
        ncfile = abidata.ref_file("tgw1_9o_DS4_SIGRES.nc")
        env = self.get_env()
        r = env.run(self.script, "abispg", ncfile, self.loglevel, self.verbose,
                    expect_stderr=self.expect_stderr)

    def test_convert(self):
        """Testing abistruct convert"""
        ncfile = abidata.ref_file("tgw1_9o_DS4_SIGRES.nc")
        env = self.get_env()
        for fmt in ["cif", "cssr", "POSCAR", "json", "mson", "abivars"]:
            r = env.run(self.script, "convert", ncfile, fmt, self.loglevel, self.verbose,
                        expect_stderr=self.expect_stderr)

    def test_supercell(self):
        """Testing abistruct supercell"""
        cif_file = abidata.cif_file("gan2.cif")
        env = self.get_env()
        r = env.run(self.script, "supercell", cif_file, "-s 2" , "-f", "abivars", self.loglevel, self.verbose,
                    expect_stderr=self.expect_stderr)

    def test_abisanitize(self):
        """Testing abistruct abisanitize"""
        ncfile = abidata.ref_file("tgw1_9o_DS4_SIGRES.nc")
        env = self.get_env()
        r = env.run(self.script, "abisanitize", ncfile, self.loglevel, self.verbose,
                   expect_stderr=self.expect_stderr)

    def test_conventional(self):
        """Testing abistruct conventional"""
        ncfile = abidata.ref_file("tgw1_9o_DS4_SIGRES.nc")
        env = self.get_env()
        r = env.run(self.script, "conventional", ncfile, self.loglevel, self.verbose,
                    expect_stderr=self.expect_stderr)

    def test_kpath(self):
        """Testing abistruct kpath"""
        env = self.get_env()
        ncfile = abidata.ref_file("si_scf_WFK.nc")
        r = env.run(self.script, "kpath", ncfile, self.loglevel, self.verbose,
                    expect_stderr=self.expect_stderr)

    def test_kmesh(self):
        """Testing abistruct kmesh"""
        env = self.get_env()
        ncfile = abidata.ref_file("tgw1_9o_DS4_SIGRES.nc")
        r = env.run(self.script, "kmesh", "--mesh", "2", "2", "2", "--is_shift", "1", "1", "1",
                    "--no-time-reversal", ncfile, expect_stderr=self.expect_stderr)

    def test_lgk(self):
        """Testing abistruct lgk"""
        env = self.get_env()
        ncfile = abidata.ref_file("tgw1_9o_DS4_SIGRES.nc")
        r = env.run(self.script, "lgk", "-k", "0", "0", "0", "--no-time-reversal", ncfile,
                    expect_stderr=self.expect_stderr)


class TestAbicomp(ScriptTest):
    script = os.path.join(script_dir, "abicomp.py")

    def test_abicomp(self):
        """Testing abicomp"""
        env = self.get_env()

        cif_paths = abidata.cif_files("al.cif", "gan.cif", "gan2.cif")
        r = env.run(self.script, "structure", cif_paths[0], cif_paths[1], cif_paths[2], self.loglevel, self.verbose,
                    expect_stderr=self.expect_stderr)

        dirpath = os.path.join(abidata.dirpath, "refs", "si_ebands")
        args = [os.path.join(dirpath, p) for p in ("si_nscf_GSR.nc", "si_scf_WFK.nc")]
        r = env.run(self.script, "ebands", args[0], args[1], self.loglevel, self.verbose,
                    expect_stderr=self.expect_stderr)

        args = [os.path.join(dirpath, p) for p in ("si_scf_GSR.nc", "si_scf_WFK.nc")]
        r = env.run(self.script, "edos", args[0], args[1], self.loglevel, self.verbose,
                    expect_stderr=self.expect_stderr)

        dirpath = os.path.join(abidata.dirpath, "refs", "znse_phonons")
        args = [os.path.join(dirpath, p) for p in ("ZnSe_hex_886.out_PHBST.nc", "ZnSe_hex_886.out_PHBST.nc")]
        r = env.run(self.script, "phbands", args[0], args[1], self.loglevel, self.verbose,
                    expect_stderr=self.expect_stderr)

        #dirpath = os.path.join(abidata.dirpath, "refs", "znse_phonons")
        args = abidata.ref_files("ZnSe_hex_886.out_PHDOS.nc", "trf2_5.out_PHDOS.nc")
        r = env.run(self.script, "phdos", args[0], args[1], self.loglevel, self.verbose,
                    expect_stderr=self.expect_stderr)

        #args = [os.path.join(dirpath, p) for p in ( "ZnSe_hex_qpt_DDB")]
        test_dir = os.path.join(os.path.dirname(__file__),  "..", 'test_files')
        args = [
            os.path.join(abidata.dirpath, "refs", "znse_phonons","ZnSe_hex_qpt_DDB"),
            os.path.join(test_dir, "AlAs_444_nobecs_DDB"),
        ]

        r = env.run(self.script, "ddb", args[0], args[1], self.loglevel, self.verbose,
                    expect_stderr=self.expect_stderr)

        args = abidata.ref_files("si_g0w0ppm_nband10_SIGRES.nc",
                                 "si_g0w0ppm_nband20_SIGRES.nc",
                                 "si_g0w0ppm_nband30_SIGRES.nc")
        r = env.run(self.script, "sigres", args[0], args[1], args[2], self.loglevel, self.verbose,
                    expect_stderr=self.expect_stderr)

        args = abidata.ref_files("si_444_MDF.nc", "si_666_MDF.nc", "si_888_MDF.nc")
        r = env.run(self.script, "mdf", args[0], args[1], args[2], self.loglevel, self.verbose,
                    expect_stderr=self.expect_stderr)

        # TODO
        # args = abidata.ref_files()
        # r = env.run(self.script, "gs_scf", *args, self.loglevel, self.verbose,
        #             expect_stderr=self.expect_stderr)

        #args = abidata.ref_files()
        #r = env.run(self.script, "dfpt2_scf", *args, self.loglevel, self.verbose,

        #args = abidata.ref_files()
        #r = env.run(self.script, "time", *args, self.loglevel, self.verbose,


class TestAbirun(ScriptTest):
    script = os.path.join(script_dir, "abirun.py")

    def test_without_flow(self):
        """Testing abirun.py commands without flow"""
        env = self.get_env()
        no_logo_colors = ["--no-logo", "--no-colors"]

        # Test doc_manager
        r = env.run(self.script, "doc_manager", self.loglevel, self.verbose, *no_logo_colors,
                expect_stderr=self.expect_stderr)
        for qtype in QueueAdapter.all_qtypes():
            r= env.run(self.script, ".", "doc_manager", qtype, self.loglevel, self.verbose, *no_logo_colors,
                       expect_stderr=self.expect_stderr)

        # Test doc_sheduler
        r = env.run(self.script, "doc_scheduler", self.loglevel, self.verbose, *no_logo_colors,
                    expect_stderr=self.expect_stderr)

        # Test abibuild
        r = env.run(self.script, "abibuild", self.loglevel, self.verbose, *no_logo_colors,
                    expect_stderr=self.expect_stderr)

    def test_with_flow(self):
        """Testing abirun.py commands with flow (no execution)"""
        env = self.get_env(check_help_version=False)
        no_logo_colors = ["--no-logo", "--no-colors"]

        # Build a flow.
        flowdir = env.base_path
        scf_input, nscf_input = make_scf_nscf_inputs()
        flow = flowtk.bandstructure_flow(flowdir, scf_input, nscf_input, manager=None)
        flow.build_and_pickle_dump()

        # Test abirun commands requiring a flow (no submission)
        for command in ["status", "debug", "deps", "inputs", "corrections", "events",
                        "history", "handlers", "cancel", "tail", "inspect"]:
            r = env.run(self.script, flowdir, command, self.loglevel, self.verbose, *no_logo_colors,
                        expect_stderr=self.expect_stderr)
            assert r.returncode == 0


class TestAbicheck(ScriptTest):
    script = os.path.join(script_dir, "abicheck.py")

    def test_abicheck(self):
        """Testing abicheck.py"""
        env = self.get_env()
        r = env.run(self.script, self.loglevel, self.verbose, expect_stderr=self.expect_stderr)

        #r = env.run(self.script, "--with-flow", self.loglevel, self.verbose, expect_stderr=self.expect_stderr)


def make_scf_nscf_inputs(paral_kgb=1):
    """Returns two input files: GS run and NSCF on a high symmetry k-mesh."""
    pseudos = abidata.pseudos("14si.pspnc")
    #pseudos = data.pseudos("Si.GGA_PBE-JTH-paw.xml")

    multi = abilab.MultiDataset(structure=abidata.cif_file("si.cif"), pseudos=pseudos, ndtset=2)
    multi.set_mnemonics(True)

    # Global variables
    ecut = 6
    global_vars = dict(ecut=ecut,
                       nband=8,
                       timopt=-1,
                       istwfk="*1",
                       nstep=15,
                       paral_kgb=paral_kgb,
                       iomode=3,
                    )

    if multi.ispaw:
        global_vars.update(pawecutdg=2*ecut)

    multi.set_vars(global_vars)

    # Dataset 1 (GS run)
    multi[0].set_kmesh(ngkpt=[8, 8, 8], shiftk=[0, 0, 0])
    multi[0].set_vars(tolvrs=1e-6)

    # Dataset 2 (NSCF run)
    kptbounds = [
        [0.5, 0.0, 0.0], # L point
        [0.0, 0.0, 0.0], # Gamma point
        [0.0, 0.5, 0.5], # X point
    ]

    multi[1].set_kpath(ndivsm=6, kptbounds=kptbounds)
    multi[1].set_vars(tolwfr=1e-12)

    # Generate two input files for the GS and the NSCF run
    scf_input, nscf_input = multi.split_datasets()
    return scf_input, nscf_input
