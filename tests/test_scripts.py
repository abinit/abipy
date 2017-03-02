# coding: utf-8
"""Test abipy command line scripts."""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import abipy.data as abidata

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

    assert not_tested == set([
        "abiGWprint.py",
        "abiGWstore.py",
        "abiGWoutput.py",
        "abiphonons.py",
        "abiGWsetup.py",
    ])


class ScriptTest(AbipyTest):
    loglevel = "--loglevel=ERROR"
    verbose = "--verbose"

    expect_stderr = True   # otherwise tests fail due to warnings and deprecation messages 

    def get_env(self):
        #import tempfile
        #env = TestFileEnvironment(tempfile.mkdtemp(suffix='', prefix='test_' + script))
        env = TestFileEnvironment()

        # Use Agg backend for plots.
        #env.writefile("matplotlibrc", "backend: Agg")
        with open(os.path.join(env.base_path, "matplotlibrc"), "wt") as fh:
            fh.write("backend: Agg\n")


        # Start with --help. If this does not work...
        env.run(self.script, "--help")

        # Script must provide a version option
        r = env.run(self.script, "--version", expect_stderr=self.expect_stderr)
        print("stderr", r.stderr)
        print("stdout", r.stdout)
        verstr = r.stderr.strip() 
        if not verstr: verstr = r.stdout.strip()  # py3k
        assert verstr == "%s version %s" % (os.path.basename(self.script), abilab.__version__)

        return env


class TestAbidoc(ScriptTest):
    script = os.path.join(script_dir, "abidoc.py")

    def test_abidoc(self):
        """Testing abidoc.py script"""
        env = self.get_env()
        env.run(self.script, "man", "ecut", self.loglevel, self.verbose, expect_stderr=self.expect_stderr)
        env.run(self.script, "apropos", "test", self.loglevel, self.verbose, expect_stderr=self.expect_stderr)
        env.run(self.script, "find", "paw", self.loglevel, self.verbose, expect_stderr=self.expect_stderr)
        env.run(self.script, "list", self.loglevel, self.verbose, expect_stderr=self.expect_stderr)
        env.run(self.script, "withdim", "natom", self.loglevel, self.verbose, expect_stderr=self.expect_stderr)


class TestAbiopen(ScriptTest):
    script = os.path.join(script_dir, "abiopen.py")

    def test_abiopen(self):
        """Testing abiopen.py script"""
        env = self.get_env()


class TestAbistruct(ScriptTest):
    script = os.path.join(script_dir, "abistruct.py")

    def test_spglib(self):
        """Testing abistruct spglib"""
        ncfile = abidata.ref_file("tgw1_9o_DS4_SIGRES.nc")
        env = self.get_env()
        env.run(self.script, "spglib", ncfile, self.loglevel, self.verbose,  
                expect_stderr=self.expect_stderr)

    def test_convert(self):
        """Testing abistruct convert"""
        ncfile = abidata.ref_file("tgw1_9o_DS4_SIGRES.nc")
        env = self.get_env()
        for fmt in ["cif", "cssr", "POSCAR", "json", "mson", "abivars"]:
            env.run(self.script, "convert", ncfile, fmt, self.loglevel, self.verbose, 
                    expect_stderr=self.expect_stderr)

    def test_abisanitize(self):
        """Testing abistruct abisanitize"""
        ncfile = abidata.ref_file("tgw1_9o_DS4_SIGRES.nc")
        env = self.get_env()
        env.run(self.script, "abisanitize", ncfile, self.loglevel, self.verbose, 
                expect_stderr=self.expect_stderr)

    def test_conventional(self):
        """Testing abistruct conventional"""
        ncfile = abidata.ref_file("tgw1_9o_DS4_SIGRES.nc")
        env = self.get_env()
        env.run(self.script, "conventional", ncfile, self.loglevel, self.verbose, 
                expect_stderr=self.expect_stderr)

    def test_kpath(self):
        """Testing abistruct kpath"""
        env = self.get_env()
        ncfile = abidata.ref_file("si_scf_WFK.nc")
        env.run(self.script, "kpath", ncfile, self.loglevel, self.verbose, 
                expect_stderr=self.expect_stderr)

    #def test_kmesh(self):
    #    """Testing abistruct kmesh"""
    #    env = self.get_env()
    #    ncfile = abidata.ref_file("tgw1_9o_DS4_SIGRES.nc")
    #    env.run(self.script, "kmesh", "--mesh=2 2 2 --shift=1 1 1 --no-time-reversal", ncfile)
    #            expect_stderr=self.expect_stderr)


class TestAbicomp(ScriptTest):
    script = os.path.join(script_dir, "abicomp.py")

    def test_abicomp(self):
        """Testing abicomp"""
        env = self.get_env()
        #env.run(self.script, "gs_scf", qtype, file1, file2, self.loglevel, self.verbose)
        #expect_stderr=self.expect_stderr)


class TestAbirun(ScriptTest):
    script = os.path.join(script_dir, "abirun.py")

    def test_manager(self):
        """Testing abirun.py manager"""
        env = self.get_env()

        no_logo_colors = ["--no-logo", "--no-colors"]

        # Test doc_manager
        env.run(self.script, ".", "doc_manager", self.loglevel, self.verbose, *no_logo_colors, 
                expect_stderr=self.expect_stderr)
        for qtype in QueueAdapter.all_qtypes():
            env.run(self.script, ".", "doc_manager", qtype, self.loglevel, self.verbose, *no_logo_colors,
                    expect_stderr=self.expect_stderr)

        # Test doc_sheduler
        env.run(self.script, ".", "doc_scheduler", self.loglevel, self.verbose, *no_logo_colors,
                expect_stderr=self.expect_stderr)


#class TestAbipsps(ScriptTest):
#    script = os.path.join(script_dir, "abipsps.py")
#
#    def test_abipsps(self):
#        """Testing abipsps.py"""
#        if not self.has_matplotlib(): return
#        si_pspnc = abidata.pseudo("14si.pspnc")
#        si_oncv = abidata.pseudo("Si.oncvpsp")
#        env = self.get_env()
#        env.run(self.script, si_pspnc.path, self.loglevel, self.verbose)
#        env.run(self.script, si_pspnc.path, si_oncv.path, self.loglevel, self.verbose,
#                expect_stderr=True)


class TestAbicheck(ScriptTest):
    script = os.path.join(script_dir, "abicheck.py")

    def test_abicheck(self):
        """Testing abicheck.py"""
        env = self.get_env()
        env.run(self.script, self.loglevel, self.verbose, expect_stderr=self.expect_stderr)


#class TestAbiinsp(ScriptTest):
#    script = os.path.join(script_dir, "abiinsp.py")
#
#    def test_abiinsp(self):
#        """Testing abiinsp.py"""
#        env = self.get_env()
