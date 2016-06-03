# coding: utf-8
"""Test abipy command line scripts."""
from __future__ import print_function, division, unicode_literals, absolute_import


import os
import abipy.data as abidata

from scripttest import TestFileEnvironment
from monty.inspect import all_subclasses
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
        "abibatch.py",
        "mrgddb.py",
        "cif2spg.py",
        "abiGWprint.py",
        "abiGWstore.py",
        "abiGWoutput.py",
        "abiphonons.py",
        "abiGWsetup.py",])


class ScriptTest(AbipyTest):
    loglevel = "--loglevel=ERROR"

    def get_env(self):
        #import tempfile
        #env = TestFileEnvironment(tempfile.mkdtemp(suffix='', prefix='test_' + script))
        env = TestFileEnvironment()

        # Use Agg backend for plots.
        env.writefile("matplotlibrc", "backend : Agg")

        # Start with --help. If this does not work...
        env.run(self.script, "--help")

        # Script must provide a version option
        r = env.run(self.script, "--version", expect_stderr=True)
        assert r.stderr.strip() == "%s version %s" % (os.path.basename(self.script), abilab.__version__)

        return env


class TestAbidoc(ScriptTest):
    script = os.path.join(script_dir, "abidoc.py")

    def test_abidoc(self):
        """Testing abidoc.py script"""
        env = self.get_env()
        env.run(self.script, self.loglevel, "man", "ecut")
        env.run(self.script, self.loglevel, "apropos", "test")
        env.run(self.script, self.loglevel, "find", "paw")
        env.run(self.script, self.loglevel, "list")


class TestAbilab(ScriptTest):
    script = os.path.join(script_dir, "abilab.py")

    def test_abidoc(self):
        """Testing abilab.py script"""
        env = TestFileEnvironment()
        env.run(self.script, "--help", expect_stderr=True)
        env.run(self.script, "--version", expect_stderr=True)


class TestAbiopen(ScriptTest):
    script = os.path.join(script_dir, "abiopen.py")

    def test_abiopen(self):
        """Testing abiopen.py script"""
        env = self.get_env()


class TestAbistruct(ScriptTest):
    script = os.path.join(script_dir, "abistruct.py")

    def test_convert(self):
        """Testing abistruct convert"""
        ncfile = abidata.ref_file("tgw1_9o_DS4_SIGRES.nc")
        env = self.get_env()
        for fmt in ["cif", "cssr", "POSCAR", "json", "mson",]:
            env.run(self.script, self.loglevel, "convert", ncfile, fmt)

    def test_bz(self):
        """Testing abistruct bz"""
        env = self.get_env()
        ncfile = abidata.ref_file("tgw1_9o_DS4_SIGRES.nc")
        #env.run(self.script, self.loglevel, "bz", ncfile)


class TestAbidiff(ScriptTest):
    script = os.path.join(script_dir, "abidiff.py")

    def test_abidiff(self):
        """Testing abidiff"""
        env = self.get_env()


class TestAbirun(ScriptTest):
    script = os.path.join(script_dir, "abirun.py")

    def test_manager(self):
        """Testing abirun.py manager"""
        env = self.get_env()
        env.run(self.script, self.loglevel, ".", "manager")
        for qtype in ['shell', 'slurm', 'pbspro', 'sge', 'moab', 'bluegene', 'torque']:
            env.run(self.script, self.loglevel, ".", "manager", qtype)


class TestAbipsps(ScriptTest):
    script = os.path.join(script_dir, "abipsps.py")

    def test_abipsps(self):
        """Testing abipsps.py"""
        if not self.has_matplotlib(): return
        si_pspnc = abidata.pseudo("14si.pspnc")
        si_oncv = abidata.pseudo("Si.oncvpsp")
        env = self.get_env()
        env.run(self.script, self.loglevel, si_pspnc.path)
        env.run(self.script, self.loglevel, si_pspnc.path, si_oncv.path, expect_stderr=True)


#class TestMrgddb(ScriptTest):
#    script = os.path.join(script_dir, "mrgdbb.py")


#class TestAbibatch(ScriptTest):
#    script = os.path.join(script_dir, "abibatch.py")


class TestAbiq(ScriptTest):
    script = os.path.join(script_dir, "abiq.py")
    def test_abiq(self):
        """Testing abiq.py"""
        env = self.get_env()


class TestAbitime(ScriptTest):
    script = os.path.join(script_dir, "abitime.py")

    def test_abitime(self):
        """Testing abitime.py"""
        env = self.get_env()


class TestAbicheck(ScriptTest):
    script = os.path.join(script_dir, "abicheck.py")

    def test_abicheck(self):
        """Testing abicheck.py"""
        env = self.get_env()
        env.run(self.script, self.loglevel)


class TestAbiinsp(ScriptTest):
    script = os.path.join(script_dir, "abiinsp.py")

    def test_abiinsp(self):
        """Testing abiinsp.py"""
        env = self.get_env()
