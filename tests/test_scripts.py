# coding: utf-8
"""
"""
from __future__ import print_function, division, unicode_literals


import os
import abipy.data as abidata  

from scripttest import TestFileEnvironment
from monty.inspect import all_subclasses
from abipy.core.testing import AbipyTest
from abipy import abilab

script_dir = os.path.abspath("../abipy/scripts/")


def test_if_all_scripts_are_tested():
    """Testing if all scripts are tested"""
    tested_scripts = set(c.script for c in all_subclasses(ScriptTest))
    all_scripts = set(f for f in os.listdir(script_dir) if f.endswith(".py"))
    not_tested = all_scripts.difference(tested_scripts)

    if not_tested:
        print("The following scripts are not tested")
        for i, s in enumerate(not_tested):
            print("[%d] %s" % (i, s))
    #assert not not_tested


class ScriptTest(AbipyTest):
    loglevel = "--loglevel=DEBUG"

    def get_env(self):
        #import tempfile
        #env = TestFileEnvironment(tempfile.mkdtemp(suffix='', prefix='test_' + script))
        env = TestFileEnvironment()

        # Use Agg backend for plots.
        env.writefile("matplotlibrc", "backend : Agg")
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
        #env.run(self.script, self.loglevel, "man", "ecut")
        #env.run(self.script, self.loglevel, "apropos", "test")
        #env.run(self.script, self.loglevel, "find", "paw")
        #env.run(self.script, self.loglevel, "list")


class TestAbiopen(ScriptTest):
    script = os.path.join(script_dir, "abiopen.py")


class TestAbistruct(ScriptTest):
    script = os.path.join(script_dir, "abistruct.py")

    def test_convert(self):
        """Testing `abistruct convert`"""
        ncfile = abidata.ref_file("tgw1_9o_DS4_SIGRES.nc")
        env = self.get_env()
        for fmt in ["cif", "cssr", "POSCAR", "json", "mson",]:
            env.run(self.script, self.loglevel, "convert", ncfile, fmt)

    def test_bz(self):
        """Testing `abistruct bz`"""
        env = self.get_env()
        ncfile = abidata.ref_file("tgw1_9o_DS4_SIGRES.nc")
        #env.run(self.script, self.loglevel, "bz", ncfile)


class TestAbidiff(ScriptTest):
    script = os.path.join(script_dir, "abidiff.py")


class TestAbidiff(ScriptTest):
    script = os.path.join(script_dir, "abirun.py")

    def test_manager(self):
        """Test `abirun.py manager`"""
        env = self.get_env()
        env.run(self.script, self.loglevel, ".", "manager")
        for qtype in ['shell', 'slurm', 'pbspro', 'sge', 'moab', 'bluegene', 'torque']:
            env.run(self.script, self.loglevel, ".", "manager", qtype)


class TestAbipsps(ScriptTest):
    script = os.path.join(script_dir, "abipsps.py")

    def test_abipsps(self):
        """Test `abipsps.py`"""
        if not self.has_matplotlib(): return
        si_pspnc = abidata.pseudo("14si.pspnc")
        si_oncv = abidata.pseudo("Si.oncvpsp")
        env = self.get_env()
        #env.run(self.script, self.loglevel, si_pspnc.path)
        #env.run(self.script, self.loglevel, si_pspnc.path, si_oncv.path)


#class TestMrgddb(ScriptTest):
#    script = os.path.join(script_dir, "mrgdbb.py")


#class TestAbibatch(ScriptTest):
#    script = os.path.join(script_dir, "abibatch.py")


class TestAbiq(ScriptTest):
    script = os.path.join(script_dir, "abiq.py")
    def test_abiq(self):
        """Test `abiq.py`"""
        env = self.get_env()


class TestAbitime(ScriptTest):
    script = os.path.join(script_dir, "abitime.py")


class TestAbicheck(ScriptTest):
    script = os.path.join(script_dir, "abicheck.py")
    def test_abicheck(self):
        env = self.get_env()
        env.run(self.script, self.loglevel)


class TestAbiinsp(ScriptTest):
    script = os.path.join(script_dir, "abiinsp.py")



if __name__ == "__main__":
    import unittest
    unittest.main()
