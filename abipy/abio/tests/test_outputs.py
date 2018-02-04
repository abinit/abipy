# coding: utf-8
"""Test for output files"""
from __future__ import unicode_literals, division, print_function

import os
import abipy.data as abidata

from abipy import abilab
from abipy.core.testing import AbipyTest
from abipy.abio.outputs import AbinitOutputFile, AbinitLogFile, AboRobot


class AbinitLogFileTest(AbipyTest):

    def test_abinit_logfile(self):
        """"Testing AbinitLogFile."""
        log_path = abidata.ref_file("refs/abinit.log")
        with AbinitLogFile(log_path) as abilog:
            repr(abilog); str(abilog)
            assert abilog.to_string(verbose=2)
            assert len(abilog.events) == 2
            if self.has_nbformat():
                abilog.write_notebook(nbpath=self.get_tmpname(text=True))


class AbinitOutputTest(AbipyTest):

    def test_gs_output(self):
        """Testing AbinitOutputFile with GS calculation."""
        abo_path = abidata.ref_file("refs/si_ebands/run.abo")
        with AbinitOutputFile(abo_path) as abo:
            repr(abo); str(abo)
            assert abo.to_string(verbose=2)
            assert abo.version == "8.0.6"
            assert abo.run_completed
            assert not abo.dryrun_mode
            assert abo.ndtset == 2
            assert abo.has_same_initial_structures
            assert abo.has_same_final_structures
            assert len(abo.initial_structures) == 2
            assert abo.initial_structure is not None
            assert abo.initial_structure.abi_spacegroup is not None
            assert abo.initial_structure == abo.final_structure
            abo.diff_datasets(1, 2, dryrun=True)

            # Test the parsing of dimension and spginfo
            dims_dataset, spginfo_dataset = abo.get_dims_spginfo_dataset()
            assert len(dims_dataset) == 2 and list(dims_dataset.keys()) == [1, 2]
            dims1 = dims_dataset[1]
            assert dims1["iscf"] == 7
            assert dims1["nfft"] == 5832
            self.assert_almost_equal(dims1["mem_per_proc_mb"], 3.045)
            self.assert_almost_equal(dims1["wfk_size_mb"], 0.717)
            self.assert_almost_equal(dims1["denpot_size_mb"], 0.046)
            assert spginfo_dataset[1]["spg_symbol"] == "Fd-3m"
            assert spginfo_dataset[1]["spg_number"] == 227
            assert spginfo_dataset[1]["bravais"] == "Bravais cF (face-center cubic)"
            dims2 = dims_dataset[2]
            assert dims2["iscf"] == -2
            assert dims2["n1xccc"] == 2501
            self.assert_almost_equal(dims2["mem_per_proc_mb"], 1.901)
            self.assert_almost_equal(dims2["wfk_size_mb"], 0.340)
            self.assert_almost_equal(dims2["denpot_size_mb"], 0.046)

            str(abo.events)
            gs_cycle = abo.next_gs_scf_cycle()
            assert gs_cycle is not None
            if self.has_matplotlib():
                assert gs_cycle.plot(show=False)
            abo.seek(0)
            assert abo.next_d2de_scf_cycle() is None

            timer = abo.get_timer()
            assert len(timer) == 1
            assert str(timer.summarize())

            if self.has_matplotlib():
                abo.compare_gs_scf_cycles([abo_path], show=False)
                timer.plot_all(show=False)
                abo.plot(show=False)

            if self.has_nbformat():
                abo.write_notebook(nbpath=self.get_tmpname(text=True))
                timer.write_notebook(nbpath=self.get_tmpname(text=True))

    def test_ph_output(self):
        """Testing AbinitOutputFile with phonon calculations."""
        abo_path = abidata.ref_file("refs/gs_dfpt.abo")
        with AbinitOutputFile(abo_path) as abo:
             repr(abo); str(abo)
             assert abo.to_string(verbose=2)

             assert abo.version == "8.3.2"
             assert abo.run_completed
             assert not abo.dryrun_mode
             assert abo.ndtset == 3
             assert abo.has_same_initial_structures
             assert abo.has_same_final_structures
             assert len(abo.initial_structures) == 3
             assert abo.initial_structure is not None
             assert abo.initial_structure.abi_spacegroup is not None
             assert abo.initial_structure == abo.final_structure

             gs_cycle = abo.next_gs_scf_cycle()
             assert gs_cycle is not None
             ph_cycle = abo.next_d2de_scf_cycle()
             assert ph_cycle is not None
             if self.has_matplotlib():
                assert ph_cycle.plot(show=False)
                assert abo.compare_d2de_scf_cycles([abo_path], show=False)
                abo.plot(show=False)

             if self.has_nbformat():
                abo.write_notebook(nbpath=self.get_tmpname(text=True))

    def test_dryrun_output(self):
        """Testing AbinitOutputFile with file produced in dry-run mode."""
        with abilab.abiopen(abidata.ref_file("refs/dryrun.abo")) as abo:
            repr(abo); str(abo)
            assert abo.to_string(verbose=2)
            assert abo.dryrun_mode
            assert abo.ndtset == 1
            assert abo.has_same_initial_structures
            assert abo.has_same_final_structures
            assert len(abo.initial_structures) == 1

            assert abo.initial_structure.abi_spacegroup is not None

            # This to test get_dims_spginfo_dataset with one dataset.
            dims_dataset, spg_dataset = abo.get_dims_spginfo_dataset()
            assert len(dims_dataset) == 1
            dims = dims_dataset[1]
            assert dims["nsppol"] == 1
            assert dims["nsym"] == 48
            assert dims["nkpt"] == 29
            self.assert_almost_equal(dims["mem_per_proc_mb"], 3.389)
            self.assert_almost_equal(dims["wfk_size_mb"], 0.717)
            self.assert_almost_equal(dims["denpot_size_mb"], 0.046)
            assert spg_dataset[1]["spg_symbol"] == "Fd-3m"
            assert spg_dataset[1]["spg_number"] == 227
            assert spg_dataset[1]["bravais"] == "Bravais cF (face-center cubic)"

    def test_abinit_output_with_ctrlm(self):
        """Testing AbinitOutputFile with file containing CTRL+M char."""
        test_dir = os.path.join(os.path.dirname(__file__), "..", "..", 'test_files')
        with abilab.abiopen(os.path.join(test_dir, "ctrlM_run.abo")) as abo:
            assert abo.version == "8.7.1"
            assert abo.run_completed
            assert abo.to_string(verbose=2)
            assert abo.ndtset == 1
            assert abo.initial_structure.abi_spacegroup is not None
            assert abo.initial_structure.abi_spacegroup.spgid == 142
            assert abo.proc0_cputime == 0.7
            assert abo.proc0_walltime == 0.7
            assert abo.overall_cputime == 0.7
            assert abo.overall_walltime == 0.7

            # Test the parsing of dimension and spginfo
            dims_dataset, spginfo_dataset = abo.get_dims_spginfo_dataset()
            dims1 = dims_dataset[1]
            assert dims1["mqgrid"] == 5580
            assert spginfo_dataset[1]["spg_symbol"] == "I4_1/acd"
            assert spginfo_dataset[1]["spg_number"] == 142

    def test_all_outputs_in_tests(self):
        """
        Try to parse all Abinit output files inside the Abinit `tests` directory.
        Requires $ABINIT_HOME_DIR env variable.
        """
        abi_homedir = os.environ.get("ABINIT_HOME_DIR")
        if abi_homedir is not None:
            #raise self.SkipTest("Environment variable `ABINIT_HOME_DIR` is required for this test.")
            abitests_dir = os.path.join(abi_homedir, "tests")
        else:
            abitests_dir = os.path.join(abidata.dirpath, "refs")

        from abipy.abio.outputs import validate_output_parser
        assert os.path.exists(abitests_dir)
        retcode = validate_output_parser(abitests_dir=abitests_dir)
        assert retcode == 0

    def test_aborobot(self):
        """Testing AboRobot."""
        abo_paths = abidata.ref_files("refs/si_ebands/run.abo", "refs/gs_dfpt.abo")
        with AboRobot.from_files(abo_paths) as robot:
            repr(robot); str(robot)
            assert robot.to_string(verbose=2)
            assert robot._repr_html_()
            dims = robot.get_dims_dataframe()
            df = robot.get_dataframe(with_geo=True)
            time_df = robot.get_time_dataframe()
            self.assert_equal(time_df["overall_walltime"].values, [4.0, 26.1])

            if self.has_nbformat():
                robot.write_notebook(nbpath=self.get_tmpname(text=True))
