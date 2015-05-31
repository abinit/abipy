# coding: utf-8
from __future__ import unicode_literals, division, print_function

import sys
import os
import mock
import yaml
import abipy.data as abidata
import abipy.abilab as abilab
import abipy.fworks.fw_tasks as fw_tasks

from abipy.abio.factories import *
from abipy.core.testing import AbipyTest
import mock_objects

from pymatgen.io.abinitio.events import Correction, DilatmxErrorHandler, DilatmxError

from fireworks import FWAction


class TestAbiFireTask(AbipyTest):

    def setUp(self):
        si = abilab.Structure.from_file(abidata.cif_file("si.cif"))
        self.si_scf_input = ebands_input(si, abidata.pseudos("14si.pspnc"), ecut=2, kppa=10).split_datasets()[0]

    def test_AbiFireTask(self):
        task = fw_tasks.AbiFireTask(self.si_scf_input)
        task.to_dict()
        self.assertPMGSONable(self.si_scf_input)
        self.assertFwSerializable(task)

    def test_ScfFireTask(self):
        task = fw_tasks.ScfFWTask(self.si_scf_input)
        task.to_dict()
        self.assertFwSerializable(task)


class TestTaskAnalysis(AbipyTest):
    def setUp(self):
        si = abilab.Structure.from_file(abidata.cif_file("si.cif"))
        self.si_scf_input = ebands_input(si, abidata.pseudos("14si.pspnc"), ecut=2, kppa=10).split_datasets()[0]

    def tearDown(self):
        os.remove(fw_tasks.HISTORY_JSON)

    @mock.patch.object(fw_tasks.AbiFireTask, 'get_event_report')
    def test_scf_unconverged(self, report):
        scf_task = fw_tasks.ScfFWTask(self.si_scf_input)
        scf_task.ftm = fw_tasks.FWTaskManager.from_user_config({})

        report.return_value = mock_objects.report_ScfConvergenceWarning()

        with mock.patch.object(fw_tasks.AbiFireTask, 'prepare_restart', return_value=(mock_objects.fake_fw, {})) as pr:
            fake_spec = {'test': 1}
            action = scf_task.task_analysis(fake_spec)
            pr.assert_called_once_with(fake_spec)
            self.assertIsInstance(action, FWAction)

            scf_task.restart_info = fw_tasks.RestartInfo(
                previous_dir='.', num_restarts=fw_tasks.FWTaskManager.fw_policy_defaults['max_restarts'])
            with self.assertRaises(fw_tasks.UnconvergedError):
                scf_task.task_analysis(fake_spec)

    @mock.patch.object(fw_tasks.AbiFireTask, 'get_event_report')
    def test_generic_error(self, report):
        scf_task = fw_tasks.ScfFWTask(self.si_scf_input)

        report.return_value = mock_objects.report_AbinitError()

        fake_spec = {'test': 1}
        with self.assertRaises(fw_tasks.AbinitRuntimeError):
            # set the returncode to avoid logging problems
            scf_task.returncode = 10
            scf_task.task_analysis(fake_spec)

    @mock.patch.object(fw_tasks.AbiFireTask, 'get_event_report')
    def test_no_report_no_err(self, report):
        scf_task = fw_tasks.ScfFWTask(self.si_scf_input)
        scf_task.set_workdir()

        report.return_value = None

        fake_spec = {'test': 1}
        with self.assertRaises(fw_tasks.AbinitRuntimeError):
            # set the returncode to avoid logging problems
            scf_task.returncode = 10
            scf_task.task_analysis(fake_spec)


class TestObjects(AbipyTest):

    def test_task_history_and_events(self):
        th = fw_tasks.TaskHistory()
        th.log_autoparal({u'time': u'12:0:0', u'ntasks': 15, u'partition': 'defq', u'nodes': 1, u'mem_per_cpu': 3000})
        th.log_finalized()
        th.log_restart(fw_tasks.RestartInfo(os.path.abspath('.'), reset=True, num_restarts=2))
        th.log_concluded()
        th.log_unconverged()
        th.log_corrections([Correction(DilatmxErrorHandler(), {}, DilatmxError('', '', '',), )])

        self.assertPMGSONable(th)

        for te in th:
            self.assertPMGSONable(te)


class TestFWTaskManager(AbipyTest):
    def tearDown(self):
        try:
            os.remove(fw_tasks.FWTaskManager.YAML_FILE)
        except OSError:
            pass

    def test_no_file(self):
        fw_tasks.FWTaskManager.from_user_config({})

    def test_ok(self):
        with open(os.path.abspath(fw_tasks.FWTaskManager.YAML_FILE), 'w') as f:
            f.write(mock_objects.MANAGER_OK)

        ftm = fw_tasks.FWTaskManager.from_user_config({'fw_policy': {'max_restarts': 30}})

        self.assertTrue(ftm.fw_policy.rerun_same_dir)
        self.assertEqual(ftm.fw_policy.max_restarts, 30)
        self.assertTrue(ftm.fw_policy.autoparal)

    def test_no_qadapter(self):
        with open(os.path.abspath(fw_tasks.FWTaskManager.YAML_FILE), 'w') as f:
            f.write(mock_objects.MANAGER_NO_QADAPTERS)

        ftm = fw_tasks.FWTaskManager.from_user_config({})

        self.assertIsNone(ftm.task_manager)

    def test_unknown_keys(self):
        with open(os.path.abspath(fw_tasks.FWTaskManager.YAML_FILE), 'w') as f:
            f.write(mock_objects.MANAGER_UNKNOWN_KEYS)

        with self.assertRaises(RuntimeError):
            fw_tasks.FWTaskManager.from_user_config({})

