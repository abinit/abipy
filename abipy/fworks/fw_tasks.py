# coding: utf-8
"""
Task classes for Fireworks.
"""
from __future__ import print_function, division, unicode_literals

#TODO check if this raises exceptions in environments without fireworks
from fireworks.core.firework import Firework, FireTaskBase, FWAction
from fireworks.utilities.fw_utilities import explicit_serialize
from fireworks.utilities.fw_serializers import serialize_fw

import abipy.data as abidata
import abipy.abilab as abilab
from pymatgen.io.abinitio.tasks import AbinitTask, ScfTask, NscfTask, RelaxTask, DdkTask, PhononTask, \
    SigmaTask, BseTask, OpticTask, AnaddbTask
from pymatgen.io.abinitio.tasks import FileNode
import numpy as np
import time
import re
import logging


class AbinitRuntimeError(Exception):
    """
    Exception raised for errors during Abinit calculation
    Initialized with a task, uses it to prepare a suitable error message
    """
    def __init__(self, task):
        self.task = task

    def __str__(self):
        return "Final status: {}. Hystory: {}".format(self.task.status, tuple(self.task.history))


@explicit_serialize
class AbiFireTaskBase(FireTaskBase):
    """
    Base class to create abinit FireTasks.
    Provides the basic methods needed to serialize and deserialize the task and the
    structure for the run_task method.
    """

    @serialize_fw
    def to_dict(self):
        """
        Explicit conversion of all the elements needed to reconstruct the Task.
        Does not rely on the recursive serialization.
        """
        workdir = self.workdir

        abinit_input = self.strategy.abinit_input
        pseudos_name_list = []
        for pseudo in abinit_input.pseudos:
            pseudos_name_list.append({'name': pseudo.name})
        abinit_input_dict = {'pseudos': pseudos_name_list}
        dtsets = []
        for ds in abinit_input:
            ds_copy = ds.deepcopy()
            for key, value in ds_copy.iteritems():
                # Convert to list to allow serialization
                if isinstance(value, np.ndarray):
                    ds_copy[key] = value.tolist()
            dtsets.append(dict(ds_copy))
        abinit_input_dict['datasets'] = dtsets

        strategy = {'abinit_input': abinit_input_dict,
                    '@module': self.strategy.__class__.__module__,
                    '@class': self.strategy.__class__.__name__}

        # All the dependencies are saved using the path as node
        deps = []
        for dep in self.deps:
            for path, ext in zip(*dep.get_filepaths_and_exts()):
                deps.append({'path': path, 'ext': ext})

        d = dict(workdir=workdir, strategy=strategy,
                 deps=deps, required_files=None,
                 status=str(self.status))

        return d

    # This are needed since during recursive de/serialization Fireworks sometimes
    # deals with the Task object as a dict (FireTaskBase is a defaultdict).
    # Do not implement __iter__, since it will mess up with the abipy code
    def __getitem__(self, item):
        return self.to_dict()[item]

    def items(self):
        return self.to_dict().items()

    @classmethod
    def _abiinput_from_dict(cls, m_dict):

        dtsets = m_dict['datasets']
        pseudos_dict = m_dict['pseudos']
        pseudos = []
        #TODO improve the way pseudos are retrieved. environment variable?
        for pseudo in pseudos_dict:
            pseudos.append(abidata.pseudo(pseudo['name']))
        abiintput = abilab.AbiInput(pseudos, ndtset=dtsets[0]['ndtset'])
        n = 0
        for ds in dtsets:
            #TODO save and set explicitly the structure object?
            abiintput.set_variables(dtset=n, **ds)
            n += 1

        return abiintput

    @classmethod
    def _basic_from_dict(cls, m_dict):
        # Load strategy
        strategy_dict = m_dict['strategy']
        strategy_cls_mod = strategy_dict['@module']
        strategy_cls_name = strategy_dict['@class']
        strategy_module = __import__(strategy_cls_mod, globals(), locals(), [strategy_cls_name])
        strategy_cls = getattr(strategy_module, strategy_cls_name)
        #TODO this will work just for StrategyWithInput.
        strategy = strategy_cls(cls._abiinput_from_dict(strategy_dict['abinit_input']))

        # Load dependencies using the FileNode as Node
        # The Task constructor requires a dict
        deps_list = m_dict['deps']
        deps = {}
        for dep in deps_list:
            deps[FileNode(dep['path'])] = dep['ext']

        return strategy, deps

    @classmethod
    def from_dict(cls, m_dict):
        """
        Explicit conversion to a Task object from the dict saved with to_dict
        Does not rely on the recursive serialization.
        """
        strategy, deps = cls._basic_from_dict(m_dict)

        task = cls(strategy, workdir=m_dict['workdir'],
                   manager=abilab.TaskManager.from_user_config(), deps=deps,
                   required_files=m_dict['required_files'])

        return task

    def _update_abinit_vars(self, fw_spec):
        """
        Internal method that retrieves the abinit variables to be updated from the spec
        """
        if 'updated_abinit_vars' in fw_spec:
            self.strategy.add_extra_abivars(fw_spec['updated_abinit_vars'])

    def _update_task_manager(self, fw_spec):
        """
        Internal method to update the task manager according to the _queueadapter spec
        """
        if '_queueadapter' in fw_spec:
            if 'ntasks' in fw_spec['_queueadapter']:
                self.manager.set_mpi_ncpus(fw_spec['_queueadapter']['ntasks'])

    def config_run(self, fw_spec):
        logging.basicConfig(filename='abipy.log', level=logging.WARNING)
        self._update_task_manager(fw_spec)
        self._update_abinit_vars(fw_spec)

    def run_abitask(self):
        """
        Runs the single task and performs the check_status check regularly.
        """
        # TODO add other checks?
        started = self.start()
        if not started:
            raise RuntimeError("Task did not start. Status: {}".format(self.status))

        #TODO get the values as external parameters. from the taskmanagerl.yml?
        self.sleep_interval = 5
        self.check_interval = 20
        check_skip = self.check_interval // self.sleep_interval
        # The check is called every check_skip intervals.
        # If the job has terminated a last check is performed
        job_terminated = False
        while not job_terminated:
            for n in range(check_skip):
                job_terminated = self.process.poll() is not None
                if job_terminated:
                    break
                time.sleep(self.sleep_interval)
            self.check_status()

    def run_task(self, fw_spec):
        """
        Implementation of Fireworks run_task. Configures, executes and analyzes
        the job execution in three separate steps
        """
        self.config_run(fw_spec)

        self.run_abitask()

        #TODO:
        # -check for errors
        # -add more info about the task
        # -remove lock file?
        # -are there other actions to be performed in abipy at the end of task?

        return self.task_analysis(fw_spec)

    #TODO find a better name. Set it to abstract?
    def task_analysis(self, fw_spec):
        """
        Analyzes the final status of the task and prepares the action.
        It is called at the end of run_task and should return a FWaction.
        Can raise an exception to fizzle the firework

        It can be overwritten by subclasses to characterize the action for each
        particular task.
        In this implementation some general information are added to the FWaction.
        """

        # In the generic case, raise an exception if the final status is not OK.
        if self.status < self.S_OK:
            raise AbinitRuntimeError(self)

        stored_data = {'history': list(self.history)}

        return FWAction(stored_data=stored_data)


#########################################################################################
# AbiFlow derived tasks
#########################################################################################

@explicit_serialize
class ScfFWTask(ScfTask, FireTaskBase):
    pass


@explicit_serialize
class NscfFWTask(NscfTask, FireTaskBase):
    pass


@explicit_serialize
class RelaxFWTask(RelaxTask, AbiFireTaskBase):
    def run_task(self, fw_spec):
        if 'updated_structure' in fw_spec:
            updated_structure = abilab.Structure.from_dict(fw_spec['updated_structure'])
            self.change_structure(updated_structure)

        return super(RelaxFWTask, self).run_task(fw_spec)

    def task_analysis(self, fw_spec):
        """
        A relax task updates forwards an updated structure for the following tasks/FWs
        """

        # Raise an exception if the final status is not OK.
        if self.status < self.S_OK:
            raise AbinitRuntimeError(self)

        stored_data = {'history': list(self.history)}

        update_spec = {'updated_structure': self.read_final_structure().as_dict()}

        return FWAction(stored_data=stored_data, update_spec=update_spec)


@explicit_serialize
class DdkFWTask(DdkTask, FireTaskBase):
    pass


@explicit_serialize
class PhononFWTask(PhononTask, FireTaskBase):
    pass


@explicit_serialize
class SigmaFWTask(SigmaTask, FireTaskBase):
    pass


@explicit_serialize
class BseFWTask(BseTask, FireTaskBase):
    pass


@explicit_serialize
class OpticFWTask(OpticTask, FireTaskBase):
    pass


@explicit_serialize
class AnaddbFWTask(AnaddbTask, FireTaskBase):
    pass

#########################################################################################
# Additional tasks
#########################################################################################


@explicit_serialize
class MultiStepRelaxFWTask(RelaxFWTask):
    """
    This class derives from RelaxFWTask, but allows a multi-FW relaxation.
    If the final status is not converged, it creates a new FW to continue the relaxation.

    This can be useful to split long relaxations in short jobs.
    """

    def __init__(self, *args, **kwargs):
        self.max_steps = kwargs.pop('max_steps', 5)
        super(MultiStepRelaxFWTask, self).__init__(*args, **kwargs)

    def to_dict(self):
        d = super(MultiStepRelaxFWTask, self).to_dict()
        d['max_steps'] = self.max_steps

        return d

    @classmethod
    def from_dict(cls, m_dict):
        strategy, deps = cls._basic_from_dict(m_dict)

        task = cls(strategy, workdir=m_dict['workdir'],
                   manager=abilab.TaskManager.from_user_config(), deps=deps,
                   required_files=m_dict['required_files'], max_steps=m_dict['max_steps'])

        return task

    # Generates the directory for the new task
    def new_workdir(self):
        if re.search(r'relaxstep\d+$', self.workdir):
            num = int(re.search(r'\d+$', self.workdir).group())
            return re.sub(r'\d+$', str(num+1), self.workdir)
        else:
            return self.workdir+'_relaxstep2'

    def task_analysis(self, fw_spec):
        """
        A relax task updates forwards an updated structure for the following tasks/FWs.

        If the status if Unconverged does raise an exception, but creates a new FW.
        """

        # Raise an exception if the final status is not Unconverged or OK.
        if self.status < self.S_UNCONVERGED or self.status == self.S_ERROR:
            raise AbinitRuntimeError(self)

        stored_data = {'history': list(self.history)}
        updated_structure = self.read_final_structure()

        if self.status == self.S_UNCONVERGED:
            if self.max_steps <= 1:
                raise AbinitRuntimeError(self)
            new_task = MultiStepRelaxFWTask(
                self.strategy.deepcopy(), workdir=self.new_workdir(),
                deps={self: "WFK"}, max_steps=self.max_steps-1)
            new_task.change_structure(updated_structure)
            new_step = Firework(new_task)
            return FWAction(stored_data=stored_data, detours=new_step)

        update_spec = {'updated_structure': updated_structure.as_dict()}
        return FWAction(stored_data=stored_data, update_spec=update_spec)


@explicit_serialize
class AutoparalFWTask(AbinitTask, AbiFireTaskBase):
    """
    Task to execute the autoparal method of abinit and apply the results to the following
    tasks in the workflow.
    """

    #TODO maybe the init should be modified to be independent of the task type
    def __init__(self, *args, **kwargs):
        self.max_ncpus = kwargs.pop('max_ncpus')
        super(AutoparalFWTask, self).__init__(*args, **kwargs)
        self.manager.policy.autoparal = 1
        if self.max_ncpus:
            self.manager.policy.max_ncpus = self.max_ncpus

    def to_dict(self):
        d = super(AutoparalFWTask, self).to_dict()
        d['max_ncpus'] = self.max_ncpus

        return d

    @classmethod
    def from_dict(cls, m_dict):
        strategy, deps = cls._basic_from_dict(m_dict)
        task = cls(strategy, workdir=m_dict['workdir'],
                   manager=abilab.TaskManager.from_user_config(), deps=deps,
                   required_files=m_dict['required_files'],
                   max_ncpus=m_dict['max_ncpus'])

        return task

    @classmethod
    def from_task(cls, task, max_ncpus=None):
        """
        Creates the AutoparalFWTask starting from the task that will be run.
        If max_ncpus is None the one provided by the local taskmanager is used.
        """
        workdir = task.workdir + '_autoparal_fake_run'
        return cls(task.strategy.deepcopy(), workdir=workdir,
                   manager=abilab.TaskManager.from_user_config(), deps=task.deps,
                   required_files=task.required_files, max_ncpus=max_ncpus)

    def run_task(self, fw_spec):
        """
        Runs the autoparal_fake_run method and updates queuadapter information
        for the following tasks.
        """
        #TODO add a try?
        confs, optimal = self.autoparal_fake_run()

        #TODO this will probably need to become a mod_spec
        update_spec = {'_queueadapter': {'ntasks': optimal.mpi_ncpus},
                       'updated_abinit_vars': optimal.vars}

        return FWAction(update_spec=update_spec)