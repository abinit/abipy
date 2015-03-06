# coding: utf-8
"""
Task classes for Fireworks.
"""
from __future__ import print_function, division, unicode_literals

try:
    from fireworks.core.firework import Firework, FireTaskBase, FWAction
    from fireworks.utilities.fw_utilities import explicit_serialize
    from fireworks.utilities.fw_serializers import serialize_fw
except ImportError:
    FireTaskBase, FWAction, Firework = 3 * [object]
    explicit_serialize = lambda x: x
    serialize_fw = lambda x: x

import abipy.data as abidata
import abipy.abilab as abilab
from pymatgen.io.abinitio.abiobjects import KSampling, RelaxationMethod
from pymatgen.io.abinitio.tasks import AbinitTask, ScfTask, NscfTask, RelaxTask, DdkTask, PhononTask, \
    SigmaTask, BseTask, OpticTask, AnaddbTask, TaskRestartError
from pymatgen.io.abinitio.nodes import Node, FileNode, Dependency
from pymatgen.io.abinitio.pseudos import PseudoParser, PseudoTable
from pymatgen.io.abinitio.strategies import ScfStrategy, NscfStrategy, RelaxStrategy, StrategyWithInput
import time
import re
import glob
import logging
import os
import inspect
from monty.json import MontyEncoder, MontyDecoder
from functools import wraps
import json

logger = logging.getLogger(__name__)


class ErrorCode(object):
    UNRECOVERABLE = 'Unrecoverable'
    UNCLASSIFIED = 'Unclassified'
    UNCONVERGED = 'Unconverged'
    INITIALIZATION = 'Initialization'


class AbinitRuntimeError(Exception):
    """
    Exception raised for errors during Abinit calculation
    Initialized with a task, uses it to prepare a suitable error message
    """
    def __init__(self, task):
        self.task = task

    def __str__(self):
        return "Final status: {}. History: {}".format(self.task.abitask.status, tuple(self.task.abitask.history))

    def to_dict(self):
        report = self.task.abitask.get_event_report()
        d = {'@class': self.__class__.__name__}
        d['num_errors'] = report.num_errors
        d['num_warnings'] = report.num_warnings
        if report.num_errors:
            errors = []
            for error in report.errors:
                errors.append({'name': error.name, 'message': error.message})
            d['errors'] = errors
        if report.num_warnings:
            warnings = []
            for warning in report.warnings:
                warnings.append({'name': warning.name, 'message': warning.message})
            d['warnings'] = warnings

        if self.task.abitask.status == self.task.abitask.S_UNCONVERGED:
            d['error_code'] = ErrorCode.UNCONVERGED
        elif report.num_errors > 0:
            d['error_code'] = ErrorCode.UNRECOVERABLE
        else:
            d['error_code'] = ErrorCode.UNCLASSIFIED

        return d

class InitializationError(Exception):
    def __init__(self, task, msg):
        self.task = task
        self.msg = msg

    def __str__(self):
        return "Error during task initilization. Status: {}. {}".format(
            self.task.abitask.status, self.msg)

    def to_dict(self):
        return {'@class': self.__class__.__name__, 'task_status': self.task.abitask.status,
                'message': self.msg, 'error_code': ErrorCode.INITIALIZATION}


@explicit_serialize
class AbiFireTask(FireTaskBase):
    """
    Base class to create abinit FireTasks.
    Provides the basic methods needed to serialize and deserialize the task and the
    structure for the run_task method.
    """

    def __init__(self, abitask):
        self.abitask = abitask

    @serialize_fw
    def to_dict(self):
        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__}
        if hasattr(self, "__init__"):
            for c in inspect.getargspec(self.__init__).args:
                if c != "self":
                    a = self.__getattribute__(c)
                    if hasattr(a, "as_dict"):
                        a = a.as_dict()
                    d[c] = a
            va = inspect.getargspec(self.__init__).varargs
            # if va:
            #     d[va] = [a.as_dict() if hasattr(a, "as_dict") else a for a in self.__getattribute__(va)]
            kwa = inspect.getargspec(self.__init__).keywords
            if kwa:
                d[kwa] = {k: v.as_dict() if hasattr(v, "as_dict") else v for k, v in self.__getattribute__(kwa).items()}

        return d

    @classmethod
    def from_dict(cls, d):
        """
        This method should return the ErrorHandler from a dict representation
        of the object given by the to_dict property.
        """
        dec = MontyDecoder()
        kwargs = {k: dec.process_decoded(v) for k, v in d.items()
                  if k in inspect.getargspec(cls.__init__).args}
        # varargs = [dec.process_decoded(c) for c in d.get(inspect.getargspec(cls.__init__).varargs, [])]
        kwargs.update({kk: dec.process_decoded(vv)
                      for kk, vv in d.get(inspect.getargspec(cls.__init__).keywords, {}).items()})
        # return cls(**args, *varargs, **kwargs)
        return cls(**kwargs)



    # @serialize_fw
    # def to_dict(self):
    #     """
    #     Explicit conversion of all the elements needed to reconstruct the Task.
    #     Does not rely on the recursive serialization.
    #     """
    #     abinit_input = self.abitask.strategy.abinit_input
    #     pseudos_name_list = []
    #     for pseudo in abinit_input.pseudos:
    #         pseudos_name_list.append({'basename': pseudo.basename, 'path': pseudo.filepath})
    #     abinit_input_dict = {'pseudos': pseudos_name_list}
    #     dtsets = []
    #     for ds in abinit_input:
    #         ds_copy = ds.deepcopy()
    #         for key, value in ds_copy.iteritems():
    #             # Convert to list to allow serialization
    #             if isinstance(value, np.ndarray):
    #                 ds_copy[key] = value.tolist()
    #         dtsets.append(dict(ds_copy))
    #     abinit_input_dict['datasets'] = dtsets
    #
    #     strategy = {'abinit_input': abinit_input_dict,
    #                 '@module': self.abitask.strategy.__class__.__module__,
    #                 '@class': self.abitask.strategy.__class__.__name__}
    #
    #     # All the dependencies are saved using the path as node
    #     deps = []
    #     for dep in self.abitask.deps:
    #         for path, ext in zip(*dep.get_filepaths_and_exts()):
    #             deps.append({'path': path, 'ext': ext})
    #
    #     #TODO remove status? does not get updated anyway
    #     d = {'abitask': {'strategy': self.abitask.strategy.as_dict(),
    #          'deps': deps,
    #          'status': str(self.abitask.status),
    #          '@module': self.abitask.__class__.__module__,
    #          '@class': self.abitask.__class__.__name__}}
    #
    #     return d

    @classmethod
    def _abiinput_from_dict(cls, m_dict):

        dtsets = m_dict['datasets']
        pseudos_dict = m_dict['pseudos']
        pseudos = []
        #TODO improve the way pseudos are retrieved. environment variable?
        for pseudo in pseudos_dict:
            # Try first to get from abidata for testing purposes
            try:
                pseudos.append(abidata.pseudo(pseudo['basename']))
            except:
                pass
            pseudo_path = os.path.expandvars(pseudo['path'])
            pseudos.append(PseudoParser().parse(pseudo_path))
        abiinput = abilab.AbiInput(pseudos, ndtset=dtsets[0]['ndtset'])
        n = 0
        for ds in dtsets:
            #TODO save and set explicitly the structure object?
            abiintput.set_vars(dtset=n, **ds)
            n += 1

        return abiinput

    @classmethod
    def _basic_from_dict(cls, m_dict):
        # Load strategy
        strategy_dict = m_dict['strategy']
        strategy_cls_mod = strategy_dict['@module']
        strategy_cls_name = strategy_dict['@class']
        strategy_module = __import__(strategy_cls_mod, globals(), locals(), [strategy_cls_name])
        strategy_cls = getattr(strategy_module, strategy_cls_name)
        #TODO this will work just for StrategyWithInput.
        #strategy = strategy_cls(cls._abiinput_from_dict(strategy_dict['abinit_input']))
        strategy = strategy_cls.from_dict(strategy_dict)

        # Load dependencies using the FileNode as Node
        # The Task constructor requires a dict
        deps_list = m_dict['deps']
        deps = {}
        for dep in deps_list:
            deps[FileNode(dep['path'])] = dep['ext']

        return strategy, deps

    # @classmethod
    # def from_dict(cls, m_dict):
    #     """
    #     Explicit conversion to a Task object from the dict saved with to_dict
    #     Does not rely on the recursive serialization.
    #     """
    #     strategy, deps = cls._basic_from_dict(m_dict['abitask'])
    #
    #     task_module = __import__(m_dict['abitask']['@module'], globals(), locals(), [m_dict['abitask']['@class']])
    #     task_cls = getattr(task_module, m_dict['abitask']['@class'])
    #     abitask = task_cls(strategy, workdir='.',
    #                        manager=abilab.TaskManager.from_user_config(), deps=deps)
    #
    #     return cls(abitask)

    def _update_abinit_vars(self, fw_spec):
        """
        Internal method that retrieves the abinit variables to be updated from the spec
        """
        if 'updated_abinit_vars' in fw_spec:
            self.abitask.strategy.add_extra_abivars(fw_spec['updated_abinit_vars'])

    def _update_task_manager(self, fw_spec):
        """
        Internal method to update the task manager according to the _queueadapter spec
        """
        if '_queueadapter' in fw_spec:
            if 'ntasks' in fw_spec['_queueadapter']:
                self.abitask.manager.set_mpi_procs(fw_spec['_queueadapter']['ntasks'])

    def _get_dep_path(self, dirname, ext):
        filepath = os.path.join(dirname, self.abitask.prefix.odata + "_" + ext)
        if os.path.isfile(filepath):
            return filepath
        else:
            # if not existent it could be a multiple file. try to retrieve it
            dirname = os.path.dirname(filepath)
            files = glob.glob(os.path.join(dirname, '*[0-9]_'+ext))
            if not files:
                return filepath
            new_path = sorted(glob.glob(os.path.join(dirname, '*[0-9]_'+ext)),
                              key=lambda k: int(re.search('(\d+)(?!.*\d)', k).group()))[-1]
            logger.warning("Path %s does not exist, used %s instead." % (filepath, new_path))
            return new_path

    def config_run(self, fw_spec):
        # Set a logger for abinitio and abipy
        log_handler = logging.FileHandler('abipy.log')
        log_handler.setFormatter(logging.Formatter(logging.BASIC_FORMAT))
        logging.getLogger('pymatgen.io.abinitio').addHandler(log_handler)
        logging.getLogger('abipy').addHandler(log_handler)

        # Check dependencies - Old
        # for dep in self.abitask.deps:
        #     if dep.status is not Node.S_OK:
        #         raise InitializationError(self, "Dependency is not satisfied: {}".format(str(dep)))

        # Resolve dependencies for structure and files and check
        abi_deps = fw_spec.get('abi_deps', {})
        for node_id, refs in self.deps.items():
            refs_list = refs.split()
            if 'structure' in refs_list:
                refs_list.remove('structure')
                if not abi_deps.get('struct_'+str(node_id), None):
                    raise InitializationError(self, "No dependency information in spec for task id: {}".format(node_id))
                updated_structure = abilab.Structure.from_dict(abi_deps['struct_'+str(node_id)])
                self.abitask.strategy.structure = updated_structure
            for ref in refs_list:
                # Check for path in fw_spec
                if not abi_deps.get('dep_'+str(node_id), None):
                    raise InitializationError(self, "No dependency information in spec for task id: {}".format(node_id))
                path = abi_deps.get('dep_'+node_id)
                path =  self._get_dep_path(path, ref)# os.path.join(path, self.abitask.prefix.odata + "_" + ref)
                dep = Dependency(FileNode(path), ref)
                # Check if dep is satisfied
                if dep.status is not Node.S_OK:
                    raise InitializationError(self, "Dependency is not satisfied: {}".format(str(dep)))
                self.abitask.add_deps(dep)

        # Update data from spec
        self._update_task_manager(fw_spec)
        self._update_abinit_vars(fw_spec)

        # set workdir to that where FW is running
        if not getattr(self.abitask, 'workdir', None):
            self.abitask.set_workdir(os.path.abspath('.'))

    def run_abitask(self, fw_spec):
        """
        Runs the single task and performs the check_status check regularly.
        """
        automatic_restart = True
        error_code = fw_spec.get('_exception_details', {}).get('error_code')
        if automatic_restart and error_code in (ErrorCode.UNCONVERGED,):
            # If previous launch set an error code try to restart the calculation
            try:
                started = self.abitask.restart()
            except TaskRestartError as e:
                # If restart fails, try to clean the run start
                logger.warning("couldn't restart task from previous run. Error: {}. Clean and restart.".format(e))
                self.abitask.reset()
                started = self.abitask.start(autoparal=False)
        else:
            # Else start the task normally
            # TODO add other checks?
            started = self.abitask.start(autoparal=False)
        if not started:
            raise InitializationError(self, "Task did not start.")

        #TODO get the values as external parameters. from the taskmanagerl.yml?
        self.sleep_interval = 5
        self.check_interval = 20
        check_skip = self.check_interval // self.sleep_interval
        # The check is called every check_skip intervals.
        # If the job has terminated a last check is performed
        job_terminated = False
        while not job_terminated:
            for n in range(check_skip):
                job_terminated = self.abitask.process.poll() is not None
                if job_terminated:
                    break
                time.sleep(self.sleep_interval)
            self.abitask.check_status()

    def run_task(self, fw_spec):
        """
        Implementation of Fireworks run_task. Configures, executes and analyzes
        the job execution in three separate steps
        """
        self.config_run(fw_spec)

        self.run_abitask(fw_spec)

        #TODO:
        # -check for errors
        # -add more info about the task
        # -remove lock file?
        # -are there other actions to be performed in abipy at the end of task?

        return self.task_analysis(fw_spec)

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
        if self.abitask.status < self.abitask.S_OK:
            raise AbinitRuntimeError(self)

        if 'create_file' in fw_spec:
            open(fw_spec['create_file'], 'a').close()

        stored_data = {'history': list(self.abitask.history)}
        report = self.abitask.get_event_report()
        stored_data['num_comments'] = report.num_comments
        if report.num_comments:
            comments = []
            for comment in report.comments:
                comments.append({'name': comment.name, 'message': comment.message})
            stored_data['comments'] = comments
        stored_data['num_warnings'] = report.num_warnings
        if report.num_warnings:
            warnings = []
            for warning in report.warnings:
                warnings.append({'name': warning.name, 'message': warning.message})
            stored_data['warnings'] = warnings

        # Forward previous dependencies and add the current one
        update_spec = dict(abi_deps=dict({'dep_'+str(self.task_id): os.getcwd()}), **fw_spec.get('abi_deps', {}))

        return FWAction(stored_data=stored_data, update_spec=update_spec)



#########################################################################################
# Additional tasks
#########################################################################################

# @explicit_serialize
# class AutoparalFWTask(AbiFireTask):
#     """
#     Task to execute the autoparal method of abinit and apply the results to the following
#     tasks in the workflow.
#     """
#
#     #TODO maybe the init should be modified to be independent of the task type
#     def __init__(self, *args, **kwargs):
#         self.max_ncpus = kwargs.pop('max_ncpus')
#         super(AutoparalFWTask, self).__init__(*args, **kwargs)
#         self.abitask.manager.policy.autoparal = 1
#         if self.max_ncpus:
#             self.abitask.manager.policy.max_ncpus = self.max_ncpus
#
#     def to_dict(self):
#         d = super(AutoparalFWTask, self).to_dict()
#         d['max_ncpus'] = self.max_ncpus
#
#         return d
#
#     @classmethod
#     def from_dict(cls, m_dict):
#         strategy, deps = cls._basic_from_dict(m_dict)
#         task = cls(strategy, workdir=m_dict['workdir'],
#                    manager=abilab.TaskManager.from_user_config(), deps=deps,
#                    max_ncpus=m_dict['max_ncpus'])
#
#         return task
#
#     @classmethod
#     def from_task(cls, task, max_ncpus=None):
#         """
#         Creates the AutoparalFWTask starting from the task that will be run.
#         If max_ncpus is None the one provided by the local taskmanager is used.
#         """
#         workdir = task.workdir + '_autoparal_fake_run'
#         return cls(task.abitask.__class__(task.strategy.deepcopy(), workdir=workdir,
#                                           manager=abilab.TaskManager.from_user_config(), deps=task.deps),
#                    max_ncpus=max_ncpus)
#
#     def run_task(self, fw_spec):
#         """
#         Runs the autoparal_fake_run method and updates queuadapter information
#         for the following tasks.
#         """
#         #TODO add a try?
#         confs, optimal = self.abitask.autoparal_fake_run()
#
#         #TODO this will probably need to become a mod_spec
#         update_spec = {'_queueadapter': {'ntasks': optimal.mpi_ncpus},
#                        'updated_abinit_vars': optimal.vars}
#
#         return FWAction(update_spec=update_spec)


@explicit_serialize
class AutoparalFireTask(AbiFireTask):
    """
    Task to execute the autoparal method of abinit and apply the results to the following
    tasks in the workflow.
    """

    def __init__(self, firetask, deps={}):
        self.firetask = firetask
        self.abitask = firetask.abitask
        self.deps = deps
        self.firetask.abitask.manager.policy.autoparal = 1
        # if self.max_ncpus:
        #     self.abitask.manager.policy.max_ncpus = self.max_ncpus

    def run_abitask(self, fw_spec):
        """
        Runs the autoparal_fake_run method and updates queuadapter information
        for the following tasks.
        """

        # Build here because it won't do that in the autoparal_run
        self.abitask.build()
        self.abitask.make_links()
        self.abitask.setup()

        #TODO add a try?
        stauts = self.abitask.autoparal_run()

        if stauts != 0:
            raise RuntimeError('Autoparall did not run')

        output = os.path.join(self.abitask.workdir, "autoparal.json")
        with open(output, "r") as f:
            self.report = json.load(f)

    def task_analysis(self, fw_spec):
        #TODO this will probably need to become a mod_spec
        update_spec = {'_queueadapter': {'ntasks': self.report['optimal_conf']['mpi_ncpus']},
                       'updated_abinit_vars': self.report['optimal_conf']['vars']}

        return FWAction(update_spec=update_spec)


def initializer(func):
    """
    Automatically assigns the parameters to variables.
    """

    varnames, varargs, varkw, defaults = inspect.getargspec(func)

    @wraps(func)
    def wrapper(self, *args, **kargs):
        l = list(zip(varnames[1:], args))
        if not varkw:
            l += list(kargs.items())
        for name, arg in l:
            setattr(self, name, arg)

        for name, default in zip(reversed(varnames), reversed(defaults)):
            if not hasattr(self, name):
                setattr(self, name, default)

        if varkw:
            setattr(self, varkw, varkw)

        func(self, *args, **kargs)

    return wrapper

#########################################################################################
# Strategy based tasks
#########################################################################################


class AbiStrategyFireTask(AbiFireTask):

    @classmethod
    def parse_deps(cls, deps):
        new_deps = {}
        for k, v in deps.items():
            if isinstance(k, basestring):
                new_deps[k] = v
            else:
                new_deps[str(k)] = v

        return new_deps


@explicit_serialize
class ScfStrategyFireTask(AbiStrategyFireTask):

    def __init__(self, structure, pseudos, ksampling=1000, accuracy="normal", spin_mode="polarized",
                 smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None, use_symmetries=True,
                 deps={}, task_id=None, **extra_abivars):

        if task_id:
            self.task_id = task_id
        else:
            self.task_id = id(self)

        extra_abivars.pop('_fw_name', None)

        if not isinstance(pseudos, PseudoTable):
            pseudos = PseudoTable(pseudos)

        pseudos_paths = []
        for p in pseudos:
            pseudos_paths.append(p.filepath)

        self.structure = structure
        self.pseudos = pseudos_paths
        self.ksampling = ksampling
        self.accuracy = accuracy
        self.spin_mode = spin_mode
        self.smearing = smearing
        self.charge = charge
        self.scf_algorithm = scf_algorithm
        self.use_symmetries = use_symmetries
        self.extra_abivars = extra_abivars
        self.deps = self.parse_deps(deps)

        if isinstance(ksampling, int):
            ksampling = KSampling.automatic_density(structure, ksampling)

        strategy = ScfStrategy(structure, pseudos, ksampling, accuracy, spin_mode, smearing, charge, scf_algorithm,
                               use_symmetries, **extra_abivars)

        self.abitask = ScfTask(strategy, workdir=None, manager=abilab.TaskManager.from_user_config(), deps=None)


@explicit_serialize
class RelaxStrategyFireTask(AbiStrategyFireTask):

    def __init__(self, structure, pseudos, ksampling=1000, relax_algo="atoms_only", accuracy="normal",
                 spin_mode="polarized", smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None,
                 deps={}, task_id=None, **extra_abivars):

        if task_id:
            self.task_id = task_id
        else:
            self.task_id = id(self)

        extra_abivars.pop('_fw_name', None)

        if not isinstance(pseudos, PseudoTable):
            pseudos = PseudoTable(pseudos)

        pseudos_paths = []
        for p in pseudos:
            pseudos_paths.append(p.filepath)

        self.structure = structure
        self.pseudos = pseudos_paths
        self.ksampling = ksampling
        self.relax_algo = relax_algo
        self.accuracy = accuracy
        self.spin_mode = spin_mode
        self.smearing = smearing
        self.charge = charge
        self.scf_algorithm = scf_algorithm
        self.extra_abivars = extra_abivars
        self.deps = self.parse_deps(deps)

        if isinstance(ksampling, int):
            ksampling = KSampling.automatic_density(structure, ksampling)

        if isinstance(relax_algo, basestring) and relax_algo in relaxation_methods:
            relax_algo = relaxation_methods[relax_algo]()

        strategy = RelaxStrategy(structure, pseudos, ksampling, relax_algo, accuracy, spin_mode, smearing, charge,
                                 scf_algorithm, **extra_abivars)

        self.abitask = RelaxTask(strategy, workdir=None, manager=abilab.TaskManager.from_user_config(), deps={})

    def task_analysis(self, fw_spec):
        """
        A relax task updates forwards an updated structure for the following tasks/FW
        """

        action = super(RelaxStrategyFireTask, self).task_analysis(fw_spec)
        # FIXME here the code assumes that 'abi_dept' is already in the update_spec
        action.update_spec['abi_deps'].update(
            {'struct_'+str(self.task_id): self.abitask.read_final_structure().as_dict()})

        return action


@explicit_serialize
class NscfStrategyFireTask(AbiStrategyFireTask):

    def __init__(self, scf_task, ksampling=1000, nscf_bands=None, nscf_algorithm=None, deps={},
                 task_id=None, **extra_abivars):

        if task_id:
            self.task_id = task_id
        else:
            self.task_id = id(self)

        extra_abivars.pop('_fw_name', None)

        self.scf_task = scf_task
        self.ksampling = ksampling
        self.nscf_bands = nscf_bands
        self.nscf_algorithm = nscf_algorithm
        self.extra_abivars = extra_abivars
        self.deps = self.parse_deps(deps)

        if isinstance(ksampling, int):
            ksampling = KSampling.automatic_density(scf_task.structure, ksampling)

        strategy = NscfStrategy(scf_task.abitask.strategy, ksampling, nscf_bands, nscf_algorithm, **extra_abivars)

        self.abitask = NscfTask(strategy, workdir=None, manager=abilab.TaskManager.from_user_config(), deps={})


@explicit_serialize
class StrategyWithInputFireTask(AbiStrategyFireTask):

    def __init__(self, abinit_input, task_type='AbinitTask', deps={}, task_id=None):

        if task_id:
            self.task_id = task_id
        else:
            self.task_id = id(self)

        self.abinit_input = abinit_input
        self.task_type = task_type
        self.deps = self.parse_deps(deps)

        strategy = StrategyWithInput(abinit_input)

        tasks_mod_name = AbinitTask.__module__
        try:
            tasks_mod = __import__(tasks_mod_name, globals(), locals(), [task_type])
            task_cls = getattr(tasks_mod, task_type)
        except (ImportError, AttributeError) as e:
            msg = "Cannot load task_type %s from module %s" % (task_type, tasks_mod_name)
            logger.error(msg)
            raise ValueError(msg)

        self.abitask = task_cls(strategy, workdir=None, manager=abilab.TaskManager.from_user_config(), deps={})



#TODO fix after refactoring
@explicit_serialize
class MultiStepRelaxStrategyFireTask(RelaxStrategyFireTask):
    """
    This class derives from RelaxFWTask, but allows a multi-FW relaxation.
    If the final status is not converged, it creates a new FW to continue the relaxation.

    This can be useful to split long relaxations in short jobs.
    """

    def __init__(self, structure, pseudos, ksampling=1000, relax_algo="atoms_only", accuracy="normal",
                 spin_mode="polarized", smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None,
                 deps={}, additional_steps=10, task_id=None, **extra_abivars):

        super(MultiStepRelaxStrategyFireTask, self).__init__(
            structure=structure, pseudos=pseudos, ksampling=ksampling, relax_algo=relax_algo, accuracy=accuracy,
            spin_mode=spin_mode, smearing=smearing, charge=charge, scf_algorithm=scf_algorithm, deps=deps,
            task_id=task_id, **extra_abivars)
        self.additional_steps = additional_steps

    @classmethod
    def new_workdir(cls, dir):
        " Generates the directory for the new task"
        if re.search(r'relaxstep\d+$', dir):
            num = int(re.search(r'\d+$', dir).group())
            return re.sub(r'\d+$', str(num+1), dir)
        else:
            return dir+'_relaxstep2'

    def task_analysis(self, fw_spec):
        """
        A relax task updates forwards an updated structure for the following tasks/FWs.
        If the status is Unconverged does not raise an exception, but creates a new FW.
        Previous dependencies are not forwarded to the new FW
        """

        # Raise an exception if the final status is not Unconverged or OK.
        if self.abitask.status < self.abitask.S_UNCONVERGED or self.abitask.status == self.abitask.S_ERROR:
            raise AbinitRuntimeError(self)

        if self.abitask.status == self.abitask.S_UNCONVERGED:
            stored_data = {'history': list(self.abitask.history)}
            if self.additional_steps <= 1:
                raise AbinitRuntimeError(self)
            new_task = self.copy()
            new_task.structure = self.abitask.read_final_structure()
            new_task.deps = self.parse_deps({self: "WFK"})
            new_task.additional_steps = self.additional_steps-1
            new_spec = {'abi_deps': {'dep_'+str(self.task_id): os.getcwd()}}
            if '_queueadapter' in fw_spec:
                new_spec['_queueadapter'] = fw_spec.get('_queueadapter')
            if '_launch_dir' in fw_spec:
                new_spec['_launch_dir'] = self.new_workdir(fw_spec['_launch_dir'])
            new_step = Firework(new_task, spec=new_spec)
            return FWAction(stored_data=stored_data, detours=new_step)
        else:
            return super(MultiStepRelaxStrategyFireTask, self).task_analysis(fw_spec)

relaxation_methods = {
    'atoms_only': RelaxationMethod.atoms_only,
    'atoms_and_cell': RelaxationMethod.atoms_and_cell
}

#########################################################################################
# AbiFlow derived tasks
#########################################################################################

@explicit_serialize
class ScfFWTask(StrategyWithInputFireTask):
    def __init__(self, *args, **kwargs):
        super(ScfFWTask, self).__init__(task_type=ScfTask, *args, **kwargs)


@explicit_serialize
class NscfFWTask(StrategyWithInputFireTask):
    def __init__(self, *args, **kwargs):
        super(NscfFWTask, self).__init__(task_type=ScfTask, *args, **kwargs)


@explicit_serialize
class RelaxFWTask(StrategyWithInputFireTask):
    def __init__(self, *args, **kwargs):
        super(RelaxFWTask, self).__init__(task_type=ScfTask, *args, **kwargs)

    def task_analysis(self, fw_spec):
        """
        A relax task updates forwards an updated structure for the following tasks/FW
        """

        action = super(RelaxStrategyFireTask, self).task_analysis(fw_spec)
        # FIXME here the code assumes that 'abi_dept' is already in the update_spec
        action.update_spec['abi_deps'].update(
            {'struct_'+str(self.task_id): self.abitask.read_final_structure().as_dict()})

        return action


@explicit_serialize
class DdkFWTask(StrategyWithInputFireTask):
    def __init__(self, *args, **kwargs):
        super(DdkFWTask, self).__init__(task_type=ScfTask, *args, **kwargs)


@explicit_serialize
class PhononFWTask(StrategyWithInputFireTask):
    def __init__(self, *args, **kwargs):
        super(PhononFWTask, self).__init__(task_type=ScfTask, *args, **kwargs)


@explicit_serialize
class SigmaFWTask(StrategyWithInputFireTask):
    def __init__(self, *args, **kwargs):
        super(SigmaFWTask, self).__init__(task_type=ScfTask, *args, **kwargs)


@explicit_serialize
class BseFWTask(StrategyWithInputFireTask):
    def __init__(self, *args, **kwargs):
        super(BseFWTask, self).__init__(task_type=ScfTask, *args, **kwargs)


@explicit_serialize
class OpticFWTask(StrategyWithInputFireTask):
    def __init__(self, *args, **kwargs):
        super(OpticFWTask, self).__init__(task_type=ScfTask, *args, **kwargs)


@explicit_serialize
class AnaddbFWTask(StrategyWithInputFireTask):
    def __init__(self, *args, **kwargs):
        super(AnaddbFWTask, self).__init__(task_type=ScfTask, *args, **kwargs)
