from __future__ import print_function, division, unicode_literals

from pymatgen.io.abinitio.events import EventReport, ScfConvergenceWarning, RelaxConvergenceWarning, AbinitError

from fireworks import Firework, FireTaskBase, FWAction

##########################
# Test reports
##########################

def report_ok():
    return EventReport('.')

def report_ScfConvergenceWarning():
    er = EventReport('.', events=[ScfConvergenceWarning(message='Fake warning', src_file=__file__, src_line=0)])
    er.set_run_completed(True,'', ' ')
    return er

def report_RelaxConvergenceWarning():
    return EventReport('.', events=[RelaxConvergenceWarning(message='Fake warning', src_file=__file__, src_line=0)])

def report_AbinitError():
    return EventReport('.', events=[AbinitError(message='Fake warning', src_file=__file__, src_line=0)])

##########################
# Fake Tasks
##########################


class FakeTask(FireTaskBase):
    def run_task(self, fw_spec):
        return FWAction()

fake_fw = Firework([FakeTask()])

##########################
# Test FWTaskManager
##########################

MANAGER_OK="""# lemaitre2 hardware: http://www.ceci-hpc.be/clusters.html#lemaitre2
hardware: &hardware
   num_nodes: 112
   sockets_per_node: 2
   cores_per_socket: 6
   mem_per_node: 48Gb

job: &job
    mpi_runner: mpirun
    shell_env:
        PATH: $HOME/local/bin:$PATH
    modules:
        - python/2.7
    # pre_run is a string in verbatim mode (note |)
    pre_run: |
        ulimit unlimited

# queues
qadapters:
  - priority: 1
    queue:
       qname: defq
       qtype: slurm
    limits:
       timelimit: 3-0:0:0
       min_cores: 1
       max_cores: 12
    hardware: *hardware
    job: *job

fw_policy:
    rerun_same_dir: True,
    max_restarts: 20
    autoparal: True
"""

MANAGER_NO_QADAPTERS="""# lemaitre2 hardware: http://www.ceci-hpc.be/clusters.html#lemaitre2
hardware: &hardware
   num_nodes: 112
   sockets_per_node: 2
   cores_per_socket: 6
   mem_per_node: 48Gb

job: &job
    mpi_runner: mpirun
    shell_env:
        PATH: $HOME/local/bin:$PATH
    modules:
        - python/2.7
    # pre_run is a string in verbatim mode (note |)
    pre_run: |
        ulimit unlimited

fw_policy:
    rerun_same_dir: True,
    max_restarts: 20
    autoparal: True
"""

MANAGER_UNKNOWN_KEYS="""# lemaitre2 hardware: http://www.ceci-hpc.be/clusters.html#lemaitre2
hardware: &hardware
   num_nodes: 112
   sockets_per_node: 2
   cores_per_socket: 6
   mem_per_node: 48Gb

job: &job
    mpi_runner: mpirun
    shell_env:
        PATH: $HOME/local/bin:$PATH
    modules:
        - python/2.7
    # pre_run is a string in verbatim mode (note |)
    pre_run: |
        ulimit unlimited

fw_policy:
    rerun_same_dir: True,
    max_restarts: 20
    autoparal: True
    wrong_key: 1
"""

