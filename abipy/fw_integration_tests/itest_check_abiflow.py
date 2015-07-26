from __future__ import print_function, division, unicode_literals

import os
import pytest
import abipy.abilab as abilab
from abipy.electrons.gsr import GsrFile
from abipy.fworks.fw_workflows import InputFWWorkflow
from abipy.fworks.fw_tasks import ScfFWTask
import abipy.fworks.fw_tasks as fw_tasks
from pymatgen.io.abinitio.utils import Directory
from abipy.core.testing import has_abinit, has_fireworks, has_mongodb

from fireworks.core.rocket_launcher import rapidfire, launch_rocket

ABINIT_VERSION = "7.11.5"

# pytestmark = [pytest.mark.skipif(not has_abinit(ABINIT_VERSION), reason="Abinit version {} is not in PATH".format(ABINIT_VERSION)),
#               pytest.mark.skipif(not has_fireworks(), reason="fireworks paackage is missing"),
#               pytest.mark.skipif(not has_mongodb(), reason="no connection to mongodb")]

pytestmark = [pytest.mark.usefixtures("cleandb"), pytest.mark.check_abiflow]


def match_results(fw, abitask):
    fw_gsr_path = Directory(os.path.join(fw.launches[-1].launch_dir, fw_tasks.OUTDIR_NAME)).has_abiext("GSR")
    with GsrFile(fw_gsr_path) as gsr1, abitask.open_gsr() as gsr2:
        if gsr1.energy - gsr2.energy > 0.0001:
            return False

    return True


class ItestCheckScf():

    def itest_scf(self, lp, fworker, fwp, tmpdir, benchmark_input_scf):
        wf = InputFWWorkflow(benchmark_input_scf, task_type=ScfFWTask)

        scf_fw_id = wf.fw.fw_id
        old_new = wf.add_to_db(lpad=lp)
        scf_fw_id = old_new[scf_fw_id]

        rapidfire(lp, fworker, m_dir=str(tmpdir))

        fw = lp.get_fw_by_id(scf_fw_id)

        assert fw.state == "COMPLETED"

        # Build the flow
        flow = abilab.Flow(fwp.workdir, manager=fwp.manager)
        work = flow.register_task(benchmark_input_scf, task_class=abilab.ScfTask)

        flow.allocate()
        flow.build_and_pickle_dump()

        # Run t0, and check status
        t0 = work[0]
        t0.start_and_wait()
        t0.check_status()
        assert t0.status == t0.S_OK

        assert match_results(fw, t0)

    def itest_scf_not_converged(self, lp, fworker, fwp, tmpdir, benchmark_input_scf):
        old_cwd = os.getcwd()
        benchmark_input_scf.set_vars(nstep=4)

        wf = InputFWWorkflow(benchmark_input_scf, task_type=ScfFWTask)

        scf_fw_id = wf.fw.fw_id
        old_new = wf.add_to_db(lpad=lp)
        scf_fw_id = old_new[scf_fw_id]

        while lp.run_exists(fworker):
            rapidfire(lp, fworker, m_dir=str(tmpdir))

        wf = lp.get_wf_by_fw_id(scf_fw_id)

        assert wf.state == "COMPLETED"

        num_restarts_fw = wf.fws[-1].tasks[0].restart_info.num_restarts

        # Build the flow
        flow = abilab.Flow(fwp.workdir, manager=fwp.manager)
        work = flow.register_task(benchmark_input_scf, task_class=abilab.ScfTask)

        flow.allocate()
        flow.build_and_pickle_dump()

        # go to the main dir (to have the abipy configuration files)
        os.chdir(old_cwd)

        # Run t0, and check status
        t0 = work[0]
        assert flow.make_scheduler().start() == 0

        num_restarts_abiflow = t0.num_restarts

        assert num_restarts_fw == num_restarts_abiflow

        assert match_results(lp.get_wf_by_fw_id(scf_fw_id).fws[-1], t0)









