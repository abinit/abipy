from __future__ import print_function, division, unicode_literals

import pytest
import abipy.abilab as abilab
from abipy.fworks.fw_workflows import InputFWWorkflow
from abipy.core.testing import has_abinit, has_fireworks, has_mongodb

from fireworks.core.rocket_launcher import rapidfire, launch_rocket

ABINIT_VERSION = "7.11.5"

# pytestmark = [pytest.mark.skipif(not has_abinit(ABINIT_VERSION), reason="Abinit version {} is not in PATH".format(ABINIT_VERSION)),
#               pytest.mark.skipif(not has_fireworks(), reason="fireworks paackage is missing"),
#               pytest.mark.skipif(not has_mongodb(), reason="no connection to mongodb")]

pytestmark = [pytest.mark.usefixtures("cleandb"), pytest.mark.check_abiflow]


def match_results(t1, t2):
    with t1.open_gsr() as gsr1, t2.open_gsr() as gsr2:
        if gsr1.energy - gsr2.energy > 0.0001:
            return False

    return True


class ItestCheck():

    def itest_scf(self, lp, fworker, fwp, tmpdir, benchmark_input_scf):
        wf = InputFWWorkflow(benchmark_input_scf)

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






