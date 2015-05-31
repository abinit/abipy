from __future__ import print_function, division, unicode_literals

import pytest
from abipy.fworks.fw_workflows import InputFWWorkflow
from abipy.fworks.fw_tasks import ScfFWTask
from abipy.core.testing import has_abinit, has_fireworks, has_mongodb

from fireworks.core.rocket_launcher import rapidfire, launch_rocket

ABINIT_VERSION = "7.11.5"

# pytestmark = [pytest.mark.skipif(not has_abinit(ABINIT_VERSION), reason="Abinit version {} is not in PATH".format(ABINIT_VERSION)),
#               pytest.mark.skipif(not has_fireworks(), reason="fireworks paackage is missing"),
#               pytest.mark.skipif(not has_mongodb(), reason="no connection to mongodb")]

pytestmark = pytest.mark.usefixtures("cleandb")

class ItestScf():

    def itest_run(self, lp, fworker, tmpdir, input_scf_si_low):
        wf = InputFWWorkflow(input_scf_si_low, task_type=ScfFWTask)

        scf_fw_id = wf.fw.fw_id
        old_new = wf.add_to_db(lpad=lp)
        scf_fw_id = old_new[scf_fw_id]

        rapidfire(lp, fworker, m_dir=str(tmpdir))

        fw = lp.get_fw_by_id(scf_fw_id)

        assert fw.state == "COMPLETED"

    def itest_not_converged(self, lp, fworker, tmpdir, input_scf_si_low):
        input_scf_si_low.set_vars(nstep=5)

        wf = InputFWWorkflow(input_scf_si_low, task_type=ScfFWTask)

        scf_fw_id = wf.fw.fw_id
        old_new = wf.add_to_db(lpad=lp)
        scf_fw_id = old_new[scf_fw_id]

        rapidfire(lp, fworker, m_dir=str(tmpdir), nlaunches=1)

        fw = lp.get_fw_by_id(scf_fw_id)

        assert fw.state == "COMPLETED"

        launch = fw.launches[-1]

        assert len(launch.action.stored_data['events']) == 1

        links = lp.get_wf_by_fw_id(scf_fw_id).links

        assert scf_fw_id in links

        fw_child_id = links[scf_fw_id][0]
        fw_child = lp.get_fw_by_id(fw_child_id)

        assert fw_child.state == "READY"

        rapidfire(lp, fworker, m_dir=str(tmpdir), nlaunches=1)

        fw_child = lp.get_fw_by_id(scf_fw_id)

        assert fw_child.state == "COMPLETED"

        # check that the new fw didn't generate further fws
        assert fw_child not in lp.get_wf_by_fw_id(scf_fw_id).links





