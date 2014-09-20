# coding: utf-8
from __future__ import print_function, division, unicode_literals

from fireworks import FireTaskBase, FWAction, FireWork, explicit_serialize
from abipy.abilab import AbinitFlow

@explicit_serialize
class FireTaskWithFlow(FireTaskBase):

    required_params = ["flow"]

    #def __init__(self, *args, **kwargs):
    #    #print("args", args, "kwargs", kwargs)
    #    #self._flow = kwargs["flow"]
    #    super(FireTaskWithFlow, self).__init__(*args, **kwargs)

    def run_task(self, fw_spec):
        # print("entering run_task: %s " % str(self))

        # Initialize the AbinitFlow from the filepath.
        # Most of the serious stuff is delegated to flow.
        # Note that FireTask converts all the objects to strings
        # hence self["flow"] is the workdir of the flow from which
        # we can reconstruct our object via pickle.
        flow = AbinitFlow.pickle_load(self["flow"])
        #print("In MyTask with flow %s" % flow)

        flow.build_and_pickle_dump()
        try:
            all_ok = fw_spec["all_ok"]
        except KeyError:
            all_ok = False

        if all_ok:
            print("all_ok will stop")
            return FWAction()
        else:
            print("launched %s" % flow.rapidfire(check_status=True))
            all_ok = flow.all_ok
            new_fw = FireWork(FireTaskWithFlow(flow=flow), {'all_ok': all_ok})
            return FWAction(stored_data={'all_ok': all_ok}, additions=new_fw)
