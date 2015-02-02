# coding: utf-8
from __future__ import print_function, division, unicode_literals

try:
    from fireworks import FireTaskBase, FWAction, Firework, explicit_serialize
except ImportError:
    FireTaskBase, FWAction, Firework = 3 * [object]
    explicit_serialize = lambda x: x

from abipy.abilab import Flow

@explicit_serialize
class FireTaskWithFlow(FireTaskBase):

    required_params = ["flow"]

    #def __init__(self, *args, **kwargs):
    #    #print("args", args, "kwargs", kwargs)
    #    #self._flow = kwargs["flow"]
    #    super(FireTaskWithFlow, self).__init__(*args, **kwargs)

    def run_task(self, fw_spec):
        # print("entering run_task: %s " % str(self))

        # Initialize the Flow from the filepath.
        # Most of the serious stuff is delegated to flow.
        # Note that FireTask converts all the objects to strings
        # hence self["flow"] is the workdir of the flow from which
        # we can reconstruct our object with pickle.
        print("In MyTask with flow %s" % self["flow"])
        workdir = self["flow"]["workdir"]
        flow = Flow.pickle_load(workdir)

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
            new_fw = Firework(FireTaskWithFlow(flow=flow), {'all_ok': all_ok})
            return FWAction(stored_data={'all_ok': all_ok}, additions=new_fw)
