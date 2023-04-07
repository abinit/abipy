#!/usr/bin/env python

import sys
import argparse

import abipy.flowtk as flowtk


def main():

    # Build the main parser.
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('flowdir', help="File or directory containing the ABINIT flow")
    parser.add_argument('--remove-lock', default=False, action="store_true",
        help="Remove the lock on the pickle file used to save the status of the flow.")
    parser.add_argument('--apply', default=False, action="store_true", help="Apply changes to the flow.")

    options = parser.parse_args()

    # Read the flow from the pickle file
    flow = flowtk.Flow.pickle_load(options.flowdir, remove_lock=options.remove_lock)

    # List with the id of the tasks that will be changed.
    nids = []

    for task in flow.iflat_tasks():

        # Select tasks according to class and status
        if task.__class__.__name__ == "PhononTask" and task.status in (task.S_UNCONVERGED, task.S_ERROR):
            nids.append(task.node_id)

            # Here we operate on the AbinitInput of the task.
            # In this particular case, we want to use toldfe instead of tolvrs.

            # Remove all the tolerance variables present in the input.
            task.input.pop_tolerances()

            # Set new stopping criterion
            task.input["toldfe"] = 1e-8
            task.input["nstep"] = 1000
            task.input["prtwf"] = 0

    if options.apply:
        print(f"Resetting {len(nids)} tasks modified by the user.")

        if nids:
            for task in flow.iflat_tasks(nids=nids):
                task.reset()

        print("Writing new pickle file.")

        flow.build_and_pickle_dump()
    else:
        print(f"Dry run mode. {len(nids)} tasks have been modified in memory.")

        if nids:
            for task in flow.iflat_tasks(nids=nids):
                print(task)

        print("Use --apply to reset these tasks and update the pickle file.")
        print("Then restart the scheduler with `nohup abirun.py FLOWDIR scheduler ...`")

    return 0


if __name__ == "__main__":
    sys.exit(main())
