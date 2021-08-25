#!/usr/bin/env python
"""
This script allows the user to submit the calculations contained in the `Flow`.
It provides a command line interface as well as a graphical interface based on wxpython.
"""
import sys
import os
import argparse
#import time
#import platform
import abipy.flowtk as flowtk

from pprint import pprint, pformat
from monty.functools import prof_main
#from monty.termcolor import cprint, colored, get_terminal_size
#from monty.string import boxed, make_banner
from abipy.core.release import __version__
from abipy.tools.printing import print_dataframe
from abipy.flowtk.worker import (WorkerClients, WorkerServer,
        print_local_workers, create_new_worker, discover_local_workers, rdiscover)

#def straceback():
#    """Returns a string with the traceback."""
#    import traceback
#    return traceback.format_exc()


def get_epilog():
    usage = """\

Usage example:

  abiw.py new_manager.py
  abiw.py clients  DONE
  abiw.py disable name
  abiw.py set_default name
  abiw.py ping [names]
  abiw.py status [names | all]
  abiw.py start name DONE
  abiw.py restart name DONE TO_BE_TESTED
  abiw.py kill name DONE
  abiw.py ladd_flows
  abiw.py send script.py DONE
  abiw.py gui
  abiw.py allgui
"""

    notes = """\

Notes:

    The standard procedure to create a new worker is as follows:

        $ abiw.py new_worker WORKER_NAME --scratch-dir=/tmp

    The new worker will create flows in `--scratch-dir`.
    Note that WORKER_NAME must be unique.

    Start the worker with:

        $ abiw.py start WORKER_NAME

    Open a new terminal and issue:

        $ abiw.py ldiscover

    to discover all the workers running on the localhost and generate
    the clients.json file in ~/.abinit/abipy/

    Finally, one can send python scripts to the worker with:

        $ abiw.py send run_si_ebands.py -w WORKER_NAME
"""

    developers = """\

############
# Developers
############

    abirun.py prof ABIRUN_ARGS               => to profile abirun.py
    abirun.py tracemalloc ABIRUN_ARGS        => to trace memory blocks allocated by Python"""

    return notes + usage # + developers


def get_parser(with_epilog=False):

    # Parent parser for commands requiring a `worker_name`
    worker_selector = argparse.ArgumentParser(add_help=False)
    worker_selector.add_argument('worker_name', type=str, help="Name of the worker server.")

    # Parent parser for commands requiring a `worker_name` that defaults to None i.e. default worker.
    worker_selector_with_default = argparse.ArgumentParser(add_help=False)
    worker_selector_with_default.add_argument("-w", '--worker-name', default=None,
                                help="Select the worker by name. If not specified, the default worker is used.")

    # Parent parser for common options.
    copts_parser = argparse.ArgumentParser(add_help=False)
    copts_parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
        help='verbose, can be supplied multiple times to increase verbosity.')

    #copts_parser.add_argument('--no-colors', default=False, action="store_true", help='Disable ASCII colors.')
    #copts_parser.add_argument('--no-logo', default=False, action="store_true", help='Disable AbiPy logo.')
    copts_parser.add_argument('--loglevel', default="ERROR", type=str,
        help="Set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG.")

    serve_parser = argparse.ArgumentParser(add_help=False)
    serve_parser.add_argument("-a", "--address", type=str, default="localhost", help="Adress")
    serve_parser.add_argument("-p", "--port", type=int, default=0,
                               help="Port. Port 0 means to select an arbitrary unused port")
    serve_parser.add_argument("-d", "--daemonize", action="store_true", default=False, help="Demonize the serve")

    # Build the main parser.
    parser = argparse.ArgumentParser(epilog=get_epilog() if with_epilog else "",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-V', '--version', action='version', version=__version__)

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help', description="Valid subcommands")

    # Subparser for clients command.
    p_clients = subparsers.add_parser("clients", parents=[copts_parser], help="List available clients.")

    # Subparser for lworkers command.
    p_servers = subparsers.add_parser("lworkers", parents=[copts_parser],
                                      help="List AbiPy workers available on the local machine")

    # Subparser for start command.
    p_start = subparsers.add_parser("start", parents=[copts_parser, worker_selector, serve_parser],
                                    help="Start worker.")
    ## FIXME ???
    #p_start.add_argument("-s", "--scratch-dir", type=str, required=True,
    #                     help="Scratch directory in which Flows will be generated")

    p_restart = subparsers.add_parser("restart", parents=[copts_parser, worker_selector, serve_parser],
                                      help="Restart worker.")

    p_new_worker = subparsers.add_parser("new_worker", parents=[copts_parser, worker_selector],
                                      help="Create new AbiPy worker.")
    p_new_worker.add_argument("-s", "--scratch-dir", type=str, required=True,
                              help="Scratch directory in which Flows will be generated")

    p_kill = subparsers.add_parser("kill", parents=[copts_parser], help="Kill worker.")
    p_kill.add_argument('worker_name', type=str, help="Name of the worker server.")

    # Subparser for send command.
    p_send = subparsers.add_parser("send", parents=[copts_parser, worker_selector_with_default],
                                   help="Send script to the worker.")
    p_send.add_argument("py_paths", nargs="+", help="Python script(s)")

    # Subparser for ladd command.
    p_lsend_flows = subparsers.add_parser("lsend_flows", parents=[copts_parser, worker_selector_with_default])
    p_lsend_flows.add_argument("flow_dirs", nargs="+", help="List of flow directories.")

    # Subparser for status command.
    p_status = subparsers.add_parser("status", parents=[copts_parser, worker_selector_with_default],
                                     help="Return status of a single worker.")

    # Subparser for .status command.
    #p_lstatus = subparsers.add_parser("lstatus", parents=[copts_parser],
    #                                  help="Return status of all the local workers.")

    # Subparser for status command.
    p_set_default = subparsers.add_parser("set_default", parents=[copts_parser, worker_selector],
                                           help="Change the default worker.")

    # Subparser for status command.
    p_all_status = subparsers.add_parser("all_status", parents=[copts_parser],
                                          help="Return status af all workers.")

    p_ldiscover = subparsers.add_parser("ldiscover", parents=[copts_parser],
                                        help="Discover AbiPy workers on localhost.")

    p_rdiscover = subparsers.add_parser("rdiscover", parents=[copts_parser],
                                        help="Discover remote AbiPy workers.")

    p_rdiscover.add_argument("hostnames", nargs="+", type=str, help="List of hostnames")

    # Subparser for gui command.
    p_gui = subparsers.add_parser("gui", parents=[copts_parser, worker_selector_with_default],
                                  help="Open GUI for the default worker.")

    # Subparser for scheduler command.
    #p_scheduler = subparsers.add_parser('scheduler', parents=[copts_parser],
    #    help="Run all tasks with a Python scheduler. Requires scheduler.yml either in $PWD or ~/.abinit/abipy.")
    #p_scheduler.add_argument('-w', '--weeks', default=0, type=int, help="Number of weeks to wait.")
    #p_scheduler.add_argument('-d', '--days', default=0, type=int, help="Number of days to wait.")
    #p_scheduler.add_argument('-hs', '--hours', default=0, type=int, help="Number of hours to wait.")
    #p_scheduler.add_argument('-m', '--minutes', default=0, type=int, help="Number of minutes to wait.")
    #p_scheduler.add_argument('-s', '--seconds', default=0, type=int, help="Number of seconds to wait.")

    # Subparser for doc_scheduler
    #p_docsched = subparsers.add_parser('doc_scheduler', parents=[copts_parser],
    #    help="Document the options available in scheduler.yml.")

    # Subparser for panel
    #p_panel = subparsers.add_parser('panel', parents=[copts_parser, flow_selector_parser],
    #                                help="Interact with the flow in the browser (requires panel package).")
    #p_panel.add_argument("-pnt", "--panel-template", default="FastList", type=str,
    #                    help="Specify template for panel dasboard." +
    #                         "Possible values are: FastList, FastGrid, Golden, Bootstrap, Material, React, Vanilla." +
    #                         "Default: FastList"
    #                    )
    #p_panel.add_argument('--no-browser', action='store_true', default=False,
    #                    help=("Start the bokeh server to serve the panel app "
    #                          "but don't open the app in the browser.\n"
    #                          "Use this option to connect remotely from localhost to the machine running the server"))
    #p_panel.add_argument("--port", default=0, type=int, help="Allows specifying a specific port when serving panel app.")

    # Subparser for manager.
    #p_manager = subparsers.add_parser('doc_manager', parents=[copts_parser], help="Document the TaskManager options.")
    #p_manager.add_argument("qtype", nargs="?", default=None, help=("Write job script to terminal if qtype='script' else "
    #    "document the qparams for the given QueueAdapter qtype e.g. slurm."))

    return parser


def serve_kwargs_from_options(options):

    print("""
Use:

    ssh -N -f -L localhost:{port}:localhost:{port} username@your_remote_cluster

for port forwarding.
""")

    import abipy.panels as mod
    assets_path = os.path.join(os.path.dirname(mod.__file__), "assets")

    d =  dict(
        debug=options.verbose > 0,
        show=False,
        port=options.port,
        static_dirs={"/assets": assets_path},
        address=options.address,
        websocket_origin="*",
    )
    print("serve_kwargs:\n", pformat(d), "\n")
    return d


def serve(worker, options):

    print(worker)
    serve_kwargs = serve_kwargs_from_options(options)

    if not options.daemonize:
        worker.serve(**serve_kwargs)
    else:
        print(f"Running server {worker.name} in demon mode...")
        import daemon
        with daemon.DaemonContext():
            worker.serve(**serve_kwargs)

    return 0


@prof_main
def main():

    def show_examples_and_exit(err_msg=None, error_code=1):
        """Display the usage of the script."""
        sys.stderr.write(get_epilog())
        if err_msg: sys.stderr.write("Fatal Error\n" + err_msg + "\n")
        sys.exit(error_code)

    parser = get_parser(with_epilog=True)

    # Parse command line.
    try:
        options = parser.parse_args()
    except Exception as exc:
        show_examples_and_exit(error_code=1)

    if not options.command:
        show_examples_and_exit(error_code=1)

    # loglevel is bound to the string value obtained from the command line argument.
    # Convert to upper case to allow the user to specify --loglevel=DEBUG or --loglevel=debug
    import logging
    numeric_level = getattr(logging, options.loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % options.loglevel)
    logging.basicConfig(level=numeric_level)

    if options.verbose > 2:
        print(options)

    if options.command == "start":
        worker = WorkerServer.init_from_config_dir(options.worker_name)
        return serve(worker, options)

    elif options.command == "restart":
        worker = WorkerServer.restart_from_config_dir(options.worker_name)
        discover_local_workers()
        return serve(worker, options)

    elif options.command == "new_worker":
        create_new_worker(options.worker_name, options.scratch_dir)
        #discover_local_workers()
        return 0

    #if os.path.basename(options.filepath) == "flows.db":
    #    from abipy.flowtk.launcher import print_flowsdb_file
    #    return print_flowsdb_file(options.filepath)

    elif options.command == "lworkers":
        print_local_workers()
        return 0

    elif options.command == "ldiscover":
        discover_local_workers()
        return 0

    elif options.command == "rdiscover":
        rdiscover(options.hostnames)
        return 0

    #elif options.command == "rupdate":
    #    rdiscover()
    #    return 0

    all_clients = WorkerClients.from_json_file()

    if options.command == "kill":
        client = all_clients.select_from_worker_name(options.worker_name)
        client.send_kill_message()

    #if options.command == "remove":
    #    client = all_clients.select_from_worker_name(options.worker_name)
    #    client.send_kill_message()

    elif options.command == "clients":
        #all_clients.print_dataframe()
        all_clients.refresh()
        #print(all_clients)
        print("\nTIP: Remember to execute `ldiscover` or `rdiscover` to discover new AbiPy workers")

    elif options.command == "send":
        client = all_clients.select_from_worker_name(options.worker_name)
        for path in options.py_paths:
            pprint(client.send_pyscript(path))

    elif options.command == "lsend_flows":
        client = all_clients.select_from_worker_name(options.worker_name)
        client.send_flow_dirs(options.flow_dirs)

    #elif options.command == "lstatus":
    elif options.command == "status":
        client = all_clients.select_from_worker_name(options.worker_name)
        json_status = client.get_json_status()

        from pandas.io.json import read_json
        json_status["dataframe"] = read_json(json_status["dataframe"])
        print_dataframe(json_status["dataframe"], title="\nWorker Status:\n")





    elif options.command == "all_status":
        for client in all_clients:
            pprint(client.get_json_state())

    elif options.command == "set_default":
        all_clients.set_default(options.worker_name)
        print(all_clients)

    elif options.command == "gui":
        client = all_clients.select_from_worker_name(options.worker_name)
        print(client)
        client.open_webgui()

    #elif options.command == "all_gui":

    return 0


if __name__ == "__main__":
    sys.exit(main())
