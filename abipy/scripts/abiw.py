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
#import tempfile
import abipy.flowtk as flowtk
#from abipy import abilab
from abipy.core.release import __version__

from pprint import pprint
#from collections import defaultdict, OrderedDict
#from socket import gethostname
#from monty import termcolor
from monty.functools import prof_main
#from monty.termcolor import cprint, colored, get_terminal_size
#from monty.string import boxed, make_banner
#from abipy.tools import duck
from abipy.flowtk.worker import WorkerClients, WorkerServer


#def straceback():
#    """Returns a string with the traceback."""
#    import traceback
#    return traceback.format_exc()


#def parse_strings(s):
#    """Parse comma separated values. Return None if s is None."""
#    return s.split(",") if s is not None else s


def get_epilog():
    usage = """\

Usage example:

  abirun.py FLOWDIR rapid                 => Keep repeating, stop when no task can be executed.
  abirun.py FLOWDIR scheduler             => Execute flow with the scheduler.


#abidb init
#db clean
#
abiw.py add name manager.yml scheduler
abiw.py list  DONE
abiw.py disable name
abiw.py set_default name
abiw.py ping [names]
abiw.py status [names | all]
abiw.py start name DONE
abiw.py restart name DONE TO_BE_TESTED
abiw.py kill name DONE
abiw.py send script.py DONE
abiw.py gui
abiw.py allgui
"""

    notes = """\

Notes:

    If FLOWDIR is not given, abirun.py automatically selects the database located within
    the working directory. An Exception is raised if multiple databases are found.

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
    worker_selector_with_default.add_argument("-n", '--worker-name', default=None,
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
    serve_parser.add_argument("-p", "--port", type=int, default=60073, help="Port")
    serve_parser.add_argument("-d", "--daemonize", action="store_true", default=False, help="Demonize the serve")

    # Build the main parser.
    parser = argparse.ArgumentParser(epilog=get_epilog() if with_epilog else "",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-V', '--version', action='version', version=__version__)

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help', description="Valid subcommands")

    # Subparser for list command.
    p_list = subparsers.add_parser("list", parents=[copts_parser], help="List workers.")

    # Subparser for start command.
    p_start = subparsers.add_parser("start", parents=[copts_parser, worker_selector, serve_parser],
                                    help="Start worker.")
    p_start.add_argument("-s", "--scratch-dir", type=str, default="/tmp",
                         help="Scratch directory in which Flows will be generated. DEFAULT: /tmp")

    p_restart = subparsers.add_parser("restart", parents=[copts_parser, worker_selector, serve_parser],
                                      help="Restart worker.")

    p_kill = subparsers.add_parser("kill", parents=[copts_parser], help="Kill worker.")
    p_kill.add_argument('worker_name', type=str, help="Name of the worker server.")

    # Subparser for send command.
    p_send = subparsers.add_parser("send", parents=[copts_parser, worker_selector_with_default],
                                   help="Send script to the worker.")
    p_send.add_argument("py_paths", nargs="+", help="Python script(s)")

    # Subparser for status command.
    p_status = subparsers.add_parser("status", parents=[copts_parser, worker_selector_with_default],
                                     help="Return status of a single worker.")

    # Subparser for status command.
    p_set_default = subparsers.add_parser("set_default", parents=[copts_parser, worker_selector],
                                           help="Change the default worker.")

    # Subparser for status command.
    p_all_status = subparsers.add_parser("all_status", parents=[copts_parser],
                                          help="Return status af all workers.")

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

    #address = "localhost"
    #if options.no_browser:
    print("""
Use:

    ssh -N -f -L localhost:{port}:localhost:{port} username@you_remote_cluster

for port forwarding.
""")

    import abipy.panels as mod
    assets_path = os.path.join(os.path.dirname(mod.__file__), "assets")

    return dict(
        debug=options.verbose > 0,
        #show=not options.no_browser,
        show=False,
        port=options.port,
        static_dirs={"/assets": assets_path},
        address=options.address,
        #websocket_origin="{address}:{port}",
    )


def serve(worker, options):

    print(worker)
    serve_kwargs = serve_kwargs_from_options(options)

    if not options.daemonize:
        worker.serve(**serve_kwargs)
    else:
        print("Running server {worker.name} in demon mode...")
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
        # Port 0 means to select an arbitrary unused port
        #host, port = "localhost", 60073
        sched_options = {"seconds": 5}
        worker = WorkerServer(options.worker_name, options.address, options.port,
                              sched_options, options.scratch_dir)

        return serve(worker, options)

    elif options.command == "restart":
        worker = WorkerServer.from_config_dir(options.worker_name)
        return serve(worker, options)

    all_clients = WorkerClients.from_json_file()

    if options.command == "kill":
        client = all_clients.select_from_name(options.worker_name)
        client.send_kill_message()

    elif options.command == "list":
        print(all_clients)

    elif options.command == "send":
        client = all_clients.select_from_name(options.worker_name)
        print(client)
        for path in options.py_paths:
            pprint(client.send_pyscript(path))

    elif options.command == "status":
        client = all_clients.select_from_name(options.worker_name)
        pprint(client.get_json_status())

    elif options.command == "set_default":
        all_clients.set_default(options.worker_name)
        print(all_clients)

    elif options.command == "all_status":
        for client in all_clients:
            pprint(client.get_json_status())

    elif options.command == "gui":
        client = all_clients.select_from_name(options.worker_name)
        client.open_webgui()

    #elif options.command == "all_gui":

    return 0


if __name__ == "__main__":
    sys.exit(main())
