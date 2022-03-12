#!/usr/bin/env python
"""
This script allows the user to submit the calculations contained in the `Flow`.
It provides a command line interface as well as a graphical interface based on wxpython.
"""
import sys
import os
import argparse

from pprint import pprint  #, pformat
from monty.functools import prof_main
from abipy.core.release import __version__
from abipy.core.config import get_config
from abipy.tools.printing import print_dataframe
from abipy.htc.base_models import MongoConnector
from abipy.htc.worker import WorkerClients, AbipyWorker, print_local_workers


def get_epilog():
    usage = """\

Usage example:

  abiw.py new WORKER_NAME        --> Create new worker
  abiw.py start WORKER_NAME      --> Start worker
  abiw.py restart WORKER_NAME    --> Restart worker

  abiw.py lscan                  --> Find workers available on the local machine
  abiw.py rscan HOST1 HOST2      --> Contact host(s) to get list of remote workers
  abiw.py status [names | all]

  abiw.py clients  --> List available clients
  abiw.py ladd_flows
  abiw.py send script.py DONE

  abiw.py gui
  abiw.py mongo_start
  abiw.py mongo_start -c COLLECTION_NAME
  abiw.py mongo_gui -c COLLECTION_NAME

  ??????
  #abiw.py kill name DONE
  #abiw.py disable name
  #abiw.py set_default name
  #abiw.py ping [names]
"""

    notes = """\

Notes:

    The standard procedure to create a new worker is as follows:

        $ abiw.py new WORKER_NAME --scratch-dir=/tmp

    The new worker will generate Flows in `--scratch-dir`.
    Note that WORKER_NAME must be unique.

    Start the worker with:

        $ abiw.py start WORKER_NAME

    Open a new terminal and issue:

        $ abiw.py lscan

    to discover all the workers running on the localhost and generate
    the clients.json file in ~/.abinit/abipy/

    Finally, one can send python scripts to the worker with:

        $ abiw.py send run_si_ebands.py -w WORKER_NAME
"""

    return notes + usage


def get_parser(with_epilog=False):

    # Parent parser for commands requiring a `worker_name`
    worker_selector = argparse.ArgumentParser(add_help=False)
    worker_selector.add_argument('worker_name', type=str, help="Name of the worker server.")

    # Parent parser for commands requiring a `worker_name` that defaults to None i.e. default worker.
    worker_selector_with_default = argparse.ArgumentParser(add_help=False)
    worker_selector_with_default.add_argument(
        "-w", '--worker-name', default=None,
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
    serve_parser.add_argument("-a", "--address", type=str, default="localhost", help="Address")
    serve_parser.add_argument("-p", "--port", type=int, default=0,
                              help="Port. Port 0 means to select an arbitrary unused port")
    serve_parser.add_argument("-d", "--daemonize", action="store_true", default=False, help="Demonize the server")

    # Build the main parser.
    parser = argparse.ArgumentParser(epilog=get_epilog() if with_epilog else "",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-V', '--version', action='version', version=__version__)

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help', description="Valid subcommands")

    # Subparser for clients command.
    p_clients = subparsers.add_parser("clients", parents=[copts_parser], help="List available clients.")

    # Subparser for lworkers command.
    p_lworkers = subparsers.add_parser("lworkers", parents=[copts_parser],
                                       help="List AbiPy workers available on the local machine")

    # Subparser for start command.
    p_start = subparsers.add_parser("start", parents=[copts_parser, worker_selector, serve_parser],
                                    help="Start worker.")

    p_restart = subparsers.add_parser("restart", parents=[copts_parser, worker_selector, serve_parser],
                                      help="Restart worker.")
    p_restart.add_argument("-f", "--force", action='store_true', default=False,
                           help="Force restart even if worker status is found to be `running`.")

    p_new = subparsers.add_parser("new", parents=[copts_parser, worker_selector], help="Create new AbiPy worker.")
    p_new.add_argument("-s", "--scratch-dir", type=str, default=None,
                       help="Scratch directory in which Flows will be generated. "
                            "Use value from ~/.abinit/abipy/config.yml if given.")
    p_new.add_argument("-c", "--collection-name", type=str, default=None,
                       help="Name of the MongDB collection used to store FlowModels.")

    p_mongo_start = subparsers.add_parser(
        "mongo_start", parents=[copts_parser, serve_parser, worker_selector_with_default],
        help="Create and start new AbiPy worker connected to a MongoDB database containing FlowModels.")
    p_mongo_start.add_argument("-c", "--collection-name", type=str,
                               default="",
                               #required=True,
                               help="Name of the MongDB collection used to store FlowModels.")
    p_mongo_start.add_argument("-s", "--scratch-dir", type=str, default=None,
                               help="Scratch directory in which Flows will be generated. "
                                    "If not given, use value from ~/.abinit/abipy/config.yml")

    p_mongo_start.add_argument("--local-db", action='store_true', default=False,
                                help="Assume MongoDB server is running on the localhost with the default port.")

    p_mongo_gui = subparsers.add_parser("mongo_gui", parents=[copts_parser, worker_selector, serve_parser],
                                        help="Start MongoGUI.")

    p_kill = subparsers.add_parser("kill", parents=[copts_parser], help="Kill AbiPy worker.")
    p_kill.add_argument('worker_name', type=str, help="Name of the AbiPy worker.")

    # Subparser for send command.
    p_send = subparsers.add_parser("send", parents=[copts_parser, worker_selector_with_default],
                                   help="Send script to the AbiPy worker.")
    p_send.add_argument("py_path", type=str, help="Python script")
    p_send.add_argument("-m", "--message", default="", type=str,
                        help="Message associated to the submission.")

    # Subparser for ladd command.
    p_lsend_flows = subparsers.add_parser("lsend_flows", parents=[copts_parser, worker_selector_with_default])
    p_lsend_flows.add_argument("flow_dir", type=str, help="Flow directory.")
    p_lsend_flows.add_argument("-m", "--message", default="", type=str,
                               help="Message associated to the submission")

    # Subparser for status command.
    p_status = subparsers.add_parser("status", parents=[copts_parser],
                                     help="Return status of a single worker.")
    p_status.add_argument("worker_names", nargs="*", help="List of worker names.")

    # Subparser for .status command.
    #p_lstatus = subparsers.add_parser("lstatus", parents=[copts_parser],
    #                                  help="Return status of all the local workers.")

    # Subparser for status command.
    p_set_default = subparsers.add_parser("set_default", parents=[copts_parser, worker_selector],
                                          help="Change the default worker.")

    # Subparser for all_status command.
    p_all_status = subparsers.add_parser("all_status", parents=[copts_parser],
                                          help="Return status af all workers.")

    # Subparser for lscan command.
    p_lscan = subparsers.add_parser("lscan", parents=[copts_parser],
                                    help="Discover AbiPy workers on localhost.")

    # Subparser for rscan command.
    p_rscan = subparsers.add_parser("rscan", parents=[copts_parser],
                                    help="Discover remote AbiPy workers.")
    p_rscan.add_argument("hostnames", nargs="*", type=str,
                             help="List of hostnames. "
                                  "If empty, the list is taken from the abipy configuration file (remote_hosts_for_workers)")

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

to activate port forwarding.
""")

    import abipy.panels as mod
    assets_path = os.path.join(os.path.dirname(mod.__file__), "assets")

    d = dict(
        debug=options.verbose > 0,
        #show=False,
        port=options.port,
        static_dirs={"/assets": assets_path},
        address=options.address,
        websocket_origin="*",
        panel_template="FastList",
    )
    #print("serve_kwargs:\n", pformat(d), "\n")
    return d


def serve(worker, options):
    s = str(worker) if options.verbose else repr(worker)
    print(s)

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

    config = get_config()

    def get_scratch_dir():
        """
        Return the scratch_dir used by the AbipyWorker either from the CLI or
        from the global Abipy configuration file.
        Raises RuntimeError if neither values are defined.
        """
        if options.scratch_dir is not None:
            return options.scratch_dir

        _scratch_dir = config.worker_scratchdir
        if _scratch_dir is None:
            raise RuntimeError(
                "The scratch directory of the AbiPy worker must be specified\n"
                "via either the command line interface or the `worker_scratchdir` variable in ~/.abinit/abipy/config.yml"
            )

        return _scratch_dir

    if options.verbose > 2:
        print(options)
        print("global abipy config options:\n", config)

    #_ = WorkerClients.lscan()

    if options.command == "new":
        mng_connector = None
        if options.collection_name:
            mng_connector = MongoConnector.from_abipy_config(collection_name=options.collection_name)
            print(mng_connector)

        scratch_dir = get_scratch_dir()
        worker = AbipyWorker.new_with_name(options.worker_name, scratch_dir, mng_connector=mng_connector)
        print(worker)
        WorkerClients.lscan()
        return 0

    elif options.command == "start":
        worker = AbipyWorker.init_from_dirname(options.worker_name)
        return serve(worker, options)

    elif options.command == "restart":
        worker = AbipyWorker.restart_from_dirname(options.worker_name, force=options.force)
        return serve(worker, options)

    elif options.command == "lscan":
        local_clients = WorkerClients.lscan()
        df = local_clients.summarize()
        print_dataframe(df, title="List of AbiPy workers available on localhost\n")
        return 0

    elif options.command == "rscan":
        hostnames = options.hostnames
        if not hostnames:
            print("Taking list of remote hosts from the AbiPy configuration file.")
            hostnames = config.remote_hosts_for_workers
            for i, rhost in enumerate(hostnames):
                print(f"\n[{i+1}] {rhost}")
            print("\n")

        remote_clients = WorkerClients.rscan(hostnames)
        df = remote_clients.summarize()
        print_dataframe(df, title="List of remote AbiPy workers\n")
        return 0

    elif options.command == "mongo_start":
        if options.local_db:
            mng_connector = MongoConnector.for_localhost(collection_name=options.collection_name)
        else:
            mng_connector = MongoConnector.from_abipy_config(collection_name=options.collection_name)
        print(mng_connector)

        if not options.collection_name:
            # If the collection is not specified, contact the MongoDB server to get the list of
            # collections and print them before exiting.
            print("Use the `-c COLLECTION_NAME` option to specify the MongoDB collection from which\n"
                  "FlowModel documents will be fetched by the AbiPy Worker.\n"
                  "Choose among the following collections:\n")
            for i, cname in enumerate(mng_connector.list_collection_names()):
                print(f"\t[{i+1}]", cname)
            return 1

        # Generate automatically the worker_name from the collection if not given.
        worker_name = options.worker_name
        if worker_name is None:
            worker_name = options.collection_name

        # Make sure there's no worker already connected to the same collection.
        errors = []
        eapp = errors.append
        for local_client in WorkerClients.lscan():
            other_connector = local_client.worker_state.mng_connector
            if other_connector is None: continue

            if options.collection_name == other_connector.collection_name:

                if local_client.worker_state.name == worker_name:
                    # Handle the case in which the user tries to start a worker that already exists.
                    status = local_client.worker_state.status

                    if status == "init":
                        # This is not critical. new_with_name will handle the problem.
                        continue
                    elif status == "dead":
                        eapp(f"Use `abiw.py restart` to restart worker: `{worker_name}`")
                    elif status == "running":
                        eapp(f"Worker `{worker_name}` is already running and connected to the collection.")
                    else:
                        raise ValueError(f"Unknown status: {status}")
            else:
                # Found a different worker connected to the same collection.
                eapp(f"There's already another AbiPy Worker connected to collection: {options.collection_name}")
                eapp(repr(local_client))

        if errors:
            for error in errors:
                print(error)
            print("Aborting now")
            return 1

        # Create the AbiPyWorker and start serving.
        scratch_dir = get_scratch_dir()
        worker = AbipyWorker.new_with_name(worker_name, scratch_dir, mng_connector=mng_connector)

        return serve(worker, options)

    #if os.path.basename(options.filepath) == "flows.db":
    #    from abipy.flowtk.launcher import print_flowsdb_file
    #    return print_flowsdb_file(options.filepath)

    elif options.command == "lworkers":
        print_local_workers()
        return 0

    #local_clients = WorkerClients.lscan()
    all_clients = WorkerClients.from_json_file()

    if options.command == "kill":
        print("Killing worker: {options.worker_name}")
        client = all_clients.select_from_worker_name(options.worker_name)
        client.send_kill_message()
        print("Calling lscan after kill")
        WorkerClients.lscan()

    #if options.command == "rm":
    #    client = all_clients.select_from_worker_name(options.worker_name)
    #    client.send_kill_message()

    elif options.command == "clients":
        all_clients.print_dataframe()
        #if options.refresh:
        #all_clients.refresh()
        #print(all_clients)
        print("\nRemember to execute lscan (rscan) to update the list of local (remote) clients...")

    elif options.command == "send":
        client = all_clients.select_from_worker_name(options.worker_name)
        pprint(client.send_pyscript(options.py_path, user_message=options.message))

    elif options.command == "lsend_flows":
        client = all_clients.select_from_worker_name(options.worker_name)
        options.flow_dir = os.path.abspath(options.flow_dir)

        data = client.send_flow_dir(options.flow_dir, user_message=options.message)
        retcode = data["returncode"]
        if retcode:
            print(data)
        else:
            print(f"Flowdir {options.flow_dir} sent to Worker: {options.worker_name}")

        return retcode

    #elif options.command == "show":  show Workername

    elif options.command == "status":

        if options.worker_names:
            from pandas.io.json import read_json
            for worker_name in options.worker_names:
                client = all_clients.select_from_worker_name(worker_name)
                json_status = client.get_json_status()
                json_status["dataframe"] = read_json(json_status["dataframe"]) #, date_format='iso')  #, date_unit="ns")
                print_dataframe(json_status["dataframe"], title="\nWorker Status:\n")

        else:
            # All workers registered in clients.json.
            for client in all_clients:
                #if client.worker_state.state != "running": continue
                try:
                    pprint(client.get_json_status())
                except Exception as exc:
                    print(f"Exception for client {repr(client)}")
                    print(exc)

    elif options.command == "set_default":
        all_clients.set_default(options.worker_name)
        print(all_clients)

    elif options.command == "gui":
        client = all_clients.select_from_worker_name(options.worker_name)
        client.open_webgui()

    elif options.command == "mongo_gui":
        client = all_clients.select_from_worker_name(options.worker_name)
        mng_connector = client.worker_state.mng_connector
        if mng_connector is None:
            raise ValueError(f"The AbiPy worker {options.worker_name} is not running in MongoDB mode!")
        print(mng_connector)
        serve_kwargs = serve_kwargs_from_options(options)
        mng_connector.open_mongoflow_gui(**serve_kwargs)

    #elif options.command == "all_gui":

    return 0


if __name__ == "__main__":
    sys.exit(main())
