#!/usr/bin/env python
"""
This script provides a command line interface to a MongoDB collection containing FlowModels
"""
from __future__ import annotations

import sys
import os
import ast
import argparse

from typing import Type
from pprint import pprint  #, pformat
from pydantic import BaseModel
from pymongo.collection import Collection
from monty.functools import prof_main
from abipy.core.release import __version__
#from abipy.core.config import get_config
from abipy.tools.printing import print_dataframe
from abipy.htc.base_models import MongoConnector
from abipy.htc.flow_models import FlowModel


MAGIC_COLLECTION_FILENAME = "_abidb_collection"


def get_epilog():
    usage = f"""\

Usage example:

  abidb.py list
  abidb.py stats -c COLLECTION_NAME


TIP: If you are using a single collection and you don't want to specify the collection name with `-c`
     you can create a `{MAGIC_COLLECTION_FILENAME}` file in the directory in which the abidb.py script is executed.
     This file should contain a single line with the name of the MongoDB collection.
"""

    return usage


def _parse_dict(string):
    obj = ast.literal_eval(string)  #; print(obj)
    if not isinstance(obj, dict):
        raise TypeError(f"Expecting dictionary, got {type(obj)}\nobj: {obj}")
    return obj


def _parse_dict_or_list(string):
    obj = ast.literal_eval(string)  #; print(obj)
    if not isinstance(obj, (dict, list)):
        raise TypeError(f"Expecting dictionary or list, got {type(obj)}\nobj: {obj}")
    return obj


def get_parser(with_epilog=False):

    # Parent parser for commands requiring a `worker_name`
    #worker_selector = argparse.ArgumentParser(add_help=False)
    #worker_selector.add_argument('worker_name', type=str, help="Name of the worker server.")

    # Parent parser for common options.
    copts_parser = argparse.ArgumentParser(add_help=False)
    copts_parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
        help='verbose, can be supplied multiple times to increase verbosity.')

    copts_parser.add_argument('--loglevel', default="ERROR", type=str,
                              help="Set the loglevel. Possible values: CRITICAL, ERROR (default), WARNING, INFO, DEBUG.")

    copts_parser.add_argument("-c", "--collection-name", type=str,
                               default="",
                               #required=True,
                               help="Name of the MongDB collection used to store FlowModels.")

    #serve_parser = argparse.ArgumentParser(add_help=False)
    #serve_parser.add_argument("-a", "--address", type=str, default="localhost", help="Address")
    #serve_parser.add_argument("-p", "--port", type=int, default=0,
    #                          help="Port. Port 0 means to select an arbitrary unused port")
    #serve_parser.add_argument("-d", "--daemonize", action="store_true", default=False, help="Demonize the server")

    # Build the main parser.
    parser = argparse.ArgumentParser(epilog=get_epilog() if with_epilog else "",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-V', '--version', action='version', version=__version__)

    # Create the parsers for the sub-commands
    subparsers = parser.add_subparsers(dest='command', help='sub-command help', description="Valid subcommands")

    p_list = subparsers.add_parser("list", parents=[copts_parser], help="List collections")

    p_status = subparsers.add_parser("status", parents=[copts_parser], help=abidb_status.__doc__)

    p_in_structure = subparsers.add_parser("in_structure", parents=[copts_parser], help=abidb_in_structure.__doc__)

    p_flowdata = subparsers.add_parser("flowdata", parents=[copts_parser], help=abidb_flowdata.__doc__)
    p_flowdata.add_argument("oid", type=str, help="Document ID")

    p_cq = subparsers.add_parser("cq", parents=[copts_parser], help=abidb_cq.__doc__)
    p_cq.add_argument("cq_inds", nargs="*",
                      help="Index or indices of query."
                           "If not given, the full list of available queries is returned."
                           "If `all` all registered queries are executed."
                            )

    p_query = subparsers.add_parser("query", parents=[copts_parser], help=abidb_query.__doc__)
    p_query.add_argument("-q", "--query", type=_parse_dict, help="Pymongo query to filter entries")
    p_query.add_argument('-p', "--projection", type=_parse_dict_or_list, default=None, help="Pymongo projection")
    p_query.add_argument("-l", "--limit", type=int, default=0,
                         help="List number of documents returned by the query, 0 means ALL docs matching the query!")

    p_schema = subparsers.add_parser("schema", parents=[copts_parser], help=abidb_schema.__doc__)

    #p_info_oid = subparsers.add_parser("info_oid", parents=[copts_parser], help=abidb_info_oid.__doc__)

    #p_mongo_start = subparsers.add_parser(
    #    "mongo_start", parents=[copts_parser, serve_parser, worker_selector_with_default],
    #    help="Create and start new AbiPy worker connected to a MongoDB database containing FlowModels.")
    #p_mongo_start.add_argument("collection_name", type=str,
    #p_mongo_start.add_argument("-c", "--collection-name", type=str,
    #                           default="",
    #                           #required=True,
    #                           help="Name of the MongDB collection used to store FlowModels.")
    #p_mongo_start.add_argument("-s", "--scratch-dir", type=str, default=None,
    #                           help="Scratch directory in which Flows will be generated. "
    #                                "If not given, use value from ~/.abinit/abipy/config.yml")

    #p_mongo_gui = subparsers.add_parser("mongo_gui", parents=[copts_parser, worker_selector, serve_parser],
    #                                    help="Start MongoGUI.")

    return parser


#def tools():
#    from pymongo import DESCENDING, ASCENDING
#    if hasattr(args, "sort") and args.sort:
#        sort = [(args.sort, ASCENDING)]
#    elif hasattr(args, "rsort") and args.rsort:
#        sort = [(args.rsort, DESCENDING)]
#    else:
#        sort = None


#def serve_kwargs_from_options(options):
#
#    print("""
#Use:
#
#    ssh -N -f -L localhost:{port}:localhost:{port} username@your_remote_cluster
#
#to activate port forwarding.
#""")
#
#    import abipy.panels as mod
#    assets_path = os.path.join(os.path.dirname(mod.__file__), "assets")
#
#    d = dict(
#        debug=options.verbose > 0,
#        #show=False,
#        port=options.port,
#        static_dirs={"/assets": assets_path},
#        address=options.address,
#        websocket_origin="*",
#        panel_template="FastList",
#    )
#    #print("serve_kwargs:\n", pformat(d), "\n")
#    return d


#def serve(worker, options):
#    s = str(worker) if options.verbose else repr(worker)
#    print(s)
#
#    serve_kwargs = serve_kwargs_from_options(options)
#
#    if not options.daemonize:
#        worker.serve(**serve_kwargs)
#    else:
#        print(f"Running server {worker.name} in demon mode...")
#        import daemon
#        with daemon.DaemonContext():
#            worker.serve(**serve_kwargs)
#
#    return 0


class Context(BaseModel):

    mongo_connector: MongoConnector

    collection: Collection

    fmodel_cls: Type[FlowModel]

    class Config:
        arbitrary_types_allowed = True


def abidb_list(options, ctx: Context):
    print(ctx.mongo_connector.list_collection_names())

    return 0


def abidb_query(options, ctx: Context):
    print("query", options.query)
    print("projection", options.projection)
    print("limit", options.limit)

    for doc in ctx.collection.find(options.query, projection=options.projection, limit=options.limit):
        print(doc)

    return 0


def abidb_status(options, ctx: Context):
    status2oids = ctx.fmodel_cls.mongo_get_status2oids(ctx.collection)
    for status, oids in status2oids.items():
        print(str(status), len(oids))

    return 0


def abidb_in_structure(options, ctx: Context):
    df = ctx.fmodel_cls.mongo_aggregate_in_structures(ctx.collection)
    print_dataframe(df)
    return 0


def abidb_flowdata(options, ctx: Context):
    flow_data = ctx.fmodel_cls.mongo_get_flowdata(options.oid, ctx.collection)
    #print(flow_data)
    print(flow_data.yaml_dump())
    return 0


def abidb_cq(options, ctx: Context):

    queries = ctx.fmodel_cls.get_common_queries()

    if not options.cq_inds:
        print("Printing list of predefined queries with their index as `cg_inds` argument is not given\n")
        for i, query in enumerate(queries):
            print(query.to_string(title=f"== Index: {i} ==="))

        print("Use e.g. `abidb.py cq 0 1` to execute queries with a given index")
        print("      or `abidb.py cq all` to select all")
        return 0

    cg_inds = set(options.cq_inds)
    if "all" in cg_inds:
        cg_inds = set(range(len(queries)))
    else:
        cg_inds = set([int(c) for c in cg_inds])

    print(cg_inds)

    for i, query in enumerate(queries):
        if i not in cg_inds: continue
        print(query)
        df = query.get_dataframe(ctx.collection)
        print(df)
        #print_dataframe(df)

    return 0


def abidb_schema(options, ctx: Context):

    #print(fmodel_cls)
    #print(fmodel_cls.schema_json(indent=2))
    schema = ctx.fmodel_cls.schema()
    import ruamel.yaml as yaml
    print(yaml.safe_dump(schema, default_flow_style=False))
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
        raise exc
        print(exc)
        show_examples_and_exit(error_code=1)

    # loglevel is bound to the string value obtained from the command line argument.
    # Convert to upper case to allow the user to specify --loglevel=DEBUG or --loglevel=debug
    import logging
    numeric_level = getattr(logging, options.loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % options.loglevel)
    logging.basicConfig(level=numeric_level)

    if not options.collection_name:
        # Take it from external file

        print(f"Taking the collection name from file: {MAGIC_COLLECTION_FILENAME} as the `-c` option is not used")
        try:
            with open(MAGIC_COLLECTION_FILENAME, "rt") as fh:
                options.collection_name = fh.read().strip()
                print(f"collection_name set to: {options.collection_name}")

        except FileNotFoundError:
            print(f"\nCannot find `{MAGIC_COLLECTION_FILENAME}`\n"
                  f"If you don't want to specify the collection name on the command line with `-c COLL_NAME`\n",
                  f"you need to create a `{MAGIC_COLLECTION_FILENAME}` file in the current working directory\n"
                  f"with one line giving the name of the MongoDB collection that should be used."
                  )
            return 1

    #if options.use_localhost:
    #mongo_connector = MongoConnector.for_localhost(collection_name=options.collection_name)

    mongo_connector = MongoConnector.from_abipy_config(collection_name=options.collection_name)
    print(mongo_connector)

    collection = mongo_connector.get_collection()
    fmodel_cls = FlowModel.get_subclass_from_collection(collection)

    ctx = Context(**locals())

    #config = get_config()
    #if options.verbose > 2:
    #    print(options)
    #    print("global abipy config options:\n", config)

    return globals()[f"abidb_{options.command}"](options, ctx)


if __name__ == "__main__":
    sys.exit(main())
