"""[summary]
"""
import sys
import os
import threading
import json
import time
import subprocess
import tempfile
import tornado.web
import tornado.escape
import requests
import panel as pn

from datetime import datetime
from pprint import pprint #, pformat
from queue import Queue, Empty
from monty import termcolor
#from monty.collections import AttrDict, dict2namedtuple
from monty.json import MSONable
from pymatgen.util.serialization import pmg_serialize
from abipy.flowtk.flows import Flow
from abipy.flowtk.launcher import MultiFlowScheduler


def yaml_safe_load_path(path):
    import ruamel.yaml as yaml
    with open(path, "rt") as fh:
        return yaml.YAML(typ='safe', pure=True).load(fh.read())


class BaseHandler(tornado.web.RequestHandler):

    def initialize(self, worker):
        self.worker = worker

    def get_json_from_body(self):
        try:
            data = tornado.escape.json_decode(self.request.body)
            return data
        except ValueError:
            # TODO: handle the error
            raise


class ActionHandler(BaseHandler):

    def post(self):
        data = self.get_json_from_body()
        action = data.get("action", None)
        if action is None: return
        reply = {}

        if action == "kill":
            reply["message"] = "Received kill action. Aborting now!"
            self.write(reply)
            # TODO: Should lock and give a change to the secondary thread to complete operations
            time.sleep(4)
            sys.exit(0)

        #elif action == "rm_errored_flows":
        #    self.worker.scheduler.rm_flows_with_status("Error")
        #    reply["message"] = f"Removed: {nflows} Errored Flows"

        else:
            reply["message"] = f"Unknown action: {action}"

        self.write(reply)


class PostFlowHandler(BaseHandler):

    def post(self):
        data = self.get_json_from_body()

        pyscript_basename = data["pyscript_basename"]
        pyscript_string = data["pyscript_string"]

        scratch_dir = self.worker.scratch_dir
        filepath = os.path.join(scratch_dir, pyscript_basename)
        with open(filepath, "wt") as fh:
            fh.write(pyscript_string)

        now = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
        workdir = tempfile.mkdtemp(prefix=f"Flow-{now}", dir=scratch_dir)

        # FIXME: problem with manager and possible collisions in filepath
        p = subprocess.Popen(f"python {filepath} -w {workdir}", stdout=subprocess.PIPE, shell=True)
        out, err = p.communicate()

        reply = dict(
            returncode=p.returncode,
            out=out.decode("utf-8") if out is not None else None,
            err=err.decode("utf-8") if err is not None else None,
            workdir=workdir,
        )

        if p.returncode == 0:
            print(f"Running flow in workdir: {workdir}")
            flow = Flow.from_file(workdir)
            self.worker.flow_scheduler.add_flow(flow)

        self.write(reply)


class _JsonHandler(BaseHandler):

    def prepare(self):
        self.set_header(name="Content-Type", value="application/json")


class JsonStateHandler(_JsonHandler):

    def get(self):
        # FIXME: This is broken now
        json_data = self.worker.flow_scheduler.to_json()
        print("json_data", json_data)
        self.write(json_data)


def find_free_port(address):
    import socket
    from contextlib import closing
    with closing(socket.socket(socket.AF_INET, socket.SOCK_STREAM)) as s:
        s.bind((address, 0))
        s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        return s.getsockname()[1]


class WorkerServer:

    def __init__(self, name: str, sched_options: dict, scratch_dir: str,
                 address="localhost", port=0):
        """
        Args:
            name: The name of the Worker. Must be unique.
            sched_options:
            scratch_dir:
            address:
                The address the server should listen on for HTTP requests.
            port:
                 Allows specifying a specific port
                 Port 0 means to select an arbitrary unused port
        """
        self.name = name
        self.address, self.port = address, port
        self.pid = os.getpid()

        # url --> callables returning panel objects.
        self.routes = {
            "/": self.serve_homepage,
            r"/flow/\d+": self.serve_panel_flow,
        }

        # url --> tornado handler
        self.extra_patterns = [
            ("/postflow", PostFlowHandler, dict(worker=self)),
            ("/json_state", JsonStateHandler, dict(worker=self)),
            ("/action", ActionHandler, dict(worker=self)),
        ]
        self.config_dir = os.path.join(os.path.expanduser("~"), ".abinit", "abipy", f"worker_{self.name}")
        if not os.path.exists(self.config_dir):
            os.mkdir(self.config_dir)

        self.scratch_dir = scratch_dir
        if not os.path.isdir(scratch_dir):
            raise ValueError(f"Scratch directory: `{scratch_dir}` does not exist!")

        sqldb_path = os.path.join(self.config_dir, "flows.db")
        self.flow_scheduler = MultiFlowScheduler(sqldb_path, **sched_options)

        state_file = os.path.join(self.config_dir, "state.json")
        if os.path.exists(state_file):
            with open(state_file, "rt") as fp:
                d = json.load(fp)
            if d["status"] == "serving":
                raise RuntimeError(f"There's already a Worker serving on this machine with pid: {d['pid']}")

        # Register function atexit
        import atexit
        atexit.register(self.write_state_file)

    @classmethod
    def _get_state_path(cls, name):
        config_dir = os.path.join(os.path.expanduser("~"), ".abinit", "abipy", f"worker_{name}")
        state_filepath = os.path.join(config_dir, "state.json")
        if not os.path.exists(state_filepath):
            raise RuntimeError(f"Cannot find state file: `{state_filepath}`")

        with open(state_filepath, "rt") as fp:
            return (json.load(fp), state_filepath)

    @classmethod
    def init_from_config_dir(cls, name):
        d, path = cls._get_state_path(name)

        if d["status"] != "init":
            raise RuntimeError(f"`status` entry in json file `{path}`\n should be `init` while it is `{d['status']}`")

        config_dir = os.path.dirname(path)
        sched_options = yaml_safe_load_path(os.path.join(config_dir, "scheduler.yml"))
        #manager = TaskManager.from_file(os.path.join(config_dir, "manager.yml")
        print("sched_options", sched_options)

        return cls(d["name"], sched_options, d["scratch_dir"],
                   address=d["address"], port=d["port"])

    @classmethod
    def restart_from_config_dir(cls, name):
        d, path = cls._get_state_path(name)

        #if d["status"] != "dead":
        if d["status"] == "serving":
            raise RuntimeError(f"There's already a worker serving on this machine with pid: {d['pid']}")

        # TODO: Problem with the default manager when creating the flow.
        #from abipy.flowtk.tasks import TaskManager
        #manager = TaskManager.from_file(os.path.join(config_dir, "manager.yml")
        config_dir = os.path.dirname(path)
        sched_options = yaml_safe_load_path(os.path.join(config_dir, "scheduler.yml"))

        return cls(d["name"], sched_options, d["scratch_dir"], address=d["address"], port=d["port"])

    def write_state_file(self, status="dead", filepath=None) -> None:
        if filepath is None:
            filepath = os.path.join(self.config_dir, "state.json")

        d = dict(
            name=self.name,
            status=status,
            pid=self.pid,
            address=self.address,
            port=self.port,
            scratch_dir=self.scratch_dir,
        )

        with open(filepath, "wt") as fp:
            json.dump(d, fp, indent=2)

    def serve(self, **serve_kwargs):
        from abipy.panels.core import abipanel
        abipanel()
        thread = threading.Thread(target=self.flow_scheduler.start, name="flow_scheduler", daemon=True)
        thread.start()
        #termcolor.enable(False)

        print("serve_kwargs:", serve_kwargs)

        self.address = serve_kwargs.pop("address", self.address)
        self.port = serve_kwargs.pop("port", self.port)
        if self.port == 0:
            self.port = find_free_port(self.address)

        # Now write state.json with the actual port.
        self.write_state_file(status="serving")

        return pn.serve(self.routes, port=self.port, address=self.address,
                        extra_patterns=self.extra_patterns, **serve_kwargs)

    def __str__(self):
        #print(f"Server running at host: {self.address}, port: {self.port}")
        #print(f"Server loop running in thread: {server_thread.name}")
        lines = []
        app = lines.append
        app("pid %d" % self.pid)
        app(str(self.flow_scheduler))

        return "\n".join(lines)

    #def remove_flows(self):
    #    """This requires locking the SQLite database."""

    def serve_homepage(self):
        d = self.flow_scheduler.groupby_status()
        md_lines = []
        if d:
            for status, values in d.items():
                if status == "error": status = "Errored"
                md_lines.append(f"## {status.capitalize()} Flows:\n")
                for row in values:
                    workdir, flow_id = row["workdir"], row["flow_id"]
                    md_lines.append(f"- [flow_id: {flow_id}, workdir: {workdir}](/flow/{flow_id})\n")
        else:
            md_lines.append("## Empty Flow list!")

        #import platform
        #from socket import gethostname
        #system, node, release, version, machine, processor = platform.uname()
        #cprint("Running on %s -- system %s -- ncpus %s -- Python %s -- %s" % (
        #      gethostname(), system, ncpus_detected, platform.python_version(), _my_name),
        #      'green', attrs=['underline'])

        from abipy.panels.viewers import JSONViewer
        return pn.Column(
                "# AbiPy Worker Homepage",
                pn.pane.Markdown("\n".join(md_lines)),
                #JSONViewer(self.flow_scheduler.to_json()),
                sizing_mode="stretch_width",
        )

    def serve_panel_flow(self):
        # URL example: /flow/123 where 123 is the flow id.
        tokens = pn.state.app_url.split("/")
        #print("In serve_panel_flow with app_url", pn.state.app_url, "tokens:", tokens)
        flow_id = int(tokens[2])
        flow, status = self.flow_scheduler.get_flow_status_by_id(flow_id)
        if flow is not None: return flow.get_panel()

        return pn.pane.Alert(f"Cannot find Flow with node ID: {flow_id}", alert_type="danger")


def create_new_worker(worker_name: str, scratch_dir: str):

    abipy_dir = os.path.join(os.path.expanduser("~"), ".abinit", "abipy")
    config_dir = os.path.join(abipy_dir, f"worker_{worker_name}")
    errors = []
    eapp = errors.append

    if os.path.exists(config_dir):
        eapp(f"Directory `{config_dir}` already exists!")

    scheduler_path = os.path.join(abipy_dir, "scheduler.yml")
    manager_path = os.path.join(abipy_dir, "manager.yml")
    if not os.path.exists(scheduler_path):
        eapp("Cannot find scheduler file at `{scheduler_path}` to initialize the worker!")
    if not os.path.exists(manager_path):
        eapp("Cannot find manager file at `{manager_path}` to initialize the worker!")
    if not os.path.isdir(scratch_dir):
        eapp("scratch_dir at `{scratch_dir}` is not a directory!")

    if errors:
        print("The following problems have been detected:")
        for i, err in enumerate(errors):
            print(f"\t [{i+1}] {err}")
        return 1

    print(f"""
Creating new worker directory `{config_dir}.`
Copying manager.yml and scheduler.yml files into it.
Now you can use:

    abiw.py start {worker_name}

to start the AbiPy worker.
Then use:

    abiw.py ldiscover

to update the list of local clients.
""")
    os.mkdir(config_dir)
    from shutil import copy
    copy(scheduler_path, config_dir)
    copy(manager_path, config_dir)

    state_dict = dict(
        name=worker_name,
        status="init",
        address="localhost",
        port=0,
        scratch_dir=scratch_dir,
    )
    with open(os.path.join(config_dir, "state.json"), "wt") as fp:
        json.dump(state_dict, fp, indent=2)

    return 0


def _get_worker_state_path_list():
    state_path_list = []

    config_dir = os.path.join(os.path.expanduser("~"), ".abinit", "abipy")
    worker_dirs = [dirname for dirname in os.listdir(config_dir) if dirname.startswith("worker_")]
    for workdir in worker_dirs:
        state_file = os.path.join(config_dir, workdir, "state.json")
        with open(state_file, "rt") as fp:
            d = json.load(fp)
            state_path_list.append((d, state_file))

    return state_path_list


def list_workers_on_localhost():

    for (state, filepath) in _get_worker_state_path_list():
        print(f"In state_file: `{filepath}`")
        pprint(state)
        print(2 * "\n")


def discover_local_workers():
    worker_states = [state for  (state, _) in _get_worker_state_path_list()]
    clients = WorkerClients.from_json_file(empty_if_not_file=True)
    clients.update_from_worker_states(worker_states)


def rdiscover(hostnames):


    # From https://docs.fabfile.org/en/2.6/api/transfer.html#fabric.transfer.Transfer
    # Most SFTP servers set the remote working directory to the connecting user’s home directory,
    # and (unlike most shells) do not expand tildes (~).
    #
    # For example, instead of saying get("~/tmp/archive.tgz"), say get("tmp/archive.tgz").

    from fabric import Connection
    from io import StringIO, BytesIO
    for host in hostnames:
        print(f"For host {host}")
        c = Connection(host)
        try:
            result = c.run("ls ~/.abinit/abipy")
            files = result.stdout.split()
            worker_dirs = [f for f in files if f.startswith("worker_")]
            #print("worker_dirs:", worker_dirs)
            #worker_dirs = [os.path.join("~/.abinit/abipy", b) for b in worker_dirs]
            worker_dirs = [os.path.join(".abinit/abipy", b) for b in worker_dirs]
        except Exception as exc:
            print(exc)
            continue

        worker_states = []
        for w in worker_dirs:
            path = os.path.join(w, "state.json")
            try:
                #strio = StringIO()
                strio = BytesIO()
                c.get(path, local=strio)
                json_bstring = strio.getvalue()
                strio.close()
                #print("json_bstring", json_bstring)
                worker_states.append(json.loads(json_bstring))

            except IOError as exc:
                print(f"Cannot find state.json file: {host}@{path}. Ignoring error.")
                print(exc)

        clients = WorkerClients.from_json_file(empty_if_not_file=True)
        clients.update_from_worker_states(worker_states)


class WorkerClient(MSONable):

    # TODO: Add server status?

    def __init__(self, server_name, server_address, server_port, default=False, timeout=None):
        self.server_name = server_name
        self.server_address = server_address
        self.server_port = server_port
        self.server_url = f"http://{server_address}:{server_port}"
        self.timeout = timeout
        self.default = bool(default)

    #@classmethod
    #def from_server_name(cls, name: str)
    #    clients = WorkerClients.from_json_file()
    #    return clients.get_client_name(name)

    #@classmethod
    #def from_dict(cls, d: dict):
    #    return cls(**d)

    @pmg_serialize
    def as_dict(self) -> dict:
        return {k: getattr(self, k) for k in [
                "server_name", "server_address", "server_port",
                "default", "timeout"
                ]}

    def __str__(self):
        return self.to_json()

    def send_pyscript(self, filepath: str):
        """
        Send a python script to the server in order to build a new Flow.
        """
        filepath = os.path.expanduser(filepath)

        with open(filepath, "rt") as fp:
            data = dict(
                pyscript_basename=os.path.basename(filepath),
                pyscript_string=fp.read(),
            )

        url=f"{self.server_url}/postflow"
        #url = "http://localhost:60073/postflow"
        print(f"Sending `{filepath}` script to url: `{url}`...\n\n")

        r = requests.post(url=url, json=data, timeout=self.timeout)
        print("ok:", r.ok, "status code:", r.status_code)
        return r.json()

    def send_kill_message(self):
        print(f"Sending kill message to server: {self.server_url}")
        data = dict(action="kill")
        r = requests.post(url=f"{self.server_url}/action", json=data, timeout=self.timeout)
        print(r.text)
        return r.json()

    def get_json_state(self):
        url = f"{self.server_url}/json_state"
        #url = "http://localhost:60073/json_state"
        print(f"Sending request to {url}")
        r = requests.get(url=url)
        print("ok:", r.ok, "status code:", r.status_code)
        print(r.text)
        return r.json()

    def open_webgui(self):
        import webbrowser
        webbrowser.open_new_tab(self.server_url)


class WorkerClients(list, MSONable):

    @classmethod
    def from_json_file(cls, empty_if_not_file=False, filepath=None):
        if filepath is None:
            filepath = os.path.join(os.path.expanduser("~"), ".abinit", "abipy", "clients.json")
        else:
            filepath = os.path.expanduser(filepath)

        new = cls()
        if empty_if_not_file and not os.path.exists(filepath):
            return new

        with open(filepath, "rt") as fp:
            d = json.load(fp)

        for w in d["clients"]:
            new.append(WorkerClient.from_dict(w))

        new._validate()
        return new

    def _validate(self):
        err_lines = []
        app = err_lines.append
        from collections import Counter

        # server_name must be unique.
        for server_name, count in Counter((w.server_name for w in self)).items():
            if count > 1:
                app(f"Server name: `{server_name}` appears: {count} times. This is forbidden!")

        # Only one default server is allowed.
        count = sum(1 if w.default == True else 0 for w in self)
        if count != 1:
            app(f"default=True appears `{count}` times. This is forbidden as only one default server_name is allowed!")

        if err_lines:
            raise RuntimeError("\n".join(err_lines))

    def __str__(self):
        lines = [str(w) for w in self]
        return (2 * "\n").join(lines)

    def as_dict(self):
        return {"clients": [w.as_dict() for w in self]}

    def write_json_file(self, filepath=None) -> None:
        if filepath is None:
            filepath = os.path.join(os.path.expanduser("~"), ".abinit", "abipy", "clients.json")
        else:
            filepath = os.path.expanduser(filepath)

        with open(filepath, "wt") as fp:
           json.dump(self.as_dict(), fp, indent=2)

    def update_from_worker_states(self, worker_states):
        print("in update clients with worker_states:", worker_states)
        for state in worker_states:
            print(state)
            server_name = state["name"]
            client = self.select_from_name(server_name, allow_none=True)
            if client is not None:
                client.server_address = state["address"]
                client.server_port = state["port"]
            else:
                new_client = WorkerClient(server_name, state["address"], state["port"])
                self.append(new_client)

        #print(self)
        self._validate()
        self.write_json_file()

    def select_from_name(self, server_name, allow_none=False):
        if server_name is None:
            for client in self:
                if client.default: return client
            raise ValueError("Cannot find default worker in list!")

        for client in self:
           if client.server_name == server_name: return client

        if allow_none: return None
        raise ValueError(f"Cannot find client with name `{server_name}` in list!")

    def set_default(self, server_name):
        the_one = self.select_from_name(server_name)
        for client in self:
            client.default = False
        the_one.default = True

        self._validate()
        self.write_json_file()
