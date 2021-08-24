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
from pprint import pprint, pformat
from queue import Queue, Empty
from socket import gethostname
from monty import termcolor
from monty.collections import AttrDict #, dict2namedtuple
from monty.string import list_strings #, marquee
from monty.json import MSONable
from monty.functools import lazy_property
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


class FlowDirsPostHandler(BaseHandler):

    def post(self):
        data = self.get_json_from_body()

        errors = []
        for flow_dir in data["flow_dirs"]:
            try:
                flow = Flow.from_file(flow_dir)
                self.worker.flow_scheduler.add_flow(flow)
            except Exception as exc:
                errors.append(str(exc))

        reply = dict(
            returncode=len(errors),
            errors=errors,
        )

        self.write(reply)


def find_free_port(address="localhost"):
    import socket
    from contextlib import closing
    with closing(socket.socket(socket.AF_INET, socket.SOCK_STREAM)) as s:
        s.bind((address, 0))
        s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        return s.getsockname()[1]


class WorkerState(AttrDict):

    @classmethod
    def new(cls, **kwargs):

        import getpass
        from socket import gethostname

        new = cls(
            name=None,
            status="init",
            pid=os.getpid(),
            address="localhost",
            port=0,
            scratch_dir=None,
            username=getpass.getuser(),
            hostname=gethostname(),
            version="0.1",
        )

        new.update(**kwargs)
        return new

    #@classmethod
    #def from_dict(cls, d):
    #    return cls(**d)

    def __str__(self):
        return self.__class__.__name__ + "\n" + pformat(self, indent=2) + "\n"

    def json_write(self, filepath):
        with open(filepath, "wt") as fp:
            json.dump(self, fp, indent=2)

    @classmethod
    def from_json_file(cls, filepath):
        with open(filepath, "rt") as fp:
            return cls(**json.load(fp))

    @classmethod
    def from_json(cls, json_string):
        return cls(**json.loads(json_string))


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
            ("/post_flow_script", PostFlowHandler, dict(worker=self)),
            ("/post_flow_dirs", FlowDirsPostHandler, dict(worker=self)),
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
                raise RuntimeError(f"There's already an AbiPy Worker serving on this machine with pid: {d['pid']}")

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

        state = WorkerState.new(
            name=self.name,
            status=status,
            pid=self.pid,
            address=self.address,
            port=self.port,
            scratch_dir=self.scratch_dir,
        )

        state.json_write(filepath)

    def serve(self, **serve_kwargs):
        from abipy.panels.core import abipanel
        abipanel()
        thread = threading.Thread(target=self.flow_scheduler.start, name="flow_scheduler", daemon=True)
        thread.start()
        #termcolor.enable(False)

        #print("serve_kwargs:", serve_kwargs)

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

    worker_state = WorkerState.new(
            name=worker_name,
            scratch_dir=scratch_dir,
    )

    with open(os.path.join(config_dir, "state.json"), "wt") as fp:
        json.dump(worker_state, fp, indent=2)

    return 0


def _get_worker_state_path_list():
    state_path_list = []

    config_dir = os.path.join(os.path.expanduser("~"), ".abinit", "abipy")
    worker_dirs = [dirname for dirname in os.listdir(config_dir) if dirname.startswith("worker_")]
    for workdir in worker_dirs:
        state_file = os.path.join(config_dir, workdir, "state.json")
        state = WorkerState.from_json_file(state_file)
        state_path_list.append((state, state_file))

    return state_path_list


def list_localhost_workers():

    for (worker_state, filepath) in _get_worker_state_path_list():
        print(f"In state_file: `{filepath}`")
        pprint(worker_state, indent=2)
        print(2 * "\n")


def discover_local_workers():
    worker_states = [state for  (state, _) in _get_worker_state_path_list()]
    clients = WorkerClients.from_json_file(empty_if_not_file=True)
    clients.update_from_worker_states(worker_states)
    return clients


def rdiscover(hostnames):
    #
    # From https://docs.fabfile.org/en/2.6/api/transfer.html#fabric.transfer.Transfer
    #
    # Most SFTP servers set the remote working directory to the connecting user’s home directory,
    # and (unlike most shells) do not expand tildes (~).
    #
    # For example, instead of saying get("~/tmp/archive.tgz"), say get("tmp/archive.tgz").

    from fabric import Connection
    from io import BytesIO
    for host in hostnames:
        print(f"Contacting host {host}. It may take some time ...")
        c = Connection(host)
        try:
            result = c.run("ls ~/.abinit/abipy")
            files = result.stdout.split()
            worker_dirs = [f for f in files if f.startswith("worker_")]
            worker_dirs = [os.path.join(".abinit/abipy", d) for d in worker_dirs]
        except Exception as exc:
            print(exc)
            continue

        worker_states = []
        for w in worker_dirs:
            state_path = os.path.join(w, "state.json")
            try:
                strio = BytesIO()
                c.get(state_path, local=strio)
                json_bstring = strio.getvalue()
                strio.close()
                worker_states.append(WorkerState.from_json(json_bstring))

            except IOError as exc:
                print(f"Cannot find state.json file: {host}@{state_path}. Ignoring error.")
                print(exc)

        clients = WorkerClients.from_json_file(empty_if_not_file=True)
        clients.update_from_worker_states(worker_states)
        return clients


class _PortForwarder:

    def __init__(self, remote_port, local_port=None, destination=None, kill_ssh=True):
        self.destination = destination
        self.pid = None

        if local_port == ":automatic:" or local_port is None:
            # Need ssh port-forwarding.
            self.local_port = find_free_port()
            self.kill_ssh = kill_ssh

            ssh_cmd = f"ssh -N -f -L localhost:{self.local_port}:localhost:{remote_port} {destination}"
            print(f"Executing ssh port-forwarding with: `{ssh_cmd}`")
            p = subprocess.run(ssh_cmd, shell=True)
            if p.returncode != 0:
                print(f"WARNING: {ssh_cmd} returned {p.returncode}")
            self.pid = p.pid
            self.url = f"http://localhost:{self.local_port}"

        else:
            # Running locally without ssh port forwarding.
            self.local_port = int(local_port)
            self.kill_ssh = False
            assert kill_ssh == False
            self.url = f"http://localhost:{self.local_port}"

    @classmethod
    def from_local_port(cls, local_port):
        return cls(remote_port=local_port, local_port=local_port, kill_ssh=False)

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        if self.pid is not None and self.kill_ssh:
            import signal
            try:
                os.kill(pid, signal.SIGTERM)
            except:
                os.kill(pid, signal.SIGKILL) # kill -9
            finally:
                pass


class WorkerClient(MSONable):

    def __init__(self, worker_state, default=False, timeout=None,
                 is_local_worker=None, ssh_destination=":automatic:", ssh_local_port=":automatic:"):

        self.worker_state = WorkerState(**worker_state.copy())

        self.default = bool(default)
        self.timeout = timeout

        if is_local_worker is None:
            # Well, this is not perfect as one may have a laptop with the same hostname as the remote cluster.
            self.is_local_worker = self.worker_state.hostname == gethostname()
        else:
            self.is_local_worker = bool(is_local_worker)

        self.ssh_destination = ssh_destination
        self.ssh_local_port = ssh_local_port

    @pmg_serialize
    def as_dict(self) -> dict:
        return dict(
                worker_state=self.worker_state,
                default=self.default,
                timeout=self.timeout,
                is_local_worker=self.is_local_worker,
                ssh_destination=self.ssh_destination,
                ssh_local_port=self.ssh_local_port,
        )

    def port_forwarder(self, kill_ssh=True):

        if not self.worker_state.status != "running":
            raise RuntimeError(f"Server status is `{self.worker_state.status} while it should be `running``")

        if self.is_local_worker:
            return _PortForwarder.from_local_port(self.worker_state.port)

        # Need ssh port-forwarding.
        #url = "http://localhost:{local_port}"

        if self.ssh_destination == ":automatic:":
            destination = "{self.worker_state.username}@{self.worker_state.remote_hostname}"
        else:
            destination = self.destination

        #if self.ssh_local_port == ":automatic:":
        #    local_port = find_free_port()
        #else:
        #    local_port = int(self.ssh_local_port)

        #ssh_cmd = f"ssh -N -f -L localhost:{local_port}:localhost:{self.worker_state.remote_port} {destination}"
        #print(f"Executing ssh port-forwarding with: `{ssh_cmd}`")
        #p = subprocess.run(ssh_cmd, shell=True)
        #if p.returncode != 0:
        #    print(f"WARNING: {ssh_cmd} returned {p.returncode}")
        # TODO: Should we kill the process?
        #p.pid

        return _PortForwarder(remote_port=self.worker_state.remote_port,
                              local_port=self.ssh_local_port,
                              destination=destination, kill_ssh=kill_ssh)

    #@lazy_property
    #def server_url(self):

    def update_state(self, worker_state):
        self.worker_state = worker_state.copy()

    def __str__(self):
        return self.to_json()

    def to_json(self) -> str:
        """
        Returns a json string representation of the MSONable object.
        """
        # TODO: monty.to_json should accept **kwargs
        from monty.json import MontyEncoder
        return json.dumps(self, cls=MontyEncoder, indent=2)

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

        with self.port_forwarder() as pf:
            url = f"{pf.url}/post_flow_script"
            print(f"Sending `{filepath}` script to url: `{url}` ...\n")

            r = requests.post(url=url, json=data, timeout=self.timeout)
            print("ok:", r.ok, "status code:", r.status_code)
            return r.json()

    def send_flow_dirs(self, flow_dirs):
        if not self.is_local_worker:
            raise RuntimeError("You cannot add a Flow to the remote worker {self.worker_state.name}")

        with self.port_forwarder() as pf:
            url = f"{pf.url}/post_flow_dirs"
            print(f"Sending `{flow_dirs}` directories to url: `{url}` ...\n")

            r = requests.post(url=url, json=dict(flow_dirs=list_strings(flow_dirs),
                              timeout=self.timeout))
            print("ok:", r.ok, "status code:", r.status_code)
            return r.json()

    def send_kill_message(self):

        with self.port_forwarder() as pf:
            print(f"Sending kill message to: {furl} ...\n")
            r = requests.post(url=f"{pf.url}/action", json=dict(action="kill"), timeout=self.timeout)
            print(r.text)
            return r.json()

    def get_json_state(self):
        with self.port_forwarder() as pf:
            url = f"{pf.url}/json_state"
            print(f"Sending request to {url} ...\n")
            r = requests.get(url=url)
            print("ok:", r.ok, "status code:", r.status_code)
            print(r.text)
            return r.json()

    def open_webgui(self):
        with self.port_forwarder(kill_ssh=False) as pf:
            print(self.worker_state)
            import webbrowser
            webbrowser.open_new_tab(pf.url)


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

        # worker_name must be unique.
        for worker_name, count in Counter((w.worker_state.name for w in self)).items():
            if count > 1:
                app(f"Server name: `{worker_name}` appears: {count} times. This is forbidden!")

        # Only zero or one default server is allowed.
        count = sum(1 if w.default == True else 0 for w in self)
        if count not in (0, 1):
            app(f"default=True appears `{count}` times. This is forbidden as only one default worker is allowed!")

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
        for state in worker_states:
            client = self.select_from_worker_name(state.name, allow_none=True)
            if client is not None:
                client.update_state(state)
            else:
                self.append(WorkerClient(state))

        self._validate()
        self.write_json_file()

    def select_from_worker_name(self, worker_name, allow_none=False):
        if worker_name is None:
            for client in self:
                if client.default: return client
            raise ValueError("Cannot find default worker in list! Use `abiw.py set_default WORKER_NAME`")

        for client in self:
           if client.worker_state.name == worker_name: return client

        if allow_none: return None
        raise ValueError(f"Cannot find client with name `{worker_name}` in list!")

    def set_default(self, worker_name):
        the_one = self.select_from_worker_name(worker_name)
        for client in self:
            client.default = False
        the_one.default = True

        self._validate()
        self.write_json_file()
