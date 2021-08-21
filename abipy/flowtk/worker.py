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
from monty.json import MSONable
from pymatgen.util.serialization import pmg_serialize
from abipy.flowtk.flows import Flow
from abipy.flowtk.launcher import MultiFlowScheduler


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
            time.sleep(4)
            sys.exit(0)
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


class JsonStatusHandler(BaseHandler):

    def prepare(self):
        self.set_header(name="Content-Type", value="application/json")

    def get(self):
        json_data = self.worker.flow_scheduler.to_json()
        print("json_data", json_data)
        self.write(json_data)



class WorkerServer: # (MSONable):

    def __init__(self, name: str, address, port, sched_options: dict, scratch_dir: str,
                 flow_scheduler=None):
        """
        Args:
            address:
                The address the server should listen on for HTTP requests.
            port:
                 Allows specifying a specific port
            sched_options:
            scratch_dir:
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
            ("/json_status", JsonStatusHandler, dict(worker=self)),
            ("/action", ActionHandler, dict(worker=self)),
        ]
        self.config_dir = os.path.join(os.path.expanduser("~"), ".abinit", "abipy", f"worker_{self.name}")
        if not os.path.exists(self.config_dir):
            os.mkdir(self.config_dir)

        self.sched_options = sched_options
        self.scratch_dir = scratch_dir
        if not os.path.isdir( scratch_dir):
            raise ValueError(f"Scratch directory: `{scratch_dir}` does not exist!")

        if flow_scheduler is None:
            sqldb_path = os.path.join(self.config_dir, "flows.db")
            self.flow_scheduler = MultiFlowScheduler(sqldb_path, **self.sched_options)
        else:
            self.flow_scheduler = flow_scheduler

        status_file = os.path.join(self.config_dir, "status.json")
        if os.path.exists(status_file):
            with open(status_file, "rt") as fp:
                d = json.load(fp)

            if d["status"] == "serving":
                raise RuntimeError(f"There's already a Worker serving on this machine with pid: {d['pid']}")

        # Register function atexit
        import atexit
        atexit.register(self.write_status_file)

    def from_config_dir(cls, name):
        config_dir = os.path.join(os.path.expanduser("~"), ".abinit", "abipy", f"worker_{self.name}")
        status_file = os.path.join(config_dir, "status.json")
        if not os.path.exists(status_file):
            raise RuntimeError(f"Cannot file status file: {status_file}")

        with open(status_file, "rt") as fp:
            d = json.load(fp)

        if d["status"] == "serving":
            raise RuntimeError(f"There's already a Worker serving on this machine with pid: {d['pid']}")

        #new = cls(name: str, address, port, sched_options: dict, scratch_dir: str,
        #          flow_scheduler=None)
        #return new

    def write_status_file(self, status="dead", filepath=None) -> None:
        if filepath is None:
            filepath = os.path.join(self.config_dir, f"status.json")

        d = dict(
            name=self.name,
            status=status,
            pid=self.pid,
            address=self.address,
            port=self.port,
            sched_options=self.sched_options,
            scratch_dir=self.scratch_dir,
        )

        with open(filepath, "wt") as fp:
            json.dump(d, fp, indent=2)

    def serve(self, **serve_kwargs):

        self.write_status_file(status="serving")

        from abipy.panels.core import abipanel
        abipanel()
        thread = threading.Thread(target=self.flow_scheduler.start, name="flow_scheduler", daemon=True)
        thread.start()
        #termcolor.enable(False)
        serve_kwargs.update(address=self.address, port=self.port)
        return pn.serve(self.routes, extra_patterns=self.extra_patterns, **serve_kwargs)

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
                md_lines.append(f"## {status.capitalize()} Flows:\n")
                for row in values:
                    workdir, flow_id = row["workdir"], row["flow_id"]
                    md_lines.append(f"- [flow_id: {flow_id}, workdir: {workdir}](/flow/{flow_id})\n")
        else:
            md_lines.append("## Empty Flow list!")

        from abipy.panels.viewers import JSONViewer
        return pn.Column(
                "# Abipy Worker Homepage",
                pn.pane.Markdown("\n".join(md_lines)),
                #JSONViewer(self.flow_scheduler.to_json()),
                sizing_mode="stretch_width",
        )

    def serve_panel_flow(self):
        # URL example: /flow/123 where 123 is the flow id.
        tokens = pn.state.app_url.split("/")
        print("In serve_panel_flow with app_url", pn.state.app_url, "tokens:", tokens)
        flow_id = int(tokens[2])
        flow, status = self.flow_scheduler.get_flow_status_by_id(flow_id)
        if flow is not None: return flow.get_panel()

        return pn.pane.Alert(f"Cannot find Flow with node ID: {flow_id}", alert_type="danger")


def list_workers_on_local_machine():

    config_dir = os.path.join(os.path.expanduser("~"), ".abinit", "abipy")
    worker_dirs = [dirname for dirname in os.listdir(config_dir) if dirname.startswith("worker_")]
    for workdir in worker_dirs:
        status_file = os.path.join(config_dir, workdir, "status.json")
        with open(status_file, "rt") as fp:
            print(f"In status_file: `{status_file}`")
            d = json.load(fp)
            pprint(d)
            print(2 * "\n")


class WorkerClient(MSONable):

    def __init__(self, server_address, server_port, server_name, priority=1, default=False, timeout=None):
        self.server_address = server_address
        self.server_port = server_port
        self.server_url = f"http://{server_address}:{server_port}"
        self.timeout = timeout
        self.server_name = server_name
        self.priority = int(priority)
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
                "server_address", "server_port", "server_name",
                "priority", "default", "timeout"
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
        url = "http://localhost:60073/postflow"
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

    def get_json_status(self):
        url = f"{self.server_url}/json_status"
        url = "http://localhost:60073/json_status"
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
    def from_json_file(cls, filepath=None):
        if filepath is None:
            filepath = os.path.join(os.path.expanduser("~"), ".abinit", "abipy", "clients.json")
        else:
            filepath = os.path.expanduser(filepath)

        with open(filepath, "rt") as fp:
            d = json.load(fp)

        new = cls()
        for w in d["workers"]:
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
        return {"workers": [w.as_dict() for w in self]}

    def write_json_file(self, filepath=None) -> None:
        if filepath is None:
            filepath = os.path.join(os.path.expanduser("~"), ".abinit", "abipy", "workers.json")
        else:
            filepath = os.path.expanduser(filepath)

        with open(filepath, "wt") as fp:
           json.dump(self.as_dict(), fp, indent=2)

    def select_from_name(self, server_name):
        if server_name is None:
            for worker in self:
                if worker.default: return worker
            raise ValueError("Cannot find default worker in list!")

        for worker in self:
           if worker.server_name == server_name: return worker
        raise ValueError(f"Cannot find worker with name `{server_name}` in list!")


    def set_default(self, server_name):
        the_one = self.select_from_name(server_name)
        for worker in self:
            worker.default = False
        the_one.default = True
        self.write_json_file()
