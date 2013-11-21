#!/usr/bin/env python
from __future__ import division, print_function

import abc
import os
import time
import subprocess
import collections
import cmd
import readline
import pprint
import cStringIO as StringIO 
import json
import yaml
import socket
import paramiko


class Cluster(object):
    __metaclass__ = abc.ABCMeta

    def __init__(self, hostname, username, password, workdir):
        self.hostname = hostname
        self.username = username
        self.password = password
        self.workdir = workdir

        # TODO
        assert password and username

        self.port = 22
        self.timeout = 10 # Timeout in seconds.

    @classmethod
    def from_qtype(cls, qtype):
        for subcls in cls.__subclasses__():
            if subcls.qtype == qtype:
                return subcls

        raise ValueError("Wrong value for qtype %s" % qtype)

    def __str__(self):
        return "%s@%s:%s" % (self.username, self.hostname, self.__class__.__name__)

    #def __del__(self):
    #    self.disconnect

    def path_inworkdir(self, basename):
        return os.path.join(self.workdir, basename)

    def ping(self):
        cmd = ["ping", "-c1", "-W100", "-t1", self.hostname]
        ping_response = subprocess.Popen(cmd, stdout=subprocess.PIPE).stdout.read()
        return ping_response

    def prefix_str(self, s):
        lines = []
        lines.extend(("[%s] %s" % (self.hostname, l) for l in s.splitlines()))

        if not lines:
            lines = ["[%s] (Empty)" % self.hostname]

        return "\n".join(lines)

    def disconnect(self):
        try:
            self.ssh.close()
        except:
            pass
        try:
            self.sftp.close()
        except:
            pass

    @property
    def has_slurm(self):
        return self.qtype == "slurm"

    @property
    def has_sge(self):
        return self.qtype == "sge"

    def _ssh_connect(self):
        ssh_client = MySSHClient()
        ssh_client.load_system_host_keys()
        ssh_client.set_missing_host_key_policy(paramiko.WarningPolicy)

        ssh_client.connect(self.hostname, port=self.port, username=self.username, password=self.password)
        return ssh_client

    @property
    def ssh(self):
        if hasattr(self, "_ssh") and self._ssh.is_connected:
            return self._ssh
        else:
            self._ssh = self._ssh_connect()
            return self._ssh

    def ssh_close(self):
        if hasattr(self, "_ssh") and self._ssh.is_connected:
            self._ssh.close()

    def ssh_exec(self, command):
        stdin, stdout, stderr = self.ssh.exec_command(command, timeout=None)
        return SSHResult(stdout, stderr)
         
    def _sftp_connect(self):
        sftp = self.ssh.open_sftp()
        sftp.__class__ = MySFTPClient
        sftp._is_connected = True
        return sftp

    @property
    def sftp(self):
        if hasattr(self, "_sftp") and self._sftp.is_connected:
            return self._sftp
        else:
            self._sftp = self._sftp_connect()
            return self._sftp

    def sftp_close(self):
        if hasattr(self, "_sftp") and self._sftp.is_connected:
            self.sftp.close()

    #def sftp_get(self, source , dest):
    #    self.sftp.get(source, dest)

    def invoke_shell(self, timeout=None):
        shell = self.ssh.invoke_shell()
        shell.__class__ = MyShell
        shell.settimeout(timeout if timeout is not None else self.timeout)
        return shell

    @abc.abstractmethod
    def get_user_jobs(self):
        """Info of the user jobs."""

    @abc.abstractmethod
    def get_qinfo(self):
        """Info of the queue."""


class SlurmCluster(Cluster):
    qtype = "slurm"

    def get_user_jobs(self):
        result = self.ssh_exec(command="squeue -u %s" % self.username)
        return self.prefix_str(result.out)

    def get_qinfo(self):
        result = self.ssh_exec(command="sinfo")
        return self.prefix_str(result.out)


class SgeCluster(Cluster):
    qtype = "sge"

    def get_user_jobs(self):
        result = self.ssh_exec(command="qstat -u %s" % self.username)
        return self.prefix_str(result.out)

    def get_qinfo(self):
        result = self.ssh_exec("qhost")
        return self.prefix_str(result.out)


class MySSHClient(paramiko.SSHClient):

    @property
    def is_connected(self):
        return self._is_connected

    def connect(self, hostname, **kwargs):
        super(MySSHClient, self).connect(hostname, **kwargs)
        self._is_connected = False

    def close(self):
        if not self.is_connected: return
        super(MySSHClient, self).close()
        self._is_connected = False


class MySFTPClient(paramiko.SFTPClient):

    @property
    def is_connected(self):
        return self._is_connected

    def close(self):
        if not self.is_connected: return
        super(MySFTPClient, self).close()
        self._is_connected = False


class MyShell(paramiko.Channel):

    def myexec(self, s):
        if not s.endswith("\n"): s += "\n"
        s = "exec " + s

        # http://stackoverflow.com/questions/5342402/can-i-get-the-exit-code-of-a-command-executed-in-a-subshell-via-ssh
        super(MyShell, self).sendall(s)

        #stdout = self.makefile()
        #stderr = self.makefile_stderr()

        # Wait for it.....
        #count = 0
        #while not self.recv_ready():
        #    time.sleep(1)
        #    count += 1
        #    if count >= 6:
        #        print('time out') #TODO: add exeption here 
        #time.sleep(10)

        #retcode = self.recv_exit_status()
        #print(retcode)
        #result = ""
        #while self.recv_ready():
        ##while self.recv_ready() and not self.exit_status_ready():
        #    o = self.recv(1024)
        #    #if not len(o): break
        #    #if o:
        #    print("o", o)
        #    result += o
        #self.settimeout(30)

        stdout, stderr = StringIO.StringIO(), StringIO.StringIO()

        #try: #TODO: add exeption here 
        # Need to read the entire buffer to caputure output
        while not self.exit_status_ready():

            if self.recv_ready():
                out_buff = self.recv(1024)
                while out_buff:
                    #print("Indside stdout", out_buff)
                    stdout.write(out_buff)
                    out_buff = self.recv(1024)


            if self.recv_stderr_ready():
                error_buff = self.recv_stderr(1024)
                while error_buff:
                    #print("Indside stderr", error_buff)
                    stderr.write(error_buff)
                    error_buff = self.recv_stderr(1024)

        #except socket.timeout:
        #  raise socket.timeout

        stdout.seek(0), stderr.seek(0)

        result = SSHResult(stdout, stderr, retcode=self.recv_exit_status())
        #self.close()
        return result

    def exit(self):
        self.sendall('exit\n')

        # Wait for it.....
        count = 0
        while not self.exit_status_ready():
            time.sleep(1)
            count += 1
            if count >= 6:
                print('time out')#TODO: add exeption here 

        return self.recv_exit_status()


class SSHResult(object):

    def __init__(self, stdout, stderr, retcode=None):
        self.stdout, self.stderr = stdout, stderr

        if retcode is None:
            self.retcode = stdout.channel.recv_exit_status()
        else:
            self.retcode = retcode

    def __str__(self):
        s  = "out: " + self.out 
        if self.err: s += "\nerr: " + self.err
        return s

    @property
    def out(self):
        try:
            self._out
        except AttributeError:
            self._out = self.stdout.read()
            return self._out

    @property
    def err(self):
        try:
            self._err
        except AttributeError:
            self._err = self.stderr.read()
            return self._err


def read_clusters():
    with open("clusters.yml", "r") as fh:
        conf = yaml.load(fh)

    global_password = conf.get("password", None)
    global_username = conf.get("username", None)
    d = conf["clusters"]

    clusters = {}
    for hostname, params in d.items():
        username = params.get("username", global_username)
        password = params.get("password", global_password)
        cls = Cluster.from_qtype(params["qtype"])
        assert hostname not in clusters

        clusters[hostname] = cls(hostname, username, password, params["workdir"])

    return clusters


_ALL_CLUSTERS = read_clusters()


def straceback():
    """Returns a string with the traceback."""
    import traceback
    return traceback.format_exc()


class RunCommand(cmd.Cmd, object):
    """ Simple shell to run a command on the host """

    prompt='abife> '

    def __init__(self):
        cmd.Cmd.__init__(self)

        # Get a copy of the dict with the clusters.
        self.clusters = _ALL_CLUSTERS.copy()

        self.flows_db = FlowsDatabase()

    def emptyline(self):
        """Disable the repetition of the last command."""

    def onecmd(self, line):
        try:
            r = super(RunCommand, self).onecmd(line)
            if r and (self.can_exit() or raw_input('Exit anyway? [Y/n]:').lower().startswith("y")):
                 return True
            return False

        except:
            print(straceback())
            return False

    def can_exit(self):
        return False

    def do_exit(self, s):
        """Exit the interpreter. You can also use the Ctrl-D shortcut."""
        return True
                                                                                   
    do_EOF = do_exit

    def do_shell(self, s):
        """Execute shell command on the local host."""
        os.system(s)
                                        
    def complete_hostnames(self, text, line, begidx, endidx):
        """
        The complete_xxx method takes four arguments:

        Args:
            text:
                the string we are matching against, all returned matches must begin with it
                line is is the current input line
            begidx: 
                the beginning index in the line of the text being matched
            endidx: 
                the end index in the line of the text being matched

        It should return a list (possibly empty) of strings representing the possible completions. 
        The arguments begidx and endidx are useful when completion depends on the position of the argument.
        """
        return [i for i in self.clusters.keys() if i.startswith(text)]

    def _select_clusters(self, s):
        """
        Returns the list of clusters corresponding to the input string
        If s is empy, the full list of clusters is returned.
        """
        if not s:
            return self.clusters.values()
        else:
            hosts = s.split()
            return filter(None, (self.clusters.get(h, None) for h in hosts))

    def do_disable_hosts(self, s):
        """
        Disable the specified list of hosts. 
        Syntax: disable_hosts host1 host2 ...
        """
        hosts = s.split()
        if not hosts:
            print(self.do_disable_hosts.__doc__)

        for h in hosts:
            self.clusters.pop(h, None)

    complete_disable_hosts = complete_hostnames

    def do_reenable_hosts(self, s):
        """
        Renable the specified list of hosts.
        Syntax: reenable_hosts host1 host2 ...
        """
        for h in hosts:
            self.clusters[h] = _ALL_CLUSTERS[h]
                                                
    def complete_reenable_hosts(self, text, line, begidx, endidx):
        return [h for h in _ALL_CLUSTERS.keys() if h not in self.clusters]

    def do_show_clusters(self, s):
        """Print the list of clusters."""
        for c in self.clusters.values():
            print(c)

    def do_ping(self, s):
        """
        Ping the hosts, e.g. ping host1 host2
        """
        for cluster in self._select_clusters(s):
            print(cluster.ping())

    complete_ping = complete_hostnames

    def do_user_jobs(self, s):
        """
        Report the jobs of the users. e.g user_jobs host1 host2
        """
        for cluster in self._select_clusters(s):
            print(cluster.get_user_jobs())

    complete_user_jobs = complete_hostnames

    def do_qinfo(self, s):
        """Report info on the queue manager."""
        for cluster in self._select_clusters(s):
            print(cluster.get_qinfo())

    complete_qinfo = complete_hostnames

    def do_run(self, s):
        """
        Execute a shell command on all hosts in the list.
        E.g. run uname -a host1
        """
        # split string into command and list of hostnames.
        tokens = s.split()
        for i, tok in enumerate(tokens):
            if tok in self.clusters.keys():
                command = " ".join(tokens[:i])
                hostnames = " ".join(tokens[i:])
                break
        else:
            command, hostnames = s, None

        #print("command", command, "hosts", hostnames)
        for cluster in self._select_clusters(hostnames):
            result = cluster.ssh_exec(command)
            print(cluster.prefix_str(result.out))

    complete_run = complete_hostnames

    def do_abicheck(self, s):
        """Test if the abinit environment is properly setup on the remote hosts."""
        print("WARNING: still under development")
        cmd = "abicheck.py"

        for cluster in self._select_clusters(s):
            print("Testing abipy environment of %s" % cluster)
            shell = cluster.invoke_shell()

            try:
                result = shell.myexec(cmd)

                print("got result", result)
                if result.retcode:
                    msg = "%s returned exited with status %d" % (cmd, result.retcode)
                    cluster.prefix_str(msg)

            except socket.timeout:
                print("socket timeout")
                
            #shell.exit()

    complete_abicheck = complete_hostnames

    def do_flow_start(self, s):
        """
        Start the flow on the remote cluster.
        Syntax: flow_start script.py hostname.
        """
        tokens = s.split()

        if len(tokens) != 2:
            print(self.do_flow_start.__doc__)
            return

        # Parse input line.
        script, hostname = os.path.basename(tokens[0]), tokens[1]
        cluster = self.clusters[hostname]

        # Build absolute paths on the remote host.
        dir_basename = os.path.basename(script).replace(".py", "")
        flow_absdir = cluster.path_inworkdir(dir_basename)
        remotepath = cluster.path_inworkdir(script)

        print("Uploading %s to %s:%s" % (script, hostname, remotepath))

        #if flow_absdir in self.flows_db:
        #    raise RuntimeError("remotepath %s is already in the database" % remotepath)

        # Upload the script and make it executable.
        cluster.sftp.put(localpath=script, remotepath=remotepath, confirm=True)
        cluster.sftp.chmod(remotepath, mode=0700)
        cluster.sftp.close()

        # Start a shell on the remote host and build the flow.
        shell = cluster.invoke_shell()
        result = shell.myexec(remotepath)

        if result.retcode:
            print("script returned %s" % result.retcode)

        # Run the scheduler with nohup.
        # This is the most delicate part: not so sure that nohup will work everywhere.
        #sched_cmd = "nohup abirun.py %s scheduler > /dev/null 2>&1 &" % flow_absdir
        #print(sched_cmd)
        #result = shell.myexec(sched_cmd)
        #print(result)

        #retcode = shell.exit()
        #print("retcode", retcode)

        # Add the flow to the local database.
        #self.flows_db.add_flow(flow_absdir, cluster)

    def complete_flow_start(self, text, line, begidx, endidx):
        tokens = line.split()

        if len(tokens) == 2:
            files = [f for f in os.listdir(os.getcwd()) if f.endswith(".py")]
            return [i for i in files if i.startswith(text)]

        elif len(tokens) == 3:
            return self.complete_hostnames(text, line, begidx, endidx)

        return []

    def do_show_flows(self, s):
        """Show the flows saved in the database."""
        print(self.flows_db)

    def do_flow_status(self, s):
        """
        Inspect the status of the flow on the remote cluster.
        Syntax: flow_status flow_workdir(s)
        """
        tokens = s.split()
        #print(tokens)

        if not tokens:
            print(self.do_flow_status.__doc__)
            return

        for workdir in tokens:
            params = self.flows_db.get(workdir, None)
            if params is None: continue

            cluster = self.clusters[params.get("hostname")]

            shell = cluster.invoke_shell()
            cmd= "abirun.py %s status" % workdir
            print(cmd)
            result = shell.myexec(cmd)
            print(result)
            #retcode = shell.exit()
            #print("retcode", retcode)

    def complete_flow_status(self, text, line, begidx, endidx):
        tokens = line.split()

        if len(tokens) == 1:
            return list(self.flows_db.keys())
                                                                              
        elif len(tokens) == 2:
            return [f for f in self.flows_db.keys() if f.startswith(text)]
                                                                              
        return []


class FlowsDatabase(collections.MutableMapping):
    JSON_FILE = "flowsdb.json"
    #flow_workdir --> {hostname, start_date, status}

    def __init__(self):
        self.filepath = os.path.join(os.getcwd(), self.JSON_FILE)

        if not os.path.exists(self.filepath):
            self.db = {}
        else:
            with open(self.filepath, "r") as fh:
                self.db = json.load(fh)

    def __str__(self):
        strio = StringIO.StringIO()
        pprint.pprint(self.db, stream=strio, indent=1)
        return strio.getvalue()

    # abc protocol.
    def __getitem__(self, key):
        return self.db[key]

    def __setitem__(self, key, value):
        self.db[key] = value

    def __delitem__(self, key):
        del self.db[key]

    def __iter__(self):
        return iter(self.db)

    def __len__(self):
        return len(self.db)
    # end abc protocol.

    def add_flow(self, flow_workdir, cluster):
        if flow_workdir in self.db:
            raise RuntimeError("Flow workdir %s is already in the database" % flow_workdir)

        d = {"hostname": cluster.hostname, "start_date": time.asctime(), "status": "init"}
        self.db[flow_workdir] = d

        self.json_dump()

    def remove_flow(self, flow_workdir):
        v = self.pop(flow_workdir, None)
        if v is not None:
            self.dump()

    def json_dump(self):
        with open(self.filepath, "w") as fh:
            json.dump(self.db, fh)


def main():
    retcode = 0

    cluster = SlurmCluster(hostname="manneback", username="gmatteo", password="e=mc^2")

    #cluster.ssh_exec("uname -a")
    #print(cluster.info)
    print(cluster.get_user_jobs())
    print(cluster.get_qinfo())

    cluster.sftp.put(localpath="run_si_ebands.py", remotepath="run_si_ebands.py", confirm=True)
    cluster.sftp.chmod("run_si_ebands.py", mode=0700)

    #shell = cluster.invoke_shell()

    # See http://stackoverflow.com/questions/2202228/how-do-i-launch-background-jobs-w-paramiko
    #result = shell.myexec("./run_si_ebands.py\n")
    return retcode


if __name__ == "__main__":
    import sys
    retcode = 0
    RunCommand().cmdloop()
    #retcode = main()
    #sys.exit(retcode)
