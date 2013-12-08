from __future__ import print_function, division

import os
import subprocess
import wx

import abipy.gui.awx as awx
import wx.lib.agw.flatnotebook as fnb
import wx.lib.newevent

from collections import OrderedDict
from pymatgen.util.io_utils import which
from abipy.gui.editor import  SimpleTextViewer, MyEditorFrame
from abipy.data.runs.abife import FlowsDatabase

# Command event used to signal that the Flows database 
# is changed and we should refresh the GUI.
DbChangedEvent, EVT_DB_CHANGED = wx.lib.newevent.NewCommandEvent()


def signal_db_changed(target):
    """Create and post the event to trigger the refresh of the GUI."""
    event = DbChangedEvent(id=-1)
    wx.PostEvent(target, event)


def straceback(color="red"):
    """
    Returns a string with the traceback.

    Use ANSII color formatting for output in terminal if color is not None.
    """
    import traceback
    s = traceback.format_exc()
    if color is not None:
        try:
            from termcolor import colored
            return colored(s, color)
        except ImportError:
            return s
    else:
        return s

def busy(func):
    """
    Decorator for functions that might become busy.
    A message window is shown while we are executing code on the remote host.

    .. note:
    
        The first argument of func is the wx parent Window.
    """
    from functools import wraps

    @wraps(func)
    def wrapper(*args, **kwargs):
        wait = wx.BusyInfo("Contacting host...")

        try:
            result = func(*args, **kwargs)
        except Exception as exc:
            result = exc
        finally:
            del wait
            if isinstance(result, Exception):
                #raise result
                awx.showErrorMessage(args[0], message=straceback(color=None))
            else:
                return result

    return wrapper


class FlowsDbViewerFrame(awx.Frame):
    """
    This frame shows the active flows and allows the user
    to open and interact with a particular `AbinitFlow`.
    """
    VERSION = "0.1"

    def __init__(self, parent, **kwargs):
        """
        Args:
            parent:
                Parent window
        """
        super(FlowsDbViewerFrame, self).__init__(parent, -1, **kwargs)

        self.statusbar = self.CreateStatusBar()

        menuBar = wx.MenuBar()

        self.ID_RUN_SCRIPT = wx.NewId()
        self.ID_CHECK_STATUS = wx.NewId()
        self.ID_USER_JOBS = wx.NewId()
        self.ID_ALL_JOBS = wx.NewId()
        self.ID_TERMINAL = wx.NewId()
        self.ID_SHOW_ABINIT_INFO = wx.NewId()
        self.ID_SHOW_ABIPY_ENV = wx.NewId()
        self.ID_RUN_COMMAND_ALL = wx.NewId()
        self.ID_RUN_COMMAND_SINGLE = wx.NewId()
        self.ID_PING = wx.NewId()
        self.ID_EDIT_TASKMANAGER = wx.NewId()
        self.ID_EDIT_SCHEDULER = wx.NewId()
        self.ID_EDIT_BASHRC = wx.NewId()

        tools_menu = wx.Menu()
        tools_menu.Append(self.ID_RUN_COMMAND_SINGLE, "Run command (selected host)", help="Execute shell command on the selected host")
        tools_menu.Append(self.ID_RUN_COMMAND_ALL, "Run command (all hosts)", help="Execute shell command on the hosts")
        tools_menu.Append(self.ID_PING, "Ping hosts", help="Ping hosts")
        menuBar.Append(tools_menu, "Tools")

        edit_menu = wx.Menu()
        edit_menu.Append(self.ID_EDIT_TASKMANAGER, 'taskmanager.yml')
        edit_menu.Append(self.ID_EDIT_SCHEDULER, 'scheduler.yml')
        edit_menu.Append(self.ID_EDIT_BASHRC, 'bashrc')
        tools_menu.AppendMenu(-1, 'E&dit', edit_menu)

        queue_menu = wx.Menu()
        queue_menu.Append(self.ID_USER_JOBS, "user jobs", help="Show the list of queued jobs belonging to the user.")
        queue_menu.Append(self.ID_ALL_JOBS, "all jobs", help="Show all the jobs in the queue.")
        menuBar.Append(queue_menu, 'Q&ueue')

        help_menu = wx.Menu()
        help_menu.Append(wx.self.ID_ABOUT, "About " + self.codename, help="Info on the application")
        menuBar.Append(help_menu, "Help")

        self.SetMenuBar(menuBar)

        # Create toolbar.
        self.toolbar = toolbar = self.CreateToolBar()

        def bitmap(path):
            return wx.Bitmap(awx.path_img(path))

        toolbar.AddSimpleTool(self.ID_RUN_SCRIPT, bitmap("script.png"), "Upload and execute the script on the remote host.")
        toolbar.AddSimpleTool(self.ID_CHECK_STATUS, bitmap("script.png"), "Check the status of the flows running on the remote host.")

        toolbar.AddSeparator()

        toolbar.AddSimpleTool(self.ID_TERMINAL, bitmap("script.png"), "Open terminal and connect to the remote host.")
        toolbar.AddSimpleTool(self.ID_SHOW_ABINIT_INFO, bitmap("script.png"), "Show the ABINIT version and the build info used on the remote host")
        toolbar.AddSimpleTool(self.ID_SHOW_ABIPY_ENV, bitmap("script.png"), "Show the abipy enviroment available on the remote host.")

        self.toolbar.Realize()
        self.Centre()

        # Associate menu/toolbar items with their handlers.
        menu_handlers = [
            (wx.ID_ABOUT, self.OnAboutBox),
            (self.ID_RUN_SCRIPT, self.OnRunScript),
            (self.ID_CHECK_STATUS, self.OnCheckStatus),
            (self.ID_USER_JOBS, self.OnUserJobs),
            (self.ID_ALL_JOBS, self.OnAllJobs),
            (self.ID_TERMINAL, self.OnTerminal),
            (self.ID_SHOW_ABINIT_INFO, self.OnShowAbinitInfo),
            (self.ID_SHOW_ABIPY_ENV, self.OnShowAbipyEnv),
            #
            (self.ID_RUN_COMMAND_ALL, self.OnRunCommandAll),
            (self.ID_RUN_COMMAND_SINGLE, self.OnRunCommandSingle),
            (self.ID_PING, self.OnPing),
            (self.ID_EDIT_TASKMANAGER, self.OnEditTaskManager),
            (self.ID_EDIT_SCHEDULER, self.OnEditScheduler),
            (self.ID_EDIT_BASHRC, self.OnEditBashrc),
        ]

        for combo in menu_handlers:
            mid, handler = combo[:2]
            self.Bind(wx.EVT_MENU, handler, id=mid)

        # Read the database with the list of flows.
        self.flows_db = FlowsDatabase.from_user_config()

        panel = wx.Panel(self, -1)
        main_sizer = wx.BoxSizer(wx.VERTICAL)

        self.notebook = FlowsDbNotebook(panel, self.flows_db)
        main_sizer.Add(self.notebook, 1, wx.EXPAND, 5)

        panel.SetSizerAndFit(main_sizer)

        # Intercept the command event associated to the database modification.
        self.Bind(EVT_DB_CHANGED, self.ReshowFlowsDb)

    @property
    def codename(self):
        """String with the code name """
        return "Flows Database"

    def OnAboutBox(self, event):
        """"Info on the application."""
        awx.makeAboutBox(codename=self.codename, version=self.VERSION,
                         description="", developers="M. Giantomassi")

    def GetSelectedCluster(self):
        """
        Return the selected cluster namely that the 
        cluster associated to the active tab. None if list is empty.
        """
        return self.notebook.GetSelectedCluster()

    def ReshowFlowsDb(self, event):
        """Refresh the GUI, called when we have modified the database."""
        #print("in refresh with event %s" % event)
        self.notebook.ReshowFlowsDb(event)

    def OnRunScript(self, event):
        """Browse all the output files produced by the selected `Workflow`."""
        cluster = self.GetSelectedCluster() 
        if cluster is None: return

        dialog = wx.FileDialog(self, wildcard="*.py")
        if dialog.ShowModal() == wx.ID_CANCEL: return 

        script = dialog.GetPath()
        #print("cluster", cluster, "script", script)

        try:
            wait = wx.BusyInfo("Contacting host...")
            self.flows_db.start_flow(script, cluster.hostname)
            self.ReshowFlowsDb(event)
            del wait

        except:
            del wait
            awx.showErrorMessage(self)

    @busy
    def OnCheckStatus(self, event):
        """Check the status of the flows."""
        cluster = self.GetSelectedCluster() 
        if cluster is None: return
                                                                               
        results, changed = self.flows_db.check_status(hostnames=cluster.hostname)

        for res in results: print(res)

        if changed:
            self.ReshowFlowsDb(event)

    @busy
    def OnUserJobs(self, event):
        """Open a new frame with the list of jobs submitted by the user."""
        cluster = self.GetSelectedCluster() 
        if cluster is None: return

        s = cluster.get_user_jobs()
        SimpleTextViewer(self, text=s, title=cluster.hostname).Show()

    @busy
    def OnAllJobs(self, event):
        """Open a new frame with the full list of jobs in the queueu."""
        cluster = self.GetSelectedCluster() 
        if cluster is None: return

        s = cluster.get_all_jobs()
        SimpleTextViewer(self, text=s, title=cluster.hostname).Show()

    def OnRunCommandAll(self, event):
        """Execute a shell command on all hosts."""
        dialog = wx.TextEntryDialog(self, message="Enter the command to execute on all hosts.") 
        if dialog.ShowModal() == wx.ID_CANCEL: return 

        self._runcmd_mode(command=dialog.GetValue(), mode="all")
        dialog.Destroy()

    def OnRunCommandSingle(self, event):
        """Execute a shell command on the selected host."""
        dialog = wx.TextEntryDialog(self, message="Enter the command to execute on the selected host.") 
        if dialog.ShowModal() == wx.ID_CANCEL: return 
                                                                               
        self._runcmd_mode(command=dialog.GetValue(), mode="single")
        dialog.Destroy()

    @busy
    def _runcmd_mode(self, command, mode="all"):
        """Execute command on the cluster(s) depending on mode."""
        if mode == "single":
            selected_cluster = self.GetSelectedCluster() 

        text = []
        for cluster in self.flows_db.clusters.values():
            # FIXME
            if cluster.hostname == "manneback": continue

            if mode == "single" and cluster != selected_cluster: continue
            result = cluster.ssh_exec(command)
            text.append(cluster.prefix_str(result.out))

        SimpleTextViewer(self, text="\n".join(text), title=command).Show()

    def OnPing(self, event):
        """Ping the cluster and show the results."""
        text = []
        for cluster in self.flows_db.clusters.values():
            text.append(cluster.ping())

        SimpleTextViewer(self, text="\n".join(text), title="Ping Statistics").Show()

    def OnTerminal(self, event):
        """Open a new terminal and ssh to the remote host."""
        cluster = self.GetSelectedCluster() 
        if cluster is None: return

        def open_terminal():
            """Try to figure which terminal is available."""
            retcode = 1

            if retcode and which("gnome-terminal") is not None:
                retcode = os.system("gnome-terminal -e 'ssh %s -X'" % cluster.hostname)

            if retcode and which("konsole") is not None:
                retcode = os.system("konsole -e 'ssh %s -X'" % cluster.hostname)

            # This does not work
            #retcode = os.system('open -a Terminal -n --args %s -X" % cluster.hostname")
                
            if retcode: # Fallback to xterm.
                retcode = os.system("xterm -e ssh %s -X" % cluster.hostname)

            return retcode

        # FIXME: Fix possible problem when the user tries to close the GUI 
        # with active terminals (maintain a list of treads and prevent the user from closing the GUI if threads?)
        try:
            thread = awx.WorkerThread(self, target=open_terminal)
            thread.start()
        except:
            awx.showErrorMessage(self)

    @busy
    def OnShowAbinitInfo(self, event):
        """Show info on the Abinit binary used on the cluster."""
        cluster = self.GetSelectedCluster() 
        if cluster is None: return

        d = cluster.get_abinit_info()
        SimpleTextViewer(self, text=str(d), title=cluster.hostname).Show()

    @busy
    def OnShowAbipyEnv(self, event):
        """Show info on the Abinit environment used on the cluster."""
        cluster = self.GetSelectedCluster() 
        if cluster is None: return

        d = dict(abinit=cluster.which("abinit"),
                 abicheck=cluster.abicheck(),
                 yaml_confs=cluster.yaml_configurations())

        SimpleTextViewer(self, text=str(d), title=cluster.hostname).Show()

    def OnEditScheduler(self, event):
        cluster = self.GetSelectedCluster() 
        self.remote_edit(cluster, remotepath=os.path.join(cluster.home, ".abinit/abipy/scheduler.yml"))

    def OnEditTaskManager(self, event):
        cluster = self.GetSelectedCluster() 
        self.remote_edit(cluster, remotepath=os.path.join(cluster.home, ".abinit/abipy/taskmanager.yml"))

    def OnEditBashrc(self, event):
        cluster = self.GetSelectedCluster() 
        self.remote_edit(cluster, remotepath=os.path.join(cluster.home, ".bashrc"))

    def remote_edit(self, cluster, remotepath):

        class LocalEditorFrame(MyEditorFrame):
            """
            Editor used to edit a file locally. The new file is 
            sfpt-transferred when the user closes the frame.
            """
            def OnClose(self, event):
                dialog = wx.MessageDialog(self, "Do you want to copy the new file to the remote host?", 
                                          style=wx.OK | wx.CANCEL | wx.CENTRE)
                                                                                                                        
                if dialog.ShowModal() == wx.ID_OK:
                    text = self.editor.getText()
                    print("will transfer text", text)

                    sftp = cluster.sftp
                    sftp.put(localpath=self.localpath, remotepath=self.remotepath, confirm=True)
                    #sftp.chmod(self.remotepath, mode=0700)

                super(MyEditorFrame, self).OnClose(event)

        # Get the file via sftp and dump it to a temporary file.
        import tempfile
        _, localpath = tempfile.mkstemp(text=True)

        with open(localpath, "w") as fo:
            cluster.sftp.getfo(remotepath, fo)

        # Open a new frame for editing and save useful variables 
        # we'll need when the user closes the editor.
        frame = LocalEditorFrame(self, localpath, title="remote: " + remotepath)
        frame.cluster, frame.localpath, frame.remotepath = cluster, localpath, remotepath
        frame.Show()


class FlowsDbNotebook(fnb.FlatNotebook):
    """
    Notebook class
    """
    def __init__(self, parent, flows_db):
        super(FlowsDbNotebook, self).__init__(parent, id=-1, 
              style=fnb.FNB_NO_X_BUTTON | fnb.FNB_NAV_BUTTONS_WHEN_NEEDED)

        self.flows_db = flows_db
        self.clusters_byname = flows_db.clusters
        self.clusters = list(self.clusters_byname.values())

        self._make_pages()

    def _make_pages(self):
        """Build the pages of the notebook."""
        for hostname, flows in self.flows_db.items():
            cluster = self.clusters_byname[hostname]
            tab = ClusterPanel(self, cluster, self.flows_db)
            self.AddPage(tab, text=hostname)

    def ReshowFlowsDb(self, event):
        """Refresh the notebook, called when we have modified the database."""
        # Save the active tab so that we can set it afterwards.
        old_selection = self.GetSelection()

        self.DeleteAllPages()
        self._make_pages()

        self.SetSelection(old_selection)

    def GetSelectedCluster(self):
        """
        Return the selected cluster namely that the cluster associated to the active tab. None is list is empy.
        """
        if self.GetPageCount() != len(self.flows_db):
            return awx.showErrorMessage(self, message="Bad user has removed pages from the notebook!")

        idx = self.GetSelection()
        if idx == -1: return None

        # Use the text of the tab to retrieve the cluster.
        hostname = self.GetPageText(idx)
        return self.clusters_byname[hostname]


class ClusterPanel(wx.Panel):
    """
    Notebook tab for a single cluster.
    """
    def __init__(self, parent, cluster, flows_db, **kwargs):
        wx.Panel.__init__(self, parent=parent, id=-1, **kwargs)

        self.cluster, self.flows_db = cluster, flows_db
        self.flows = self.flows_db[cluster.hostname]

        self.BuildUi()

    def BuildUi(self):

        main_sizer = wx.BoxSizer(wx.VERTICAL)

        # List Control with info on the flows.
        self.flows_listctrl = FlowsListCtrl(self, self.flows, self.cluster, self.flows_db)
        main_sizer.Add(self.flows_listctrl, 1, wx.EXPAND, 5)

        self.SetSizerAndFit(main_sizer)


class FlowsListCtrl(wx.ListCtrl):
    """
    ListCtrl with the list of flows being executed on a cluster.
    """
    def __init__(self, parent, flows, cluster, flows_db, **kwargs):
        """
        Args:
            parent:
                Parent window.
            flows:
                List of flows info
            cluster:
                The cluster where the flows are executed.
            flows_db:
                The database of flows.
        """
        super(FlowsListCtrl, self).__init__(parent, id=-1, style=wx.LC_REPORT | wx.BORDER_SUNKEN, **kwargs)

        self.flows, self.cluster, self.flows_db = flows, cluster, flows_db

        columns = ["Workdir", "Status", "Start Date"]

        for (index, col) in enumerate(columns):
            self.InsertColumn(index, col)

        # Used to store the Max width in pixels for the data in the column.
        column_widths = [awx.get_width_height(self, s)[0] for s in columns]

        for flow in flows:
            entry = map(str, [flow.workdir, flow.status, flow.start_date])
            self.Append(entry)

            w = [awx.get_width_height(self, s)[0] for s in entry]
            column_widths = map(max, zip(w, column_widths))

        # Set the width in pixel for each column.
        for (index, col) in enumerate(columns):
            self.SetColumnWidth(index, column_widths[index])

        self.Bind(wx.EVT_LIST_ITEM_RIGHT_CLICK, self.OnRightClick)
        self.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.OnItemActivated)

    def OnRightClick(self, event):
        currentItem = event.m_itemIndex
        if currentItem == -1: return

        # Open the popup menu then destroy it to avoid mem leak.
        flow = self.flows[currentItem]
        menu = FlowPopupMenu(self, self.cluster, flow, self.flows_db)
        self.PopupMenu(menu, event.GetPoint())
        menu.Destroy()

    def OnItemActivated(self, event):
        idx = int(self.GetFirstSelected())
        flow = self.flows[idx]
        print(flow)


# Callbacks for the PopupMenu.
@busy
def flow_show_status(parent, cluster, flow, flows_db):
    """Show the status of the flow."""
    results, changed = flows_db.check_status(cluster.hostname, flow.workdir)
    SimpleTextViewer(parent, text=str(results[0]), title=cluster.hostname).Show()

    if changed: # Generate the event to refresh the GUI.
        signal_db_changed(parent)


def flow_cancel(parent, cluster, flow, flows_db):
    """Cancel the flow i.e. remove all the jobs that are in the queue."""
    dialog = wx.MessageDialog(parent, "Do you want to cancel all the queued jobs belonging to this flow?", 
                              style=wx.OK | wx.CANCEL | wx.CENTRE)
    if dialog.ShowModal() == wx.ID_CANCEL: return
    flows_db.cancel_flow(cluster.hostname, flow.workdir)
    # Generate event to refresh the GUI.
    signal_db_changed(parent) 


def flow_remove(parent, cluster, flow, flows_db):
    """Remove the working directory of the flow."""
    dialog = wx.MessageDialog(parent, "Do you want to remove the workdir of the flow on the remote host?", 
                              style=wx.OK | wx.CANCEL | wx.CENTRE)
    if dialog.ShowModal() == wx.ID_CANCEL: return
    cluster.rmdir(flow.workdir)
    flows_db.remove_flow(flow)
    # Generate event to refresh the GUI.
    signal_db_changed(parent) 


@busy
def flow_sched_log(parent, cluster, flow, flows_db):
    """Show the content of the scheduler log file."""
    sched_log = os.path.join(flow.workdir, "sched.log")
    s = cluster.read_file(sched_log)
    SimpleTextViewer(parent, text=s, title=sched_log).Show()


class FlowPopupMenu(wx.Menu):
    """
    A `FlowPopupMenu` has a list of callback functions indexed by the menu title. 
    The signature of the callback function is func(parent, cluster, flow) where parent is 
    the wx Window that will become the parent of the new frame created by the callback.
    and flow is the dictionary with info on the `AbinitFlow`.
    """
    MENU_TITLES = OrderedDict([
        ("show_status", flow_show_status),
        ("cancel", flow_cancel),
        ("remove", flow_remove),
        ("sched_log", flow_sched_log),
    ])

    def __init__(self, parent, cluster, flow, flows_db):
        super(FlowPopupMenu, self).__init__()
        self.parent, self.cluster, self.flow = parent, cluster, flow
        self.flows_db = flows_db

        self._make_menu()

    def _make_menu(self):
        """Build the menu"""
        self.menu_title_by_id = OrderedDict()

        for title in self.MENU_TITLES:
            self.menu_title_by_id[wx.NewId()] = title

        for (id, title) in self.menu_title_by_id.items():
            # Register menu handlers with EVT_MENU, on the menu.
            self.Append(id, title)
            wx.EVT_MENU(self, id, self.OnMenuSelection)

    def _get_callback(self, title):
        return self.MENU_TITLES[title]

    def OnMenuSelection(self, event):
        title = self.menu_title_by_id[event.GetId()]
        callback = self._get_callback(title)

        #print("Calling callback %s with cluster %s and flow %s" % (callback, self.cluster, self.flow))
        try:
            callback(self.parent, self.cluster, self.flow, self.flows_db)
        except:
            awx.showErrorMessage(parent=self.parent)


def wxapp_flowsdb_viewer():
    """Standalone application for `FlowsDbViewerFrame"""
    app = awx.App()
    FlowsDbViewerFrame(None).Show()
    return app

#class JobsPanel(wx.Panel):
#    def __init__(self, parent, cluster, **kwargs):
#        super(JobsPanel, self).__init__(parent, -1, **kwargs)
#        self.cluster = cluster
#
#        main_sizer = wx.BoxSizer(wx.VERTICAL)
#
#        #user_jobs_button = wx.Button(self, -1, label='User Jobs')
#        #all_jobs_button = wx.Button(self, -1, label='All Jobs')
#        #self.Bind(wx.EVT_BUTTON, self.ShowAllJobs, all_jobs_button)
#        #self.Bind(wx.EVT_BUTTON, self.ShowUserJobs, user_jobs_button)
#
#        #hbox = wx.BoxSizer(wx.HORIZONTAL)
#        #hbox = wx.BoxSizer(wx.HORIZONTAL)
#        #hbox.Add(user_jobs_button)
#        #hbox.Add(all_jobs_button, flag=wx.LEFT, border=5)
#        #main_sizer.Add(hbox, proportion=0, flag=wx.ALIGN_CENTER | wx.TOP | wx.BOTTOM, border=10)
#
#        #self.text_ctrl = wx.TextCtrl(self, -1, value="", style=wx.TE_MULTILINE|wx.TE_LEFT|wx.TE_READONLY)
#        #main_sizer.Add(self.text_ctrl, 1, wx.ALIGN_CENTER_HORIZONTAL, 5)
#
#        self.SetSizerAndFit(main_sizer)

if __name__ == "__main__":
    wxapp_flowsdb_viewer().MainLoop()
