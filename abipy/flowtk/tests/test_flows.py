# coding: utf-8
import os
import tempfile
import shutil
import abipy.data as abidata

from monty.functools import lazy_property
from pymatgen.core.lattice import Lattice
from abipy.core.structure import Structure
from abipy.flowtk.launcher import BatchLauncher
from abipy.flowtk.flows import *
from abipy.flowtk.works import *
from abipy.flowtk.tasks import *
from abipy.core.testing import AbipyTest
from abipy import abilab
from abipy import flowtk


class FakeAbinitInput(object):
    """Emulate an Abinit input."""
    @lazy_property
    def pseudos(self):
        return [abidata.pseudo("14si.pspnc")]

    @lazy_property
    def structure(self):
        coords = []
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        lattice = Lattice([[3.8401979337, 0.00, 0.00],
                          [1.9200989668, 3.3257101909, 0.00],
                          [0.00, -2.2171384943, 3.1355090603]])
        return Structure(lattice, ["Si", "Si"], coords)

    def get(self, key, default=None):
        """The real AbinitInput is a dict-like object."""
        if default is not None: return default
        return key


class FlowUnitTest(AbipyTest):
    """Provides helper function for testing Abinit flows."""
    MANAGER = """\
policy:
    autoparal: 1
qadapters:
    - &batch
      priority: 1
      queue:
        qtype: slurm
        qname: Oban
        qparams:
            mail_user: nobody@nowhere
      limits:
        timelimit: 0:20:00
        min_cores: 4
        max_cores: 12
        #condition: {"$eq": {omp_threads: 2}}
      hardware:
        num_nodes: 10
        sockets_per_node: 1
        cores_per_socket: 2
        mem_per_node: 4 Gb
      job:
        modules:
            - intel/compilerpro/13.0.1.117
            - fftw3/intel/3.3
        shell_env:
            PATH: /home/user/tmp_intel13/src/98_main/:/home/user//NAPS/intel13/bin:$PATH
            LD_LIBRARY_PATH: /home/user/NAPS/intel13/lib:$LD_LIBRARY_PATH
        mpi_runner: mpirun

# Connection to the MongoDb database (optional)
db_connector:
    database: abinit
    collection: test
    #host: 0.0.0.0
    #port: 8080
    #user: gmatteo
    #password: helloworld

batch_adapter: *batch
"""
    def setUp(self):
        """Initialization phase."""
        super().setUp()

        # Temporary directory for the flow.
        self.workdir = tempfile.mkdtemp()

        # Create the TaskManager.
        self.manager = TaskManager.from_string(self.MANAGER)

        # Fake input file
        self.fake_input = FakeAbinitInput()

    def tearDown(self):
        """Delete workdir"""
        shutil.rmtree(self.workdir)


class FlowTest(FlowUnitTest):

    def test_base(self):
        """Testing Flow..."""
        aequal = self.assertEqual
        flow = Flow(workdir=self.workdir, manager=self.manager)
        assert flow.isinstance(Flow)
        assert not flow.isinstance(None)
        assert not flow.has_scheduler

        # Build a work with a task
        work = flow.register_task(self.fake_input)
        assert work.is_work
        assert len(work.color_hex) == 7
        assert work.color_hex.startswith("#")
        task0_w0 = work[0]
        assert task0_w0.is_task
        str(task0_w0.status.colored)
        assert len(flow) == 1
        assert flow.num_tasks == 1
        assert flow.has_db

        #print(task0_w0.input_structure)
        str(task0_w0.make_input)

        # Task history
        assert len(task0_w0.history) == 0
        task0_w0.history.info("Hello %s", "world")
        assert len(task0_w0.history) == 1
        str(task0_w0.history)
        record = task0_w0.history.pop()
        print(record, repr(record))
        assert record.get_message(asctime=False) == "Hello world"
        assert len(task0_w0.history) == 0
        assert flow.select_tasks(nids=task0_w0.node_id)[0] == task0_w0
        assert flow.select_tasks(wslice=slice(0,1,1)) == [task0_w0]
        assert flow.select_tasks(task_class="DfptTask") == []
        assert flow.get_task_scfcycles() == []

        # Build a workflow containing two tasks depending on task0_w0
        work = Work()
        assert work.is_work
        work.register(self.fake_input)
        work.register(self.fake_input)
        assert len(work) == 2

        flow.register_work(work, deps={task0_w0: "WFK"})
        assert flow.is_flow
        assert len(flow) == 2

        # Add another work without dependencies.
        task0_w2 = flow.register_task(self.fake_input)[0]
        assert len(flow) == 3
        assert not flow.is_work

        # Allocate internal tables
        flow.allocate()

        # Check dependecies.
        #task0_w1 = flow[1][0]
        assert flow[1].depends_on(task0_w0)
        assert flow[1][0].depends_on(task0_w0)
        assert flow[1][0] in task0_w0.get_children()
        assert task0_w0 in flow[1][0].get_parents()
        assert flow[1][0].find_parent_with_ext("WFK") == task0_w0
        assert flow[1][0].find_parent_with_ext("FOOBAR") is None
        assert not flow[2][0].depends_on(task0_w0)
        assert not flow[2][0] in task0_w0.get_children()
        assert not task0_w0 in flow[2][0].get_parents()
        assert flow[1].pos == 1
        assert flow[1][0].pos == (1, 0)
        assert flow[2][0].pos == (2, 0)

        assert not flow.all_ok
        assert flow.num_tasks == 4
        assert flow.ncores_used == 0

        # API for iterations
        aequal(len(list(flow.iflat_tasks(status="Initialized"))), sum(len(work) for work in flow))
        aequal(list(flow.iflat_tasks(nids=task0_w0.node_id)), [task0_w0])

        aequal([task0_w0], flow.tasks_from_nids(task0_w0.node_id))
        aequal([(0, 0)], flow.wti_from_nids(task0_w0.node_id))
        aequal([task0_w2], flow.tasks_from_nids([task0_w2.node_id]))
        aequal([(2, 0)], flow.wti_from_nids([task0_w2.node_id]))

        # Check for deadlocks
        flow.check_dependencies()

        # Save the flow in pickle format.
        flow.build_and_pickle_dump()

        data = {"foo": 1, "structure": flow[0][0].input.structure}
        flow[0][0].write_json_in_outdir("foo.json", data)

        # Find the pickle file in workdir and recreate the flow.
        same_flow = Flow.pickle_load(self.workdir)
        assert same_flow == flow

        # to/from string
        # FIXME This does not work with py3k
        #s = flow.pickle_dumps(protocol=0)
        #same_flow = Flow.pickle_loads(s)
        #aequal(same_flow, flow)

        self.assertMSONable(flow)

        flow.show_info()
        flow.show_summary()
        flow.show_inputs()
        flow.show_inputs(varnames="znucl")

        df_vars = flow.get_vars_dataframe("ecut", "acell")
        assert "ecut" in df_vars

        # Test show_status
        flow.show_status()
        flow.show_tricky_tasks()
        flow.show_event_handlers()

        df = flow.compare_abivars(varnames=["ecut", "natom"], printout=True, with_colors=True)
        assert "ecut" in df

        dfs = flow.compare_structures(with_spglib=False, verbose=2, printout=True, with_colors=True)
        assert "alpha" in dfs.lattice

        dfs, ebands_plotter = flow.compare_ebands(verbose=0)
        assert not dfs and not ebands_plotter

        dfs, hist_plotter = flow.compare_hist(with_spglib=False, verbose=2, printout=True, with_colors=True)
        assert not dfs and not hist_plotter

        if self.has_networkx():
            assert flow.plot_networkx(mode="network", with_edge_labels=False, arrows=False,
                      node_size="num_cores", node_label="name_class", layout_type="spring", show=False)
            assert flow.plot_networkx(mode="status", with_edge_labels=True, arrows=True,
                      node_size="num_cores", node_label="name_class", layout_type="spring", show=False)

        if self.has_python_graphviz():
            assert flow.get_graphviz(engine="automatic", graph_attr=None, node_attr=None, edge_attr=None)
            assert flow.graphviz_imshow(ax=None, figsize=None, dpi=300, fmt="png", show=False)

        if self.has_panel():
            assert hasattr(flow.get_panel(), "show")

    def test_workdir(self):
        """Testing if one can use workdir=None in flow.__init__ and then flow.allocate(workdir)."""
        flow = Flow(workdir=None, manager=self.manager)
        flow.register_task(self.fake_input)
        #flow.register_work(work)
        work = Work()
        work.register_scf_task(self.fake_input)
        flow.register_work(work)

        # If flow.workdir is None, we should used flow.allocate(workdir)
        with self.assertRaises(RuntimeError): flow.allocate()

        tmpdir = tempfile.mkdtemp()
        flow.allocate(workdir=tmpdir)

        str(flow)
        assert len(flow) == 2
        flow.build()

        for i, work in enumerate(flow):
            assert work.workdir == os.path.join(tmpdir, "w%d" % i)
            for t, task in enumerate(work):
                assert task.workdir == os.path.join(work.workdir, "t%d" % t)

    def test_nscf_flow_with_append(self):
        """Test creation of NSCF tasks from flow with append = True"""

        def make_scf_nscf_inputs():
            """Build ands return the input files for the GS-SCF and the GS-NSCF tasks."""

            multi = abilab.MultiDataset(structure=abidata.cif_file("si.cif"),
                                        pseudos=abidata.pseudos("14si.pspnc"), ndtset=2)

            # Set global variables (dataset1 and dataset2)
            multi.set_vars(ecut=6, nband=8)

            # Dataset 1 (GS-SCF run)
            multi[0].set_kmesh(ngkpt=[8, 8, 8], shiftk=[0, 0, 0])
            multi[0].set_vars(tolvrs=1e-6)

            # Dataset 2 (GS-NSCF run on a k-path)
            kptbounds = [
                [0.5, 0.0, 0.0], # L point
                [0.0, 0.0, 0.0], # Gamma point
                [0.0, 0.5, 0.5], # X point
            ]

            multi[1].set_kpath(ndivsm=6, kptbounds=kptbounds)
            multi[1].set_vars(tolwfr=1e-12)

            # Return two input files for the GS and the NSCF run
            scf_input, nscf_input = multi.split_datasets()
            return scf_input, nscf_input

        scf_input, nscf_input = make_scf_nscf_inputs()
        hello_flow = flowtk.Flow(workdir=self.mkdtemp())
        hello_flow.register_scf_task(scf_input, append=True)
        hello_flow.register_nscf_task(nscf_input, deps={hello_flow[0][0]: "DEN"}, append=True)
        #flow[0].get_graphviz_dirtree()
        #abilab.print_doc(flowtk.PhononWork)

        hello_flow = flowtk.Flow(workdir=self.mkdtemp())
        hello_flow.register_scf_task(scf_input, append=True)
        assert len(hello_flow) == 1
        hello_flow.register_nscf_task(nscf_input, deps={hello_flow[0][0]: "DEN"}, append=False)
        assert len(hello_flow) == 2
        empty_work = hello_flow.new_work()
        assert len(empty_work) == 0
        assert len(hello_flow) == 3



class TestFlowInSpectatorMode(FlowUnitTest):

    def test_spectator(self):
        flow = Flow(workdir=self.workdir, manager=self.manager)

        work0 = Work()
        gs_task = work0.register_scf_task(self.fake_input)
        assert gs_task.isinstance(ScfTask)
        assert gs_task.isinstance("ScfTask")
        task = work0.register_scf_task(self.fake_input)
        assert task.is_abinit_task
        assert not task.is_optic_task
        assert not task.is_anaddb_task

        work1 = Work()
        work1.register_scf_task(self.fake_input)

        flow.register_work(work0)
        flow.register_work(work1)

        flow.disconnect_signals()
        flow.disconnect_signals()

        flow.connect_signals()
        flow.connect_signals()

        for mode in [False, True]:
            flow.set_spectator_mode(mode=mode)
            assert flow.in_spectator_mode == mode
            for node in flow.iflat_nodes():
                assert node.in_spectator_mode == mode

        assert len(list(flow.iflat_nodes())) == 1 + len(flow.works) + sum(len(work) for work in flow)
        assert flow.node_from_nid(flow.node_id) == flow

        flow.set_spectator_mode(mode=False)
        flow.build_and_pickle_dump()

        # picke load always returns a flow in spectator mode.
        flow = Flow.pickle_load(flow.workdir)
        assert flow.in_spectator_mode

        #with self.assertRaises(flow.SpectatorError): flow.pickle_dump()
        #with self.assertRaises(flow.SpectatorError): flow.make_scheduler().start()

        work = flow[0]
        assert work.send_signal(work.S_OK) is None
        #with self.assertRaises(work.SpectatorError): work.on_ok()
        #with self.assertRaises(work.SpectatorError): work.on_all_ok()

        task = work[0]
        assert task.send_signal(task.S_OK) is None
        #with self.assertRaises(task.SpectatorError): task._on_done()
        #with self.assertRaises(task.SpectatorError): task.on_ok()
        #with self.assertRaises(task.SpectatorError): task._on_ok()


class TestBatchLauncher(FlowUnitTest):

    def test_batchlauncher(self):
        """Testing BatchLauncher methods."""
        # Create the TaskManager.
        manager = TaskManager.from_string(self.MANAGER)
        print("batch_adapter", manager.batch_adapter)
        assert manager.batch_adapter is not None

        def build_flow_with_name(name):
            """Build a flow with workdir None and the given name."""
            flow = Flow(workdir=None, manager=self.manager)
            flow.set_name(name)

            flow.register_task(self.fake_input)
            work = Work()
            work.register_scf_task(self.fake_input)
            flow.register_work(work)

            return flow

        tmpdir = tempfile.mkdtemp()
        batch = BatchLauncher(workdir=tmpdir, manager=manager)
        str(batch)

        flow0 = build_flow_with_name("flow0")
        flow1 = build_flow_with_name("flow1")
        flow2_same_name = build_flow_with_name("flow1")

        batch.add_flow(flow0)

        # Cannot add the same flow twice.
        with self.assertRaises(batch.Error):
            batch.add_flow(flow0)

        batch.add_flow(flow1)

        # Cannot add two flows with the same name.
        with self.assertRaises(batch.Error):
            batch.add_flow(flow2_same_name)

        batch.submit(dry_run=True)

        for i, flow in enumerate([flow0, flow1]):
            assert flow.workdir == os.path.join(batch.workdir, "flow%d" % i)

        batch.pickle_dump()
        batch_from_pickle = BatchLauncher.pickle_load(batch.workdir)
        assert all(f1 == f2 for f1, f2 in zip(batch.flows, batch_from_pickle.flows))
