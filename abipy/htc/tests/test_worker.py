import os

from monty.shutil import remove
from abipy.core.testing import AbipyTest
from abipy.flowtk import TaskManager
from abipy.htc.worker import AbipyWorker, WorkerState, WorkerClient, WorkerClients, ABIPY_DIRPATH


class TestWorker(AbipyTest):

    def test_worker_api(self):
        """Testing AbipyWorker."""
        worker_name = "__abipy_worker_unittest__"
        worker_dir = os.path.join(ABIPY_DIRPATH, "worker_" + worker_name)
        if os.path.isdir(worker_dir):
            remove(worker_dir)

        worker = AbipyWorker.new_with_name(worker_name=worker_name,
                                           scratch_dir="/tmp",
                                           scheduler_path=None,
                                           manager_path=None,
                                           mongo_connector=None,
                                           verbose=1)
        assert repr(worker)
        assert str(worker)
        assert isinstance(worker.manager, TaskManager)

        filepath = worker.write_state_file(status="dead", filepath=None)
        state = WorkerState.from_json_file(filepath)
        assert state.name == worker_name
        assert state.status == "dead"

        clients = WorkerClients.lscan(dirpath=None)
        assert len(clients) > 0
        assert repr(clients)
        assert str(clients)
        assert all(isinstance(c, WorkerClient) for c in clients)
        clients.print_dataframe()
        assert clients.get_dataframe() is not None
        d = clients.as_dict()
        same_clients = WorkerClients.from_dict(d)
        assert type(clients) is type(same_clients)
        print("same_clients\n", same_clients)
        # FIXME
        #assert len(clients) == len(same_clients)
        #self.assertMSONable(clients)

        c = clients.select_from_worker_name(worker_name)
        assert repr(c)
        assert str(c)
        assert c.worker_state.name == worker_name
        assert c.is_local_worker and not c.is_remote_worker
        d = c.as_dict()
        same_c = WorkerClient.from_dict(d)
        assert same_c.worker_state.name == c.worker_state.name
        #c.check_server_url()

        #c.get_json_status()
        #assert c.ssh_destination == ""

        if os.path.isdir(worker_dir):
            remove(worker_dir)

        clients = WorkerClients.lscan(dirpath=None)
        with self.assertRaises(ValueError)
            clients.select_from_worker_name(worker_name)
