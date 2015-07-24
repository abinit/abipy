# from __future__ import print_function, division, unicode_literals
# from pymatgen.io.abinitio.tasks import ParalHints
#
# try:
#     from fireworks import Workflow
# except ImportError:
#     Workflow = object
#
# import traceback
# from abipy.fworks.fw_tasks import FWTaskManager
# import pymatgen.io.abinitio.qutils as qu
#
# SHORT_SINGLE_CORE_SPEC = {'_queueadapter': {'ntasks': 1, 'time': '00:10:00'}, 'mpi_ncpus': 1}
#
#
# def parse_workflow(fws, links_dict):
#     new_list = []
#     for fw in fws:
#         if isinstance(fw, Workflow):
#             new_list.extend(fw.fws)
#         else:
#             new_list.append(fw)
#
#     new_links_dict = {}
#     for parent, children in links_dict.items():
#         if isinstance(parent, Workflow):
#             new_links_dict.update(parent.links)
#             for leaf_fw_id in parent.leaf_fw_ids:
#                 new_links_dict[leaf_fw_id] = children
#         else:
#             new_links_dict[parent] = children
#
#     # dict since the collection will be updated
#     for parent, children in dict(new_links_dict).items():
#         final_childrens = []
#         for child in children:
#             if isinstance(child, Workflow):
#                 new_links_dict.update(child.links)
#                 final_childrens.extend(child.root_fw_ids)
#             else:
#                 final_childrens.append(child)
#         new_links_dict[parent] = final_childrens
#
#     return new_list, new_links_dict
#
#
# def append_fw_to_wf(new_fw, wf):
#
#     if new_fw.fw_id in wf.id_fw:
#         raise ValueError('FW ids must be unique!')
#
#     for leaf_id in wf.leaf_fw_ids:
#         wf.links[leaf_id].append(new_fw.fw_id)
#
#     wf.links[new_fw.fw_id] = []
#     wf.id_fw[new_fw.fw_id] = new_fw
#
#     # add depends on
#     for pfw in new_fw.parents:
#         if pfw.fw_id not in wf.links:
#             raise ValueError(
#                 "FW_id: {} defines a dependent link to FW_id: {}, but the latter was not added to the workflow!".format(
#                     new_fw.fw_id, pfw.fw_id))
#         wf.links[pfw.fw_id].append(new_fw.fw_id)
#
#     # sanity: make sure the set of nodes from the links_dict is equal to
#     # the set of nodes from id_fw
#     if set(wf.links.nodes) != set(map(int, wf.id_fw.keys())):
#         raise ValueError("Specified links don't match given FW")
#
#     wf.fw_states[new_fw.fw_id] = new_fw.state
#
#
# def get_short_single_core_spec(fw_manager_path=None):
#     if fw_manager_path:
#         ftm = FWTaskManager.from_file(fw_manager_path)
#     else:
#         ftm = FWTaskManager.from_user_config()
#
#     if ftm.has_task_manager():
#         #TODO add mem_per_cpu?
#         pconf = ParalHints({}, [{'tot_ncpus': 1, 'mpi_ncpus': 1, 'efficiency': 1}])
#         try:
#             tm = ftm.task_manager
#             tm.select_qadapter(pconf)
#             #TODO make a FW_task_manager parameter
#             tm.qadapter.set_timelimit(600)
#             qadapter_spec = tm.qadapter.get_subs_dict()
#             return qadapter_spec
#         except RuntimeError as e:
#             traceback.print_exc()
#
#     # No taskmanger or no queue available
#     #FIXME return something else? exception?
#     return {}
#
#
#
#
#
#
#
#
#
