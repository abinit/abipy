# # coding: utf-8
# """
# Utility tasks for Fireworks.
# """
#
# from __future__ import print_function, division, unicode_literals
#
# try:
#     from fireworks.core.firework import Firework, FireTaskBase, FWAction
#     from fireworks.core.launchpad import LaunchPad
#     from fireworks.utilities.fw_utilities import explicit_serialize
#     from fireworks.utilities.fw_serializers import serialize_fw
# except ImportError:
#     FireTaskBase, FWAction, Firework, LaunchPad = 4 * [object]
#     explicit_serialize = lambda x: x
#     serialize_fw = lambda x: x
#
# import os
# import shutil
# import logging
# import traceback
# import importlib
# from abipy.fworks.fw_tasks import INDIR_NAME, OUTDIR_NAME, TMPDIR_NAME
# from abipy.fworks.fw_databases import MongoDatabase
# from monty.serialization import loadfn
#
#
# logger = logging.getLogger(__name__)
#
#
# @explicit_serialize
# class FinalCleanUpTask(FireTaskBase):
#
#     def __init__(self, out_exts=["WFK", "1WF"]):
#         self.out_exts = out_exts
#
#     @serialize_fw
#     def to_dict(self):
#         return dict(out_exts=self.out_exts)
#
#     @classmethod
#     def from_dict(cls, m_dict):
#         return cls(out_exts=m_dict['out_exts'])
#
#     @staticmethod
#     def delete_files(d, exts=None):
#         deleted_files = []
#         if os.path.isdir(d):
#             for f in os.listdir(d):
#                 if exts is None or "*" in exts or any(ext in f for ext in exts):
#                     fp = os.path.join(d, f)
#                     try:
#                         if os.path.isfile(fp):
#                             os.unlink(fp)
#                         elif os.path.isdir(fp):
#                             shutil.rmtree(fp)
#                         deleted_files.append(fp)
#                     except:
#                         logger.warning("Couldn't delete {}: {}".format(fp, traceback.format_exc()))
#
#         return deleted_files
#
#     def run_task(self, fw_spec):
#         # the FW.json/yaml file is mandatory to get the fw_id
#         # no need to deserialize the whole FW
#         try:
#             fw_dict = loadfn('FW.json')
#         except IOError:
#             try:
#                 fw_dict = loadfn('FW.yaml')
#             except IOError:
#                 raise RuntimeError("No FW.json nor FW.yaml file present: impossible to determine fw_id")
#
#         fw_id = fw_dict['fw_id']
#         lp = LaunchPad.auto_load()
#         wf = lp.get_wf_by_fw_id_lzyfw(fw_id)
#
#         deleted_files = []
#         # iterate over all the fws and launches
#         for fw_id, fw in wf.id_fw.items():
#             for l in fw.launches:
#                 l_dir = l.launch_dir
#
#                 deleted_files.extend(self.delete_files(os.path.join(l_dir, TMPDIR_NAME)))
#                 deleted_files.extend(self.delete_files(os.path.join(l_dir, INDIR_NAME)))
#                 deleted_files.extend(self.delete_files(os.path.join(l_dir, OUTDIR_NAME), self.out_exts))
#
#         logging.info("Deleted files:\n {}".format("\n".join(deleted_files)))
#
#         return FWAction(stored_data={'deleted_files': deleted_files})
#
