#!/usr/bin/env python
"""Tests for cube module"""
from __future__ import print_function, division, unicode_literals, absolute_import

import tempfile
import unittest
import numpy as np
import abipy.data as data

from abipy.core.testing import AbipyTest
from abipy.core.fields import Density, core_density_from_file
from abipy.core.mesh3d import Mesh3D

#filepath = data.ref_file("si_DEN-etsf.nc")
#density = Density.from_file(filepath=filepath)

#mesh = density.mesh

# density.structure.get_sites_in_sphere([0.0, 0.0, 0.0], 10.0)
# print(mesh.i_closest_gridpoints(points=[density.structure[1].coords]))
# print(mesh.dist_gridpoints_in_spheres(points=[density.structure[0].coords], radius=2.0))

# print(np.sum(density.datar)* density.mesh.dv)
#
# density.export_to_cube(filename='valence.cube')
#
# density2 = Density.from_cube(filename='valence.cube')
#
#
# for maxr in [0.01, 0.1, 0.3, 0.5, 1.0, 8.0]:
#
#
# class TestCubeUtils(AbipyTest):
#
#     @unittest.skip("Si.in.rhoc file is missing!")
#     def test_aecore_density(self):
#         """Testing ae_core_density_on_mesh."""
#         maxr = 8.0
#         rhoc = {'Si': core_density_from_file('Si.in.rhoc')}
#         ae_density_new = Density.ae_core_density_on_mesh(valence_density=density, structure=density.structure,
#                                                          rhoc=rhoc, maxr=maxr,
#                                                          method='mesh3d_dist_gridpoints', small_dist_factor=5.0,
#                                                          small_dist_mesh=(20, 20, 20))
#
#         ae_density = Density.ae_core_density_on_mesh(valence_density=density, structure=density.structure,
#                                                      rhoc=rhoc, maxr=maxr,
#                                                      method='get_sites_in_sphere', small_dist_factor=5.0,
#                                                      small_dist_mesh=(20, 20, 20))
#
#         print(ae_density_new.nelect_updown)
#         print(ae_density.nelect_updown)
#
# for maxr in [0.9]:
#     print('maxr', maxr)
#     ae_density_new = Density.ae_core_density_on_mesh(valence_density=density, structure=density.structure,
#                                                      rhoc_files={'Si': 'Si.in.rhoc'}, maxr=maxr,
#                                                      method='mesh3d_dist_gridpoints')
#     print(ae_density_new.nelect_updown)
#     ae_density = Density.ae_core_density_on_mesh(valence_density=density, structure=density.structure,
#                                                  rhoc_files={'Si': 'Si.in.rhoc'}, maxr=maxr, method='get_sites_in_sphere')
#     print(ae_density.nelect_updown)

#print(ae_density.datar/ae_density_new.datar)
# print(np.allclose(ae_density.datar/ae_density_new.datar,
#                   ae_density.datar[0]/ae_density_new.datar[0]*np.ones_like(ae_density.datar)))
#
# print(np.allclose(ae_density.datar, ae_density_new.datar))
#
# ae_density.export_to_cube(filename='core.cube')
#
# total_density = density+ae_density
#
# total_density.export_to_cube(filename='total.cube')
#
# # class TestCubeUtils(AbipyTest):
# #     """Unit tests for Density."""
# #
# #
# #     def test_cube_export(self):
# #         """Testing density in the cube format."""
# #         filepath = data.ref_file("si_DEN-etsf.nc")
# #         density = Density.from_file(filepath=filepath)
# #         density.export_to_cube(filename='test.cube')
