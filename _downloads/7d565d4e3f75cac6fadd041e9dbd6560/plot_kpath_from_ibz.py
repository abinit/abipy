#!/usr/bin/env python
r"""
K-path from IBZ
===============

This example shows how to extract energies along a k-path
from a calculation done with a (relatively dense) IBZ sampling.
"""
from abipy.abilab import abiopen
import abipy.data as abidata

# Open the file with energies computed with a homogeneous sampling of the BZ
# and extract the band structure object.
with abiopen(abidata.ref_file("si_scf_GSR.nc")) as gs_file:
    ebands_ibz = gs_file.ebands

# This is a GS calculation done with a 8x8x8 k-mesh.
print(ebands_ibz.kpoints)

# Build new ebands with energies along G-X-L-G path.
# Smooth bands require dense meshes.
r = ebands_ibz.with_points_along_path(knames=["G", "X", "L", "G"])

print(r.ebands)
r.ebands.plot()
