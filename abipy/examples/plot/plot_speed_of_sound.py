#!/usr/bin/env python
r"""
Speed of Sound
==============

This example shows how to compute the speed of sound by fitting phonon frequencies
along selected directions by linear least-squares fit.
For a command line interface use:

.. code-block:: bash

    abiview.py ddb_vs DDB_FILE
"""

#%%
# Initialize object from DDB file.

import os
import abipy.data as abidata

from abipy import abilab
from abipy.dfpt.vsound import SoundVelocity

ddb_path = os.path.join(abidata.dirpath, "refs", "si_sound_vel", "Si_DDB")
sv = SoundVelocity.from_ddb(ddb_path)

#%%
# Get pandas dataframe with results.

df = sv.get_dataframe()
abilab.print_dataframe(df)

#%%
# Plot fit with matplotlib

sv.plot()

#%%
# Plot fit with plotly

sv.plotly(template="plotly_dark")
