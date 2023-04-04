#!/usr/bin/env python

from abipy.electrons.effmass_analyzer import EffMassAnalyzer
import sys

from monty.json import json, MontyEncoder
import numpy as np

#carr = np.array([1j, 2j])
#print(json.dumps(carr, cls=MontyEncoder))

emana = EffMassAnalyzer.from_file(sys.argv[1])
print(emana)
#print(emana.kpoints)

#emana.select_kpoint_band((0, 0, 0), band=7, etol_ev=0.1)
#emana.set_kpoint_band((0, 0, 0), band=3)
#emana.select_band_edges()
#emana.select_cbm()
emana.select_vbm(etol_ev=1e-1)

emana.summarize()
emana.plot_all_segments()
emana.plot_emass(acc=4)

#for segment in emana.segments:
#    segment.get_effmass_line(acc=2)

segment = emana.segments[0]
segment.plot_emass()
print(segment.get_dataframe_with_accuracies())
