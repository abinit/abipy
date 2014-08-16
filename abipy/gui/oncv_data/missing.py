#!/usr/bin/env python
import os

from pymatgen import periodic_table

input_files = [f for f in os.listdir(".") if f.endswith(".dat")]

#for f in input_files:
#    os.rename(f, f.replace(".dat", ".in"))

input_files = [f.split(".dat")[0] for f in os.listdir(".") if f.endswith(".dat")]
symbols_done = [f.split("-")[0] for f in input_files]
print("done:\n", symbols_done, "\nend done")

for s in periodic_table.all_symbols():
  if s not in symbols_done:
    print(s)
