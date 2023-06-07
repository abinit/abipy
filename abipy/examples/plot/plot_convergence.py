#!/usr/bin/env python
r"""
Plot convergence
================

This example shows how to plot the convergence of
scalar quantities wrt to one independent variable.
For instance, the convergence of the energy per atom
with the cutoff energy ecut.
"""

#%%
# Let's assume we have the total energy per atom
# and the pressure in GPa as a function of ecut
# stored in three python lists:

ecut_list = [35., 40., 45., 50., 55., 60., 65.]

energy_per_atom_ev = [
-444.29611271, -444.3113239, -444.3124781,
-444.31251066, -444.31254606, -444.31256315, -444.31258609,
]

pressure_gpa =  [
-17.56316746, -16.00634717, -15.80295881,
-15.79975677, -15.79858973, -15.79811464, -15.79929941,
]

#%%
# To plot energy vs ecut with a convergence window of 1e-3
# centered around our best estimate i.e. the value obtained
# with the highest ecut, use the high-level API `from_xy_label_vals`:

from abipy.tools.plotting import ConvergenceAnalyzer

ca = ConvergenceAnalyzer.from_xy_label_vals("ecut (Ha)", ecut_list,
                                            "E/natom (eV)", energy_per_atom_ev, tols=1e-3)

print(ca)
ca.plot()

#%%
# To analyze the convergence of two or more quantities, use the low-level API
# that expects two dictionaries mapping the y-label used in the plot to
# (1) the list of values and (2) the associated absolute tolerance as in:

yvals_dict = {
    "E/natom (eV)": energy_per_atom_ev,
    "P (GPa)": pressure_gpa,
}

ytols_dict = {
    "E/natom (eV)": 1e-3,
    "P (GPa)": 1e-1,
}

ca = ConvergenceAnalyzer("ecut (Ha)", ecut_list, yvals_dict, ytols_dict)

#print(ca)
#ca.plot()

#%%
# To specify multiple convergence criteria, replace scalars with tuples as in:

ytols_dict = {
    "E/natom (eV)": (1e-3, 1e-4),
    "P (GPa)": (1e-1, 1e-2),
}

ca = ConvergenceAnalyzer("ecut (Ha)", ecut_list, yvals_dict, ytols_dict)
#print(ca)
ca.plot()


#%%
# Finally, if you have a (xlsf or csv) file with columns named e.g.: `ecut`, `energy_per_atom`, `pressure`
# and you want to plot the results with ConvergenceAnalyzer, use the `from_file` method
# with a `ytols_dict` whose keys define the columns in tabular data:

ytols_dict = {
    "ecut": (1e-3, 1e-4),
    "pressure": (1e-1, 1e-2),
}
#ca = ConvergenceAnalyzer.from_file(filepath, ytols_dict)

# To change the x/y-labels in the plot, use:
#ca.set_label("ecut", "ecut (Ha)")
#ca.set_label("energy_per_atom", "E/natom (eV)")
#ca.set_label("pressure", "P (GPa)")

#print(ca)
#ca.plot()
