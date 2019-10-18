#!/usr/bin/env python
r"""
Electrons and Phonons from the materials project website
========================================================

This example shows how to dowload the electronic band structure
and the DDB file using the mp identifier and use the AbiPy API
to generate a matplotlib grid with electrons + phonons.

IMPORTANT: Electrons and Phonons have been obtained with different codes
and different computational settings! Of course, one can always
initialize ElectronBands and PhononBands from local netcdf files
obtained with Abinit
"""
from abipy import abilab

# List of mp ids for Si, Diamond
mpids = ["mp-149", "mp-66"]

# Get list of AbiPy ebands from mpids
ebands_list = [abilab.ElectronBands.from_mpid(mpid) for mpid in mpids]

# Get list of DDB files from the MP website and run anaddb to get the phonon bands.
phbands_list = []
for i, mpid in enumerate(mpids):
    print("Downloading DDB for mpid %s (%s) ..." % (mpid, ebands_list[i].structure.formula))
    ddb = abilab.DdbFile.from_mpid(mpid)
    if ddb is None:
        raise RuntimeError("%d does not provide DDB" % mpid)
    print("Invoking anaddb to compute phonon bands...")
    phbst, _ = ddb.anaget_phbst_and_phdos_files(nqsmall=0)
    phbands_list.append(phbst.phbands)
    phbst.close()

# The figure has [len(mpids), 2] subplots
# The i-th row contains electrons and phonons for the i-th mp identifier.
nrows, ncols = len(mpids), 2
ax_mat, fig, plt = abilab.get_axarray_fig_plt(None, nrows=nrows, ncols=ncols,
                                              sharex=False, sharey=False, squeeze=False)

# Use the `ax` keyword argument to select the matplotlib Axes used to plot the object.
# In the band structure plot, we show the fundamental/direct gap as well as the possible
# phonon-absorption (-emission) processes allowed by energy-conservation.
# (This is a qualitative analysis of e-ph scattering, quasi-momentum and ph dispersion are not taken into account).
for i, (ebands, phbands) in enumerate(zip(ebands_list, phbands_list)):
    ebands.plot(ax=ax_mat[i, 0], with_gaps=True, ylims=(-5, 10), max_phfreq=phbands.maxfreq, show=False)
    phbands.plot(ax=ax_mat[i, 1], show=False)

    # Hide xlabel if not last row.
    if i != len(ebands_list) - 1:
        for ax in ax_mat[i]:
            ax.xaxis.label.set_visible(False)

plt.show()
