""""Tests for analyzer module."""
import abipy.data as abidata

from abipy.core.testing import AbipyTest
from abipy.dynamics.analyzer import MdAnalyzer, MultiMdAnalyzer, ArrheniusPlotter


class AnalyzerTest(AbipyTest):
    pass

    #def test_md_analyzer(self):
    #def test_multi_md_analyzer(self):

    #def test_arrhenius_plotter(self):
    #    key_path = {
    #        "matgl-MD":  "diffusion_cLLZO-matgl.csv",
    #        "m3gnet-MD": "diffusion_cLLZO-m3gnet.csv",
    #    }

    #    mpl_style_key = {
    #        "matgl-MD" : dict(c='blue'),
    #        "m3gnet-MD": dict(c='purple'),
    #    }

    #    plotter = ArrheniusPlotter()
    #    for key, path in key_path.items():
    #        plotter.add_entry_from_file(path, key, mpl_style=mpl_style_key[key])

    #    # temperature grid refined
    #    thinvt_arange = (0.6, 2.5, 0.01)
    #    xlims, ylims = (0.5, 2.0), (-7, -4)

    #    plotter.plot(thinvt_arange=thinvt_arange,
    #                 xlims=xlims, ylims=ylims, text='LLZO cubic', savefig=None)
