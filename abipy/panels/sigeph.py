"""Panels to interact with the SIGEPH file."""
import param
import panel as pn
import panel.widgets as pnw

from .core import PanelWithElectronBands, mpl, ply, depends_on_btn_click #, PanelWithEbandsRobot


class SigEPhFilePanel(PanelWithElectronBands):
    """
    Panel with widgets to interact with a |SigEphFile|.
    """
    def __init__(self, sigeph, **params):
        PanelWithElectronBands.__init__(self, ebands=sigeph.ebands, **params)
        self.sigeph = sigeph

        self.sigma_spin_select = pnw.Select(name="Spin index", options=list(range(sigeph.nsppol)))
        self.sigma_kpoint_select = pnw.Select(name="Kpoint in Sigma_nk",
                options={"[%d]: %s" % (ik, repr(k)): ik for ik, k in enumerate(sigeph.sigma_kpoints)})

        self.plot_qpsolution_btn = pnw.Button(name="Plot Sigma_nk", button_type='primary')
        #sigma_band_select  = param.ObjectSelector(default=0, objects=[0], doc="Band index in sigma_nk")

    def plot_lws(self):
        # Insert results in grid.
        gspec = pn.GridSpec(sizing_mode='scale_width')
        for i, rta_type in enumerate(("serta", "mrta")):
            fig = self.sigeph.plot_lws_vs_e0(rta_type=rta_type, **self.mpl_kwargs)
            gspec[i, 0] = fig
            fig = self.sigeph.plot_tau_vtau(rta_type=rta_type, **self.mpl_kwargs)
            gspec[i, 1] = fig

        return gspec

    def plot_qpgaps(self):
        # Insert results in grid.
        gspec = pn.GridSpec(sizing_mode='scale_width')

        for i, qp_kpt in enumerate(self.sigeph.sigma_kpoints):
            fig = self.sigeph.plot_qpgaps_t(qp_kpoints=qp_kpt, qp_type="qpz0", **self.mpl_kwargs)
            gspec[i, 0] = fig
            fig = self.sigeph.plot_qpgaps_t(qp_kpoints=qp_kpt, qp_type="otms", **self.mpl_kwargs)
            gspec[i, 1] = fig

        return gspec

    def plot_qps_vs_e0(self):
        return mpl(self.sigeph.plot_qps_vs_e0(**self.mpl_kwargs))

    @depends_on_btn_click('plot_qpsolution_btn')
    def on_plot_qpsolution_sk(self):
        fig = self.sigeph.plot_qpsolution_sk(self.sigma_spin_select.value, self.sigma_kpoint_select.value,
                                             itemp=0, with_int_aw=True,
                                             ax_list=None, xlims=None, fontsize=8, **self.mpl_kwargs)
        return mpl(fig)

    def get_panel(self, as_dict=False, **kwargs):
        """Return tabs with widgets to interact with the SIGEPH file."""
        d = {}
        d["Summary"] = self.get_summary_view_for_abiobj(self.sigeph)
        #app(("e-Bands", pn.Row(self.get_plot_ebands_widgets(), self.on_plot_ebands_btn)))

        # Build different tabs depending on the calculation type.
        # TODO
        #if self.sigeph.imag_only:
        #    d["LWS"] = self.plot_lws

        #else:
        #    d["QP-gaps"] = self.plot_qpgaps
        #    d["QP_vs_e0"] = self.plot_qps_vs_e0
        #    d["QPSolution"] = pn.Row(
        #        pn.Column("## Quasiparticle solution",
        #                  self.sigma_spin_select,
        #                  self.sigma_kpoint_select,
        #                  self.plot_qpsolution_btn),
        #        self.on_plot_qpsolution_sk
        #    )

        #d["Structure"] = self.get_structure_view()
        d["NcFile"] = self.sigeph.get_ncfile_view()

        if as_dict: return d
        return self.get_template_from_tabs(d, template=kwargs.get("template", None))
