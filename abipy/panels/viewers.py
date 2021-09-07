"""Viewer objects."""

import os
#import param
import panel as pn
import panel.widgets as pnw

from panel.viewable import Viewer


class AceViewer(Viewer):

    def __init__(self, filepath, theme="terminal", height=1200, **params):
        self.filepath = filepath
        super().__init__(**params)

        basename = os.path.basename(filepath)
        self.open_btn = pnw.Button(name=f"Open {basename}", button_type='primary')
        self.open_btn.on_click(self.open_ace_editor)

        self.ace = pnw.Ace(language='text', readonly=True, theme=theme,
                           sizing_mode='stretch_width', height=height, visible=False)

        self.controls = pn.Card(self.ace.param.height, self.ace.param.theme, self.ace.param.visible,
                                title="ACE controls", collapsed=True)

        self.layout = pn.Column(f"## File: <small>{filepath}</small>",
                                pn.Row(self.open_btn, self.controls),
                                self.ace,
                                pn.layout.Divider(),
                                sizing_mode="stretch_width")

    def open_ace_editor(self, event):
        self.ace.visible = True
        self.ace.value = open(self.filepath, "rt").read()
        self.open_btn.name = "Reopen %s" % os.path.basename(self.filepath)

    def __panel__(self):
        return self.layout


class JSONViewer(Viewer):

    def __init__(self, dictionary, theme="dark", hover_preview=True, depth=2, **params):

        super().__init__(**params)

        self.json = pn.pane.JSON(dictionary,
                                 depth=depth, # -1 for full expansion
                                 hover_preview=hover_preview,
                                 theme=theme,
                                 sizing_mode="stretch_width",
                                 )

        self.controls = pn.Card(self.json.param.depth,
                                self.json.param.hover_preview,
                                self.json.param.theme,
                                self.json.param.visible,
                                title="JSON controls",
                                collapsed=True,
                                #header_color="black",
                                #header_background="CornflowerBlue",
                                )

        #self.controls = pn.Card(self.json.controls(jslink=True), title="JSON controls", collapsed=True)

        self.layout = pn.Column(self.controls,
                                self.json,
                                sizing_mode="stretch_width")

    def __panel__(self):
        return self.layout
