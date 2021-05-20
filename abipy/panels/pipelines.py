"""
Web applications based on panel pipelines: https://panel.holoviz.org/user_guide/Pipelines.html
"""

import param

import panel as pn
import panel.widgets as pnw
import bokeh.models.widgets as bkw

#from abipy.panels.core import AbipyParameterized
from abipy.abilab import abiopen, extcls_supporting_panel


class UploadFile(param.Parameterized):

    file_input = pnw.FileInput()

    abipath = param.Filename()
    #ready = param.Boolean(default=False, precedence=-1)

    @param.output("abipath")
    def output(self):
        #print("filename", self.file_input.filename)
        if self.file_input.value is None: return None
        #print("value", self.file_input.value)
        import tempfile
        workdir = tempfile.mkdtemp()

        fd, tmp_path = tempfile.mkstemp(suffix=self.file_input.filename)
        print(tmp_path)
        with open(tmp_path, "wb") as fh:
            fh.write(self.file_input.value)

        #print("abipath", self.abipath)
        self.abipath = tmp_path
        return tmp_path

    @param.depends("file_input.filename")
    def view(self):
        if self.file_input.value is None: return None

        abipath = self.output()
        #with abiopen(abipath) as abifile:
        abifile = abiopen(abipath)
        row = pn.Row(bkw.PreText(text=abifile.to_string())) #, sizing_mode="scale_both"))
        if hasattr(abifile, "close"): abifile.close()
        return row

    def panel(self):

        help_str = """
This web app exposes some of the post-processing capabilities of AbiPy.

Use the **Choose File** to upload one of the files supported by the app, then
click the **Next** button to analyze the file.
"""

        #table_str = "\n\n" + extcls_supporting_panel(tablefmt="simple") + "\n\n"
        table_str = extcls_supporting_panel(tablefmt="html")

        # Add accordion after the button with warning and help taken from the docstring of the callback
        col = pn.Column(); ca = col.append
        acc = pn.Accordion(
                ("Help", pn.pane.Markdown(help_str, name="help")),
                ("Supported Extensions", pn.pane.HTML(table_str, name="extensions")),
        )
        ca(pn.layout.Divider())
        ca(acc)

        return pn.Row(
                pn.Column(#"## Upload file:",
                          self.file_input,
                          col,
                         ),
                self.view)


#class Stage2(param.Parameterized):
class Analyze(UploadFile):

    #ready = param.Boolean(default=True)

    @param.depends("abipath")
    def view(self):
        print("stage2 abifile", self.abipath)
        #with abiopen(self.abipath) as abifile:
        #    return abifile.get_panel()

        abifile = abiopen(self.abipath)
        return abifile.get_panel()

    def panel(self):
        print("in panel")
        return pn.Row(self.view)


def analyze_file_app(**kwargs):

    pn.extension("plotly")
    pipeline = pn.pipeline.Pipeline(**kwargs)

    #pipeline.add_stage('Stage 1', stage1)
    #stage2 = Analyze(abipath=stage1.output())
    #pipeline.add_stage('Stage 2', stage2)
    pipeline.add_stage('Upload File', UploadFile)
    pipeline.add_stage('Analyze File', Analyze)
    print(pipeline)
    pipeline.show()



#class UploadFiles(param.Parameterized):
#
#    file_inputs = pnw.FileInput(multiple=True)
#
#    #abipath = param.Filename()
#
#    @param.output("abipath")
#    def output(self):
#        #print("filename", self.file_inputs.filename)
#        if self.file_inputs.value is None: return None
#        #print("value", self.file_inputs.value)
#        import tempfile
#        workdir = tempfile.mkdtemp()
#
#        for filename in self.file_inputs.value:
#
#            fd, tmp_path = tempfile.mkstemp(suffix=filename)
#            print(tmp_path)
#            with open(tmp_path, "wb") as fh:
#                fh.write(self.file_inputs.value)
#
#        #print("abipath", self.abipath)
#        self.abipath = tmp_path
#        return tmp_path
#
#    @param.depends("file_inputs.filename")
#    def view(self):
#        if self.file_inputs.value is None: return None
#
#        abipath = self.output()
#        #with abiopen(abipath) as abifile:
#        abifile = abiopen(abipath)
#        row = pn.Row(bkw.PreText(text=abifile.to_string())) #, sizing_mode="scale_both"))
#        if hasattr(abifile, "close"): abifile.close()
#        return row
#
#    def panel(self):
#
#        help_str = """
#This web app exposes some of the post-processing capabilities of AbiPy.
#
#Use the **Choose File** to upload one of the files supported by the app, then
#click the **Next** button to analyze the file.
#"""
#
#        #table_str = "\n\n" + extcls_supporting_panel(tablefmt="simple") + "\n\n"
#        table_str = extcls_supporting_panel(tablefmt="html")
#
#        # Add accordion after the button with warning and help taken from the docstring of the callback
#        col = pn.Column(); ca = col.append
#        acc = pn.Accordion(
#                ("Help", pn.pane.Markdown(help_str, name="help")),
#                ("Supported Extensions", pn.pane.HTML(table_str, name="extensions")),
#        )
#        ca(pn.layout.Divider())
#        ca(acc)
#
#        return pn.Row(
#                pn.Column(#"## Upload file:",
#                          self.file_inputs,
#                          col,
#                         ),
#                self.view)
#
#
#def analyze_files_app(**kwargs):
#
#    pn.extension("plotly")
#    pipeline = pn.pipeline.Pipeline(**kwargs)
#
#    #pipeline.add_stage('Stage 1', stage1)
#    #stage2 = Analyze(abipath=stage1.output())
#    #pipeline.add_stage('Stage 2', stage2)
#    pipeline.add_stage('Upload File', UploadFiles)
#    #pipeline.add_stage('Analyze File', Analyze)
#    #print(pipeline)
#    pipeline.show()
