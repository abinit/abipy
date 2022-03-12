import os
import tempfile

from abipy.core.testing import AbipyTest
from abipy.flowtk.abiinspect import *
import abipy.data as abidata


class YamlTokenizerTest(AbipyTest):
    """Test YamlTokenizer."""

    def test_base(self):
        string = """---
none: [~, null]
bool: [true, false, on, off]
int: 42
float: 3.14159
list: [LITE, RES_ACID, SUS_DEXT]
dict: {hp: 13, sp: 5}
...

this is not a YAML document!
and the tokenizer will ignore it

--- !Monster
name: Cave spider
hp: [2,6]    # 2d6
ac: 16
attacks: [BITE, HURT]
...

This is not a proper document since it does not start with ---
the end tag below is ignored
...
--- !Monster
name: Dragon
hp: [2,6]    # 2d6
ac: 32
attacks: [BITE, HURT]
...
            """
        # for i, line in enumerate(string.splitlines()): print(i, line)
        fd, filename = tempfile.mkstemp(text=True)

        with open(filename, "w") as fh:
            fh.write(string)

        doc_tags = [None, "!Monster", "!Monster"]
        doc_linenos = [1, 13, 23]

        with YamlTokenizer(filename) as r:
            # Iterate the docs
            n = 0
            for i, doc in enumerate(r):
                n += 1
                str(doc)
                #print("doc", doc)
                assert doc.tag == doc_tags[i]
                assert doc.lineno == doc_linenos[i]

            assert n == len(doc_tags)

            # Read all docs present in the file.
            r.seek(0)
            all_docs = r.all_yaml_docs()
            # print(all_docs)
            assert len(all_docs) == 3

            # We should be at the begining at the file.
            assert all_docs == r.all_yaml_docs()

            # Find documents by tag.
            r.seek(0)
            monster = r.next_doc_with_tag("!Monster")
            # print("monster",monster)
            assert monster == all_docs[1]

            monster = r.next_doc_with_tag("!Monster")
            assert monster == all_docs[2]

            # this should raise StopIteration
            with self.assertRaises(StopIteration):
                monster = r.next_doc_with_tag("!Monster")


class AbinitInpectTest(AbipyTest):
    def test_scfcycle(self):
        """Testing ScfCycle."""
        filepath = os.path.join(abidata.dirpath, "refs", "text_files", "mgb2_scf.abo")
        cycle = GroundStateScfCycle.from_file(filepath)
        str(cycle)
        cycle.to_string(verbose=2)

        assert cycle.num_iterations == 6
        last = cycle.last_iteration

        assert last["Etot(hartree)"] == -7.1476241568657 and last["vres2"] == 3.879e-08
        assert list(cycle["vres2"]) == [
            1.769e02,
            7.920e-01,
            1.570e-01,
            4.259e-03,
            4.150e-05,
            3.879e-08,
        ]

        # Testing CyclesPlotter.
        p = CyclesPlotter()
        p.add_label_cycle("mgb2 SCF", cycle)
        p.add_label_cycle("same SCF", cycle)

        if self.has_matplotlib():
            assert cycle.plot(show=False)
            assert p.combiplot(show=False)
            #p.slideshow()

        if self.has_plotly():
            assert cycle.plotly(show=False)

    def test_relaxation(self):
        """Testing Relaxation object."""
        filepath = os.path.join(abidata.dirpath, "refs", "text_files", "sic_relax.abo")
        relaxation = Relaxation.from_file(filepath)
        str(relaxation)
        assert len(relaxation) == 4

        assert relaxation[0]["Etot(hartree)"][-1] == -8.8077409200473
        assert relaxation[-1]["Etot(hartree)"][-1] == -8.8234906607147

        for scf_step in relaxation:
            str(scf_step.num_iterations)

        if self.has_matplotlib():
           assert relaxation.plot(show=False)
           #assert relaxation.slideshow(show=False)

        #if self.has_plotly():
        #   assert relaxation.plotly(show=False)
