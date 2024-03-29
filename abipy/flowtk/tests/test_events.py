# coding: utf-8
import os
import datetime

from abipy.core.testing import AbipyTest
from abipy.flowtk import events

_test_dir = os.path.join(os.path.dirname(__file__), "..", "..", 'test_files')

def ref_file(filename):
    return os.path.join(_test_dir, filename)


def ref_files(*filenames):
    return list(map(ref_file, filenames))


class EventsParserTest(AbipyTest):

    def test_mgb2_outputs(self):
        """Testing MgB2 output files."""
        # Analyze scf log
        parser = events.EventsParser()
        report = parser.parse(ref_file("mgb2_scf.log"), verbose=1)
        self.assert_msonable(report)

        assert str(report)
        assert (report.num_errors, report.num_warnings, report.num_comments) == (0, 0, 0)
        assert report.run_completed
        fmt = "%a %b %d %H:%M:%S %Y"
        assert report.start_datetime ==  datetime.datetime.strptime("Fri Mar 13 20:08:51 2015", fmt)
        assert report.end_datetime ==  datetime.datetime.strptime("Fri Mar 13 20:08:57 2015", fmt)

        # Analyze nscf log
        report = events.EventsParser().parse(ref_file("mgb2_nscf.log"), verbose=0)
        assert (report.num_errors, report.num_warnings, report.num_comments) == (0, 2, 0)
        print(report)
        self.assert_msonable(report)

        #d = report.as_dict()
        #print(d)
        #assert 0

        for i, warning in enumerate(report.warnings):
            print(warning)
            assert warning == report[i]
            # Msonable is conflict with YAMLObject
            #self.assert_msonable(warning, check_inst=False)

        report = parser.report_exception(ref_file("mgb2_scf.log"), "exception")
        assert len(report.errors) == 1

    def test_parse_bad_yaml_doc(self):
        """Parsing Abinit log file with wrong YAML document."""
        parser = events.EventsParser()
        report = parser.parse(ref_file("badyaml.log"), verbose=1)
        print(report)
        assert not report.run_completed
        assert (report.num_errors, report.num_warnings, report.num_comments) == (1, 1, 0)

        # The event parser should have registered a AbinitYamlWarning and an AbinitYamlError
        assert len(report.get_events_of_type(events.AbinitYamlWarning)) == 1
        assert len(report.get_events_of_type(events.AbinitYamlError)) == 1


class EventHandlersTest(AbipyTest):
    def test_events(self):
        # Autodoc
        events.autodoc_event_handlers()

        for cls in events.get_event_handler_classes():
            # Test pickle
            handler = cls()
            self.serialize_with_pickle(handler, test_eq=False)
            self.assert_msonable(handler)

        assert events.as_event_class(events.AbinitWarning) == events.AbinitWarning
        assert events.as_event_class('!WARNING') == events.AbinitWarning
