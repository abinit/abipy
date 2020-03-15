# coding: utf-8
from abipy.core.testing import AbipyTest
from abipy.flowtk.utils import *


# FIXME
#class FilePathFixerTest(AbipyTest):
#
#    def test_base(self):
#        fixer = FilepathFixer()
#        #assert fixer.fix_paths('/foo/out_1WF17') == {'/foo/out_1WF17': '/foo/out_1WF'}
#        #assert fixer.fix_paths('/foo/out_1WF5.nc') == {'/foo/out_1WF5.nc': '/foo/out_1WF.nc'}
#        assert fixer.fix_paths('/foo/out_1DEN17') == {'/foo/out_1DEN17': '/foo/out_1DEN'}
#        assert fixer.fix_paths('/foo/out_1DEN5.nc') == {'/foo/out_1DEN5.nc': '/foo/out_1DEN.nc'}

class RpnTest(AbipyTest):

    def test_mongodb_like_conditions(self):
        class Foo(object):
            one = 1.0
            two = 2.0
            three = 3.0
            four = 4.0

        map_res = [
            ( {"one": 1.0}, True),
            ( {"one": {"$eq": 1.0}}, True),
            ( {"one": {"$eq": "one"}}, True),
            ( {"one": {"$ne": "two"}}, True),
            ( {"one": {"$ne": 1.0}}, False),
            ( {"four": {"$divisible": 2.0}}, True),
            ( {"four": {"$divisible": 3.0}}, False),
            ( {"two": {"$gt": "one"}}, True ),
            ( {"$and": [ {"one": 1.0}, {"two": {"$lt": 3}}]}, True),
            ( {"$and": [{"one": {"$ge": 0.8}}, {"two": {"$le": 6.0}}]}, True),
            ( {"$or": [ {"$not": {"one": 1.0}}, {"two": {"$lt": 20}}]}, True),
            ( {"$not": {"$and": [ {"$not": {"one": 1.0}}, {"two": {"$lt": 3}}] }}, True),
        ]

        for map, res in map_res:
            print("map", map)
            rpn = map2rpn(map, obj=Foo)
            print("rpn", rpn)
            self.assertEqual(res, evaluate_rpn(rpn), msg="map %s" % map)


class ConditionTest(AbipyTest):
    def test_condition(self):
        c = Condition({})
        assert not c
        print(c)

        class A(object):
            def __init__(self):
                self.one = 1.0

        aobj = A()
        assert Condition({"one": 1.0})(aobj)
        assert not Condition({"one": 2.0})(aobj)


class SparseHistogramTest(AbipyTest):
    def test_sparse(self):
        items = [1, 2, 2.9, 4]
        hist = SparseHistogram(items, step=1)
        assert hist.binvals == [1.0, 2.0, 3.0]
        assert hist.values == [[1], [2, 2.9], [4]]
        #hist.plot()

        hist = SparseHistogram([iv for iv in enumerate(items)], key=lambda t: t[1], step=1)
        assert hist.binvals == [1.0, 2.0, 3.0]
        assert hist.values == [[(0, 1)], [(1, 2), (2, 2.9)], [(3, 4)]]
