import os
from abipy.flowtk import PseudoTable

_root = os.path.dirname(__file__)
_paths = [f for f in os.listdir(_root) if f.endswith("hgh")]
# Need one pseudo for element
d = {}
for p in _paths:
    i = p.index(".")
    head, tail = p[:i], p[i:]
    d[head] = tail

_paths = [k + v for k, v in d.items()]
_paths = [os.path.join(_root, f) for f in _paths]
del d

HGH_TABLE = PseudoTable(_paths)

# Add fake hints.
for pseudo in HGH_TABLE:
    pseudo.dojo_report = {}
    pseudo.dojo_report["hints"] = {}
    for accuracy in ["low", "normal", "high"]:
        pseudo.dojo_report["hints"][accuracy] = {"ecut": 50}
    #assert pseudo.has_hints
