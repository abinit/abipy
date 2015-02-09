"""Directory with abipy + abinit tutorials."""

def help():
    """List the available tutorials with a brief description."""
    import os
    from importlib import import_module
    
    docs = []
    for pyfile in os.listdir(__file__):
        if not (pyfile.startswith("lesson") and pyfile.endswith(".py"): continue
        path = os.path.join(os.path.dirname(__file__), pyfile)
        mod = import_module(path)
        docs.appen(mod.__doc__)

    return "\n".join(docs)
