"""
Deployment file to facilitate AbiPy releases.
Use invoke --list to get list of tasks
"""

import os

from invoke import task
from monty.os import cd

#from abipy.core.release import __version__ as CURRENT_VER
#NEW_VER = datetime.datetime.today().strftime("%Y.%-m.%-d")

ABIPY_ROOTDIR = os.path.dirname(__file__)
DOCS_DIR = os.path.join(ABIPY_ROOTDIR, "docs")

@task
def make_doc(ctx):
    with cd(DOCS_DIR):
        ctx.run("make clean")
        ctx.run("make", env=dict(READTHEDOCS="1"), pty=True)
        open_doc(ctx)


@task
def push_doc(ctx):
    make_doc(ctx)
    with cd(DOCS_DIR):
        ctx.run("./ghp_import.py _build/html/ -n -p")


@task
def open_doc(ctx):
    import webbrowser
    webbrowser.open_new_tab("file://" + os.path.join(ABIPY_ROOTDIR, "docs/_build/html/index.html"))


@task
def twine(ctx):
    with cd(ABIPY_ROOTDIR):
        ctx.run("rm dist/*.*", warn=True)
        ctx.run("python setup.py register sdist bdist_wheel")
        ctx.run("twine upload dist/*")


@task
def pytest(ctx):
    pytest_cmd = r"""\
pytest -n 2 --cov-config=.coveragerc --cov=abipy -v --doctest-modules abipy \
    --ignore=abipy/integration_tests --ignore=abipy/data/refs --ignore=abipy/scripts/ \
    --ignore=abipy/examples/plot --ignore=abipy/examples/flows --ignore=abipy/gui
"""
    with cd(ABIPY_ROOTDIR):
        ctx.run(pytest_cmd, pty=True)


@task
def plots(ctx):
    with cd(os.path.join(ABIPY_ROOTDIR, "abipy", "examples")):
        ctx.run("_runplots.py", pty=True)

@task
def flows(ctx):
    with cd(os.path.join(ABIPY_ROOTDIR, "abipy", "examples")):
        ctx.run("_runflows.py", pty=True)

#@task
#def move_to_master(ctx):
#    ctx.run("git tag -a v%s -m \"v%s release\"" % (NEW_VER, NEW_VER))
#    ctx.run("git push --tags")
#    ctx.run("git checkout master")
#    ctx.run("git pull")
#    ctx.run("git merge develop")
#    ctx.run("git push")
#    ctx.run("git checkout develop")


#@task
#def update_changelog(ctx):
#
#    output = subprocess.check_output(["git", "log", "--pretty=format:%s",
#                                      "v%s..HEAD" % CURRENT_VER])
#    lines = ["* " + l for l in output.decode("utf-8").strip().split("\n")]
#    with open("CHANGES.rst") as f:
#        contents = f.read()
#    l = "=========="
#    toks = contents.split(l)
#    head = "\n\nv%s\n" % NEW_VER + "-" * (len(NEW_VER) + 1) + "\n"
#    toks.insert(-1, head + "\n".join(lines))
#    with open("CHANGES.rst", "w") as f:
#        f.write(toks[0] + l + "".join(toks[1:]))


#@task
#def release(ctx, run_tests=True):
#    ctx.run("rm -r dist build abipy.egg-info", warn=True)
#    set_ver(ctx)
#    if run_tests: pytest(ctx)
#    publish(ctx)
#    log_ver(ctx)
#    update_doc(ctx)
#    merge_stable(ctx)
#    release_github(ctx)
