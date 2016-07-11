#!/usr/bin/env python
"""Deployment file to facilitate releases of abipy."""
from __future__ import division, print_function, unicode_literals

from fabric.api import local, lcd


def unittests():
    """Run the unit tests."""
    local("nosetests")


def commit():
    """Commit changes."""
    local("git commit -a")


def push():
    """Push to github."""
    local("git push")


def deploy():
    """Unit tests + git commit + git push."""
    unittests()
    commit()
    push()


def publish():
    local("python setup.py release")


#def setver():
#    local("sed s/version=.*,/version=\\\"{}\\\",/ setup.py > newsetup"
#          .format(ver))
#    local("mv newsetup setup.py")


def make_doc():
    """Build docs."""

    # Generate html pages from the notebooks.
    with lcd("abipy/examples/notebooks"):
       local("make html") # "ipython nbconvert --to html *.ipynb"
       #local("mv *.html ../../docs/_static")

    # Run sphynx.
    with lcd("docs"):
        local("make clean")
        local("make html")
        #local("cp _static/* _build/html/_static")


def makedoc_check():
    """Check the docstrings and the links."""
    with lcd("docs"):
        #local("make doctest")
        local("make linkcheck")


#def update_doc():
#    with lcd("docs/_build/html/"):
#        local("git pull")
#    make_doc()
#    with lcd("docs/_build/html/"):
#        local("git add .")
#        local("git commit -a -m \"Update dev docs\"")
#        local("git push origin gh-pages")


#def update_coverage():
#    with lcd("docs/_build/html/"):
#        local("git pull")
#    local("nosetests --config=nose.cfg --cover-html --cover-html-dir=docs/_build/html/coverage")
#    update_doc()


#def release():
#    setver()
#    test()
#    publish()
#    log_ver()
#    update_doc()


#def merge_stable():
#    local("git commit -a -m \"v%s release\"" % ver)
#    local("git push")
#    local("git checkout stable")
#    local("git pull")
#    local("git merge master")
#    local("git push")
#    local("git checkout master")

#def release_github():
#    with open("CHANGES.rst") as f:
#        contents = f.read()
#    toks = re.split("\-+", contents)
#    desc = toks[1].strip()
#    toks = desc.split("\n")
#    desc = "\n".join(toks[:-1]).strip()
#    payload = {
#        "tag_name": "v" + ver,
#        "target_commitish": "master",
#        "name": "v" + ver,
#        "body": desc,
#        "draft": False,
#        "prerelease": False
#    }
#    response = requests.post(
#        "https://api.github.com/repos/materialsproject/pymatgen/releases",
#        data=json.dumps(payload),
#        headers={"Authorization": "token " + os.environ["GITHUB_RELEASES_TOKEN"]})
#    print response.text
#
#
#def update_changelog():
#    output = subprocess.check_output(["git", "log", "--pretty=format:%s",
#                                      "v%s..HEAD" % ver])
#    lines = ["* " + l for l in output.strip().split("\n")]
#    with open("CHANGES.rst") as f:
#        contents = f.read()
#    l = "=========="
#    toks = contents.split(l)
#    toks.insert(-1, "\n\nvXXXX\n--------\n" + "\n".join(lines))
#    with open("CHANGES.rst", "w") as f:
#        f.write(toks[0] + l + "".join(toks[1:]))
#
#
#def log_ver():
#    filepath = os.path.join(os.environ["HOME"], "Dropbox", "Public",
#                            "pymatgen", ver)
#    with open(filepath, "w") as f:
#        f.write("Release")
#
#
#def release(skip_test=False):
#    setver()
#    if not skip_test:
#        local("nosetests")
#    publish()
#    log_ver()
#    update_doc()
#    merge_stable()
#    release_github()
#
#
#def open_doc():
#    pth = os.path.abspath("docs/_build/html/index.html")
#    webbrowser.open("file://" + pth)
