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
def style(ctx):
    with cd(ABIPY_ROOTDIR):
        ctx.run("pycodestyle 2>&1 | tee style.log", pty=True)
        ctx.run("flake8 --count --show-source --statistics | tee -a style.log", pty=True)
        #ctx.run("pydocstyle abipy | tee -a style.log", pty=True)


@task
def plots(ctx):
    with cd(os.path.join(ABIPY_ROOTDIR, "abipy", "examples")):
        ctx.run("_runplots.py", pty=True)


@task
def flows(ctx):
    with cd(os.path.join(ABIPY_ROOTDIR, "abipy", "examples")):
        ctx.run("_runflows.py", pty=True)


@task
def pygrep(ctx, pattern):
    """
    Grep for `pattern` in all py files contained in
    """
    # grep -r -i --include \*.h
    # Syntax notes:
    #    -r - search recursively
    #    -i - case-insensitive search
    #    --include=\*.${file_extension} - search files that match the extension(s) or file pattern only
    with cd(os.path.join(ABIPY_ROOTDIR, "abipy",)):
        cmd = 'grep -r -i --color --include "*.py" "%s" .' % pattern
        print("Executing:", cmd)
        ctx.run(cmd, pty=True)


@task
def update_vars(ctx, abinit_repo_path):
    abinit_repo_path = os.path.abspath(abinit_repo_path)

    dir_with_pyfiles = os.path.join(ABIPY_ROOTDIR, "abipy", "abio", "abivar_database")

    local_files = [f for f in os.listdir(dir_with_pyfiles) if f.startswith("variables_")]
    for local_file in local_files:
        # "vimdiff $ABINIT_REPOPATH/abimkdocs/variables_abinit.py variables_abinit.py
        source = os.path.join(abinit_repo_path, "abimkdocs", local_file)
        cmd = f"vimdiff {source} {local_file}"
        print(f"Executing: {cmd}")
        os.system(cmd)


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


#@task
#def watchdog(ctx, jobs="auto", sleep_time=5):
#    """
#    Start watchdog service to watch F90 files and execute `make` when changes are detected.
#    """
#    from monty.termcolor import cprint
#    cprint("Starting watchdog service to watch F90 files and execute `make` when changes are detected", "green")
#    cprint("Enter <CTRL + C> in the terminal to kill the service.", "green")
#
#    cprint(f"Start watching py files with sleep_time {sleep_time} s ....", "green")
#    top = find_top_build_tree(".", with_abinit=True)
#    jobs = max(1, number_of_cpus() // 2) if jobs == "auto" else int(jobs)
#
#    # http://thepythoncorner.com/dev/how-to-create-a-watchdog-in-python-to-look-for-filesystem-changes/
#    # https://stackoverflow.com/questions/19991033/generating-multiple-observers-with-python-watchdog
#    import time
#    from watchdog.observers import Observer
#    from watchdog.events import PatternMatchingEventHandler
#    event_handler = PatternMatchingEventHandler(patterns="*.py", ignore_patterns="",
#                                                   ignore_directories=False, case_sensitive=True)
#
#    def on_created(event):
#        print(f"hey, {event.src_path} has been created!")
#
#    def on_deleted(event):
#        print(f"what the f**k! Someone deleted {event.src_path}!")
#
#    def on_modified(event):
#        print(f"hey buddy, {event.src_path} has been modified")
#        cmd = "abicheck.py"
#        cprint("Executing: %s" % cmd, "yellow")
#        with cd(top):
#            try:
#                result = ctx.run(cmd, pty=True)
#                if result.ok:
#                    cprint("Command completed successfully", "green")
#                    cprint("Watching for changes ...", "green")
#            except Exception:
#                cprint(f"Command returned non-zero exit status", "red")
#                cprint(f"Keep on watching for changes hoping you get it right ...", "red")
#
#    def on_moved(event):
#        print(f"ok ok ok, someone moved {event.src_path} to {event.dest_path}")
#
#    event_handler.on_created = on_created
#    event_handler.on_deleted = on_deleted
#    event_handler.on_modified = on_modified
#    event_handler.on_moved = on_moved
#
#    path = ABIPY_ROOTDIR
#
#    observer = Observer()
#    observer.schedule(event_handler, path, recursive=True)
#    observer.start()
#
#    try:
#        while True:
#            time.sleep(sleep_time)
#    except KeyboardInterrupt:
#        observer.stop()
#        observer.join()
