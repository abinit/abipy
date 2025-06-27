"""
Tools and helper functions to build the command line interface of the AbiPy scripts.
"""
from __future__ import annotations

import argparse
import sys
import os

from functools import wraps
from pprint import pformat


def user_wants_to_abort() -> bool:
    """Interactive prompt, return False if user entered `n` or `no`."""
    try:
        answer = input("\nDo you want to continue [Y/n]")
    except EOFError:
        return False

    return answer.lower().strip() in ["n", "no"]


def set_loglevel(loglevel: str) -> None:
    # loglevel is bound to the string value obtained from the command line argument.
    # Convert to upper case to allow the user to specify --loglevel=DEBUG or --loglevel=debug
    import logging
    numeric_level = getattr(logging, loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % loglevel)
    logging.basicConfig(level=numeric_level)


def pn_serve_parser(**kwargs) -> argparse.ArgumentParser:
    """
    Parent parser implementing cli options for panel.serve
    """
    p = argparse.ArgumentParser(add_help=False)

    p.add_argument("--port", default=0, type=int, help="Port to listen on.")
    p.add_argument("--address", default=None,
                              help="The address the server should listen on for HTTP requests.")
    #p.add_argument("--show", default=True, action="store_true", help="Open app in web browser")
    p.add_argument("--num_procs", default=1, type=int,
                              help="Number of worker processes for the app. Defaults to 1")
    p.add_argument('--panel-template', "-pnt", default="FastList",
                  help="Specify template for panel dashboard." +
                       "Possible values are: FastList, FastGrid, Golden, Bootstrap, Material, React, Vanilla." +
                       "Default: FastList")
    p.add_argument('--has-remote-server', default=False, action="store_true",
                  help="True if we are running on the ABINIT server. " +
                       "This flag activates limitations on what the user can do." +
                       "Default: False")
    p.add_argument("--websocket-origin", default=None, type=str,
            help="Public hostnames which may connect to the Bokeh websocket.\n Syntax: " +
                  "HOST[:PORT] or *. Default: None")
    p.add_argument('--max_size_mb', default=150, type=int,
                help="Maximum message size in MB allowed by Bokeh and Tornado. Default: 150")

    p.add_argument('--no-browser', action='store_true', default=False,
                   help=("Start the jupyter server to serve the notebook "
                         "but don't open the notebook in the browser.\n"
                         "Use this option to connect remotely from localhost to the machine running the kernel"))

    return p


def get_pn_serve_kwargs(options) -> dict:
    """
    Return dict with the arguments to be passed to pn.serve.
    """
    import abipy.panels as mod
    assets_path = os.path.join(os.path.dirname(mod.__file__), "assets")

    serve_kwargs = dict(
        address=options.address,
        port=options.port,
        #dev=True,
        #start=True,
        #show=options.show,
        show=not options.no_browser,
        debug=options.verbose > 0,
        #title=app_title,
        num_procs=options.num_procs,
        static_dirs={"/assets": assets_path},
        websocket_origin=options.websocket_origin,
        #
        # Increase the maximum websocket message size allowed by Bokeh
        # https://panel.holoviz.org/reference/widgets/FileInput.html
        websocket_max_message_size=options.max_size_mb * 1024**2,
        # Increase the maximum buffer size allowed by Tornado
        http_server_kwargs={'max_buffer_size': options.max_size_mb * 1024**2},
    )

    if getattr(options, "verbose"):
        print("Calling pn.serve with serve_kwargs:\n", pformat(serve_kwargs), "\n")

    if options.no_browser:
        print("""
Use:

    ssh -N -f -L localhost:{port}:localhost:{port} username@your_remote_cluster

for port forwarding.
""")


    return serve_kwargs


def customize_mpl(options) -> None:

    if options.mpl_backend is not None:
        # Set matplotlib backend
        import matplotlib
        matplotlib.use(options.mpl_backend)

    if options.seaborn:
        # Use seaborn settings.
        import seaborn as sns
        sns.set(context=options.seaborn, style='darkgrid', palette='deep',
                font='sans-serif', font_scale=1, color_codes=False, rc=None)


def add_expose_options_to_parser(parser, with_mpl_options=True) -> None:
    """
    Add Expose options to the parser.
    """

    parser.add_argument('-e', '--expose', action='store_true', default=False,
        help="Open file and generate matplotlib figures automatically by calling expose method.")
    parser.add_argument("-s", "--slide-mode", default=False, action="store_true",
        help="Iterate over figures. Expose all figures at once if not given on the CLI.")
    parser.add_argument("-t", "--slide-timeout", type=int, default=None,
        help="Close figure after slide-timeout seconds (only if slide-mode). Block if not specified.")
    parser.add_argument("-ew", "--expose-web", default=False, action="store_true",
            help='Generate matplotlib plots in $BROWSER instead of X-server. WARNING: Not all the features are supported.')
    parser.add_argument("-ply", "--plotly", default=False, action="store_true",
            help='Generate plotly plots in $BROWSER instead of matplotlib. WARNING: Not all the features are supported.')
    parser.add_argument("-cs", "--chart-studio", default=False, action="store_true",
            help="Push figure to plotly chart studio ." +
                 "Requires --plotly option and user account at https://chart-studio.plotly.com.")

    if with_mpl_options:
        parser.add_argument('-sns', "--seaborn", const="paper", default=None, action='store', nargs='?', type=str,
            help='Use seaborn settings. Accept value defining context in ("paper", "notebook", "talk", "poster"). Default: paper')
        parser.add_argument('-mpl', "--mpl-backend", default=None,
            help=("Set matplotlib interactive backend. "
                  "Possible values: GTKAgg, GTK3Agg, GTK, GTKCairo, GTK3Cairo, WXAgg, WX, TkAgg, Qt4Agg, Qt5Agg, macosx."
                  "See also: https://matplotlib.org/faq/usage_faq.html#what-is-a-backend."))


class EnumAction(argparse.Action):
    """
    Argparse action for handling Enums

    Usage:

        class Do(enum.Enum):
            Foo = "foo"
            Bar = "bar"

        parser = argparse.ArgumentParser()
        parser.add_argument('do', type=Do, action=EnumAction)

    Taken from https://stackoverflow.com/questions/43968006/support-for-enum-arguments-in-argparse
    """
    def __init__(self, **kwargs):
        # Pop off the type value
        enum_type = kwargs.pop("type", None)

        # Ensure an Enum subclass is provided
        if enum_type is None:
            raise ValueError("type must be assigned an Enum when using EnumAction")
        if not issubclass(enum_type, enum.Enum):
            raise TypeError("type must be an Enum when using EnumAction")

        # Generate choices from the Enum
        kwargs.setdefault("choices", tuple(e.value for e in enum_type))

        super().__init__(**kwargs)

        self._enum = enum_type

    def __call__(self, parser, namespace, values, option_string=None):
        # Convert value back into an Enum
        value = self._enum(values)
        setattr(namespace, self.dest, value)


def fix_omp_num_threads() -> int:
    """
    Set OMP_NUM_THREADS to 1 if env var is not defined. Return num_threads.
    """
    num_threads = os.getenv("OMP_NUM_THREADS", default=None)
    if num_threads is None:
        num_threads = 1
        os.environ["OMP_NUM_THREADS"] = str(num_threads)

    return num_threads


def range_from_str(string: str) -> range:
    """
    Convert string into a range object.
    """
    if string is None: return None

    tokens = string.split(":")
    start, stop, step = 0, None, 1
    if len(tokens) == 1:
        stop = int(tokens[0])
    elif len(tokens) == 2:
        start, stop = map(int, tokens)
    elif len(tokens) == 3:
        start, stop, step = map(int, tokens)
    else:
        raise ValueError(f"Cannot interpret {string=} as range object.")

    return range(start, stop, step)


def prof_main(main):
    """
    Decorator for profiling main programs.

    Profiling is activated by prepending the command line options
    supported by the original main program with the keyword `prof`.

    Examples:

            $ script.py arg --foo=1

        becomes

            $ script.py prof arg --foo=1

        The decorated main accepts two new arguments:

            prof_file: Name of the output file with profiling data
                If not given, a temporary file is created.
            sortby: Profiling data are sorted according to this value.
                default is "time". See sort_stats.
    """

    @wraps(main)
    def wrapper(*args, **kwargs):
        arg0 = None
        try:
            do_prof = sys.argv[1] in ("prof", "tuna", "snakeviz")
            if do_prof:
                print("Entering profiling mode...")
                arg0 = sys.argv.pop(1)
        except Exception:
            do_prof = False

        if not do_prof:
            # No profiling. Run script as usual.
            sys.exit(main())

        import cProfile
        import pstats
        import tempfile
        prof_file = kwargs.get("prof_file", None)
        if prof_file is None:
            _, prof_file = tempfile.mkstemp()
            print(f"Profiling data stored in: {prof_file}")

        sortby = kwargs.get("sortby", "time")
        cProfile.runctx("main()", globals(), locals(), prof_file)
        s = pstats.Stats(prof_file)
        s.strip_dirs().sort_stats(sortby).print_stats()

        if arg0 in ("tuna", "snakeviz"):
            cmd = f"{arg0} {prof_file}"
            print(f"Executing {cmd=}")
            os.system(cmd)

        if "retval" not in kwargs:
            sys.exit(0)
        else:
            return kwargs["retval"]

    return wrapper
