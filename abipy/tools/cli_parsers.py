"""
"""
import argparse
import os

from pprint import pformat


def pn_serve_parser(**kwargs):
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
                  help="Specify template for panel dasboard." +
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
                help="Maximum message size in Mb allowed by Bokeh and Tornado. Default: 150")

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


def customize_mpl(options):

    if options.mpl_backend is not None:
        # Set matplotlib backend
        import matplotlib
        matplotlib.use(options.mpl_backend)

    if options.seaborn:
        # Use seaborn settings.
        import seaborn as sns
        sns.set(context=options.seaborn, style='darkgrid', palette='deep',
                font='sans-serif', font_scale=1, color_codes=False, rc=None)
