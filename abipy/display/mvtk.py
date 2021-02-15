# coding: utf-8
"""
mayavi_ toolkit.

WARNING: This code is still under development.
"""
import itertools
import numpy as np

DEFAULT_FIGURE_KWARGS = dict(size=(1024, 768), bgcolor=(1, 1, 1), fgcolor=(0, 0, 0))


def get_fig_mlab(figure=None, **kwargs):  # pragma: no cover
    try:
        from mayavi import mlab
    except ImportError as exc:
        print("mayavi is not installed. Use `conda install mayavi` or `pip install mayavi`")
        raise exc

    # To use the full envisage application
    #mlab.options.backend = "envisage"
    #mlab.options.backend = "test"
    #mlab.options.offscreen = True

    if figure is None:
        # Add defaults
        for k, v in DEFAULT_FIGURE_KWARGS.items():
            if k not in kwargs: kwargs[k] = v
        figure = mlab.figure(**kwargs)

    return figure, mlab


def plot_wigner_seitz(lattice, figure=None, **kwargs):  # pragma: no cover
    """
    Adds the skeleton of the Wigner-Seitz cell of the lattice to a mayavi_ figure

    Args:
        lattice: Reciprocal-space |Lattice| object
        figure: mayavi figure, None to plot on the curretn figure
        kwargs: kwargs passed to the mayavi function ``plot3d``. Color defaults to black
            and line_width to 1.

    Returns: mayavi figure
    """
    figure, mlab = get_fig_mlab(figure=figure)

    if "color" not in kwargs:
        kwargs["color"] = (0, 0, 0)
    if "line_width" not in kwargs:
        kwargs["line_width"] = 1
    if "tube_radius" not in kwargs:
        kwargs["tube_radius"] = None

    bz = lattice.get_wigner_seitz_cell()
    for iface in range(len(bz)):
        for line in itertools.combinations(bz[iface], 2):
            for jface in range(len(bz)):
                if iface < jface and any(np.all(line[0] == x) for x in bz[jface])\
                        and any(np.all(line[1] == x) for x in bz[jface]):
                    #do_plot = True
                    #if in_unit_cell:
                    #    kred0 = lattice.get_fractional_coords(line[0])
                    #    kred1 = lattice.get_fractional_coords(line[1])
                    #    do_plot = np.alltrue((kred0 >= 0) & (kred0 <= 0.5) &
                    #                         (kred1 >= 0) & (kred1 <= 0.5))
                    #    print(kred0, kred1, do_plot)
                    #if not do_plot: continue
                    mlab.plot3d(*zip(line[0], line[1]), figure=figure, **kwargs)

    return figure


def plot_unit_cell(lattice, figure=None, **kwargs):  # pragma: no cover
    """
    Adds the unit cell of the lattice to a mayavi_ figure.

    Args:
        lattice: Lattice object
        figure: mayavi figure, None to plot on the curretn figure
        kwargs: kwargs passed to the mayavi function ``plot3d``. Color defaults to black
            and line_width to 1.

    Returns: mayavi figure
    """
    figure, mlab = get_fig_mlab(figure=figure)

    if "color" not in kwargs:
        kwargs["color"] = (0, 0, 0)
    if "line_width" not in kwargs:
        kwargs["line_width"] = 1
    if "tube_radius" not in kwargs:
        kwargs["tube_radius"] = None

    v = 8 * [None]
    v[0] = lattice.get_cartesian_coords([0.0, 0.0, 0.0])
    v[1] = lattice.get_cartesian_coords([1.0, 0.0, 0.0])
    v[2] = lattice.get_cartesian_coords([1.0, 1.0, 0.0])
    v[3] = lattice.get_cartesian_coords([0.0, 1.0, 0.0])
    v[4] = lattice.get_cartesian_coords([0.0, 1.0, 1.0])
    v[5] = lattice.get_cartesian_coords([1.0, 1.0, 1.0])
    v[6] = lattice.get_cartesian_coords([1.0, 0.0, 1.0])
    v[7] = lattice.get_cartesian_coords([0.0, 0.0, 1.0])

    for i, j in ((0, 1), (1, 2), (2, 3), (0, 3), (3, 4), (4, 5), (5, 6),
                 (6, 7), (7, 4), (0, 7), (1, 6), (2, 5), (3, 4)):
        mlab.plot3d(*zip(v[i], v[j]), figure=figure, **kwargs)

    #mlab.xlabel("x-axis")
    #mlab.ylabel("y-axis")
    #mlab.zlabel("z-axis")

    return figure


def plot_lattice_vectors(lattice, figure=None, **kwargs): # pragma: no cover
    """
    Adds the basis vectors of the lattice provided to a mayavi_ figure.

    Args:
        lattice: |Lattice| object.
        figure: mayavi figure, None if a new figure should be created.
        kwargs: kwargs passed to the mayavi function ``plot3d``. Color defaults to black
            and line_width to 1.

    Returns: mayavi figure
    """
    figure, mlab = get_fig_mlab(figure=figure)

    if "color" not in kwargs:
        kwargs["color"] = (0, 0, 0)
    if "line_width" not in kwargs:
        kwargs["line_width"] = 1
    if "tube_radius" not in kwargs:
        kwargs["tube_radius"] = None

    vertex1 = lattice.get_cartesian_coords([0.0, 0.0, 0.0])
    vertex2 = lattice.get_cartesian_coords([1.0, 0.0, 0.0])
    mlab.plot3d(*zip(vertex1, vertex2), figure=figure, **kwargs)
    vertex2 = lattice.get_cartesian_coords([0.0, 1.0, 0.0])
    mlab.plot3d(*zip(vertex1, vertex2), figure=figure, **kwargs)
    vertex2 = lattice.get_cartesian_coords([0.0, 0.0, 1.0])
    mlab.plot3d(*zip(vertex1, vertex2), figure=figure, **kwargs)

    return figure


def plot_structure(structure, frac_coords=False, to_unit_cell=False, style="points+labels",
                   unit_cell_color=(0, 0, 0), color_scheme="VESTA", figure=None, show=False, **kwargs):  # pragma: no cover
    """
    Plot structure with mayavi.

    Args:
        structure: |Structure| object
        frac_coords:
        to_unit_cell: True if sites should be wrapped into the first unit cell.
        style: "points+labels" to show atoms sites with labels.
        unit_cell_color:
        color_scheme: color scheme for atom types. Allowed values in ("Jmol", "VESTA")
        figure:
        kwargs:

    Returns: mayavi figure
    """
    figure, mlab = get_fig_mlab(figure=figure)

    #if not frac_coords:
    plot_unit_cell(structure.lattice, color=unit_cell_color, figure=figure)
    from pymatgen.analysis.molecule_structure_comparator import CovalentRadius
    from pymatgen.vis.structure_vtk import EL_COLORS

    for site in structure:
        symbol = site.specie.symbol
        color = tuple(i / 255 for i in EL_COLORS[color_scheme][symbol])
        radius = CovalentRadius.radius[symbol]
        if to_unit_cell and hasattr(site, "to_unit_cell"): site = site.to_unit_cell
        x, y, z = site.frac_coords if frac_coords else site.coords

        if "points" in style:
            mlab.points3d(x, y, z, figure=figure, scale_factor=radius,
                          resolution=20, color=color, scale_mode='none', **kwargs)
        if "labels" in style:
            mlab.text3d(x, y, z, symbol, figure=figure, color=(0, 0, 0), scale=0.2)

    if show: mlab.show()
    return figure


def plot_labels(labels, lattice=None, coords_are_cartesian=False, figure=None, **kwargs):  # pragma: no cover
    """
    Adds labels to a mayavi_ figure.

    Args:
        labels: dict containing the label as a key and the coordinates as value.
        lattice: |Lattice| object used to convert from reciprocal to cartesian coordinates
        coords_are_cartesian: Set to True if you are providing.
            coordinates in cartesian coordinates. Defaults to False.
            Requires lattice if False.
        figure: mayavi figure, None to plot on the curretn figure
        kwargs: kwargs passed to the mayavi function `text3d`. Color defaults to blue and size to 25.

    Returns: mayavi figure
    """
    figure, mlab = get_fig_mlab(figure=figure)

    #if "color" not in kwargs:
    #    kwargs["color"] = "b"
    #if "size" not in kwargs:
    #    kwargs["size"] = 25
    #if "width" not in kwargs:
    #    kwargs["width"] = 0.8
    if "scale" not in kwargs:
        kwargs["scale"] = 0.1

    for k, coords in labels.items():
        label = k
        if k.startswith("\\") or k.find("_") != -1:
            label = "$" + k + "$"
        off = 0.01
        if coords_are_cartesian:
            coords = np.array(coords)
        else:
            if lattice is None:
                raise ValueError("coords_are_cartesian False requires the lattice")
            coords = lattice.get_cartesian_coords(coords)
        x, y, z = coords + off
        mlab.text3d(x, y, z, label, figure=figure, **kwargs)

    return figure


class MayaviFieldAnimator(object): # pragma: no cover

    def __init__(self, filepaths):
        self.filepaths = filepaths
        self.num_files = len(filepaths)

    def volume_animate(self):
        from abipy import abilab
        with abilab.abiopen(self.filepaths[0]) as nc:
            nsppol, nspden, nspinor = nc.nsppol, nc.nspden, nc.nspinor
            structure = nc.structure
            datar = nc.field.datar
            # [nspden, nx, ny, nz] array
            nx, ny, nz = datar.shape[1:]
            s = datar[0]
            print(s.dtype, s.shape)

        #cart_coords = np.empty((nx*ny*nz, 3))
        #cnt = 0
        #for i in range(nx):
        #    for j in range(ny):
        #        for k in range(nz):
        #            cart_coords[ctn, :] = (i/nx, j/ny, k/nz)
        #            cnt += 1
        #cart_coords = structure.lattice.get_cartesian_coords(cart_coords)
        # We reorder the points, scalars and vectors so this is as per VTK's
        # requirement of x first, y next and z last.
        #pts = pts.transpose(2, 1, 0, 3).copy()
        #pts.shape = pts.size / 3, 3
        #scalars = scalars.T.copy()
        #vectors = vectors.transpose(2, 1, 0, 3).copy()
        #vectors.shape = vectors.size / 3, 3

        #from tvtk.api import tvtk
        #sgrid = tvtk.StructuredGrid(dimensions=(dims[1], dims[0], dims[2]))
        #sgrid.points = pts
        #s = random.random((dims[0]*dims[1]*dims[2]))
        #sgrid.point_data.scalars = ravel(s.copy())
        #sgrid.point_data.scalars.name = 'scalars'

        figure, mlab = get_fig_mlab(figure=None)
        source = mlab.pipeline.scalar_field(s)
        data_min, data_max = s.min(), s.max()
        print(data_min, data_max)
        #mlab.pipeline.volume(source)
        #                     #vmin=data_min + 0.65 * (data_max - data_min),
        #                     #vmax=data_min + 0.9 * (data_max - data_min))
        #mlab.pipeline.iso_surface(source)
        mlab.pipeline.image_plane_widget(source, plane_orientation='x_axes', slice_index=0)
        mlab.pipeline.image_plane_widget(source, plane_orientation='y_axes', slice_index=0)
        mlab.pipeline.image_plane_widget(source, plane_orientation='z_axes', slice_index=0)

        @mlab.show
        @mlab.animate(delay=1000, ui=True)
        def anim():
            """Animate."""
            t = 1
            while True:
                #vmin, vmax = .1 * np.max(data[t]), .2 * np.max(data[t])
                #print 'animation t = ',tax[t],', max = ',np.max(data[t])
                with abilab.abiopen(self.filepaths[t]) as nc:
                    print("Animation step", t, "from file:", self.filepaths[t])
                    #nsppol, nspden, nspinor = nc.nsppol, nc.nspden, nc.nspinor
                    datar = nc.field.datar
                    # [nspden, nx, ny, nz] array
                    #nx, ny, nz = datar.shape[1:]
                    scalars = datar[0]

                #data_min, data_max = scalars.min(), scalars.max(),
                #mlab.pipeline.volume(source, vmin=data_min + 0.65 * (data_max - data_min),
                #                     vmax=data_min + 0.9 * (data_max - data_min))
                source.mlab_source.scalars = scalars

                t = (t + 1) % self.num_files
                yield

        anim()
