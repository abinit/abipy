---
name: abipy-plotting
description: Add, extend, and test AbiPy Matplotlib and Plotly visualization APIs, including plot methods, caller-supplied axes, figure decorators, Robot convergence plots, yield_figs integration, and headless examples. Use for scientific plotting code; not for Panel application layout.
---

# Implement AbiPy plotting

Keep data preparation separable from rendering, follow the plotting conventions of the nearest domain object, and
return figure objects that callers can compose, save, test, or display. Use `abipy-reference-data` when a real file
is needed and `abipy-run-tests` for focused validation.

## Match the established API

Matplotlib methods normally:

- use a `plot_*` name;
- accept `ax=None` for a single axes or an appropriate axes-array argument for composed plots;
- accept domain options plus `**kwargs`;
- use `@abipy.tools.plotting.add_fig_kwargs`;
- return a Matplotlib `Figure`, or `None` only when the documented data/mode cannot be plotted.

Plotly methods normally use a `plotly_*` name, `get_fig_plotly`/`get_figs_plotly`, and
`@add_plotly_fig_kwargs`. Do not implement Plotly by mechanically converting Matplotlib when a native interactive
figure needs different traces, hover data, or subplot behavior.

Inspect a mature nearby object before adding a new signature. Preserve established parameter names such as
`fontsize`, `xlims`, `ylims`, `sortby`, `hue`, and `with_legend` when their semantics match.

## Matplotlib figure ownership

Use AbiPy helpers instead of calling `plt.subplots` directly when the method supports composition:

```python
ax, fig, plt = get_ax_fig_plt(ax=ax)
ax_array, fig, plt = get_axarray_fig_plt(
    ax_array, nrows=nrows, ncols=ncols, squeeze=False,
)
```

When the caller supplies axes, draw into them and return their existing figure. Do not clear unrelated artists,
change global Matplotlib state, call `plt.show`, or close the figure inside the core function.

Let `@add_fig_kwargs` handle `show`, `savefig`, `title`, `size_kwargs`, `tight_layout`, grids, annotations, closure,
and optional Matplotlib-to-Plotly conversion. Nested calls to another decorated plotting method must pass
`show=False` so a composite plot does not display intermediate figures.

For arrays of axes, normalize to a predictable shape, fill axes deterministically, and hide unused panels. Do not
depend on Matplotlib's squeezed scalar/one-dimensional return shape unless the helper call explicitly requests it.

## Scientific presentation

Labels must state physical quantities and units. Apply conversions before plotting and keep the converted data
available for tests where practical. Be explicit about energy references, normalization, spin/channel selection,
band or atom indexing, and whether coordinates are reduced or Cartesian.

Choose defaults that reveal the physical signal without silently discarding data. Filtering non-finite values,
clipping ranges, wrapping angles, smoothing, interpolation, or logarithmic scaling must be documented and should be
optional when it changes interpretation. Do not hide outliers merely to improve appearance.

Use `set_grid_legend`, `set_axlims`, `set_ax_xylabels`, shared-axis helpers, and existing color/marker conventions.
Honor caller style kwargs and avoid hard-coding a global seaborn style, backend, or rcParams in a library method.

Legends and colorbars should identify all encodings without obscuring data. For convergence plots, make the sort
variable and hue grouping explicit and deterministic. Preserve input order when no sort is requested.

## Separate computation from rendering

Move substantial extraction, aggregation, fitting, or statistical logic into a method that returns arrays,
dataframes, or a small result object. The plotting method should orchestrate these results and axes. This makes the
scientific transformation testable without inspecting pixels and lets Matplotlib and Plotly share the same data.

Do not reread a large NetCDF variable once per subplot when it can be loaded or reduced once. Conversely, avoid
eagerly materializing an entire large dataset when the plot needs a small slice.

## File objects, Robots, and automatic figures

Put a single-file scientific view on the file/domain object. Put cross-file comparison and convergence plots on its
`Robot`. Robot methods should use the existing label/path mapping, `sortby` and `hue` conventions, and dataframe or
parameter helpers rather than reopen files independently.

Add a plot to `yield_figs` only when it is useful as a default automatic view for `abiopen -e`. Yield figures with
display disabled and avoid expensive combinatorial plots. An object may expose useful plotting methods without
including all of them in `yield_figs`.

Add notebook cells or Plotly/Panel integration only when the corresponding object already supports that layer. A
Matplotlib plot does not automatically require a new notebook or dashboard API.

## Plotly behavior

Construct traces with meaningful names, legend groups, axes labels, hover text, and units. Respect a caller-provided
figure and the requested subplot row/column. Use the global Plotly show setting and `@add_plotly_fig_kwargs` rather
than calling `fig.show()` unconditionally.

Keep Matplotlib and Plotly numerically consistent when they represent the same method, but allow presentation to
differ where native interaction benefits from it. Test shared data generation separately from backend-specific
trace construction.

## Examples and documentation

Add an example under `abipy/examples/plot` when the plot introduces a user-facing workflow not already demonstrated.
Use bundled data through `abipy.data`, keep runtime modest, and write the example so it works with a noninteractive
backend. Do not commit generated PDFs, PNGs, notebooks, or text exports from running the example unless they are
intentional maintained artifacts.

Examples may call plotting methods with their normal display defaults for gallery use. Automated checks must select
the `Agg` backend or pass `show=False`.

## Tests

Place tests beside the domain object. Guard optional backends with `AbipyTest.has_matplotlib()` or
`has_plotly()` as established in the suite.

For Matplotlib, call the method with `show=False` and assert the returned figure plus meaningful invariants such as:

- expected axes count and visibility;
- line, image, collection, or trace count;
- axis labels and units;
- plotted numerical data for a representative slice;
- support for a caller-provided axes;
- important option branches and documented `None` cases.

Close figures created by tests when a loop or large suite could accumulate them. Do not use interactive windows,
browsers, clipboard integration, or `show=True` in automated tests. Prefer semantic assertions to image snapshots;
use an image comparison only for a genuine layout/rendering regression that structural checks cannot capture.

For Plotly, assert figure type, trace/subplot count, names, axes labels, and representative trace data with display
disabled. Do not require Chart Studio or a browser.

Report the new public methods, backend coverage, data transformations, default `yield_figs` impact, example added,
and the headless test commands used.
