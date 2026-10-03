---
name: abipy-panels
description: Build, extend, and test AbiPy Panel dashboards for files, structures, Robots, and FlowTK objects using the repository's Parameterized classes, lazy callbacks, panes, and templates. Use for interactive Panel UI behavior; not for the underlying scientific computation or standalone plotting API.
---

# Develop AbiPy Panel interfaces

Keep the scientific operation on its domain object and make the Panel layer a thin, lazy interactive adapter. Inspect
`abipy/panels/core.py` and the nearest domain panel before choosing widgets, callback structure, or layout. Use
`abipy-plotting` when the underlying figure API must change, `abipy-robots` for collection semantics, and
`abipy-flowtk` for workflow behavior rather than dashboard presentation.

## Locate the integration point

Panel implementations live in `abipy/panels`. File and domain objects usually expose `get_panel`, which imports the
panel class locally so importing the scientific module does not require Panel. Reuse the most specific base class:

- `AbipyParameterized` for common verbosity, execution, plotting-template, and server-mode behavior;
- `PanelWithStructure` for structure viewers, analysis, and structure-derived inputs;
- `PanelWithElectronBands` for electronic-band and DOS controls;
- `BaseRobotPanel` or `PanelWithEbandsRobot` for multi-file interfaces;
- existing node, task, work, or flow panels for FlowTK dashboards.

Do not copy a large inherited section into a specialized panel. Add shared behavior to `core.py` only when multiple
panels genuinely need the same semantics.

## Define parameters and widgets

Declare user-controlled state with `param` descriptors on the class. Give selectors explicit objects, numerical
parameters meaningful bounds, and controls clear documentation. Keep transient computation results out of shared
class-level mutable state.

Build widgets through Param or existing AbiPy helpers when possible. Reuse `mpl`, `ply`, and `dfc` for Matplotlib,
Plotly, and dataframe presentation instead of embedding backend-specific display logic. Preserve the repository's
sizing and template conventions; do not set global Panel, Bokeh, Matplotlib, or Plotly configuration from a domain
panel.

`has_remote_server` is a functional safety mode, not decoration. Disable or constrain operations that assume a local
desktop, unrestricted paths, arbitrary commands, excessive MPI resources, or trusted uploaded content. Do not
weaken those restrictions to make a local-only feature appear in a served application.

## Make expensive work lazy

Calculations, filesystem scans, external programs, and substantial plots should run only after an explicit user
action. Follow the existing button pattern:

```python
run_btn = param.Action(lambda self: self.param.trigger("run_btn"), label="Run")

@depends_on_btn_click("run_btn")
def on_run_btn(self):
    """Explain what will be computed when the button is pressed."""
    return mpl(self.obj.plot_result(show=False))
```

`depends_on_btn_click` uses the button's click count, shows the callback docstring before the first click, manages
the busy state, and normally converts exceptions into an inspectable pane. Use it instead of an eager `pn.bind` or
`param.depends` callback when the operation is costly or has side effects.

Do not perform work in `get_panel` merely because a tab is constructed. Cache only results whose inputs and
invalidation rules are clear. If widgets are shared across tabs, retain the shared-widget warning or arrange for
dependent output to recompute; stale output must not silently look current.

## Assemble the view

Follow the prevalent `get_panel(as_dict=False, **kwargs)` contract. Construct a deterministic mapping of descriptive
tab names to Panel objects, then return that mapping when `as_dict=True` or the appropriate tabs/template when
false. Preserve any established signature on the associated domain object's `get_panel` method.

Keep layout code readable: group controls separately from results, put advanced or expensive operations in their own
tabs, and avoid deeply nested anonymous rows and columns when a helper has clearer meaning. A callback should return
a Panel-compatible object consistently across successful branches, including a useful message for an empty result.

Use existing helpers for Markdown, dataframes, JSON, terminal output, loading indicators, clipboard content, and
figures. Escape or render untrusted text as text; do not interpolate uploaded or external content into raw HTML.

## Separate execution and UI concerns

A panel may invoke an existing calculation, reader, Robot, or Flow method, but it should not become the only home of
scientific logic. First implement reusable computation that returns a domain object, arrays, a dataframe, or a
figure; then adapt it to panes.

Treat buttons that run ABINIT, anaddb, schedulers, shell commands, file removal, or other mutations according to the
underlying API's execution and authorization constraints. Rendering a dashboard must never start jobs. Use temporary
paths for generated files unless the user deliberately selected a destination.

## Tests

Panel is optional, so use the repository's optional-dependency guard and place focused assertions beside the domain
object or relevant Panel module. Prefer tests that instantiate the panel without starting a server or browser and
verify:

- the expected panel/template object and tab names;
- `as_dict=True` behavior where supported;
- parameter bounds, selector contents, and remote-server restrictions;
- pre-click callbacks do not execute expensive work;
- a simulated button trigger returns the expected pane or figure wrapper;
- exception and empty-result paths remain visible to the user;
- the associated object's `get_panel` performs its local import correctly.

Patch or use small reference objects for expensive executables and network operations. Do not launch `panel serve`,
open a browser, contact Materials Project, or submit FlowTK jobs in unit tests. For a regression that depends on live
server callbacks, add the smallest integration check and document the extra requirement.

Report the domain object and panel class changed, new controls and callback triggers, any remote-server restrictions,
and the focused headless tests used.
