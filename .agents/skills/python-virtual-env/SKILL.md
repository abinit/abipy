---
name: python-virtual-env
description: Select and activate the correct existing Python environment before running Python tools, tests, or build scripts. Prefer Conda environments under ~/miniconda3; ask when repository evidence does not identify a reliable environment.
---

# Select a Python environment

Choose and activate the correct existing Python environment before running project Python commands. Do not assume
that `python`, `pip`, or `pytest` on the default `PATH` belongs to the repository.

Use this skill for Python tools, tests, documentation builders, import diagnostics, and dependency checks. Skip
discovery only when the user supplied an absolute interpreter or the active interpreter has been verified in the
current command session.

## Selection policy

Use evidence in this order:

1. An environment or interpreter explicitly named by the user.
2. Repository instructions and environment files.
3. An active environment whose interpreter and relevant imports have been verified.
4. Names and package contents of existing Conda environments.

Never select an environment only because its Python version seems plausible. If multiple environments remain
plausible and the choice could affect results, ask the user. Do not create an environment or install packages merely
to avoid that question.

## Discover and activate Conda environments

The user normally installs Miniconda under `$HOME/miniconda3`. In a non-interactive shell:

```bash
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda env list
```

If that installation does not exist, use `command -v conda` or ask where Conda is installed. Do not run
`conda init`; it changes shell configuration and is unnecessary for one command.

Activation and the dependent command must run in the same shell invocation:

```bash
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate ENVIRONMENT_NAME
python --version
python -c 'import sys; print(sys.executable)'
python -c 'import abipy; print(abipy.__file__)'
```

`conda run -n ENVIRONMENT_NAME python ...` is an alternative when activation is impractical. Do not silently switch
away from an explicitly requested environment if activation fails.

## Verify suitability

Before a long or mutating workflow, verify the interpreter and important project import:

```bash
which python
python --version
python -c 'import sys; print(sys.executable)'
python -c 'import abipy; print(abipy.__file__)'
```

For this repository, the `abipy` path should resolve to the working tree unless an installed-package test was
explicitly requested. Prefer `python -m pytest` so pytest belongs to the selected interpreter.

## Package changes

Install or update packages only when requested or when it is an authorized, necessary part of the task. Prefer the
project's established package manager and use `python -m pip` instead of bare `pip`. Never install into Conda `base`
by default, and do not perform bulk environment updates as routine troubleshooting.

Environment creation, export, and bulk updates are separate mutations and require explicit user intent.

When environment selection affected the task, report the environment name and interpreter path.
