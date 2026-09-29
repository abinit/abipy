# Integration tests for Abinit/AbiPy

## Requirements

- pytest (framework for unit tests)
- pytest-timeout (each test is killed after the `timeout` given in pytest.ini)
- abinit
- abipy with configured manager.yaml


## How to run the tests

Change the configuration options for the task manager reported in the file `manager.yaml`

Execute pytest in the current directory.
Results are produced in the directory `_integration_tests_`.
Use `pytest -v` to print the TaskManager configurations in the header.

By default, the abipy autoparal is disabled and the tasks are executed with a fixed number
of MPI processes so that the results do not depend on the machine.
Inputs with `paral_kgb 1` and without `np_spkpt`, `npband`, `npfft` get the Abinit variable
`autoparal 1` so that Abinit selects the distribution at runtime.
The following environment variables change this behaviour:

- `ABIPY_ITEST_MPI_PROCS`: number of MPI processes (default: 2).
- `ABIPY_ITEST_AUTOPARAL=1`: run the tests also with the abipy autoparal.

Hanging tasks are detected by the scheduler via the `frozen_timeout` policy (2 minutes).
When a test fails, the report includes the flow workdir and the last lines of the files
produced by the tasks that failed.

## How to add a new integration test

The file pytest.ini contains the configuration options passed to py.test

Test functions should start with the prefix `itest_` where `i` stands for integration.

Each test function receives the fixture arguments `fwp` and `tvars` defined in conftest.py.
`fwp` contains the parameters used to generate the `AbinitFlow`, whereas `tvars` is a dictionary
with the Abinit variables used to generate the input files.
pytest will generate multiple test for the each `itest_` function
that receives `fwp` and `tvars` in input. In pseudo-code:

```python
for fwp in fwp_list:
    for tvars in tvars_list:
        test = build_test_with(fwp, tvars)
        test.run()
```
