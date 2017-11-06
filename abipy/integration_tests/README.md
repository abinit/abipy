# Integration tests for Abinit/AbiPy

## Requirements

- abinit >= 8.6
- abipy, pymatgen and corresponding dependencies
- pytest (framework for unit tests)

## How to run the tests

Change the configuration options for the task manager reported in the file `manager.yaml`

Execute pytest in the current directory. 
Results are produced in the directory in `_mytests`

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
