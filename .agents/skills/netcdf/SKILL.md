---
name: netcdf
description: Diagnose NetCDF schemas and values across AbiPy's Python readers and ABINIT's Fortran writers, including dimensions, slices, units, complex values, and Fortran/C axis ordering. Use for .nc parsing or data discrepancies.
---

# Diagnose NetCDF data

Trace a value from the ABINIT writer to the file schema and then to the AbiPy reader. Determine whether a discrepancy
comes from the value constructed by ABINIT, serialization, units, axis order, or AbiPy interpretation.

## Start from both ends

In AbiPy, locate the reader with searches for the variable name, `read_value`, `read_variable`, `ETSF_Reader`, and the
relevant `*File` class. Check transformations after reading: transpose/reshape operations, unit conversion, complex
reconstruction, spin conventions, and index normalization.

In an ABINIT checkout, locate the corresponding `nf90_def_var` and `nf90_put_var` calls and inspect the exact array or
slice written. If its `abisrc.py` supports `nc_explain`, prefer:

```bash
./abisrc.py nc_explain PROCEDURE --follow 2
./abisrc.py nc_explain PROCEDURE --check FILE.nc
```

Otherwise search the Fortran source directly. Do not assume the file was produced by the currently checked-out
ABINIT version.

## Axis ordering

The NetCDF Fortran interface maps Fortran array order to the declared file dimension order. C-model tools such as
Python `netCDF4`, xarray, `ncdump`, and `h5dump` commonly show dimensions in reverse order relative to the Fortran
declaration. For example, a Fortran `x(n1,n2,n3)` can appear as `(n3,n2,n1)` in Python.

This is not by itself evidence of a bad transpose. Confirm the dimension definitions and write call. Never reverse
axes mechanically when the call uses `start`, `count`, an array slice, or an intermediate buffer.

## Investigation checklist

1. Inspect the file header and variable attributes with `ncdump -h FILE.nc` or Python.
2. Record dimension names, not just sizes, plus units, fill values, groups, and dtype.
3. Locate the exact AbiPy read path and every transformation applied afterward.
4. Locate the exact Fortran declaration, construction, and write call or slice.
5. Account for unit conversion, packed complex numbers (often a final axis of size 2), one-based indices, symmetry,
   parallel decomposition, and unwritten/fill entries.
6. Compare a few representative physical indices rather than raw memory layouts.

Minimal Python inspection:

```python
from netCDF4 import Dataset

with Dataset("file.nc") as nc:
    var = nc.variables["variable_name"]
    print(var.dimensions, var.shape, var.dtype)
    print(var.__dict__)
    print(var[:])
```

When changing a reader, add or update a focused test beside its module and use a small file from `abipy/data/refs`
when one already exercises the schema. Do not add large binary fixtures casually; first check whether an existing
reference file covers the case.
