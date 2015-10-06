This directory aims at collecting integration tests for the Zenobe infrastructure

Don't run these jobs on your local machine because they are quite heavy !

Description of the files :
- test_abipy_limits.py contains a GS, a NSCF and a BSE task that will fail at first with the manager because of the too small walltime, then because of memory problems... This should trigger automatic fixers for time, cpu and memory !
