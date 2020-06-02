To update the python modules using the more recent version available in the Abinit repository use vimdiff:

    vimdiff $ABINIT_REPO/abimkdocs/variables_abinit.py variables_abinit.py

and patch the file manually!
Do not change the initial part of the module since it's needed by AbiPy.
