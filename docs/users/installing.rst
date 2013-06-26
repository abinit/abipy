.. The source of this document is INSTALL. During the doc build process,
.. this file is copied over to doc/users/installing.rst.
.. Therefore, you must edit INSTALL, *not* doc/users/installing.rst!

**********
Installing
**********

.. _install_requirements:

Build requirements
==================

These are external packages which you will need to install before installing abipy:

#. python 2.7 (or later but not python3)

#. pymatgen 2.7.3 (or later)

#. numpy 1.5 (or later) array support for python 
   (`download <http://sourceforge.net/project/showfiles.php?group_id=1369&package_id=175103>`__)

#. scipy 0.10 (or later).

#. matplotlib 1.1 (or later)

.. _install_from_source:

Installing from source
======================

Once you have satisfied the requirements detailed above, you can build abipy::

  cd abipy
  python setup.py build
  python setup.py install


