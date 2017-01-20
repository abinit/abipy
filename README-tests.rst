Build/install process use Python distutils
==========================================

Tests use `nose <http://somethingaboutorange.com/mrl/projects/nose>`_ framework.

Running tests
=============

To run unittest::

    nosetests

If nosetest fails with error::

   `Error reading config file 'setup.cfg': no such option 'exclude-dir`

install the nose-exclude plugin with::

   `pip install nose-exclude`

and scriptest with::

   `pip install scripttest`


