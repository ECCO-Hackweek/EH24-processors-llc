Introduction
============
This package offers routines to help MITgcm users format ungridded datasets into objects readable by ``pkg/profiles`` and ``pkg/obsfit`` using ``xarray``. In particular, ``obsprep`` automates certain routines related to gridding and interpolation parameters.

Requirements
^^^^^^^^^^^^

obsprep is compatible with python 3 (>= version 3.11).

Installation from GitHub
^^^^^^^^^^^^^^^^^^^^^^^^

obsprep is under active development. To obtain the latest development version, you may clone the `source repository <https://github.com/ECCO-Hackweek/EH24-processors-llc/>`_ and install it::

    git clone https://github.com/ECCO-Hackweek/EH24-processors-llc.git
    cd EH24-processors-llc/obsprep
    python setup.py install


