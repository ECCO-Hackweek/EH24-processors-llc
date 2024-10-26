Introduction
============
This package offers routines to help MITgcm users format ungridded datasets into objects readable by ``pkg/profiles`` and ``pkg/obsfit`` using ``xarray``. In particular, ``obsprep`` automates certain routines related to gridding and interpolation parameters.

.. code-block:: python

    from obsprep import Prep

profiles
~~~~~~~~
Here is how to use the `Prep` with the 'profiles' package:

.. code-block:: python

     # grid_ds has fields XC and YC and was generated without blank tiles
     OP = Prep('profiles', grid_ds)

obsfit
~~~~~~
Here is how to use the `Prep` with the 'obsfit' package:

.. code-block:: python

    OP = Prep('obsfit', grid_ds)

