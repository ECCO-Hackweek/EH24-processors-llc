Introduction
============
This package offers routines to help MITgcm users format ungridded datasets into objects readable by ``pkg/profiles`` and ``pkg/obsfit`` using ``xarray``. In particular, ``obsprep`` automates certain routines related to gridding and interpolation parameters.

.. code-block:: python

    from obsprep.preproc import UngriddedObsPreprocessor

profiles
~~~~~~~~
Here is how to use the `UngriddedObsPreprocessor` with the 'profiles' package:

.. code-block:: python

     # grid_noblank_ds has fields XC and YC and was generated without blank tiles
     UOP = UngriddedObsPreprocessor('profiles', grid_noblank_ds)

obsfit
~~~~~~
Here is how to use the `UngriddedObsPreprocessor` with the 'obsfit' package:

.. code-block:: python

    UOP = UngriddedObsPreprocessor('obsfit', grid_noblank_ds)

