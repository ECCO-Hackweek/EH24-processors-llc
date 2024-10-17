Introduction
=====================
What does this look like?

.. code-block:: python

    from preproc.preproc import InSituPreprocessor

profiles
~~~~~~~~
Here is how to use the `InSituPreprocessor` with the 'profiles' package:

.. code-block:: python

     # grid_ds has fields XC and YC
     ISP = InSituPreprocessor('profiles', grid_ds)

obsfit
~~~~~~
Here is how to use the `InSituPreprocessor` with the 'obsfit' package:

.. code-block:: python

    ISP = InSituPreprocessor('obsfit', grid_ds)

