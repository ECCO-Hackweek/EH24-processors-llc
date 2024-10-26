Temporal metadata
-----------------
Users' input datasets can include a datetime field used to derive the time ``YYYYMMDD``, ``HHMMSS``, and ``date`` (Julian date). 

.. code-block:: python

    >>> print(ds_in)

    <xarray.Dataset> Size: 2MB
    Dimensions:       (iOBS: 243157)
    Dimensions without coordinates: iOBS
    Data variables:
        obs_datetime  (iOBS) datetime64[ns] 2MB ...

    >>> OP = op.Prep('obsfit', ds_in)
    >>> OP.get_obs_datetime(time_vars = 'obs_datetime')
    >>> ds_out = OP.ds
    >>> print(OP.ds)

    <xarray.Dataset> Size: 8MB
    Dimensions:       (iOBS: 243157)
    Dimensions without coordinates: iOBS
    Data variables:
        obs_datetime  (iOBS) datetime64[ns] 2MB ...
        obs_date      (iOBS) float64 2MB 7.393e+05 7.393e+05 ... 7.393e+05 7.393e+05
        obs_YYYYMMDD  (iOBS) int64 2MB 20240104 20240104 ... 20240104 20240104
        obs_HHMMSS    (iOBS) int64 2MB 190633 190633 190633 ... 193939 193939 193939
