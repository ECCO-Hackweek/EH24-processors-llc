Grid types
----------
You can also specify the grid type from the options ``sphericalpolar`` (default), ``llc``, or ``cubedsphere``. The latter two  options perform an additional routine generating ``_interp_`` fields so that MITgcm can recognize the processor on which a given observation resides

.. code-block:: python

    >>> from MITpreprobs.utils import generate_random_points
    >>> num_obs = 10
    >>> ungridded_lons, ungridded_lats = generate_random_points(num_obs)
    >>> UOP = UngriddedObsPreprocessor('profiles')
    >>> UOP.get_obs_point(ungridded_lons,
    >>>                   ungridded_lats,
    >>>                   grid_type = 'llc',
    >>>                   grid_noblank_ds = grid_noblank_ds,
    >>>                   num_interp_points = 4)
    >>> print(UOP.ungridded_obs_ds)

    <xarray.Dataset> Size: 3kB
    Dimensions:              (iPROF: 10, iINTERP: 4)
    Dimensions without coordinates: iPROF, iINTERP
    Data variables:
        prof_lon             (iPROF) float64 80B 87.36 56.06 97.61 ... 77.79 36.23
        prof_lat             (iPROF) float64 80B 14.06 59.9 -42.31 ... 84.97 47.71
        prof_point           (iPROF, iINTERP) uint32 160B 66365 66005 ... 81074
        prof_interp_XC11     (iPROF, iINTERP) float64 320B 82.5 82.5 ... 22.5 22.5
        prof_interp_YC11     (iPROF, iINTERP) float64 320B 10.46 10.46 ... 37.59
        prof_interp_XCNINJ   (iPROF, iINTERP) float64 320B 111.5 111.5 ... 51.5 51.5
        prof_interp_YCNINJ   (iPROF, iINTERP) float64 320B 36.8 36.8 ... 56.74 56.74
        prof_interp_i        (iPROF, iINTERP) float64 320B 6.0 6.0 5.0 ... 15.0 15.0
        prof_interp_j        (iPROF, iINTERP) float64 320B 5.0 4.0 5.0 ... 14.0 16.0
        prof_interp_weights  (iPROF, iINTERP) float64 320B 0.5348 0.2633 ... 0.1309

Alternatively, the user can provide an incomplete ungridded observation dataset and get back the gridding fields of interest

.. code-block:: python

    >>> print(ungridded_obs_ds_in)

    <xarray.Dataset> Size: 8kB
    Dimensions:       (iPROF: 10, iDEPTH: 50)
    Dimensions without coordinates: iPROF, iDEPTH
    Data variables:
        prof_lon      (iPROF) float64 80B -1.162 100.2 -9.685 ... 149.5 132.9 120.0
        prof_lat      (iPROF) float64 80B 55.35 81.67 27.57 ... -9.436 -31.64 -34.37
        prof_T        (iPROF, iDEPTH) float64 4kB 28.72 28.61 27.95 ... 4.05 3.513
        prof_Tweight  (iPROF, iDEPTH) float64 4kB 0.5819 0.6558 ... 2.741 2.835

    >>> UOP = UngriddedObsPreprocessor('profiles', ungridded_obs_ds = ungridded_obs_ds_in)
    >>> UOP.get_obs_point(grid_type = 'llc',
    >>>                   grid_noblank_ds = grid_noblank_ds,
    >>>                   num_interp_points = 4)
    >>> ungridded_obs_ds_out = UOP.ungridded_obs_ds
    >>> print(ungridded_obs_ds_out)

    <xarray.Dataset> Size: 11kB
    Dimensions:              (iPROF: 10, iDEPTH: 50, iINTERP: 4)
    Dimensions without coordinates: iPROF, iDEPTH, iINTERP
    Data variables:
        prof_lon             (iPROF) float64 80B -1.162 100.2 -9.685 ... 132.9 120.0
        prof_lat             (iPROF) float64 80B 55.35 81.67 27.57 ... -31.64 -34.37
        prof_T               (iPROF, iDEPTH) float64 4kB 28.72 28.61 ... 4.05 3.513
        prof_Tweight         (iPROF, iDEPTH) float64 4kB 0.5819 0.6558 ... 2.835
        prof_point           (iPROF, iINTERP) uint32 160B 84996 85356 ... 44078
        prof_interp_XC11     (iPROF, iINTERP) float64 320B -7.5 -7.5 ... 112.5 112.5
        prof_interp_YC11     (iPROF, iINTERP) float64 320B 37.59 37.59 ... -36.8
        prof_interp_XCNINJ   (iPROF, iINTERP) float64 320B 21.5 21.5 ... 141.5 141.5
        prof_interp_YCNINJ   (iPROF, iINTERP) float64 320B 56.74 56.74 ... -10.46
        prof_interp_i        (iPROF, iINTERP) float64 320B 7.0 7.0 8.0 ... 8.0 9.0
        prof_interp_j        (iPROF, iINTERP) float64 320B 27.0 28.0 ... 3.0 3.0
        prof_interp_weights  (iPROF, iINTERP) float64 320B 0.3112 0.271 ... 0.1552

