Grid types
----------
You can also specify the grid type from the options ``sphericalpolar`` (default), ``llc``, or ``cubedsphere``. The latter two  options perform an additional routine generating ``_interp_`` fields so that MITgcm can recognize the processor on which a given observation resides

.. code-block:: python

    ISP = InSituPreprocessor('profiles')
    ISP.get_obs_point(ungridded_lons,
                      ungridded_lats,
                      grid_type = 'llc',
                      grid_noblank_ds = grid_noblank_ds,
                      num_interp_points = 4)


The resulting dataset ``ISP.ds`` will have the fields::

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