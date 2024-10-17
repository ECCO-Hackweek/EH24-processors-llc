Grid types
----------
You can also specify the grid type from the options ``sphericalpolar`` (default), ``llc``, or ``cubedsphere``. The latter two  options perform an additional routine generating ``_interp_`` fields so that MITgcm can recognize the processor on which a given observation resides

.. code-block:: python

    ISP = InSituPreprocessor('obsfit', grid_ds, grid='llc')

    def generate_random_points(nobs, lon_range=(-180, 180), lat_range=(60, 90)):
         lons = np.random.uniform(low=lon_range[0], high=lon_range[1], size=nobs)
         lats = np.random.uniform(low=lat_range[0], high=lat_range[1], size=nobs)
         return lons, lats

    nobs = 10
    ungridded_lons, ungridded_lats = generate_random_points(nobs)
    ISP = InSituPreprocessor('obsfit', grid_ds, grid='llc')
    ISP.get_obs_point(ungridded_lons, ungridded_lats)
    ISP.get_sample_interp_info()


The resulting dataset ``ISP.ds`` will have the fields::

    <xarray.Dataset> Size: 720B
    Dimensions:            (iOBS: 10, iINTERP: 1)
    Dimensions without coordinates: iOBS, iINTERP
    Data variables:
        obs_point          (iOBS) int64 80B 88733 114170 90578 ... 124240 112345
        obs_lon            (iOBS) float64 80B 135.2 122.5 -179.2 ... -165.4 -76.1
        obs_lat            (iOBS) float64 80B 60.06 88.17 63.53 ... 78.58 83.08
        obs_interp_XC11    (iOBS, iINTERP) float64 80B 112.5 52.0 ... 143.9 -39.94
        obs_interp_YC11    (iOBS, iINTERP) float64 80B 57.28 82.38 ... 82.11 82.11
        obs_interp_XCNINJ  (iOBS, iINTERP) float64 80B 141.8 -128.0 ... -101.1
        obs_interp_YCNINJ  (iOBS, iINTERP) float64 80B 67.47 82.38 ... 71.6 71.6
        obs_interp_i       (iOBS, iINTERP) float64 80B 24.0 18.0 19.0 ... 16.0 13.0
        obs_interp_j       (iOBS, iINTERP) float64 80B 7.0 10.0 9.0 ... 20.0 5.0
