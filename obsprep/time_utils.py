import numpy as np
import pandas as pd
from datetime import datetime as dt

def datenum(d):
    return 366 + d.toordinal() + (d - dt.fromordinal(d.toordinal())).total_seconds()/(24*60*60)

def get_obs_datetime(ds, time_var, pkg_str, dims_obs, remove_original_time_field=False):
    """
    Get obsfit datetime fields

    Parameters
    ----------
    ds : xarray.Dataset 
        the dataset containing observational data
    time_var : str
        name of datetime dimension in ds
    pkg_str : str
        prefix for datetime fields (e.g. 'obs' or 'prof')
    dims_obs : list
        dimensions of datetime fields (e.g. ['iOBS'] or ['iPROF'])
    remove_original_time_field : bool
        option to remove time_var field from ds

    Returns
    -------
    xarray.Dataset
        ds with added ObsFit datetime variables
    """

    # generate time fields
    df = ds[time_var].to_dataframe()
    date = pd.to_datetime(df[time_var]).apply(lambda x:datenum(x)).values
    YYYYMMDD = df[time_var].apply(lambda x: int(x.strftime('%Y%m%d'))).values
    HHMMSS = df[time_var].apply(lambda x: int(x.strftime('%H%M%S'))).values

    # store time info in dictionary
    time_vars =  [f'{pkg_str}_{x}' for x in ['HHMMSS', 'YYYYMMDD', 'date']]
    time_flds = [HHMMSS, YYYYMMDD, date]
    time_dict = dict(zip(time_vars, time_flds))

    # add fields to dataset
    for key, val in time_dict.items():
        ds[key] = (dims_obs, val)

    # reorder time fields to first
    variables = list(ds.data_vars.keys())
    [variables.remove(x) for x in time_vars]
    [variables.insert(1, x) for x in time_vars]

    # remove time_var from dataset
    if remove_original_time_field:
        variables.remove(time_var)
    
    return ds[variables]
