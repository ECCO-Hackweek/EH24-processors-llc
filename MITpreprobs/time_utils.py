import numpy as np
import pandas as pd
from datetime import datetime as dt

def datenum(d):
    return 366 + d.toordinal() + (d - dt.fromordinal(d.toordinal())).total_seconds()/(24*60*60)

def get_obs_datetime(ds, time_var, remove_original_time_field=False):
    """
    Get obsfit datetime fields

    Parameters
    ----------
    ds: xarray.Dataset 
    time_var: str
        name of datetime dimension in ds
    
    Returns
    -------
    xarray.Dataset
        ds with added ObsFit datetime variables
    """

    df = ds[time_var].to_dataframe()

    obs_date = pd.to_datetime(df[time_var]).apply(lambda x:datenum(x)).values
    obs_YYYYMMDD = df[time_var].apply(lambda x: int(x.strftime('%Y%m%d'))).values
    obs_HHMMSS = df[time_var].apply(lambda x: int(x.strftime('%H%M%S'))).values

    ds['obs_date'] = ("iOBS", obs_date)
    ds['obs_YYYYMMDD'] = ("iOBS", obs_YYYYMMDD)
    ds['obs_HHMMSS'] = ("iOBS", obs_HHMMSS)

    variables = list(ds.data_vars.keys())

    variables.remove('obs_HHMMSS')
    variables.remove('obs_YYYYMMDD')
    variables.remove('obs_date')
    variables.insert(1, 'obs_HHMMSS')
    variables.insert(1, 'obs_YYYYMMDD')
    variables.insert(1, 'obs_date')

    if remove_original_time_field:
        variables.remove('obs_datetime')
    
    return ds[variables]