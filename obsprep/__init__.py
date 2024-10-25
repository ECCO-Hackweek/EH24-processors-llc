_all__ = [
    "utils","preproc","interp",
]

## Add to __init__.py
# from .prep import Prep
#
#
## Then usage looks like:
#
# import obsprep as op
#
# swot = xr.open_dataset('/path/to/swot.nc')
# OP = op.Prep(swot)
# OP.get_obs_point()
# ...
# 
# and if the user wants to access individual functions from the library they can run
# op.interp.get_interp_points

