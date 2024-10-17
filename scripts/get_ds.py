import pyresample
import numpy as np
from os.path import expanduser,join,isdir
import sys
user_home_dir = expanduser('~')
sys.path.append('/home/jovyan/ECCO-v4-Python-Tutorial/ecco_access/')


import ecco_v4_py as ecco
import ecco_access as ea

def get_ds():
    access_mode = 's3_open_fsspec'
    
    
    # ECCO_dir specifies parent directory of all ECCOv4r4 downloads
    # ECCO_dir = None downloads to default path ~/Downloads/ECCO_V4r4_PODAAC/
    ECCO_dir = join('/efs_ecco','ECCO_V4r4_PODAAC')
    
    # for access_mode = 's3_open_fsspec', need to specify the root directory 
    # containing the jsons
    jsons_root_dir = join('/efs_ecco','mzz-jsons')
    
    
    ShortNames_list = ["ECCO_L4_TEMP_SALINITY_LLC0090GRID_MONTHLY_V4R4"]
    
    # retrieve files
    StartDate = '2010-01'
    EndDate = '2010-12'
    ds = ea.ecco_podaac_to_xrdataset(ShortNames_list,\
                                     StartDate=StartDate,EndDate=EndDate,\
                                     mode=access_mode,\
                                     download_root_dir=ECCO_dir,\
                                     max_avail_frac=0.5,\
                                     jsons_root_dir=jsons_root_dir)
    return ds
