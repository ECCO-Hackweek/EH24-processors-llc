import xarray as xr
import numpy as np
from scipy.spatial import KDTree
from ecco_v4_py.llc_array_conversion import llc_tiles_to_faces, llc_tiles_to_compact
from .utils import patchface3D_5f_to_wrld, compact2worldmap
import copy

class InSituPreprocessor:
    def __init__(self, pkg, ds=None, grid='sphericalpolar', sNx=30, sNy=30):
        """
        Initialize the InSituPreprocessor class.

        Parameters:
        - pkg: The package type (either 'profiles' or 'obsfit').
        - grid_noblank_ds: xarray dataset containing the grid information (XC, YC).
        - grid: The type of grid being used (default 'sphericalpolar').
        - sNx, sNy: MPI-partition tile sizes.
        """

        # check inputs
        if pkg not in ['profiles', 'obsfit']:
            raise ValueError(f"Invalid pkg '{pkg}'. Must be either 'profiles' or 'obsfit'.")

        self.pkg = pkg.lower()
        self.grid = grid
#        self.msk = grid_noblank_ds.mskC.where(grid_noblank_ds.mskC).isel(k=0).values
        self.sNx = sNx
        self.sNy = sNy

        if ds is not None:
            self.ds = xr.Dataset()
            
        self.get_pkg_fields()

    def get_pkg_fields(self):
        """
        Set package-specific attributes based on the 'pkg' input.
        """
        # Define dimensions for in-situ, interpolation, and depth based on the package 
        self.pkg_str = 'prof' if self.pkg == 'profiles' else 'obs'
        self.iPKG = f'i{self.pkg_str.upper()}'
        self.dims_insitu = [self.iPKG]
        self.dims_interp = self.dims_insitu + ['iINTERP'] * (self.pkg_str == 'prof')
        self.dims_depth = ['iDEPTH' if self.pkg_str == 'prof' else '']
        
        # Combine dimensions for spatial fields, including depth if applicable
        self.dims_spatial = self.dims_insitu + self.dims_depth * (len(self.dims_depth[0]) > 0)
        self.lon_str = 'prof_lon' if self.pkg_str == 'prof' else 'sample_lon'
        self.lat_str = 'prof_lat' if self.pkg_str == 'prof' else 'sample_lat'
   
    def get_obs_point(self, grid_noblank_ds, ungridded_lons, ungridded_lats):
        """
        Find the nearest grid point for given ungridded longitude and latitude coordinates.
        
        Parameters:
        - ungridded_lons: List of longitudes for observation points.
        - ungridded_lats: List of latitudes for observation points.

        TODO: Avoid land, investigate pyresample
        """
        
        if 'XC' not in list(grid_noblank_ds.coords) or 'YC' not in list(grid_noblank_ds.coords):
            raise ValueError(f"grid_noblank_ds must have fields XC and YC.")
        if len(ungridded_lons) == 0:
            ungridded_lons = [ungridded_lons]
        if len(ungridded_lats) == 0:
            ungridded_lats = [ungridded_lats]

        
        self.xc = grid_noblank_ds.XC.values
        self.yc = grid_noblank_ds.YC.values

        obs_point = _get_nearest_grid_index(self.xc, self.yc, ungridded_lons, ungridded_lats)
        self.obs_point_str = f'{self.pkg_str}_point'
        self.interp_str = 'prof' if self.pkg_str == 'prof' else 'sample'
        setattr(self, self.obs_point_str, obs_point)
        
        # add fields to dataset
        self.ds[self.obs_point_str] = (self.dims_insitu, obs_point)
        self.ds[self.lon_str] = (self.dims_insitu, ungridded_lons)
        self.ds[self.lat_str] = (self.dims_insitu, ungridded_lats)

        if self.grid in ['llc', 'cubedsphere']:
            self.get_sample_interp_info()

    def get_sample_interp_info(self):
        """
        Interpolate sample grid information and store in the dataset.
        """

        def empty_faces(xgrid):
            return {face: np.zeros_like(xgrid[face]) for face in range(1, 6)}
        
        # transform xc and yc from worldmap back to 5faces
        xc_5faces = llc_tiles_to_faces(self.xc, less_output=True)
        yc_5faces = llc_tiles_to_faces(self.yc, less_output=True)
        
        xc_5faces.copy()
        yc_5faces.copy()
        
        # Create dictionary of fields using deepcopy
        tile_keys = ['XC11', 'YC11', 'XCNINJ', 'YCNINJ', 'i', 'j']
        tile_fields = {key: copy.deepcopy(empty_faces(xc_5faces)) for key in tile_keys}
        XC11, YC11, XCNINJ, YCNINJ, iTile, jTile = (tile_fields[key] for key in tile_keys)
        
        tile_count = 0
        # Iterate through each face and divide into tiles        
        for iF in range(1, len(xc_5faces)+1):
            face_XC = xc_5faces[iF]
            face_YC = yc_5faces[iF]
            for ii in range(face_XC.shape[0] // self.sNx):
                for jj in range(face_XC.shape[1] // self.sNy):
                    tile_count += 1
                    tmp_i = slice(self.sNx * ii, self.sNx * (ii + 1))
                    tmp_j = slice(self.sNy * jj, self.sNy * (jj + 1))
                    tmp_XC = face_XC[tmp_i, tmp_j]
                    tmp_YC = face_YC[tmp_i, tmp_j]
                    XC11[iF][tmp_i, tmp_j] = tmp_XC[0, 0]
                    YC11[iF][tmp_i, tmp_j] = tmp_YC[0, 0]
                    XCNINJ[iF][tmp_i, tmp_j] = tmp_XC[-1, -1]
                    YCNINJ[iF][tmp_i, tmp_j] = tmp_YC[-1, -1]
                    iTile[iF][tmp_i, tmp_j] = np.outer(np.ones(self.sNx), np.arange(1, self.sNy + 1))
                    jTile[iF][tmp_i, tmp_j] = np.outer(np.arange(1, self.sNx + 1), np.ones(self.sNy))

        # Grab values of these fields at [prof/obs]_point    
        for tile_key, tile_field in tile_fields.items():
            tile_field = {iF: tile_field[iF][None, :] for iF in tile_field}

            # TODO: 
            #     Replace with ecco.faces_to_tiles
            #     then tiles_to_world
            tile_field_wm = patchface3D_5f_to_wrld(tile_field)
            tile_field_at_obs_point = tile_field_wm.ravel()[getattr(self, self.obs_point_str)]

            # profiles default: make singleton iINTERP dimension
            #                   set sample_interp_weight np.ones
            #                   this is encoded in len(self.dims_interp)
            tile_field_at_obs_point = tile_field_at_obs_point.reshape(tile_field_at_obs_point.shape +\
                                      (1,) * (len(self.dims_interp) - tile_field_at_obs_point.ndim))

            # set interp field in dataset
            self.ds[f'{self.interp_str}_interp_{tile_key}'] = (self.dims_interp, tile_field_at_obs_point)
            

def _get_nearest_grid_index(xc, yc, ungridded_lons, ungridded_lats):
    """
    Find the nearest grid index for given longitude and latitude points using KDTree.

    Parameters:
    - xc, yc: Gridded longitude and latitude arrays.
    - ungridded_lons, ungridded_lats: Ungridded longitude and latitude values for observations.
    """
    # turn from tiles to worldmap
    nx = len(xc[0, 0, :]) # (last two dimensions of xc with shape (ntile, nx, nx)
    xc_wm = compact2worldmap(llc_tiles_to_compact(xc, less_output=True), nx, 1)[0, :, :]    
    yc_wm = compact2worldmap(llc_tiles_to_compact(yc, less_output=True), nx, 1)[0, :, :]    
    grid_shape = xc_wm.shape

    # set up nearest neighbors KDtree
    gridded_coords = np.c_[yc_wm.ravel(), xc_wm.ravel()]
    ungridded_coords = np.c_[ungridded_lats, ungridded_lons]

    kd_tree = KDTree(gridded_coords)
    distance, nearest_grid_idx = kd_tree.query(ungridded_coords, k=1)
    # Ensure all indices are valid (i.e. land in the worldmap-sized array)
    assert((nearest_grid_idx>np.prod(grid_shape)).sum()==0)
    return nearest_grid_idx
