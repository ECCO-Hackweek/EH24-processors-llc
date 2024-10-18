import copy
import xarray as xr
import numpy as np
from ecco_v4_py.llc_array_conversion import llc_tiles_to_faces, llc_tiles_to_compact
from MITpreprobs.utils import patchface3D_5f_to_wrld, compact2worldmap
from MITpreprobs.interp import *

class UngriddedObsPreprocessor:
    """
    A class for preprocessing ungridded observational data for profiles or obsfit.

    This class handles the initialization and methods required for
    preprocessing ungridded data, including finding the nearest grid points
    and performing interpolation.
    """    
    def __init__(self, pkg, ungridded_obs_ds=None, sNx=30, sNy=30):
        """
        Initialize the UngriddedObsPreprocessor class.

        Parameters
        ----------
        pkg: str
            The package type (either 'profiles' or 'obsfit').
        ungridded_obs_ds: xarray.Dataset, optional
            The dataset containing in-situ data. If None, an empty dataset is created.
        sNx: int, optional
            The size of the MPI-partition tiles in the x-direction (default is 30).
        sNy: int, optional
            The size of the MPI-partition tiles in the y-direction (default is 30).
        """

        # check inputs
        if pkg not in ['profiles', 'obsfit']:
            raise ValueError(f"Invalid pkg '{pkg}'. Must be either 'profiles' or 'obsfit'.")

        self.pkg = pkg.lower()
#        self.msk = grid_noblank_ds.mskC.where(grid_noblank_ds.mskC).isel(k=0).values
        self.sNx = sNx
        self.sNy = sNy

        self.ungridded_obs_ds = xr.Dataset() if ungridded_obs_ds is None else ungridded_obs_ds
            
        self.get_pkg_fields()

    def get_pkg_fields(self):
        """
        Set package-specific attributes based on the 'pkg' input.
        """
        # Define dimensions for in-situ, interpolation, and depth based on the package 
        self.pkg_str = 'prof' if self.pkg == 'profiles' else 'obs'
        self.iPKG = f'i{self.pkg_str.upper()}'
        self.dims_obs = [self.iPKG]
        self.dims_depth = ['iDEPTH' if self.pkg_str == 'prof' else '']
        
        # Combine dimensions for spatial fields, including depth if applicable
        self.dims_spatial = self.dims_obs + self.dims_depth * (len(self.dims_depth[0]) > 0)
        self.lon_str = 'prof_lon' if self.pkg_str == 'prof' else 'sample_lon'
        self.lat_str = 'prof_lat' if self.pkg_str == 'prof' else 'sample_lat'
   
    def get_obs_point(self, ungridded_lons=None, ungridded_lats=None, 
                                grid_type='sphericalpolar', grid_noblank_ds=None,
                                num_interp_points=1,
                               ):
        """
        Find the nearest grid point for given ungridded longitude and latitude coordinates.
    
        Parameters
        ----------
        ungridded_lons : list of float, optional
            List of longitudes for observation points.
        ungridded_lats : list of float, optional
            List of latitudes for observation points.
        grid_type : str, optional
            The type of grid being used ('sphericalpolar', 'llc', or 'cubedsphere').
        grid_noblank_ds : xarray.Dataset, optional
            Dataset containing grid information without blanks (must have fields XC and YC).
        num_interp_points : int, optional
            Number of interpolation points to consider (default is 1).
    
        Raises
        ------
        ValueError
            If the grid type is invalid or if required data is not provided.
        """
        ds_has_lon = f'{self.lon_str}' in self.ungridded_obs_ds.keys()
        ds_has_lat = f'{self.lat_str}' in self.ungridded_obs_ds.keys()

        if (ds_has_lon) & (ds_has_lat):
            if grid_type == 'sphericalpolar':
                print(f'Ungridded dataset already has fields {self.lon_str} and {self.lat_str}. Leaving get_obs_point.')
                return
            
        if (ungridded_lons is None) & (not ds_has_lon):
            raise ValueError(f"Ungridded longitudes not provided")
        if (ungridded_lats is None) & (not ds_has_lat):
            raise ValueError(f"Ungridded latitudes not provided")
            
        if not ds_has_lon:
            if len(ungridded_lons) == 0:
                ungridded_lons = [ungridded_lons]
            self.ungridded_obs_ds[self.lon_str] = (self.dims_obs, ungridded_lons)           
        
        if not ds_has_lat:
            if len(ungridded_lats) == 0:
                ungridded_lats = [ungridded_lats]
            self.ungridded_obs_ds[self.lat_str] = (self.dims_obs, ungridded_lats)

        if grid_type == 'sphericalpolar':
            return
            
        if grid_type not in {'llc', 'cubedsphere'}:
            raise ValueError(f"Invalid grid_type: '{grid_type}'. Options supported are 'sphericalpolar', 'llc', or 'cubedsphere'.")

        # grid_type is llc or curvilinear
        if grid_noblank_ds is not None:
            if 'XC' not in list(grid_noblank_ds.coords) or 'YC' not in list(grid_noblank_ds.coords):
                raise ValueError(f"grid_noblank_ds must have fields XC and YC.")
        else:
                raise ValueError("provide grid_noblank_ds")
            
        # add attributes relevant to interpolation field creation
        self.xc = grid_noblank_ds.XC.values
        self.yc = grid_noblank_ds.YC.values
        self.mask = grid_noblank_ds.hFacC.where(grid_noblank_ds.hFacC).isel(k=0).values

        # same xc, yc in worldmap form for later        
        nx = len(self.xc[0, 0, :]) # (last two dimensions of xc with shape (ntile, nx, nx)
        self.xc_wm = compact2worldmap(llc_tiles_to_compact(self.xc, less_output=True), nx, 1)[0, :, :]
        self.yc_wm = compact2worldmap(llc_tiles_to_compact(self.yc, less_output=True), nx, 1)[0, :, :]
        self.mask_wm = compact2worldmap(llc_tiles_to_compact(self.mask, less_output=True), nx, 1)[0, :, :]
        
        self.max_target_grid_radius = np.max(np.abs([grid_noblank_ds.dxG.values, grid_noblank_ds.dyG.values]))

        self.num_interp_points = num_interp_points        
        self.dims_interp = self.dims_obs + ['iINTERP']

        # run interpolation routine
        self.obs_points = self.interp()

        # add fields to ungridded_obs_ds
        self.interp_str = 'prof' if self.pkg_str == 'prof' else 'sample'
        self.obs_point_str = f'{self.interp_str}_point'

        if self.obs_points.ndim == 1:
            self.obs_points = self.obs_points[:, None]
        self.ungridded_obs_ds[self.obs_point_str] = (self.dims_interp, self.obs_points)

        # add grid interp fields
        self.get_sample_interp_info()

    def interp(self):
        """
        Perform interpolation from observation points to grid points using pyresample.
    
        Returns
        -------
        index_array : ndarray
            Array containing indices of the nearest grid points for the given observation points.
        """
        # turn from tiles to worldmap                
        valid_input_index, valid_output_index, index_array, distance_array =\
            get_interp_points(
                self.ungridded_obs_ds[self.lon_str],
                self.ungridded_obs_ds[self.lat_str],
                self.xc_wm.ravel(),
                self.yc_wm.ravel(),
                nneighbours=self.num_interp_points,
               
 max_target_grid_radius=self.max_target_grid_radius
            )
        self.interp_distance = distance_array
        
        return index_array
        
        
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
        # save tiles in worldmap views for interp later
        # probably a bad idea memory-wise for high res grids
        tile_fields_wm = dict.fromkeys(tile_keys)

        for tile_key, tile_field in tile_fields.items():
            tile_field = {iF: tile_field[iF][None, :] for iF in tile_field}

            # TODO: 
            #     Replace with ecco.faces_to_tiles
            #     then tiles_to_world
            tile_field_wm = patchface3D_5f_to_wrld(tile_field)
            tile_field_at_obs_point = tile_field_wm.ravel()[self.obs_points]

            tile_fields_wm[tile_key] = tile_field_wm

            # profiles default: make singleton iINTERP dimension
            #                   set sample_interp_weight np.ones
            #                   this is encoded in len(self.dims_interp)
            tile_field_at_obs_point = tile_field_at_obs_point.reshape(tile_field_at_obs_point.shape +\
                                      (1,) * (len(self.dims_interp) - tile_field_at_obs_point.ndim))

            # set interp field in dataset
            self.ungridded_obs_ds[f'{self.interp_str}_interp_{tile_key}'] = (self.dims_interp, tile_field_at_obs_point)

            # save tile_fields
            self.tile_fields_wm = tile_fields_wm

        self.get_interp_weights()

    def get_interp_weights(self):
        """ 
        Set interpolation weights
        """
        print('Warning: currently only supported inverse distance weighting')
        inv_dist = 1 / self.interp_distance

        interp_weights = inv_dist / np.sum(inv_dist, axis=1, keepdims=True)
        # if self.num_interp_points == 4:
        # self.method = 'bilinear'
        
        # compute interp_weights from obs_points
        # convert from lat-lon to xyz
        # compute bilinear interp weights
        
        self.ungridded_obs_ds[f'{self.interp_str}_interp_weights'] = (self.dims_interp, interp_weights)
        
