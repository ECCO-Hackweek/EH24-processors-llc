import copy
import xarray as xr
import numpy as np
from ecco_v4_py.llc_array_conversion import llc_tiles_to_faces, llc_tiles_to_compact
from obsprep.utils import patchface3D_5f_to_wrld, compact2worldmap, get_sample_type
from obsprep.interp import *
from obsprep.time_utils import get_obs_datetime

class Prep:
    """
    A class for preprocessing ungridded observational data for profiles or obsfit.

    This class handles the initialization and methods required for
    preprocessing observational data, including generating interpolation and 
    temporal metadata.
    """    
    def __init__(self, pkg, ds=None):
        """
        Initialize the Prep class.

        Parameters
        ----------
        pkg: str
            The package type (either 'profiles' or 'obsfit').
        ds: xarray.Dataset, optional
            The dataset containing observational data. If None, an empty dataset is created.
        """

        # check inputs
        if pkg not in ['profiles', 'obsfit']:
            raise ValueError(f"Invalid pkg '{pkg}'. Must be either 'profiles' or 'obsfit'.")

        self.pkg = pkg.lower()

        self.ds = xr.Dataset() if ds is None else ds
            
        self.get_pkg_fields()

    def get_pkg_fields(self):
        """
        Set package-specific attributes based on the 'pkg' input.
        """
        # Define dimensions for in-situ, interpolation, and depth based on the package 
        self.pkg_str, self.dims_obs, self.dims_depth = ('prof', ['iPROF'], ['iDEPTH']) if self.pkg == 'profiles' else ('obs', ['iSAMPLE'], [''])
        
        # Combine dimensions for spatial fields, including depth if applicable
        self.dims_spatial = self.dims_obs + self.dims_depth * (len(self.dims_depth[0]) > 0)
        self.lon_str, self.lat_str, self.depth_str = tuple([f'{self.dims_obs[0][1:].lower()}_{x}' for x in ['lon', 'lat', 'depth']])


    def get_obs_point(self, lons=None, lats=None, depths=None,
                                grid_type='sphericalpolar', grid_ds=None,
                                num_interp_points=1, sNx=30, sNy=30, max_target_grid_radius=None,
                               ):
        """
        Find the nearest grid point for given ungridded longitude and latitude coordinates.
    
        Parameters
        ----------
        lons : list of float, optional
            List of longitudes for observation points.
        lats : list of float, optional
            List of latitudes for observation points.
        grid_type : str, optional
            The type of grid being used ('sphericalpolar', 'llc', or 'cubedsphere').
        grid_ds : xarray.Dataset, optional
            Dataset containing grid information without blanks (must have fields XC and YC).
        num_interp_points : int, optional
            Number of interpolation points to consider (default is 1; only 1, 4, and 8 are supported).
        sNx: int, optional
            The size of the MPI-partition tiles in the x-direction (default is 30).
        sNy: int, optional
            The size of the MPI-partition tiles in the y-direction (default is 30).
        max_target_grid_radius: float, optional
            The nearest neighbours search radius (default None)
        Raises
        ------
        ValueError
            If the grid type is invalid or if required data is not provided.
        """
        self.sNx = sNx
        self.sNy = sNy
        ds_has_lon = f'{self.lon_str}' in self.ds.keys()
        ds_has_lat = f'{self.lat_str}' in self.ds.keys()
        ds_has_depth = f'{self.depth_str}' in self.ds.keys()

        if (ds_has_lon) & (ds_has_lat):
            if grid_type == 'sphericalpolar':
                print(f'Dataset already has fields {self.lon_str} and {self.lat_str}. Leaving get_obs_point.')
                return
            
        if (lons is None) & (not ds_has_lon):
            raise ValueError(f"longitudes not provided")
        if (lats is None) & (not ds_has_lat):
            raise ValueError(f"latitudes not provided")
        if (depths is None) & (not ds_has_depth) & (self.pkg_str != 'obs') & (num_interp_points == 8):
            raise ValueError(f"depths not provided")
            
        if not ds_has_lon:
            if len(lons) == 0:
                lons = [lons]
            self.ds[self.lon_str] = (self.dims_obs, lons)           
        
        if not ds_has_lat:
            if len(lats) == 0:
                lats = [lats]
            self.ds[self.lat_str] = (self.dims_obs, lats)

        if (not ds_has_depth):
            if depths is not None:
              if len(depths) == 0:
                  depths = [depths]
            else:
                print('Depth not provided. Assuming surface data.')
                depths = [0] * len(lons)
            self.ds[self.depth_str] = (self.dims_obs, depths)

        if grid_type == 'sphericalpolar':
            return
            
        if grid_type not in {'llc', 'cubedsphere'}:
            raise ValueError(f"Invalid grid_type: '{grid_type}'. Options supported are 'sphericalpolar', 'llc', or 'cubedsphere'.")

        # grid_type is llc or curvilinear
        if grid_ds is not None:
            if 'XC' not in list(grid_ds.coords) or 'YC' not in list(grid_ds.coords):
                raise ValueError(f"grid_ds must have fields XC and YC.")
        else:
                raise ValueError("provide grid_ds")

        if (num_interp_points not in [1, 4]) & (self.pkg_str == 'prof'):
            raise ValueError(f"Invalid num_interp_points '{num_interp_points}'. Must be either 1 or 4 for pkg/profiles.")
        elif (num_interp_points not in [1, 4, 8]) & (self.pkg_str == 'obs'):
            raise ValueError(f"Invalid num_interp_points '{num_interp_points}'. Must be either 1, 4, or 8 for pkg/obsfit.")
            
        # add attributes relevant to interpolation field creation
        self.xc = grid_ds.XC.values
        self.yc = grid_ds.YC.values
        self.rc = grid_ds.Z.values
        self.maskc = grid_ds.maskC.values
        self.maskw = grid_ds.maskW.values
        self.masks = grid_ds.maskS.values

        # same xc, yc in worldmap form for later        
        nx = len(self.xc[0, 0, :]) # (last two dimensions of xc with shape (ntile, nx, nx)
        self.nz = len(self.rc)
        self.xc_wm = compact2worldmap(llc_tiles_to_compact(self.xc, less_output=True), nx, 1)[0, :, :]
        self.yc_wm = compact2worldmap(llc_tiles_to_compact(self.yc, less_output=True), nx, 1)[0, :, :]
        self.maskc_wm = compact2worldmap(llc_tiles_to_compact(self.maskc, less_output=True), nx, self.nz)
        self.maskw_wm = compact2worldmap(llc_tiles_to_compact(self.maskw, less_output=True), nx, self.nz)
        self.masks_wm = compact2worldmap(llc_tiles_to_compact(self.masks, less_output=True), nx, self.nz)
        
       # self.mask_wm = compact2worldmap(llc_tiles_to_compact(self.mask, less_output=True), nx, 1)[0, :, :]
        
        # set nearest neighbours search radius
        if max_target_grid_radius == None:
            max_target_grid_radius = np.max(np.abs([grid_ds.dxG.values, grid_ds.dyG.values]))
        self.max_target_grid_radius = max_target_grid_radius

        self.num_interp_points = num_interp_points
        self.nneighbours_horizontal = min(self.num_interp_points, 4)
        
        self.dims_interp = self.dims_obs + ['iINTERP']

        # run interpolation routine
        self.obs_points = self.interp() # indices in WM view 

        # add fields to ds
        self.interp_str = 'prof' if self.pkg_str == 'prof' else 'sample'
        self.interp_sfx = 'weights' if self.pkg_str == 'prof' else 'frac'
        self.obs_point_str = f'{self.interp_str}_point'

        if self.obs_points.ndim == 1:
            self.obs_points = self.obs_points[:, None]
        
        self.ds[self.obs_point_str] = (self.dims_interp, self.obs_points)

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
        valid_input_index, valid_output_index, index_array, distance_array, max_target_grid_radius =\
            get_interp_points(
                self.ds[self.lon_str],
                self.ds[self.lat_str],
                self.xc_wm.ravel(),
                self.yc_wm.ravel(),
                #nneighbours=self.num_interp_points,
                nneighbours=self.nneighbours_horizontal,
                max_target_grid_radius=self.max_target_grid_radius
            )
        
        if self.num_interp_points == 8:
            valid_input_index = np.tile(valid_input_index,2)
            valid_output_index = np.tile(valid_output_index, 2)
            index_array = np.tile(index_array, 2) # nearest grid index
            distance_array = np.tile(distance_array, 2)
            
        self.interp_distance = distance_array
        self.max_target_grid_radius = max_target_grid_radius
        
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
            self.ds[f'{self.interp_str}_interp_{tile_key}'] = (self.dims_interp, tile_field_at_obs_point)

            # save tile_fields
            self.tile_fields_wm = tile_fields_wm

        # set interp_k field in dataset
        if self.pkg_str == 'obs': 
            sample_interp_k = get_sample_interp_k(self.rc, self.ds[self.depth_str], self.num_interp_points)
            self.ds[f'{self.interp_str}_interp_k'] = (self.dims_interp, sample_interp_k)

        self.get_sample_type 
        self.get_interp_weights()

    def get_interp_weights(self):
        """ 
        Set interpolation weights
        """
        
        # find land points
        # is_land = np.zeros((self.nz, self.nneighbours_horizontal, len(self.ds.iSAMPLE)))
        # if sample_type == 3:
        #     mask = self.maskw
        # elif sample_type == 4:
        #     mask = self.masks
        # else:
        #     mask = self.maskc
        # for k in range(self.nz):
        #     is_land[k,:,:] = mask[k,:,:].ravel()[self.obs_points]
        
            
        # has size [num_interp by num_samples]
    #if (is_land.sum(axis=0) < num_interp): # for a given sample, did any of the nearest points end up on land?
        # do weighting routine with as many points as we have left
        
    
        
        
        print('Warning: currently only supported inverse distance weighting')
        inv_dist = 1 / self.interp_distance
        if inv_dist.ndim == 1:
            inv_dist = np.expand_dims(inv_dist, axis=1)

        interp_weights = inv_dist / np.sum(inv_dist, axis=1, keepdims=True)

        
        # if self.num_interp_points == 1:
        #     interp_weights = np.ones_like(self.interp_distance)
            
        # elif self.num_interp_points == 4:
        #     # bilinear interpolation
        #     lon_cur = self.ds[self.lon_str]
        #     lat_cur = self.ds[self.lat_str]
            
        #     lon_fac=(lon_cur-lon_1)/(lon_2-lon_1)
        #     lat_fac=(lat_cur-lat_1)/(lat_2-lat_1)
            
        #     interp_weights = lon_fac * lat_fac
            
        # elif self.num_interp_points == 8:
        #     # trilinear interpolation
        #     interp_weights = lon_fac * lat_fac * depth_fac
        
        # self.method = 'bilinear'
        
        # compute interp_weights from obs_points
        # convert from lat-lon to xyz
        # compute bilinear interp weights
        
        self.ds[f'{self.interp_str}_interp_{self.interp_sfx}'] = (self.dims_interp, interp_weights)

    def get_obs_datetime(self, *args, **kwargs):
        """
        Set temporal metadata
        """
        self.ds =  get_obs_datetime(self.ds,
                               *args,
                               pkg_str=self.pkg_str,
                               dims_obs=self.dims_obs,
                               **kwargs)

    def get_sample_type(self, *args, **kwargs):
        """
        Set field
        """
        self.ds =  get_sample_type(self.ds,
                               *args,
                               **kwargs)
