import numpy as np
from xmitgcm.llcreader.llcmodel import faces_dataset_to_latlon 

def generate_random_points(nobs, lon_range=(-180, 180), lat_range=(-90, 90), depth_range=(0, 500)):
    lons = np.random.uniform(low=lon_range[0], high=lon_range[1], size=nobs)
    lats = np.random.uniform(low=lat_range[0], high=lat_range[1], size=nobs)
    depths = np.random.uniform(low=depth_range[0], high=depth_range[1], size=nobs)
    return lons, lats, depths

def faces_dataset_to_latlon_pole(ds):
    """
    Convert a 13-faces llc dataset to a latitude-longitude dataset format and handle polar transformations.


    Parameters
    ----------
    ds : xarray.Dataset
        The input dataset contains fields with 13-face horizontal spatial dimensions.
        This dataset is assumed to have coordinate dimensions 'i', 'j', 'i_g', and 'j_g'

    Returns
    -------
    xarray.Dataset
        A new dataset with coordinates transformed to a latitude-longitude coordinate system. The polar face is
        rotated and placed above Europe.

    """

    nx = len(ds.i)

    # use xmitgcm's faces_dataset_to_latlon, which does everything except format the pole
    ds_latlon = faces_dataset_to_latlon(ds)

    # pad above so that the fields in the dataset are shape (4*nx, 4*nx)
    ds_latlon = ds_latlon.pad(j=(0, nx), j_g=(0, nx), constant_values=0).sortby(["j", "j_g"])
    ds_latlon = ds_latlon.assign_coords(j=np.arange(4*nx), j_g=np.arange(4*nx))
    target_slice = ds_latlon.isel(j=slice(3*nx, None), i=slice(None, nx),
                                 j_g=slice(3*nx, None), i_g=slice(None, nx))

    # Get pole subset of dataset
    pole = ds.isel(face=6)
    
    # Reassign horizontal coordinates
    pole = pole.assign_coords({
        'i': np.arange(nx),
        'j': np.arange(3*nx, 4*nx),
        'i_g': np.arange(nx),
        'j_g': np.arange(3*nx, 4*nx),
    })

    # loop through coordinates, adding pole and rotating
    dims = pole.dims
    
    for coord in pole.coords:
        if coord in pole.dims:
            continue
        coord_dims = pole[coord].dims
    
        if 'i' in pole[coord].dims or 'j' in pole[coord].dims or \
           'i_g' in pole[coord].dims or 'j_g' in pole[coord].dims:

            i_dim = next((dim for dim in coord_dims if 'i' in dim.lower()), None)
            j_dim = next((dim for dim in coord_dims if 'j' in dim.lower()), None)
            extra_dims = [dim for dim in dims if dim not in {j_dim, i_dim}]

            # reassign coords
            target_slice = target_slice.assign_coords({coord: pole[coord]})
    
            # rotate the pole: transpose, then flip along j_dim, leaving other dimensions intact
            target_slice = target_slice.transpose(*extra_dims, i_dim, j_dim).isel({j_dim: slice(None, None, -1)})
            target_slice = target_slice.rename({j_dim: i_dim, i_dim: j_dim})
            target_slice = target_slice.assign_coords({i_dim: pole.i, j_dim: pole.j})
        else:
            continue

        # store pole coord in ds_latlon[coord]
        ds_latlon[coord].loc[dict({j_dim: slice(3*nx, None), i_dim:slice(None, nx-1)})] = target_slice[coord]

    return ds_latlon


def patchface3D_5f_to_wrld(faces):
    # Extract dimensions from one of the faces
    nz, _, nx = faces[1].shape

    array_out = np.zeros((nz, 4*nx, 4*nx))

    # Face 1
    array_out[:, :3*nx, :nx] = faces[1]

    # Face 2
    array_out[:, :3*nx, nx:2*nx] = faces[2]

    # Face 4
    face4 = np.transpose(np.flip(faces[4], 2), (0,2,1))
    array_out[:, :3*nx, 2*nx:3*nx] = face4

    # Face 5
    face5 = np.transpose(np.flip(faces[5], 2), (0,2,1))
    array_out[:, :3*nx, 3*nx:4*nx] = face5

    # Face 3
    face3 = np.rot90(faces[3][0,:,:], 3)
    array_out[:, 3*nx:4*nx, :nx] = face3

    return array_out

def patchface3D_wrld_to_5f(array_in):
    nz, ny, nx = np.shape(array_in)
    nx = int(nx/4)
    faces = dict()

    # face 1
    faces[1] = array_in[:,:3*nx,:nx]

    # face 2
    faces[2] = array_in[:,:3*nx,nx:2*nx]

    # face 4
    face4 = array_in[:,:3*nx,2*nx:3*nx]
    faces[4] = sym_g_mod(face4, 5)

    # face 5
    face5 = array_in[:,:3*nx,3*nx:4*nx]
    faces[5] = sym_g_mod(face5, 5)

    # face 3
    face3 = array_in[:,3*nx:4*nx,:nx]
    faces[3] = sym_g_mod(face3, 7)

    return faces

def sym_g_mod(field_in, sym_in):
    field_out = field_in
    for icur in range(sym_in-4):
        field_out = np.flip(np.transpose(field_out, (0, 2, 1)), axis=2)

    return field_out

def compact2worldmap(fldin,nx,nz):
    #add a new dimension in case it's only 2d field:
    if nz == 1:
        fldin=fldin[np.newaxis, :, :]
    #defining a big face:
    a=np.zeros((nz,4*nx,4*nx))       #(50,270,360)
    #face1
    tmp=fldin[:,0:3*nx,0:nx]        #(50,270,90)
    a[:,0:3*nx,0:nx]=tmp
    #face2
    tmp=fldin[:,(3*nx):(6*nx),0:nx] #(50, 270,90)
    a[:,0:3*nx,nx:2*nx]=tmp
    #face3
    tmp=fldin[:,(6*nx):(7*nx),0:nx] #(50, 90, 90)
    tmp=np.transpose(tmp, (1,2,0))  #(90, 90, 50)
    ##syntax to rotate ccw:
    tmp1=list(zip(*tmp[::-1]))
    tmp1=np.asarray(tmp1)
    tmp1=np.transpose(tmp1,[2,0,1]) #(50, 90, 90)
    a[:,3*nx:4*nx,0:nx]=tmp1
    #face4
    tmp=np.reshape(fldin[:,7*nx:10*nx,0:nx],[nz,nx,3*nx]) #(50,90,270)
    tmp=np.transpose(tmp, (1,2,0))
    #syntax to rotate cw:
    tmp1=list(zip(*tmp))[::-1]      #type is <class 'list'>
    tmp1=np.asarray(tmp1)           #type <class 'numpy.ndarray'>, shape (270,90,50)
    tmp1=np.transpose(tmp1,[2,0,1]) #(50,270,90)
    a[:,0:3*nx,2*nx:3*nx]=tmp1
    #face5
    tmp=np.reshape(fldin[:,10*nx:13*nx,0:nx],[nz,nx,3*nx]) #(50,90,270)
    tmp=np.transpose(tmp, (1,2,0))                         #(90,270,50)
    tmp1=list(zip(*tmp))[::-1]      #type is <class 'zip'> --> <class 'list'>
    tmp1=np.asarray(tmp1)           #type <class 'numpy.ndarray'>, shape (270,90,50)
    tmp1=np.transpose(tmp1,[2,0,1]) #(50,270,90)
    a[:,0:3*nx,3*nx:4*nx]=tmp1
    return a

def patchface3D(fldin,nx,nz):

    print(nz)
    #add a new dimension in case it's only 2d field:
    if nz == 1:
        fldin=fldin[np.newaxis, :, :]

    #defining a big face:
    a=np.zeros((nz,4*nx,4*nx))       #(50,270,360)

    #face1
    tmp=fldin[:,0:3*nx,0:nx]        #(50,270,90)
    a[:,0:3*nx,0:nx]=tmp

    #face2
    tmp=fldin[:,(3*nx):(6*nx),0:nx] #(50, 270,90)
    a[:,0:3*nx,nx:2*nx]=tmp

    #face3
    tmp=fldin[:,(6*nx):(7*nx),0:nx] #(50, 90, 90)
    tmp=np.transpose(tmp, (1,2,0))  #(90, 90, 50)
    ##syntax to rotate ccw:
    tmp1=list(zip(*tmp[::-1]))
    tmp1=np.asarray(tmp1)
    tmp1=np.transpose(tmp1,[2,0,1]) #(50, 90, 90)
    a[:,3*nx:4*nx,0:nx]=tmp1

    #face4
    tmp=np.reshape(fldin[:,7*nx:10*nx,0:nx],[nz,nx,3*nx]) #(50,90,270)
    tmp=np.transpose(tmp, (1,2,0))
    print(tmp.shape)                                      #(90,270,50)
    #syntax to rotate cw:
    tmp1=list(zip(*tmp))[::-1]      #type is <class 'list'>
    tmp1=np.asarray(tmp1)           #type <class 'numpy.ndarray'>, shape (270,90,50)
    tmp1=np.transpose(tmp1,[2,0,1]) #(50,270,90)
    a[:,0:3*nx,2*nx:3*nx]=tmp1

    #face5
    tmp=np.reshape(fldin[:,10*nx:13*nx,0:nx],[nz,nx,3*nx]) #(50,90,270)
    tmp=np.transpose(tmp, (1,2,0))                         #(90,270,50)
    tmp1=list(zip(*tmp))[::-1]      #type is <class 'zip'> --> <class 'list'>
    tmp1=np.asarray(tmp1)           #type <class 'numpy.ndarray'>, shape (270,90,50)
    tmp1=np.transpose(tmp1,[2,0,1]) #(50,270,90)
    a[:,0:3*nx,3*nx:4*nx]=tmp1

    return a

def get_sample_type(ds, fld):

    # Get the factor using the dictionary, default to 0 if fld is not found
    fac = field_map.get(fld, 0)

    # Create the sample_type variable
    sample_type = (fac * np.ones_like(ds.sample_lat.values)).astype(int)
    ds['sample_type'] = ("iSAMPLE", sample_type) 

    return ds


grid_map = {
    'C': {'grid_lon': 'XC', 'grid_lat': 'YC', 'mask': 'maskC'},
    'W': {'grid_lon': 'XG', 'grid_lat': 'YC', 'mask': 'maskW'},
    'S': {'grid_lon': 'XC', 'grid_lat': 'YG', 'mask': 'maskS'},
}

field_data = {
    'T': {'sample_type': 1, 'mask_sfx': 'C'},
    'S': {'sample_type': 2, 'mask_sfx': 'C'},
    'U': {'sample_type': 3, 'mask_sfx': 'W'},
    'V': {'sample_type': 4, 'mask_sfx': 'S'},
    'SSH': {'sample_type': 5, 'mask_sfx': 'C'}
}

field_map = {
    field: {
        'sample_type': data['sample_type'],
        'mask': grid_map[data['mask_sfx']]['mask'],
        'grid_lon': grid_map[data['mask_sfx']]['grid_lon'],
        'grid_lat': grid_map[data['mask_sfx']]['grid_lat']
    }
    for field, data in field_data.items()
}

def get_mask_sfx_from_sample_type(sample_types):
    sample_type_to_mask_sfx = {data['sample_type']: data['mask_sfx'] for data in field_data.values()}

    # Retrieve corresponding sample codes
    mask_sfxs = [sample_type_to_mask_sfx[stype] for stype in sample_types]
    return mask_sfxs