# Ensure the input DataArray has dimensions 'k', 'tile', 'j', 'i'
assert 'tile' in da.dims, "The DataArray must have a 'tile' dimension."

# Handle j/i dimension names, check if 'j_g' or 'i_g' are used
j_dim = 'j' if 'j' in da.dims else 'j_g'
i_dim = 'i' if 'i' in da.dims else 'i_g'

nx = da.sizes[i_dim]
other_dims = [dim for dim in da.dims if dim not in ['tile', j_dim, i_dim]]

# Initialize the output DataArray as a 4x4 patch, using the same coordinates as the input
coords = {dim: da.coords[dim] for dim in other_dims}
coords[j_dim] = np.arange(4 * nx)
coords[i_dim] = np.arange(4 * nx)

da_wm = xr.DataArray(np.nan, dims=other_dims + [j_dim, i_dim], coords=coords)

# face 1
face1_stack = xr.concat(
    [da.sel(tile=i).assign_coords(j=np.arange(i * nx, (i + 1) * nx)) for i in [0, 1, 2]],
    dim='j'
)
da_wm.loc[{j_dim: slice(0, 3 * nx - 1), i_dim: slice(0, nx - 1)}] = face1_stack

# face 2
face2_stack = xr.concat(
    [da.sel(tile=i).assign_coords(**{j_dim: np.arange((i - 3) * nx, (i - 2) * nx)},
                                  **{i_dim: np.arange(nx, 2 * nx)}) for i in [3, 4, 5]],
    dim=j_dim
)
da_wm.loc[{j_dim: slice(0, 3 * nx - 1), i_dim: slice(nx, 2 * nx - 1)}] = face2_stack

# face 4
# Rotate 
face4_stack = xr.concat(
    [da.sel(tile=i).transpose(*other_dims, i_dim, j_dim).isel({i_dim: slice(None, None, -1)})\
                   .assign_coords(
#                      **{j_dim: np.arange((i - 6) * nx, (i - 5) * nx),
#                      i_dim: np.arange(2 * nx, 3 * nx)}
                      **{i_dim: np.arange((i - 7) * nx, (i - 6) * nx)[::-1],
                         j_dim: np.arange(2 * nx, 3 * nx)}
                   ) for i in [9, 8, 7]],
    dim = i_dim
)

#face4_stack = face4_stack.transpose(..., i_dim, j_dim).isel({i_dim:slice(None, None, -1)})

face4_stack = face4_stack.rename({'i': 'temp_i', 'j': 'temp_j'})

face4_stack = face4_stack.isel({'temp_i':slice(None, None, -1)})
face4_stack = face4_stack.assign_coords({'temp_i':np.arange(0, 3*nx), 'temp_j':np.arange(2*nx, 3*nx)})
face4_stack = face4_stack.rename({'temp_i': 'i', 'temp_j': 'j'})



#face4_stack[0][0].plot()