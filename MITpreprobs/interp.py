import pyresample as pr

def get_interp_points(ungridded_lons, ungridded_lats,
                      grid_lon_wm_flat, grid_lat_wm_flat, nneighbours=4,
                      max_target_grid_radius = int(15e4),
                     ):
    """
    Get interpolation points for ungridded longitude and latitude based on a source grid.

    Parameters
    ----------
    ungridded_lons : array-like
        Longitudes of ungridded points.
    ungridded_lats : array-like
        Latitudes of ungridded points.
    grid_lon_wm_flat : array-like
        Longitudes of the source grid (flattened).
    grid_lat_wm_flat : array-like
        Latitudes of the source grid (flattened).
    nneighbours : int, optional
        Number of nearest neighbours to find (default is 4).
    max_target_grid_radius : int, optional
        Maximum radius of influence for neighbours (default is 150000).

    Returns
    -------
    tuple
        A tuple of (valid_input_index, valid_output_index, index_array, distance_array).
    """
    target_grid = pr.geometry.SwathDefinition(
        lats=ungridded_lats, lons=ungridded_lons
    )
    
    source_grid = pr.geometry.SwathDefinition(
        lons=grid_lon_wm_flat, lats=grid_lat_wm_flat
    )
    
    valid_input_index, valid_output_index, index_array, distance_array = pr.kd_tree.get_neighbour_info(
        source_grid,
        target_grid,
        radius_of_influence=int(max_target_grid_radius),
        neighbours=nneighbours,
    )

    # edge case: kdtree can return grid_lon_wm_flat.size as a valid index
    index_array[index_array == grid_lon_wm_flat.size] -= 1
    return (valid_input_index, valid_output_index, index_array, distance_array)
    
