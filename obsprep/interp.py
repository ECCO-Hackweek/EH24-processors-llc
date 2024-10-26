import pyresample as pr
import numpy as np

def get_interp_points(lons, lats,
                      grid_lon_wm_flat, grid_lat_wm_flat, nneighbours=4,
                      max_target_grid_radius = int(15e4), max_attempts=10, radius_factor=1.5,
                     ):
    """
    Get interpolation points for ungridded longitude and latitude based on a source grid.

    Parameters
    ----------
    lons : array-like
        Longitudes of ungridded points.
    lats : array-like
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
        lats=lats, lons=lons
    )
    
    source_grid = pr.geometry.SwathDefinition(
        lons=grid_lon_wm_flat, lats=grid_lat_wm_flat
    )

    invalid_search = True
    attempt = 0

    while (invalid_search) & (attempt < max_attempts):
      valid_input_index, valid_output_index, index_array, distance_array = pr.kd_tree.get_neighbour_info(
          source_grid,
          target_grid,
          radius_of_influence=int(max_target_grid_radius),
          neighbours=nneighbours,
      )

      # invalid interp points show up as having inf distance
      invalid_search = np.any(~np.isfinite(distance_array))
      max_target_grid_radius = max_target_grid_radius * radius_factor
      print(f'Found invalid interp point. Increasing max_target_grid_radius by factor {radius_factor} and retrying kdtree')

    if attempt == max_attempts:
        raise ValueError(f"Couldn't find valid interp points in {max_attempts} attempts")

    # edge case: kdtree can return grid_lon_wm_flat.size as a valid index
    index_array[index_array == grid_lon_wm_flat.size] -= 1

    return (valid_input_index, valid_output_index, index_array, distance_array, max_target_grid_radius)
    


def get_depth_indices(rC, depth_cur):

    """
    Get vertical interpolation levels for ungridded depths based on a source grid.

    Parameters
    ----------
    rC : array-like
        depths of vertical layers in model (negative values)
    depth_cur : array-like
        depths of ungridded points (positive values)
  
    Returns
    -------
    tuple
        A tuple of (sample_k1, sample_k2, depth_fac)
    """

    sample_k1 = np.zeros_like(depth_cur)#, dtype=int)
    sample_k2 = np.zeros_like(depth_cur)#, dtype=int)
    depth_fac = np.zeros_like(depth_cur)#, dtype=float)
    
    # Case 1: above first depth level
    mask_above_first = -rC[0] > depth_cur
    sample_k1[mask_above_first] = 1
    sample_k2[mask_above_first] = 1
    depth_fac[mask_above_first] = 1.0
    
    # Case 2: below last depth level
    mask_below_last = -rC[Nr-1] <= depth_cur
    sample_k1[mask_below_last] = Nr
    sample_k2[mask_below_last] = Nr
    depth_fac[mask_below_last] = 1.0
    
    # Case 3: between two depth levels
    mask_between = ~mask_above_first & ~mask_below_last
    
    for k in range(Nr-1):
        mask_between_k = (-rC[k] <= depth_cur) & (-rC[k+1] > depth_cur) & mask_between
        sample_k1[mask_between_k] = k+1
        sample_k2[mask_between_k] = k+2
        depth_1 = -rC[k]
        depth_2 = -rC[k+1]
        depth_fac[mask_between_k] = (depth_cur[mask_between_k] - depth_1) / (depth_2 - depth_1)

    return (sample_k1, sample_k2, depth_fac)
