import pytest
import numpy as np
import xarray as xr
from obsprep.prep import Prep

@pytest.fixture
def mock_dataset():
    """Create a mock xarray dataset for testing."""
    dataset_llc = xr.Dataset(
        {
            'data_variable': (('tile','j', 'i'), np.random.rand(13, 90, 90))
        },
        coords={
            'XC': (('tile','j', 'i'), np.random.rand(13, 90, 90)),
            'YC': (('tile','j', 'i'), np.random.rand(13, 90, 90)),
            'Z': (('k'), np.random.rand(50)),
            'dxG': (('tile','j', 'i'), np.random.rand(13, 90, 90)),
            'dyG': (('tile','j', 'i'), np.random.rand(13, 90, 90)),
            'hFacC': (('k', 'tile','j', 'i'), np.random.rand(50, 13, 90, 90)),
        }
    )
    return dataset_llc

def test_init_invalid_pkg():
    """Test initialization with an invalid pkg parameter."""
    with pytest.raises(ValueError, match="Invalid pkg 'invalid_pkg'"):
        Prep(pkg='invalid_pkg')

def test_init_valid_pkg(mock_dataset):
    """Test initialization with a valid pkg parameter."""
    prep = Prep(pkg='profiles', ds=mock_dataset)
    assert prep.pkg == 'profiles'

def test_get_pkg_fields_profiles():
    """Test package-specific fields for 'profiles'."""
    prep = Prep(pkg='profiles')
    prep.get_pkg_fields()
    assert prep.dims_obs == ['iPROF']
    assert prep.lon_str == 'prof_lon'
    assert prep.lat_str == 'prof_lat'

def test_get_pkg_fields_obsfit():
    """Test package-specific fields for 'obsfit'."""
    prep = Prep(pkg='obsfit')
    prep.get_pkg_fields()
    assert prep.dims_obs == ['iSAMPLE']
    assert prep.lon_str == 'sample_lon'
    assert prep.lat_str == 'sample_lat'

def test_get_obs_point_with_valid_coords(mock_dataset):
    """Test get_obs_point with valid coordinates."""
    prep = Prep(pkg='profiles')
    prep.get_obs_point([0.5], [0.5], grid_ds=mock_dataset)

    # Check if expected fields are created in dataset
    assert 'prof_lon' in prep.ds
    assert 'prof_lat' in prep.ds

def test_get_obs_point_llc_coords(mock_dataset):
    prep = Prep(pkg='profiles')
    prep.get_obs_point([0.5], [0.5], grid_type='llc', grid_ds=mock_dataset)

    # Check if expected fields are created in dataset
    assert prep.ds['prof_point'].values.shape == (1,1)
    assert 'prof_interp_XC11' in prep.ds
    assert 'prof_interp_YC11' in prep.ds
