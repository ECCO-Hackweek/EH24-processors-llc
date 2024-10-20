import pytest
import numpy as np
import xarray as xr
from MITpreprobs.preproc import UngriddedObsPreprocessor

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
        }
    )
    return dataset_llc

def test_init_invalid_pkg():
    """Test initialization with an invalid pkg parameter."""
    with pytest.raises(ValueError, match="Invalid pkg 'invalid_pkg'"):
        UngriddedObsPreprocessor(pkg='invalid_pkg')

def test_init_valid_pkg(mock_dataset):
    """Test initialization with a valid pkg parameter."""
    preprocessor = UngriddedObsPreprocessor(pkg='profiles', ds=mock_dataset)
    assert preprocessor.pkg == 'profiles'
    assert preprocessor.grid == 'sphericalpolar'
    assert preprocessor.sNx == 30
    assert preprocessor.sNy == 30

def test_get_pkg_fields_profiles():
    """Test package-specific fields for 'profiles'."""
    preprocessor = UngriddedObsPreprocessor(pkg='profiles')
    preprocessor.get_pkg_fields()
    assert preprocessor.dims_obs == ['iPROF']
    assert preprocessor.lon_str == 'prof_lon'
    assert preprocessor.lat_str == 'prof_lat'

def test_get_pkg_fields_obsfit():
    """Test package-specific fields for 'obsfit'."""
    preprocessor = UngriddedObsPreprocessor(pkg='obsfit')
    preprocessor.get_pkg_fields()
    assert preprocessor.dims_obs == ['iOBS']
    assert preprocessor.lon_str == 'sample_lon'
    assert preprocessor.lat_str == 'sample_lat'

def test_get_obs_point_with_valid_coords(mock_dataset):
    """Test get_obs_point with valid coordinates."""
    preprocessor = UngriddedObsPreprocessor(pkg='profiles', ds=mock_dataset)
    preprocessor.get_obs_point(mock_dataset, [0.5], [0.5])
    assert preprocessor.ds['prof_point'].values.shape == (1,)
    assert 'prof_lon' in preprocessor.ds
    assert 'prof_lat' in preprocessor.ds

def test_get_sample_interp_info(mock_dataset):
    """Test get_sample_interp_info method."""
    preprocessor = UngriddedObsPreprocessor(pkg='profiles', ds=mock_dataset, grid='llc')
    preprocessor.get_obs_point(mock_dataset, [0, 1], [0, 1])
    preprocessor.get_sample_interp_info()

    # Check if expected fields are created in dataset
    assert 'prof_interp_XC11' in preprocessor.ds
    assert 'prof_interp_YC11' in preprocessor.ds


