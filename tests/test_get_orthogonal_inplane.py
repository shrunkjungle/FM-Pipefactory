import numpy as np
import pytest
from pipefactory import get_orthogonal_inplane

def test_orthogonal_vector_xy_plane():
    d = np.array([1, 0, 0])
    expected_result = np.array([0, -1, 0])
    result = get_orthogonal_inplane(d)
    np.testing.assert_almost_equal(result, expected_result)

def test_orthogonal_vector_to_z_axis():
    d = np.array([0, 1, 0])
    expected_result = np.array([1, 0, 0])
    result = get_orthogonal_inplane(d)
    np.testing.assert_almost_equal(result, expected_result)

def test_vector_parallel_to_z_axis():
    d = np.array([0, 0, 1])
    with pytest.raises(ValueError):
        get_orthogonal_inplane(d)

def test_non_3d_vector():
    d = np.array([1, 0])
    with pytest.raises(ValueError):
        get_orthogonal_inplane(d)

