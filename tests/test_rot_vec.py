import numpy as np
import pytest as py
from ..pipefactory import rot_vec

def test_rotation_90_degrees():
    vec = np.array([1, 0, 0])
    theta = 90
    expected_result = np.array([0, -1., 0])
    result = rot_vec(vec, theta)
    np.testing.assert_almost_equal(result, expected_result)

def test_rotation_negative_90_degrees():
    vec = np.array([1., 0, 0])
    theta = -90
    expected_result = np.array([0, 1., 0])
    result = rot_vec(vec, theta)
    np.testing.assert_almost_equal(result, expected_result)

def test_rotation_negative_90_degrees_y_axis():
    vec = np.array([1., 0., 0.])
    expected_results = np.array([0.0, 0.0, 1.0])
    results = rot_vec(vec, -90., "y")
    np.testing.assert_almost_equal(expected_results, results)

def test_rotation_180_degrees():
    vec = np.array([1, 0, 0])
    theta = 180
    expected_result = np.array([-1, 0, 0])
    result = rot_vec(vec, theta)
    np.testing.assert_almost_equal(result, expected_result)

def test_rotation_negative_180_degrees():
    vec = np.array([1, 0, 0])
    theta = -180
    expected_result = np.array([-1, 0, 0])  # Rotation by -180 degrees should give the same result as 180 degrees
    result = rot_vec(vec, theta)
    np.testing.assert_almost_equal(result, expected_result)

def test_invalid_vector_length():
    vec = np.array([1, 0])
    theta = 90
    with py.raises(Exception):
        rot_vec(vec, theta)
