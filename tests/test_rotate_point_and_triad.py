import numpy as np

from ..pipefactory import rotate_point_about_point, rotate_triad

def test_rotate_point_and_triad():
    # Define the original triad
    triad = [np.array([1., 0, 0]), 
             np.array([0, 1., 0]), 
             np.array([0, 0, 1.])]
    
    # Define the rotation axis and angle
    axis = np.array([0., 0., 1.])  # Example rotation axis
    axis = axis / np.linalg.norm(axis)  # Normalize the axis
    angle = -90.0  # Example rotation angle in degrees
    
    # Define the point of origin
    p = np.array([0.0, -10., 0.0])  # Example point of origin
    x0 = np.array([0.0, 0.0, 0.0])
    
    # First rotate about [0.0, 0.0, 1.0], 90 degrees left
    x1 = rotate_point_about_point(x0, axis, -90.0, p)
    rotated_triad = rotate_triad(triad, axis, -90.0)
    
    assert np.allclose(x1, np.array([10., -10., 0.]))
    assert np.allclose(rotated_triad[0], np.array([0.0, -1.0, 0.0]))
    
    # Second rotate about [0.0, 0.0, 1.0], 90 degrees right. Triad should be back to where we were
    x2 = rotate_point_about_point(x1, axis, 90.0, x1 + np.array([10.0, 0.0, 0.0]))
    rotated_triad_next = rotate_triad(rotated_triad, axis, 90.0)
    
    assert np.allclose(x2, np.array([20.0, -20.0, 0.0]))
    assert np.allclose(rotated_triad_next[0], np.array([1.0, 0.0, 0.0]))
    assert np.allclose(rotated_triad_next[1], np.array([0.0, 1.0, 0.0]))
    assert np.allclose(rotated_triad_next[2], np.array([0.0, 0.0, 1.0]))