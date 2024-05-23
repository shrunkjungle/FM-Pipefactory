import numpy as np


def rot_vec(vec: np.array, 
            theta: float, 
            axis: str = "z", 
            angle_type: str = "deg"):
    """
        Rotate a 3D vector about a specified axis.

        Args:
        vec (np.array): A 3D vector to be rotated.
        theta (float): The angle of rotation.
        axis (str, optional): Axis of rotation ('x', 'y', or 'z'). Defaults to 'z'.
        angle_type (str, optional): Type of angle measurement ('deg' for degrees, 'rad' for radians). Defaults to 'deg'.

        Raises:
        Exception: If the input vector is not a 3x1 array.

        Returns:
        np.array: The rotated vector.
    """
    # Check if the input vector is a 3x1 array
    if vec.shape != (3,):
        raise Exception("Input vector must be a 3x1 array.")

    # Convert angle from degrees to radians if specified
    if angle_type == "deg":
        theta = np.deg2rad(theta)

    # Initialize a 3x3 identity matrix
    R = np.eye(3)


    if axis == "x":
        # Update the rotation matrix for x-axis rotation
        R[1][1] = np.cos(theta)
        R[1][2] = -np.sin(theta)
        R[2][1] = np.sin(theta)
        R[2][2] = np.cos(theta)

    elif axis == "y" or axis == "dr":
        # Update the rotation matrix for x-axis rotation
        R[0][0] = np.cos(theta)
        R[0][2] = np.sin(theta)
        R[2][0] = -np.sin(theta)
        R[2][2] = np.cos(theta)

    elif axis == "z" or axis == "dn":
        # Update the rotation matrix for z-axis rotation
        R[0][0] = np.cos(theta)
        R[0][1] = np.sin(theta)
        R[1][0] = -np.sin(theta)
        R[1][1] = np.cos(theta)
    else:
        # Raise an exception if axis is not recognised
        raise Exception("Rotation axis not recognised, should be x, y or z.")

    # Return the rotated vector
    return np.dot(R, vec)


def get_orthogonal_inplane(d : np.array,
                           v : np.array = np.array([0.0, 0.0, 1.0])):
    """
    Get a unit vector that is orthogonal to the given vector 'd' and another vector 'v'.

    Args:
    d (np.array): A 3D vector.

    Raises:
    ValueError: If the input vector 'd' is not a 3-dimensional vector.

    Returns:
    np.array: A unit vector that is orthogonal to 'd' and lies in the XY plane.
    """
    # Check if the input vector is a 3-dimensional vector
    if (d.shape != (3,) or v.shape != (3,)):
        raise ValueError("Input vector 'd' must be a 3-dimensional numpy array.")

    # Calculate the cross product of 'd' and the unit vector along the z-axis
    y = np.cross(d, v)

    # Check if the norm of the vector is zero (to avoid division by zero)
    norm_y = np.linalg.norm(y)
    if norm_y == 0:
        raise ValueError("The vector 'd' is parallel to the z-axis, resulting in a zero orthogonal vector.")

    # Normalize the resulting vector to get a unit vector
    return y / norm_y


def find_section(section_ends: list, x: float):
    """
    Find the index of the section where the given value 'x' falls based on the section end points.

    Args:
    section_ends (list): A list of end points of the sections, expected to be in sorted order.
    x (float): The value to locate within the sections.

    Returns:
    int: The index of the section where 'x' falls. If 'x' does not fall in any section, returns None.
    """
    for i, element in enumerate(section_ends):
        if x <= element + 1e-4:
            return i
    return None

def rotate_vector(v: np.array, axis: np.array, angle: float) -> np.array:
    """
    Rotate a vector about an arbitrary axis using Rodrigues' rotation formula.

    Parameters:
        v (np.array): The vector to be rotated.
        axis (np.array): The axis of rotation (unit vector).
        angle (float): The angle of rotation in degrees.

    Returns:
        np.array: The rotated vector.
    """
    # Convert angle from degrees to radians
    angle = np.deg2rad(angle)
    
    # Calculate cosine and sine of the angle
    cos_theta = np.cos(angle)
    sin_theta = np.sin(angle)
    
    # Apply Rodrigues' rotation formula
    v_rotated = (v * cos_theta +
                 np.cross(axis, v) * sin_theta +
                 axis * np.dot(axis, v) * (1 - cos_theta))

    return v_rotated



def rotate_triad(triad: list[np.array], 
                 axis: np.array, 
                 angle: float):
    """
    Rotate a triad of vectors around a given axis by a specified angle while considering a given point as the new origin.

    Parameters:
        triad (list[np.array]): A list of numpy arrays representing the original triad.
        axis (np.array): The axis of rotation.
        angle (float): The rotation angle in degrees.

    Returns:
        list[np.array]: The rotated triad.
    """
    
    # Rotate each vector in the translated triad using the rotate_vector function
    rotated_triad = [rotate_vector(v, axis, angle) for v in triad]

    return rotated_triad


def rotate_point_about_point(x0: np.array, 
                             v: np.array, 
                             angle: float, 
                             p: np.array) -> np.array:
    """
    Rotate a point around another point using a given axis and angle.

    Parameters:
        x0 (np.array): The point to be rotated.
        v (np.array): The rotation axis.
        angle (float): The rotation angle in degrees.
        p (np.array): The point about which the rotation occurs.

    Returns:
        np.array: The rotated point.
    """
    # Translate point x0 to be relative to point p
    translated_point = x0 - p

    # Rotate the translated point around axis v
    cos_theta = np.cos(np.radians(angle))
    sin_theta = np.sin(np.radians(angle))
    v_unit = v / np.linalg.norm(v)
    rotation_matrix = np.eye(3) * cos_theta + \
                      np.outer(v_unit, v_unit) * (1 - cos_theta) + \
                      np.array([[0, -v_unit[2], v_unit[1]],
                                [v_unit[2], 0, -v_unit[0]],
                                [-v_unit[1], v_unit[0], 0]]) * sin_theta
    rotated_point = np.dot(rotation_matrix, translated_point)

    # Translate the rotated point back to its original position relative to point p
    rotated_point += p

    return rotated_point


def cylinder_geodesic_distance(s0 : float, 
                               theta0 : float, 
                               s1 : float,
                               theta1 : float,
                               radius : float):
    
        """
        Calculate the geodesic distance on the surface of a cylinder between two points.

        The surface of the cylinder is assumed to be unwrapped into a 2D plane where one axis (s) corresponds 
        to the axial direction of the cylinder, and the other axis (theta) corresponds to the angular position 
        around the cylinder's circumference. The geodesic distance is the shortest path between two points on 
        this surface, taking into account the periodic boundary conditions in the angular direction.

        Parameters:
        - s0 (float): The axial coordinate of the first point.
        - theta0 (float): The angular coordinate (in radians) of the first point, where 0 <= theta0 < 2*pi.
        - s1 (float): The axial coordinate of the second point.
        - theta1 (float): The angular coordinate (in radians) of the second point, where 0 <= theta1 < 2*pi.
        - radius (float): The radius of the cylinder.

        Returns:
        - float: The geodesic distance between the two points on the cylinder's surface.

        
        """
    
        # Calculate geodesic distance
        delta_x = s0 - s1
        delta_theta = theta0 - theta1
        
        # Account for periodicity of angles
        delta_theta = min(abs(delta_theta), 2*np.pi - abs(delta_theta))
        
        return np.sqrt(delta_x**2 + radius**2 * delta_theta**2)




