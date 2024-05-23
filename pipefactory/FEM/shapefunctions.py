import numpy as np


class ShapeFunction():
    """
        Base Class for shape functions
    """
    def __init__(self, 
                 type : str, 
                 dof : int, 
                 dim : int):

        self.type = type
        self.dof = dof
        self.dim = dim

    def N(self, local_coords : np.array):
        return NotImplementedError
    
    def dNdX(self, local_coords : np.array):
        return NotImplementedError
    

class LinearQuad(ShapeFunction):

    """
        Implementation of Shape Function
    """

    def __init__(self):

        super().__init__("Quad4", 4 , 2) 

    def N(self, local_coords : np.array):
        """
        Compute linear shape functions for a 4-node quadrilateral element.

        :param xi: Natural coordinate in the xi-direction.
        :param eta: Natural coordinate in the eta-direction.
        :return: Array of shape functions.
        """

        xi, eta, zeta = local_coords
        N = np.zeros(4)
        N[0] = 0.25 * (1 - xi) * (1 - eta)  # N1
        N[1] = 0.25 * (1 + xi) * (1 - eta)  # N2
        N[2] = 0.25 * (1 + xi) * (1 + eta)  # N3
        N[3] = 0.25 * (1 - xi) * (1 + eta)  # N4
        return N

    def dNdX(self, local_coords : np.array):
        """
        Compute derivatives of linear shape functions for a 4-node quadrilateral element.

        :param xi: Natural coordinate in the xi-direction.
        :param eta: Natural coordinate in the eta-direction.
        :return: Matrix of shape function derivatives. Each row corresponds to a node, 
                the first column is dN/dxi, and the second column is dN/deta.
        """

        xi, eta, zeta = local_coords
        dN = np.zeros((4, 3))
        dN[0, :] = [-0.25 * (1 - eta), -0.25 * (1 - xi), 0.0]  # dN1/dxi, dN1/deta
        dN[1, :] = [0.25 * (1 - eta), -0.25 * (1 + xi), 0.0]   # dN2/dxi, dN2/deta
        dN[2, :] = [0.25 * (1 + eta), 0.25 * (1 + xi), 0.0]    # dN3/dxi, dN3/deta
        dN[3, :] = [-0.25 * (1 + eta), 0.25 * (1 - xi), 0.0]   # dN4/dxi, dN4/deta
        
        return dN

class LinearHex(ShapeFunction):
    """
    Implementation of Shape Function for an 8-node hexahedral element.
    """

    def __init__(self):
        super().__init__("Hex8", 8, 3)

    def N(self, local_coords : np.array):
        """
        Compute linear shape functions for an 8-node hexahedral element.

        :param xi: Natural coordinate in the xi-direction.
        :param eta: Natural coordinate in the eta-direction.
        :param zeta: Natural coordinate in the zeta-direction.
        :return: Array of shape functions.
        """
        xi, eta, zeta = local_coords

        N = np.zeros(8)
        N[0] = 0.125 * (1 - xi) * (1 - eta) * (1 - zeta)  # N1
        N[1] = 0.125 * (1 + xi) * (1 - eta) * (1 - zeta)  # N2
        N[2] = 0.125 * (1 + xi) * (1 + eta) * (1 - zeta)  # N3
        N[3] = 0.125 * (1 - xi) * (1 + eta) * (1 - zeta)  # N4
        N[4] = 0.125 * (1 - xi) * (1 - eta) * (1 + zeta)  # N5
        N[5] = 0.125 * (1 + xi) * (1 - eta) * (1 + zeta)  # N6
        N[6] = 0.125 * (1 + xi) * (1 + eta) * (1 + zeta)  # N7
        N[7] = 0.125 * (1 - xi) * (1 + eta) * (1 + zeta)  # N8
        return N
    

    def dNdX(self, local_coords : np.array):
        """
        Compute derivatives of linear shape functions for an 8-node hexahedral element.

        :param xi: Natural coordinate in the xi-direction.
        :param eta: Natural coordinate in the eta-direction.
        :param zeta: Natural coordinate in the zeta-direction.
        :return: Matrix of shape function derivatives. Each row corresponds to a node,
                 and columns correspond to derivatives with respect to xi, eta, and zeta.
        """

        xi, eta, zeta = local_coords

        dN = np.zeros((8, 3))
        dN[0, :] = [-0.125 * (1 - eta) * (1 - zeta), -0.125 * (1 - xi) * (1 - zeta), -0.125 * (1 - xi) * (1 - eta)] # Checked
        dN[1, :] = [0.125 * (1 - eta) * (1 - zeta), -0.125 * (1 + xi) * (1 - zeta), -0.125 * (1 + xi) * (1 - eta)] # Checked
        dN[2, :] = [0.125 * (1 + eta) * (1 - zeta), 0.125 * (1 + xi) * (1 - zeta), -0.125 * (1 + xi) * (1 + eta)] # Checked
        dN[3, :] = [-0.125 * (1 + eta) * (1 - zeta), 0.125 * (1 - xi) * (1 - zeta), -0.125 * (1 - xi) * (1 + eta)] # Checked
        dN[4, :] = [-0.125 * (1 - eta) * (1 + zeta), -0.125 * (1 - xi) * (1 + zeta), 0.125 * (1 - xi) * (1 - eta)] # Checked
        dN[5, :] = [0.125 * (1 - eta) * (1 + zeta), -0.125 * (1 + xi) * (1 + zeta), 0.125 * (1 + xi) * (1 - eta)] # Checked
        dN[6, :] = [0.125 * (1 + eta) * (1 + zeta), 0.125 * (1 + xi) * (1 + zeta), 0.125 * (1 + xi) * (1 + eta)] # Checked
        dN[7, :] = [-0.125 * (1 + eta) * (1 + zeta),0.125 * (1 - xi) * (1 + zeta), 0.125 * (1 - xi) * (1 + eta)] # Checked

        return dN
