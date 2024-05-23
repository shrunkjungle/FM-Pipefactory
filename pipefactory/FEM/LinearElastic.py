import numpy as np
from pipefactory import LinearHex, GaussQuadrature, GaussQuadrature2D


class LinearElasticity:


    def __init__(self, reduced : bool = False):


        self.LinearHex = LinearHex()

        self.dofel = 24

        self.index_u = 3 * np.arange(8)
        self.index_v = 3 * np.arange(8) + 1
        self.index_w = 3 * np.arange(8) + 2

        self.gauss = GaussQuadrature(reduced)

        self.gauss_2d = GaussQuadrature2D(reduced)

    def getMaterialTensor(self, 
                          E : float,
                          nu : float):
        """
        Calculate the isotropic elasticity tensor for a material based on its Young's modulus and Poisson's ratio.

        This function constructs and returns the elasticity tensor (often denoted as 'D') for an isotropic material. 
        The tensor is constructed using the material's Young's modulus (E) and Poisson's ratio (nu). The elasticity 
        tensor is essential in the field of continuum mechanics and material science, particularly for stress-strain 
        analysis in isotropic materials.

        Parameters:
        - E (float): Young's modulus of the material. This is a measure of the material's stiffness, defined as the ratio 
        of stress (force per unit area) to strain (proportional deformation) in the linear elasticity regime of a 
        uniaxial deformation.
        - nu (float): Poisson's ratio of the material. This dimensionless coefficient describes the ratio of transverse 
        strain to axial strain in a material subjected to uniaxial stress.

        Returns:
        - numpy.ndarray: A 6x6 matrix representing the isotropic elasticity tensor. The tensor is constructed based on 
        Lamé's first parameter (lambda) and the shear modulus (mu), both of which are derived from the input Young's 
        modulus and Poisson's ratio.

        Note:
        - The function assumes isotropic material behavior, where material properties are identical in all directions.
        - The elasticity tensor is used in the context of linear elasticity, where stress is proportional to strain.
    
        """

        lam = E * nu / ((1 + nu) * (1 - 2 * nu))  # Lamé's first parameter
        mu = E / (2. * (1. + nu))  # Shear modulus
        D = np.zeros((6, 6))  # Initialize the elasticity tensor
        D[0:3, 0:3] = lam * np.ones((3, 3))  # Fill the major diagonal blocks with Lamé's first parameter
        for i in range(6):
            if i < 3:
                D[i, i] += 2.0 * mu  # Add 2*mu to the first three diagonal elements
            else:
                D[i, i] += mu  # Add mu to the last three diagonal elements

        return D


    def computeStiffnessMatrix(self, 
                               x : np.ndarray, 
                               E : float, 
                               nu : float):
        """
        Compute the stiffness matrix for an element in a finite element analysis.

        This function calculates the stiffness matrix (Ke) for a finite element based on its nodal coordinates, Young's modulus, 
        and Poisson's ratio. The stiffness matrix is a fundamental component in finite element analysis (FEA), representing the 
        relationship between nodal displacements and applied forces within an element. The computation involves integrating the 
        product of the transpose of the strain-displacement matrix (B), the material elasticity matrix (D), and B over the volume 
        of the element, weighted by the determinant of the Jacobian (J) and the integration weight (ip.wgt) at each Gauss point.

        Parameters:
        - x (numpy.ndarray): An array of nodal coordinates for the element. The shape of this array should correspond to the 
        number of nodes and the spatial dimensionality of the problem (e.g., 2D or 3D).
        - E (float): Young's modulus of the material. This is a measure of the stiffness of the material, defining the ratio 
        of stress to strain under uniaxial deformation within the elastic limit.
        - nu (float): Poisson's ratio of the material. This dimensionless coefficient represents the ratio of transverse 
        strain to axial strain under uniaxial stress conditions.

        Returns:
        - numpy.ndarray: The stiffness matrix (Ke) for the element, with dimensions corresponding to the degrees of freedom 
        of the element (self.dofel x self.dofel). This matrix quantifies the element's resistance to deformation under 
        applied loads.

        Notes:
        - The function assumes linear isotropic elasticity for the material behavior, where material properties are identical 
        in all directions and stress is proportional to strain.
        - Integration over the element's volume is performed using the Gauss quadrature method, where 'self.gauss.points' 
        contains the Gauss points and their corresponding weights (ip.wgt) and derivatives of shape functions (ip.dN).
        """
        Ke = np.zeros((self.dofel, self.dofel))  # Initialize the stiffness matrix

        D = self.getMaterialTensor(E, nu)  # Compute elasticity Tensor

        for ip in self.gauss.points:

            invJ, detJ = self.computeInvJacobian(x, ip)  # Compute the Jacobian matrix at the Gauss point
            dNdX = np.matmul(ip.dN, invJ)  # Derivative of shape functions w.r.t global coordinates

            B = self.computeBMatrix(dNdX)  # Compute the strain-displacement matrix

            # Integrate the contribution of the Gauss point to the stiffness matrix
            Ke += np.matmul(np.transpose(B), np.matmul(D, B)) * detJ * ip.wgt

        return Ke

    
    def computeInvJacobian(self, x, ip):
            
            J = np.matmul(x, ip.dN)

            return np.linalg.inv(J), np.linalg.det(J)

    def computeBMatrix(self, 
                       dNdX : np.ndarray):
        """
        Construct the strain-displacement matrix (B-matrix) for an element in a finite element analysis (FEA) model.

        The B-matrix relates the nodal displacements of an element to the strains within that element. This function 
        computes the B-matrix based on the spatial derivatives of the shape functions (dNdX) with respect to the 
        coordinates. The computed B-matrix is essential for forming the stiffness matrix of the element in the context 
        of FEA, particularly in structural and solid mechanics problems.

        Parameters:
        - dNdX (numpy.ndarray): A matrix of shape function derivatives with respect to the coordinates. Each row 
        corresponds to a shape function, and the columns correspond to the spatial derivatives in each direction 
        (x, y, z).

        Returns:
        - numpy.ndarray: A 6x24 matrix representing the strain-displacement relationships for the element. The B-matrix 
        is structured to facilitate the computation of strains (e_11, e_22, e_33, e_23, e_13, e_12) from the nodal 
        displacements (u, v, w) of the element.

        Note:
        - This function assumes a three-dimensional analysis where each node has three degrees of freedom (u, v, w), 
        corresponding to displacements in the x, y, and z directions, respectively.
        - The `index_u`, `index_v`, and `index_w` attributes used in the function should be predefined and represent the 
        indices in the B-matrix where the displacements u, v, and w contribute to the strain components.
        - The strains are computed using the engineering strain definitions, where e_11, e_22, and e_33 are normal 
        strains, and e_23, e_13, and e_12 are shear strains.
        """

        B = np.zeros((6, 24))

        B[0,self.index_u] = dNdX[:,0]; # e_11 = u_1,1
        B[1,self.index_v] = dNdX[:,1]; # e_22 = u_2,2
        B[2,self.index_w] = dNdX[:,2]; # e_33 = u_3,3

        B[3,self.index_v] = dNdX[:,2];	B[3,self.index_w] = dNdX[:,1];	# e_23 = u_2,3 + u_3,2
        B[4,self.index_u] = dNdX[:,2];	B[4,self.index_w] = dNdX[:,0];	# e_13 = u_1,3 + u_3,1
        B[5,self.index_u] = dNdX[:,1];	B[5,self.index_v] = dNdX[:,0];	# e_12 = u_1,2 + u_2,1

        return B
    
    def computeMassMatrix(self, 
                          x : np.array, 
                          rho : float, 
                          diagonal: bool = False):
        """
        Compute the mass matrix for an element in a finite element analysis.

        This function calculates the mass matrix (Me) for a finite element based on its nodal coordinates and the material density. 
        The mass matrix is a key component in dynamic finite element analysis, representing the distribution of mass within an element. 
        The computation involves integrating the product of the shape functions over the volume of the element, weighted by the material 
        density (rho).

        Parameters:
        - x (numpy.ndarray): An array of nodal coordinates for the element. The shape of this array should correspond to the number of 
        nodes and the spatial dimensionality of the problem (e.g., 2D or 3D).
        - rho (float): The density of the material. This represents the mass per unit volume of the material.
        - diagonal (bool) : Flag to return a diagonal (or lumped) mass matrix

        Returns:
        - numpy.ndarray: The mass matrix (Me) for the element, with dimensions corresponding to the degrees of freedom of the element 
        (self.dofel x self.dofel). This matrix quantifies the distribution of mass within the element under the assumption of a 
        consistent mass formulation.

        Notes:
        - The function assumes a consistent mass formulation, where the mass matrix is derived from the integral of the shape function 
        products, providing a more accurate representation of the mass distribution than a lumped mass matrix - which can also be returned if requested using "digaonal" option
        - Integration over the element's volume is performed using the Gauss quadrature method, where 'self.gauss.points' contains the 
        Gauss points and their corresponding weights (ip.wgt) for numerical integration.
        - The shape functions and their derivatives should be defined according to the element type and geometry (e.g., linear, quadratic).
        """
        Me = np.zeros((self.dofel, self.dofel))  # Initialize the mass matrix

        for ip in self.gauss.points:
            J = np.matmul(x, ip.dN)  # Compute the Jacobian matrix at the Gauss point
            N = ip.N  # Shape function evaluated at the Gauss point

            # Compute the outer product of the shape functions and integrate over the element's volume
            Me += np.outer(N, N) * rho * np.linalg.det(J) * ip.wgt

        if diagonal:
            # Sum the rows of the consistent mass matrix to obtain the diagonal (lumped) mass matrix
            Me = np.diag(np.sum(Me, axis=1))

        return Me
    
    def computeLoadVec(self,
                       x : np.array,
                       v : np.array,
                       f ):
        
        assert x.shape == (4, 3), "Coordinates should be in a 4x3 array."
        assert v.shape == (4, 3), "Normals should be in a 4x3 array."
        
        fe = np.zeros((4, 3))
        
        for ip in self.gauss_2d.points:
            
            J = np.matmul(x, ip.dN) # Compute the Jacobian Matrix at the Gauss Point

            xip = np.matmul(x.transpose(), ip.N)

            vip = np.matmul(v.transpose(), ip.N)

            fe += f(xip) * np.outer(vip, ip.N) * np.linalg.det(J) * ip.wgt

        return fe



       