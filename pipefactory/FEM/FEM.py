import numpy as np

from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh

# from ...pipefactory import LinearElasticity, CollarLoad
from .LinearElastic import LinearElasticity
from .Loads import CollarLoad

class FEM:

    def __init__(self, 
                 mesh,
                 f, 
                 reduced : bool = False):

        self.mesh = mesh

        self.f = f

        self.tdof = 3 * mesh.nnodes

        self.le = LinearElasticity(reduced)

        self.assemble()

    def getX(self,e):
        ln = e.list_of_nodes
        X = []
        for n in ln:
            X.append(self.mesh.nodes[n].coords)
        return np.array(X).transpose()

    def assemble(self):

        from concurrent.futures import ThreadPoolExecutor

        def process_element(e):
            return self.le.computeStiffnessMatrix(self.getX(e), E = 1.0, nu = 0.2)
        
        def process_element_M(e):
            return self.le.computeMassMatrix(self.getX(e), rho = 1.0, diagonal = True)
        
        def process_e2g(e):
            return e.e2g
        
        def process_load(e):
            return self.le.computeLoadVec(self.getX(e), self.f)


        with ThreadPoolExecutor() as executor:
            self.Kes = list(executor.map(process_element, self.mesh.elements))
            self.Mes = list(executor.map(process_element_M, self.mesh.elements))
            self.e2g = list(executor.map(process_e2g, self.mesh.elements))
            self.fes = list(executor.amp(process_load, self.mesh.outer_face_elements))
            self.e2g_face = list(executor.map(process_e2g, self.mesh.outer_face_elements))
        

        # Initialize lists to hold COO format data
        values = []
        values_M = []
        indices_row = []
        indices_col = []

        # Example loop over elements (assuming Kes and e2g are your element stiffness matrices and mappings)
        for Ke, mapping in zip(self.Kes, self.e2g):
            local_dof = Ke.shape[0]
            for i in range(local_dof):
                for j in range(local_dof):
                    # Get the global DOF indices
                    global_i = mapping[i]
                    global_j = mapping[j]

                    # Add the value and its indices to the lists
                    values.append(Ke[i, j])
                    indices_row.append(global_i)
                    indices_col.append(global_j)

        for Me, mapping in zip(self.Mes, self.e2g):
            local_dof = Me.shape[0]
            for i in range(local_dof):
                for j in range(local_dof):
                    values_M.append(Me[i, j])    

        # Create the sparse global stiffness matrix

        self.K_global = coo_matrix((values, (indices_row, indices_col)), shape=(self.tdof, self.tdof))
        self.M_global = coo_matrix((values_M, (indices_row, indices_col)), shape=(self.tdof, self.tdof))
        
    def assemble_load(self, load : CollarLoad):
        """
            This function accepts a CollarLoad and evaluates it over the surface and assembles the load vector.
        """

        self.element_loads = []
        for i, e in enumerate(self.mesh.outer_face_elements):
            id = e.midline_indx
            s_midpoint = 0.5 * (self.mesh.midline[id] + self.mesh.midline[id + 1])
            phi = 0.0
            for n in e.list_of_nodes:
                phi += self.mesh.nodes[n].phi / 4.
            self.element_loads.append(load(s_midpoint, phi))

        self.nodal_loads = np.zeros((self.mesh.nnodes, 3))
        for i, e in enumerate(self.mesh.outer_face_elements):
            for j in e.list_of_nodes:
                self.nodal_loads[j,:] += 0.25 * self.element_loads[i] * (self.mesh.nodes[j].v / np.linalg.norm(self.mesh.nodes[j].v))

   
    def eig_decomposition(self, numEigs = 10, maxIter = 1000):

        #lam, v = eigsh(self.K_global, k = numEigs, which ='SM', maxiter=maxIter)

        lam, v = np.linalg.eigh(self.K_global.todense())

        return lam, v
