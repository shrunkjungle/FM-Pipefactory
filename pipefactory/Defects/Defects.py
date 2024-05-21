from pipefactory import Element, Node, cylinder_geodesic_distance
import numpy as np

class Defect():

    def __init__(self):

        return NotImplementedError
    
    def __call__(self):

        return NotImplementedError


class Dimple(Defect):

    def __init__(self,
                 s0 : float,
                 phi0 : float,
                 A : float,
                 ell = float):
        
        super().__init__()
        
        self.s0 = s0
        self.phi0 = np.deg2rad(phi0)
        self.A = A
        self.ell = ell


    def __call__(self,
                 xm : np.array,
                 s : float, 
                 node : Node):
    
        d = cylinder_geodesic_distance(s, node.phi, self.s0, self.phi0, np.linalg.norm(node.coords - xm))

        dr = -self.A * np.exp(-(d / self.ell)**2) # Square Exponential Hump
        
        return dr


class Hole(Defect):

    def __init__(self,
                 s0 : float,
                 phi0 : float,
                 radius : float):
        
        super().__init__()
        
        self.s0 = s0
        self.phi0 = np.deg2rad(phi0)
        self.radius = radius

    def __call__(self,
                 xm : np.array,
                 s : float, 
                 e : Element):
    
        d = cylinder_geodesic_distance(s, e.midpoint_phi, self.s0, self.phi0, np.linalg.norm(e.midpoint - xm))

        if (d  < self.radius):
            return True
        else:
            return False
        


class Radial_Slit(Defect):

    def __init__(self,
                 s0 : float,
                 phi0 : float,
                 phi1 : float,
                 slit_width : float):
        
        super().__init__()
        
        self.s0 = s0

        self.phi0 = np.deg2rad(phi0)
        self.phi1 = np.deg2rad(phi1)
        self.slit_width = slit_width
        

    def __call__(self,
                 xm : np.array,
                 s : float, 
                 e : Element):
        
        if (e.midpoint_phi > self.phi0) and (e.midpoint_phi < self.phi1) and (np.abs(s - self.s0) < 0.5 * self.slit_width): # Need to account for periodicity of angle.
            return True
        else:
            return False


        
    
        d = cylinder_geodesic_distance(s, e.midpoint_phi, self.s0, self.phi0, np.linalg.norm(e.midpoint - xm))

        if (d  < self.radius):
            return True
        else:
            return False