from pipefactory import Element, Node, cylinder_geodesic_distance
import numpy as np

class Defect():
    """
        Defect represents the base class for a defect, defining the API which must be implemented.
    """
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
    

class Weld(Defect):

    def __init__(self,
                 s0 : float,
                 Aout :  float, 
                 Ain : float,
                 ell_out : float,
                 ell_in : float,
                 outer_radius : float,
                 thickness : float):
        
        self.s0 = s0
        self.Aout = Aout
        self.Ain = Ain
        self.ell_out = ell_out
        self.ell_in = ell_in
        self.radius = outer_radius - thickness/2.
        self.thickness = thickness
    
    def __call__(self,
                 xm: np.array,
                 s : float,
                 node : Node):
        
        z = 2. * (np.linalg.norm(node.coords - xm) - self.radius) / self.thickness
        
        d = cylinder_geodesic_distance(s, 0.0, self.s0, 0.0, np.linalg.norm(node.coords - xm))

        alpha = 0.5*(z + 1)
        ell =  alpha * self.ell_out + (1 - alpha) * self.ell_in 

        if ( z >= 0.0 ):
            A = self.Aout + 0.5 * self.thickness
        else:
            A = self.Ain + 0.5 * self.thickness

        return z * A * np.exp(-(d / ell)**2)


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
                 slit_width : float,
                 outer_radius : float,
                 thickness: float,
                 partial: bool = False,
                 profile = None):
        
        super().__init__()
        
        self.s0 = s0

        self.phi0 = np.deg2rad(phi0)
        self.phi1 = np.deg2rad(phi1)
        self.slit_width = slit_width
        self.r = outer_radius - thickness/2.
        self.thickness = thickness
        self.partial = partial
        self.profile = profile

        if(self.phi0 > self.phi1):
            self.dphi = 2.*np.pi - self.phi0 + self.phi1
        else:
            self.dphi = self.phi1 - self.phi0

    def where_in_slit(self,
                      phi : float):
        """
        This function returns the factor through the slit a point is, if it is outside then return "None"
        """
        if (self.phi0 > self.phi1): # Goes through the north pole
            if (phi > self.phi0) and (phi <= 2.0*np.pi):
                return (phi - self.phi0) / self.dphi
            elif (phi < self.phi1) and (phi > 0.0):
                return (2*np.pi - self.phi0 + phi) / self.dphi
            else:
                return None
        else: # phi1 > phi0
            if (phi >= self.phi0) and (phi <= self.phi1):
                return (phi - self.phi0) / self.dphi
            else:
                return None
            
    def distance_to_slit(self,
                         s):
        """
            This function computes the distance of the midpoint of the element to the slit.
        """
        return np.abs(s - self.s0) < 0.5 * self.slit_width

    def __call__(self,
                 xm : np.array,
                 s : float, 
                 e : Element):
        
        ds_slit = self.where_in_slit(e.midpoint_phi)
        
        if (np.abs(s - self.s0) < 0.5 * self.slit_width) and (ds_slit is not None): # Is the distance to the slit less than half the slit_width.
            if(self.partial) and (self.profile is not None):
                z = 2. * (np.linalg.norm(e.midpoint - xm) - self.r) / self.thickness
            
                return self.profile(ds_slit, z)
            else:
                return True
        else:
            return False

class RadialCrack(Defect):

    def __init__(self,
                 s0 : float,
                 phi0 : float,
                 phi1 : float,
                 crack_width : float,
                 crack_depth : float,
                 outer_radius : float,
                 thickness: float):
        
        super().__init__()
        
        self.s0 = s0

        self.phi0 = np.deg2rad(phi0)
        self.phi1 = np.deg2rad(phi1)
        self.w = crack_width
        self.d = crack_depth
        self.r = outer_radius
        self.thickness = thickness

        if(self.phi0 > self.phi1):
            self.dphi = 2.*np.pi - self.phi0 + self.phi1
        else:
            self.dphi = self.phi1 - self.phi0

    def closest_idx(self, midline : np.ndarray):
        crack_mid_idx = np.argmin(np.abs(midline-self.s0))
        return crack_mid_idx
    
    def is_in_crack(self, n : Node):
        phi = n.phi
        if np.linalg.norm(n.v) > self.r - self.d:
            if (self.phi0 > self.phi1): # Goes through the north pole
                if (phi > self.phi0) and (phi <= 2.0*np.pi):
                    return True
                elif (phi < self.phi1) and (phi > 0.0):
                    return True
                else:
                    return False
            else: # phi1 > phi0
                if (phi > self.phi0) and (phi < self.phi1):
                    return True
                else:
                    return False
        else:
            return False
    
    def __call__(self, n : Node, N : int):

        factor = (np.linalg.norm(n.v) - self.r + self.d)/ self.d

        if n.global_id >= N:
            return self.w * factor/2
        else:
            return -self.w * factor/2


