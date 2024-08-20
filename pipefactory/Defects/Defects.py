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
        
class Cuboid(Defect):

    def __init__(self,
                 s0 : float,
                 phi0 : float,
                 phi1 : float,
                 length : float,
                 height : float,
                 outer_radius : float,
                 thickness: float,
                 el_thru_thick : int):
        
        super().__init__()
        
        self.s0 = s0
        self.phi0 = np.deg2rad(phi0+0.001) # Ensure there is no issue with rounding due to deg2rad when angle is equal to n.phi in mesh.
        self.phi1 = np.deg2rad(phi1-0.001) #* Breaks for >~ 360000 el around circum (very unlikely)
        
        if(self.phi0 > self.phi1):
            self.dphi = 2.*np.pi - self.phi0 + self.phi1
        else:
            self.dphi = self.phi1 - self.phi0

        self.r = outer_radius
        self.l = length
        self.h = height
        self.thickness = thickness
        self.ett = el_thru_thick

        self.r_el, self.el_dr = self.elements_deep()

    def elements_deep(self):
        dx = self.thickness/self.ett

        if self.h < 1.5*dx:
            r_el = 1
            el_dr = self.h
        else:
            r_el = int(np.floor(self.h/dx+0.5))
            el_dr = self.h/r_el
        return r_el, el_dr
    
    def where_in_wedge(self,
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

    def affected_el_idx(self, midline : np.ndarray):
        el_midline = np.array([(midline[i]+midline[i+1])*0.5 for i in range(len(midline)-1)])
        # cuboid_mid_idx = np.argmin(np.abs(el_midline-self.s0))
        all_idx = np.nonzero(np.abs(el_midline-self.s0)<=self.l/2)[0].tolist()
        return all_idx
    
    def __call__(self, xm : np.array, e : Element):
        ds_wedge = self.where_in_wedge(e.midpoint_phi)
        if ds_wedge is not None:
            if (self.r - np.linalg.norm(e.midpoint - xm)) < self.thickness/self.ett:
                return True, ds_wedge
            else:
                return False, ds_wedge
        else:
            return False, ds_wedge


    


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

        self.phi0 = np.deg2rad(phi0+0.001) # Ensure there is no issue with rounding due to deg2rad when angle is equal to n.phi in mesh.
        self.phi1 = np.deg2rad(phi1-0.001) #* Breaks for >~ 360000 el around circum (very unlikely)
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
                 thickness : float,
                 smoothing_dist : float,
                 el_thru_thick : int):
        
        super().__init__()
        
        self.s0 = s0

        self.phi0 = np.deg2rad(phi0+0.001) # Ensure there is no issue with rounding due to deg2rad when angle is equal to n.phi in mesh.
        self.phi1 = np.deg2rad(phi1-0.001) #* Breaks for >~ 360000 el around circum (very unlikely)
        self.w = crack_width
        self.d = crack_depth
        self.r = outer_radius
        self.thickness = thickness
        self.sd = smoothing_dist
        self.ett = el_thru_thick

        if(self.phi0 > self.phi1):
            self.dphi = 2.*np.pi - self.phi0 + self.phi1
        else:
            self.dphi = self.phi1 - self.phi0

        self.stepped_depth()

    def affected_idx(self, midline : np.ndarray):
        crack_mid_idx = np.argmin(np.abs(midline-self.s0))
        all_idx = np.nonzero(np.abs(midline-self.s0)<self.sd)[0].tolist()

        if len(all_idx) == 0:
            left_idx = []
            right_idx = []
        else:
            left_idx = np.arange(all_idx[0],crack_mid_idx).tolist()
            right_idx = np.arange(all_idx[-1],crack_mid_idx, -1).tolist()

        return crack_mid_idx, left_idx, right_idx
    
    def stepped_depth(self):
        dx = self.thickness/self.ett

        if self.ett == 1:
            if self.d > self.thickness:
                self.degen_cutoff = self.r - self.thickness - dx*0.01
            else:
                self.degen_cutoff = self.r - 0.5*dx
        
        else:
            if self.d < 1.5*dx:
                self.n_above = 0
                self.n_below = self.ett - 2
                self.degen_cutoff = self.r - 0.5*dx
            elif (self.d >= (self.ett-1.5)*dx) and (self.d < self.thickness):
                self.n_above = self.ett - 2
                self.n_below = 0
                self.degen_cutoff = self.r - (self.ett-1.5)*dx
            elif self.d == self.thickness:
                self.degen_cutoff = self.r - (self.ett-0.5)*dx
            elif self.d > self.thickness:
                self.degen_cutoff = self.r - self.thickness - dx*0.01 # hundreth of an element
            else:
                self.n_above = np.floor(self.d/dx+0.5) - 1
                self.n_below = self.ett - 2 - (np.floor(self.d/dx+0.5) - 1)
                self.degen_cutoff = self.r - (np.floor(self.d/dx+0.5)-0.5)*dx

            self.apex = self.degen_cutoff - 0.5*dx

        
    def is_in_wedge(self, n : Node):
        phi = n.phi
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
    
    def is_in_crack(self, n : Node):
        if np.linalg.norm(n.v) > self.degen_cutoff:
            return self.is_in_wedge(n)
        else:
            return False
    
    def __call__(self, n : Node, zfactor : float, N : int):

        dr = 0.
        dz = 0.

        if self.d >= self.thickness or self.ett == 1:
            if self.is_in_crack(n):
                factor = (np.linalg.norm(n.v) - self.r + self.d)/ self.d
                if n.global_id >= N:
                    dz = self.w * factor/2
                else:
                    dz = -self.w * factor/2 * zfactor
        
        else:
            dx = self.thickness/self.ett
            dr0 =  self.r - self.d - self.apex
            if np.abs(np.linalg.norm(n.v) - self.apex) < dx*0.25:
                dr = dr0
            else:
                for i in range(self.n_above):
                    if np.abs(np.linalg.norm(n.v) - self.apex-(i+1)*dx) < dx*0.25:
                        dist_left = self.r - self.apex - dr0
                        new_dx = dist_left/(self.n_above+1)
                        dr =  dr0 + (i+1)*(new_dx - dx)
                for i in range(self.n_below):
                    if np.abs(np.linalg.norm(n.v) - self.apex+(i+1)*dx) < dx*0.25:
                        dist_left = self.apex + dr0 - (self.r - self.thickness)
                        new_dx = dist_left/(self.n_below+1)
                        dr =  dr0 - (i+1)*(new_dx - dx)

            if self.is_in_crack(n):
                factor = (np.linalg.norm(n.v) + dr - self.r + self.d)/ self.d
                if n.global_id >= N:
                    dz = self.w * factor/2
                else:
                    dz = -self.w * factor/2 * zfactor

            dr = dr*np.abs(zfactor)

        return dz, dr


