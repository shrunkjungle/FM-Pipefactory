import numpy as np

from pipefactory import cylinder_geodesic_distance

class CollarLoad:

    def __init__(self,
                 N : int,
                 s0 : float,
                 ell : float,
                 radius : float,
                 dt : float,
                 freq : float = 10.,
                 cycles : int = 5,):
        
        self.N = N
        self.s0 = s0
        self.ell = ell
        self.radius = radius

        self.dt = dt
        self.sampling_rate = 1. / dt
        self.cycles = cycles
        self.freq = freq

      
        samples_per_cycle = int(self.sampling_rate / self.freq)
        self.hamming_window = np.hamming(samples_per_cycle * cycles)
        self.tend = self.cycles / self.freq
        
        # Set up N loads around
        dphi = 2. * np.pi / N
        self.phi = []
        for i in range(N):
            self.phi.append(i*dphi)

    def __call__(self, s, phi):

        """
            This the function call for Collar Loads which evaluates to the magnitude of the load at given s and phi (radial angle)
        
        """
        f = 0.
        for phi0 in self.phi:
            d = cylinder_geodesic_distance(s, phi, self.s0, phi0, self.radius)
            f += np.exp(-(d/self.ell)**2)
        return f
    
    def amplitude(self, t):

        if t > self.tend:
            return 0.0
        else:
            n = np.floor(t / self.dt)
            return self.hamming_window[n] * np.sin(2. * np.pi * self.freq * t)
