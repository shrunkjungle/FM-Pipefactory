import numpy as np

from pipefactory import LinearHex, LinearQuad

class IntPoint():

    def __init__(self,
                 local_coords : np.array,
                 wgt : float,
                 shapefunction):
        
        self.local_coords = local_coords
        self.N = shapefunction.N(local_coords)
        self.dN = shapefunction.dNdX(local_coords)
        self.wgt = wgt


class GaussQuadrature():

    def __init__(self, reduced = False):

        self.points = []
        if (reduced):
            self.points.append(IntPoint(np.array([0.0, 0.0, 0.0]), 8.0, LinearHex()))
        else:
            gauss_points_1d = [-1.0 / np.sqrt(3), 1.0 / np.sqrt(3)]
            for x in gauss_points_1d:
                for y in gauss_points_1d:
                    for z in gauss_points_1d:
                        self.points.append(IntPoint(np.array([x, y, z]), 1.0, LinearHex()))

        self.nip = len(self.points)

class GaussQuadrature2D():

    def __init__(self, reduced = False):

        self.points = []
        if (reduced):
            self.points.append(IntPoint(np.array([0.0, 0.0, 0.0]), 4.0, LinearQuad()))
        else:
            gauss_points_1d = [-1.0 / np.sqrt(3), 1.0 / np.sqrt(3)]
            for x in gauss_points_1d:
                for y in gauss_points_1d:
                    self.points.append(IntPoint(np.array([x, y, 0.0]), 1.0, LinearQuad()))

        self.nip = len(self.points)


        

