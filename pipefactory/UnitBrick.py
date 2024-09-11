import numpy as np
from ..pipefactory import Node, Element

import meshio


class UnitBrick:

    def __init__(self):

        self.nodes = []

        x = []
        x.append([0.0, 0.0, 0.0])
        x.append([1.0, 0.0, 0.0])
        x.append([1.0, 1.0, 0.0])
        x.append([0.0, 1.0, 0.0])
        x.append([0.0, 0.0, 1.0])
        x.append([1.0, 0.0, 1.0])
        x.append([1.0, 1.0, 1.0])
        x.append([0.0, 1.0, 1.0])

        for i, xi in enumerate(x):
            self.nodes.append(Node(np.array(xi), i))

        self.elements = []

        self.elements.append(Element([0,1,2,3,4,5,6,7], 'hex', 0, 0))

        self.nnodes = len(self.nodes)
        self.nel = len(self.elements)

    def export(self, 
               filename: str = "foo.vtk",
               point_data : dict = {},
               cell_data : dict = {}
               ):
        
        points = [] 
        for n in self.nodes:
            points.append(n.coords)
            
        connectivity = []
        for e in self.elements:
            if(e.active):
                connectivity.append(e.list_of_nodes)

        elem_type = "hexahedron"

        cells = [
            (elem_type, connectivity),
        ]

        # Alternative with the same options
        meshio.write_points_cells(filename, points, cells, point_data=point_data, cell_data = cell_data)




class TwoBricks:

    def __init__(self):

        self.nodes = []

        x = []
        x.append([0.0, 0.0, 0.0])
        x.append([1.0, 0.0, 0.0])
        x.append([1.0, 1.0, 0.0])
        x.append([0.0, 1.0, 0.0])
        x.append([0.0, 0.0, 1.0])
        x.append([1.0, 0.0, 1.0])
        x.append([1.0, 1.0, 1.0])
        x.append([0.0, 1.0, 1.0])
        x.append([2.0, 0.0, 0.0])
        x.append([2.0, 1.0, 0.0])
        x.append([2.0, 0.0, 1.0])
        x.append([2.0, 1.0, 1.0])

        for i, xi in enumerate(x):
            self.nodes.append(Node(np.array(xi), i))

        self.elements = []

        self.elements.append(Element([0, 1, 2, 3, 4, 5, 6, 7], 'hex', 0, 0))
        self.elements.append(Element([1, 8, 9, 2, 5, 10, 11, 6], 'hex', 0, 0))

        self.nnodes = len(self.nodes)
        self.nel = len(self.elements)
    
    def export(self, 
               filename: str = "foo.vtk",
               point_data : dict = {},
               cell_data : dict = {}
               ):
        
        points = [] 
        for n in self.nodes:
            points.append(n.coords)
            
        connectivity = []
        for e in self.elements:
            if(e.active):
                connectivity.append(e.list_of_nodes)

        
        elem_type = "hexahedron"

        cells = [
            (elem_type, connectivity),
        ]

        # Alternative with the same options
        meshio.write_points_cells(filename, points, cells, point_data=point_data, cell_data = cell_data)
