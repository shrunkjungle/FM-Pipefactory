import numpy as np

class Node():

    def __init__(self, 
                 coords : np.array, 
                 global_id: int,
                 local_id: int = None, 
                 proc: int = None,
                 v: np.array = None,
                 phi: float = None):
        
        self.coords = coords
        self.global_id = global_id
        if local_id == None:
            self.local_id = global_id
        else:
            self.local_id = local_id
        self.proc = proc
        self.elements = []
        self.midline_indx = None
        self.section_indx = None
        self.phi = phi

        self.v = v

    def for_xml(self):
        return self.global_id, self.coords[0], self.coords[1], self.coords[2]

    def set_global_id(self, i):
        self.set_global_id = i

    def add_element(self, ie):
        self.elements.append(ie)

    def transform(self, u : np.array):
        self.coords += u

    def set_midline_indx(self, id):
        self.midline_indx = id

    def set_section_indx(self, id):
        self.section_indx = id


class Element():

    def __init__(self,
                 list_of_nodes : list[int],
                 elem_type: str,
                 global_id: int,
                 midline_indx: int,
                 active : bool = True,
                 local_id: int = None,
                 proc: int = None,
                 edge_list: list[int] = None):
        
        self.list_of_nodes = list_of_nodes
        self.elem_type = type
        self.nodel = len(list_of_nodes)
        self.global_id = global_id
        self.edge_list = edge_list
        self.active = active
        self.midline_indx = midline_indx
        self.midpoint = None
        self.midpoint_phi = None

        if self.elem_type == "hex":
            self.e2g = np.zeros((24,), dtype=np.int32)
            for i in range(3):
                for j in range(8):
                    self.e2g[j*3 + i] = 3 * list_of_nodes[j] + i

        if self.elem_type == "quad":
            self.e2g = np.zeros((12,), dtype=np.int32)
            for i in range(3):
                for j in range(4):
                    self.e2g[j*3 + i] = 3 * list_of_nodes[j] + i

    def change_node(self, old, new):
        loc = np.where(np.array(self.list_of_nodes) == old)[0][0]
        self.list_of_nodes[loc] = new
    
    def for_xml(self):
        return self.global_id, self.list_of_nodes
    
    def set_xi(self, xi):
        self.xi = xi

    def add_edges(self, edge_list):
        self.edge_list = edge_list







            




                
        





            
        
