from pipefactory import Node, Element, find_section, rot_vec, get_orthogonal_inplane, rotate_point_about_point, rotate_triad, RadialCrack, Cuboid

import meshio
import numpy as np
from scipy.sparse import coo_matrix

class Pipe():

    

    def __init__(self,
                 outer_radius : float,
                 thickness: float,
                 section_list : list,
                 elem_type: tuple,
                 element_size: float,
                 element_around_circum: int,
                 elements_through_thickness: int = 1,
                 partition: bool = False,
                 nparts: int = None,
                 size_overlap: int = None,
                 mesh_refinement = None):
    

        self.section_list = section_list
        self.thickness = thickness
        self.element_size = element_size
        self.element_around_circum = element_around_circum
        self.elements_through_thickness = elements_through_thickness

        self.v = {
        "z": np.array([0.0, 0.0, 1.0]),
        "y": np.array([0.0, -1.0, 0.0]),
        "x": np.array([1.0, 0.0, 0.0])
        }

        self.midline_to_section = []

        self.process_section_list()

        self.radius = outer_radius - thickness/2


        self.elem_type = elem_type[0]
        self.higher_order = elem_type[1]

        self.midline = self.make_midline(0.0, self.section_ends, element_size, self.higher_order, MR=mesh_refinement)

        self.nodes, self.elements, self.outer_face_elements = self.build()

        self.nnodes = len(self.nodes)
        self.nel = len(self.elements)

        self.transform(partition, nparts, size_overlap)

        self.find_outside_face()


    def process_section_list(self):

        self.num_sections = len(self.section_list)
    
        self.length = 0.0
        self.old_length = 0.0
        self.section_lengths = []

        self.init_dir = None

        # def cleanup_bend_new(Bdict):
        #     v1 = Bdict['param']['dir1']
        #     v2 = Bdict['param']['dir2']

        #     v1 = v1/np.linalg.norm(v1)
        #     v2 = v2/np.linalg.norm(v2)

        #     angle = np.arccos(np.dot(v1,v2))

        #     Bdict['param']['angle'] = angle
        #     Bdict['param']['dir1'] = v1
        #     Bdict['param']['dir2'] = v2

        # def cleanup_straight_new(Sdict):
        #     v1 = Sdict['dir']

        #     Sdict['dir'] = v1/np.linalg.norm(v1)

        for idx, s in enumerate(self.section_list):

            if self.init_dir is None:
                if s['type'].lower() == 'straight_new':
                    self.init_dir = s['dir']
                elif s['type'].lower() == 'bend_new':
                    self.init_dir = s['param']['dir1']
                else:
                    self.init_dir = np.array([1.,0.,0.])

            if s['type'].lower() == 'straight':
                self.length += s['length']
            elif s['type'].lower() == 'straight_new':
                # cleanup_straight_new(self.section_list[idx])
                self.length += s['length']
            elif s['type'].lower() == 'bend':
                self.length += np.abs(s['param']['radius'] * s['param']['angle'])
            elif s['type'].lower() == 'bend_new':
                # cleanup_bend_new(self.section_list[idx])
                self.length += np.abs(s['param']['radius'] * s['param']['angle'])
            else:
                raise Exception("Section can only be of type straight or bend.")
            self.section_lengths.append(self.length - self.old_length)
            self.old_length = self.length

        self.section_ends = np.cumsum(self.section_lengths)


    def make_midline(self, 
                     start: float, 
                     section_ends: list, 
                     element_size: float,
                     higher_order: bool = False,
                     MR=None):
        """
        Create a NumPy array of equally spaced points within specified sections.

        Each section is defined by its start point (inclusive) and end point (inclusive),
        and the points within each section are spaced according to the specified element size.

        Parameters:
        start (float): The starting point of the first section.
        section_ends (list of float): A list of end points for each section. Each element in the list
                                    represents the end point of a section, starting from the given start point.
        element_size (float): The distance between each point in a section.
        higher_order (bool): If the mesh is higher order then we need to make more points.

        Returns:
        np.array: A NumPy array of type np.float64 containing the midline points. Each point is spaced at 
                'element_size' distance, starting from 'start' and ending at the last point of the final 
                section in 'section_ends'.
        """
        midline = [start]

        if MR is not None:
            def el_func(x):
                return MR(x)*(MR.dx - element_size) + element_size

        for end in section_ends:
            if MR is None:
                num_elements = int(np.ceil((end - start) / element_size))
                num_nodes = num_elements + 1
                if(higher_order):
                    num_nodes += num_elements
                tmp = np.linspace(start, end , num_nodes).tolist()
                tmp.pop(0)

            else:
                cur_dist = start
                rtmp=[]
                while el_func(cur_dist)+cur_dist < end:
                    rtmp.append(el_func(cur_dist)+cur_dist)
                    cur_dist = el_func(cur_dist)+cur_dist
                last = rtmp[-1]
                tmp = [start+(x-start)*(end-start)/(last-start) for x in rtmp]

            midline += tmp
            start = end

        return np.array(midline, dtype=np.float64)

    def build(self):
        
        # Elements around circumference
        elements_per_radius = int(self.element_around_circum)
        if (self.higher_order):
            elements_per_radius *= 2
        theta = np.linspace(0.0, 2. * np.pi, elements_per_radius+1).tolist()
        theta.pop(-1)

        # Elements through thickness

        k = 0
        nr = len(theta)
        ns = len(self.midline)

        self.nr = nr
        self.ns = ns

        if (self.elem_type == "hex"):
            if self.higher_order:
                z = np.linspace(-0.5*self.thickness, 0.5*self.thickness, 2 * self.elements_through_thickness + 1).tolist()
            else:
                z = np.linspace(-0.5*self.thickness, 0.5*self.thickness, self.elements_through_thickness + 1).tolist()
        else:
            z = [0.0]

        # Add nodes
        nodes = []
        for l, zi in enumerate(z):
            for i, s in enumerate(self.midline):
                for j, th in enumerate(theta):
                    x = np.array([s, (self.radius + zi) * np.sin(th), (self.radius + zi) * np.cos(th)])
                    v = np.array([0.0, (self.radius + zi) * np.sin(th), (self.radius + zi) * np.cos(th)])
                    nodes.append(Node(x, k, v=v, phi = th))
                    nodes[-1].set_midline_indx(i)
                    k += 1 # Increment Global Counter


        # Build Connectivity
        ie = 0
        elements = []

        i_face_e = 0
        outside_face = []

        if(self.elem_type == "quad"):
            if(self.higher_order == False):
                for i in range(self.ns - 1):
                    for j in range(nr):
                        if j < nr-1:
                            list_of_nodes = [i*nr + j, (i+1)*nr + j, (i+1)*nr+j+1, i*nr+j+1]
                        else:
                            list_of_nodes = [i*nr + j, (i+1)*nr + j, (i+1)*nr, i*nr]
                        elements.append(Element(list_of_nodes, "quad4", ie))
                        ie += 1
    
            else: # This is a higher order element
                
                for i in range(0, ns - 2, 2):
                    for j in range(0,nr,2):
                        if j < nr - 2:
                            list_of_nodes = [ i*nr +j, (i+2)*nr + j, (i+2)*nr + j + 2, i * nr + j + 2,
                                            (i+1)*nr + j, (i+2)*nr + j + 1, (i+1)*nr + j + 2, i * nr + j + 1,
                                            (i+1)*nr + j + 1]
                        else:
                            list_of_nodes = [i*nr + j, (i+2)*nr + j, (i+2)*nr, i*nr,
                                             (i+1)*nr + j, (i+2)*nr+j+1, (i+1)*nr, i*nr+j+1,
                                             (i+1)*nr + j + 1]
                        elements.append(Element(list_of_nodes,"quad9", ie))
                        ie += 1



        elif(self.elem_type == "tri"):
            if(self.higher_order == False):
                for i in range(ns - 1):
                    for j in range(nr):
                            if j < nr-1:
                                list_of_nodes = [i*nr + j, (i+1)*nr + j + 1, i*nr+j+1]
                                elements.append(Element(list_of_nodes, "tri", ie))
                                ie += 1

                                list_of_nodes = [i*nr + j, (i+1)*nr+j, (i+1)*nr+j+1]
                                elements.append(Element(list_of_nodes, "tri", ie))
                                ie += 1
                            else:
                                list_of_nodes = [i*nr + j, (i+1)*nr, i*nr]
                                elements.append(Element(list_of_nodes, "tri", ie))
                                ie += 1

                                list_of_nodes = [i*nr + j, (i+1)*nr + j, (i+1)*nr]
                                elements.append(Element(list_of_nodes, "tri", ie))
                                ie += 1

            else: # This is a higher order element

                for i in range(0, ns - 2, 2):
                    for j in range(0,nr,2):
                        if j < nr - 2:

                            list_of_nodes = [ i*nr + j, (i+2)*nr + j+2, i*nr + j+2,
                                            (i+1)*nr + j+1, (i+1)*nr + j+2, i*nr + j + 1]
                            elements.append(Element(list_of_nodes, "tri6", ie))
                            ie += 1

                            list_of_nodes = [ i*nr + j, (i+2)*nr + j, (i+2)*nr + j+2,
                                            (i+1)*nr + j, (i+2)*nr + j + 1, (i+1)*nr + j + 1]
                            elements.append(Element(list_of_nodes, "tri6", ie))
                            ie += 1
                            
                            
                        else:
                            list_of_nodes = [i*nr + j, (i+2)*nr, i*nr, 
                                             (i+1)*nr + j+1, (i+1)*nr, i*nr + j+1]
                            elements.append(Element(list_of_nodes, "tri6", ie))
                            ie += 1

                            list_of_nodes = [i*nr+j, (i+2)*nr + j, (i+2)*nr, (i+1)*nr+j, (i+2)*nr+j+1, (i+1)*nr+j+1]
                            elements.append(Element(list_of_nodes, "tri6", ie))
                            ie += 1

        elif (self.elem_type == "hex"):
            nz = self.ns * nr
            if(self.higher_order == False):
                for k in range(len(z)-1):
                    for i in range(self.ns - 1):
                        for j in range(nr):
                            if j < nr-1:
                                list_of_nodes = [i*nr + j + k*nz, (i+1)*nr + j + k*nz, (i+1)*nr+j+1 + k*nz, i*nr+j+1 + k*nz,
                                                i*nr + j + (k+1)*nz, (i+1)*nr + j + (k+1)*nz, (i+1)*nr+j+1 + (k+1)*nz, i*nr+j+1 + (k+1)*nz]
                                
                                if(k == len(z) - 2):
                                    list_of_nodes_face = [i*nr + j + (k+1)*nz, (i+1)*nr + j + (k+1)*nz, (i+1)*nr+j+1 + (k+1)*nz, i*nr+j+1 + (k+1)*nz]

                            else:
                                list_of_nodes = [i*nr + j + k*nz, (i+1)*nr + j + k*nz, (i+1)*nr + k*nz, i*nr + k*nz,
                                                i*nr + j + (k+1)*nz, (i+1)*nr + j + (k+1)*nz, (i+1)*nr + (k+1)*nz, i*nr + (k+1)*nz]
                                if (k == len(z) - 2):
                                    list_of_nodes_face = [i*nr + j + (k+1)*nz, (i+1)*nr + j + (k+1)*nz, (i+1)*nr + (k+1)*nz, i*nr + (k+1)*nz]
                            
                            elements.append(Element(list_of_nodes, "hex8", ie, midline_indx=i))

                            if(k == len(z) - 2):
                                outside_face.append(Element(list_of_nodes_face, "quad", i_face_e, midline_indx=i))
                                i_face_e += 1

                            ie += 1
            
            else:
                raise("This element has yet to be tested, Work In Progress.")
                """
                for k in range(0,len(z)-1,2):
                    for i in range(0,self.ns-2,2):
                        for j in range(0,nr,2):
                            if j < nr - 2:
                                list_of_nodes = [i*nr + j + k*nz, (i+1)*nr + j + k*nz, (i+1)*nr+j+1 + k*nz, i*nr+j+1 + k*nz,
                                     i*nr + j + (k+2)*nz, (i+1)*nr + j + (k+2)*nz, (i+1)*nr+j+1 + (k+2)*nz, i*nr+j+1 + (k+2)*nz,
                                     (i+1)*nr + j + k*nz, (i+1)*nr + j+1 + k*nz, (i+1)*nr + j+2 + k*nz, i*nr + j+1 +k*nz,
                                     (i+1)*nr + j + (k+2)*nz, (i+1)*nr + j+1 + (k+2)*nz, (i+1)*nr + j+2 + (k+2)*nz, i*nr + j+1 +(k+2)*nz,
                                     i*nr + j + (k+1)*nz, (i+2)*nr + j + 1 + (k+1)*nz, (i+2)*nr + j + 2 + (k+1)*nz, i*nr + j + 2 + (k+1)*nz,
                                     i*nr + j+1 + (k+1)*nz, (i+2)*nr + j+1 + (k+1)*nz, (i+1)*nr + j + (k+1)*nz, (i+1)*nr + j+2 + (k+1)*nz,
                                     (i+1)*nr + j + 1 + k*nz, (i+1)*nr + j + 1 + (k+2)*nz, (i+1)*nr + j+1 + (k+1)*nz]
                """

        else:
            
            raise Exception("We have currently only implemented quad and tri elements")
        
        for i, e in enumerate(outside_face):
            mid = np.array([0.0, 0.0, 0.0])
            for id_n in e.list_of_nodes:
                try:
                    mid += nodes[id_n].coords 
                except:
                    print(id_n)
            mid /= 4
            e.midpoint = mid
    
        
        return nodes, elements, outside_face
    

    def find_outside_face(self):

        self.outer_face = []
        for i, n in enumerate(self.nodes):
            if (np.abs(np.linalg.norm(self.midline_x[n.midline_indx] - n.coords) - (self.radius + 0.5 * self.thickness))) < 1e-4:
                self.outer_face.append(i)
    
    def transform(self, 
                  partition : bool = False,
                  nparts : int = None, 
                  size_overlap : int = None):
        
        self.map_to_section()
        self.transform_ends()

        self.build_midline()
        if(partition):
            self.partition_pipe(nparts, size_overlap)
        
        self.transform_mesh()
    
    def map_to_section(self):
        """
            This function assigns each point on the midline to a section, using the utility function find_section()

            It takes no inputs or returns nothing. It updates each update each "node" with the index of the section.
        """
        self.midline_to_section = []
        for si in self.midline:
            self.midline_to_section.append(find_section(self.section_ends, si))

        for n in self.nodes:
            id = self.midline_to_section[n.midline_indx]
            n.set_section_indx(id)

    def transform_ends(self):

        self.tri = [np.array([1., 0., 0.]), 
                         np.array([0., 1., 0.]), 
                         np.array([0., 0., 1.])]

        init_angle = np.arccos(self.init_dir[0]) # np.dot([1,0,0], init_dir)

        if init_angle > 0.0:
            init_axis = np.cross(np.array([1.,0.,0.]),self.init_dir)
            init_axis = init_axis/np.linalg.norm(init_axis)
            self.tri = rotate_triad(self.tri, init_axis, init_angle)

        self.triads = [ self.tri ] * ( len(self.section_ends) + 1 )
        
        self.x_sec_ends = [np.array([0., 0., 0.])] * ( len(self.section_ends) + 1 ) # origin is [0,0,0]
        
        for i in range(len(self.section_ends)): # for each of the sections ends
            x0 = self.x_sec_ends[i] # This is the coordinates of the last end

            # Straight Section
            if self.section_list[i]['type'].lower() == 'straight':
                l = self.section_list[i]['length']
                self.triads[i+1] = self.triads[i].copy()
                self.x_sec_ends[i + 1] = self.get_coords_straight(l, x0, self.triads[i][0])

            #############NEW

            elif self.section_list[i]['type'].lower() == 'straight_new':
                l = self.section_list[i]['length']
                self.triads[i+1] = self.triads[i].copy()
                self.x_sec_ends[i + 1] = self.get_coords_straight(l, x0, self.triads[i][0])

            # Bend Section
            elif self.section_list[i]['type'].lower() == 'bend': # If it is a bend

                angle = self.section_list[i]['param']['angle']

                axis, v = self.which_rotation_axis(self.section_list[i]['param']['axis'], i)

                # Need a point to rotate about
                roc = self.section_list[i]['param']['radius']
                p = x0 + np.sign(angle) * roc * v

                self.x_sec_ends[i+1] = rotate_point_about_point(x0, axis, angle, p)
                self.triads[i+1] = rotate_triad(self.triads[i], axis, angle)

            ############# NEW Bend Section
            elif self.section_list[i]['type'].lower() == 'bend_new': # If it is a bend

                angle = self.section_list[i]['param']['angle']

                axis, v = self.which_rotation_axis_new(self.section_list[i]['param']['dir1'],self.section_list[i]['param']['dir2'], i)

                # Need a point to rotate about
                roc = self.section_list[i]['param']['radius']
                p = x0 + roc * v

                self.x_sec_ends[i+1] = rotate_point_about_point(x0, axis, angle, p)
                self.triads[i+1] = rotate_triad(self.triads[i], axis, angle)

            else:
                raise Exception("Unknown section type " + self.section_list[i]['type'])
            
    def which_rotation_axis(self,direction, i):
        
        if (direction == "left_right"):
            axis = self.triads[i][2]
            v = self.triads[i][1]
        elif (direction == "up_down"):
            axis = -self.triads[i][1]
            v = self.triads[i][2]
        else:
            raise Exception("Unknown section type " + direction)

        return axis, v
    
    def which_rotation_axis_new(self,dir1,dir2, i):
        
        axis = np.cross(dir1,dir2)
        axis = axis/np.linalg.norm(axis)

        v = np.cross(axis,dir1)

        return axis, v
    
    def add_defect_displacement(self, defect_instance):

        for n in self.nodes: # loop through nodes 
            idx = n.midline_indx 
            x0 = self.midline_x[idx]
            x = n.coords # coor
            r = np.linalg.norm(x - x0)
            dr = defect_instance(self.midline_x[idx], 
                                 self.midline[idx], n)
            n.coords = x0 + (r + dr) * ((x - x0) / np.linalg.norm(x - x0))

    def add_elements(self, cuboid_instance : Cuboid):

        if self.elem_type != "hex":
            raise Exception("Only Hex elements supported currently.")
        
        el_mid_idxs = cuboid_instance.affected_el_idx(self.midline)
        nx = len(el_mid_idxs)

        nr, dr = cuboid_instance.elements_deep()

        def correct_angle(x):
            if x < np.pi:
                return x+2*np.pi
            else:
                return x

        el_factor = []
        for e in self.elements: # find elements under cuboid
            if e.midline_indx in el_mid_idxs:
                nodes = e.list_of_nodes
                m = np.array([0., 0., 0.])
                x0 = np.array([0., 0., 0.])
                for i in nodes:
                    n = self.nodes[i]
                    m += n.coords
                    x0 += self.midline_x[n.midline_indx]
                
                n_phi_list = [self.nodes[i].phi for i in nodes]
                if (np.max(n_phi_list) - np.min(n_phi_list))>np.pi:
                    n_phi_list = list(map(correct_angle, n_phi_list))
                m_phi = sum(n_phi_list)

                e.midpoint = m / len(nodes)
                e.midpoint_phi = m_phi / len(nodes)

                if cuboid_instance(x0 / len(nodes), e)[0]:
                    el_factor.append((e.global_id, cuboid_instance(x0 / len(nodes), e)[1]))

        assert len(el_factor)/nx == float(len(el_factor)//nx) # check divides into
        nth = int(len(el_factor)/nx)

        for i in range(nx): # sort elements under cuboid based on position
            a = el_factor[i*nth:(i+1)*nth]
            a.sort(key = lambda x: x[1])
            el_factor[i*nth:(i+1)*nth] = [entry[0] for entry in a]

        surface_nodes = []
        new_ie = 0
        for i in range(nx):
            for j in range(nth):
                surface_nodes.append(self.nodes[self.elements[el_factor[new_ie]].list_of_nodes[4]])
                if j == nth-1:
                    surface_nodes.append(self.nodes[self.elements[el_factor[new_ie]].list_of_nodes[7]])
                new_ie += 1

        new_ie -= nth
        for j in range(nth):
            surface_nodes.append(self.nodes[self.elements[el_factor[new_ie]].list_of_nodes[5]])
            if j == nth-1:
                surface_nodes.append(self.nodes[self.elements[el_factor[new_ie]].list_of_nodes[6]])
            new_ie += 1

        from copy import deepcopy

        new_nodes = []
        nn = self.nnodes
        k=0
        nz = (nx+1)*(nth+1)

        for l in range(nr):
            for i in range(nx+1):
                for j in range(nth+1):
                    new_node = deepcopy(surface_nodes[k%nz])
                    idx = new_node.midline_indx 
                    x0 = self.midline_x[idx]
                    x = new_node.coords
                    new_node.coords += (l+1)*dr*((x - x0) / np.linalg.norm(x - x0))
                    new_node.global_id = nn + k

                    new_nodes.append(new_node.global_id)
                    self.nodes.append(new_node)

                    k += 1 # Increment Global Counter

        ie = self.nel
        cuboid_nodes = [node.global_id for node in surface_nodes]
        cuboid_nodes = cuboid_nodes + new_nodes

        for k in range(nr): 
            for i in range(nx): 
                for j in range(nth): 
                    list_of_new_nodes = [i*(nth+1) + j + k*nz, (i+1)*(nth+1) + j + k*nz, (i+1)*(nth+1)+j+1 + k*nz, i*(nth+1)+j+1 + k*nz,
                                    i*(nth+1) + j + (k+1)*nz, (i+1)*(nth+1) + j + (k+1)*nz, (i+1)*(nth+1)+j+1 + (k+1)*nz, i*(nth+1)+j+1 + (k+1)*nz]
                    list_of_nodes = [cuboid_nodes[idx] for idx in list_of_new_nodes]

                    self.elements.append(Element(list_of_nodes, "hex8", ie, midline_indx=i+el_mid_idxs[0]))
                    ie += 1

        self.nnodes = len(self.nodes)
        self.nel = len(self.elements)


    def remove_elements(self, hole_instance):

        for e in self.elements: # loop through all elements

            nodes = e.list_of_nodes
            m = np.array([0., 0., 0.])
            m_s = 0.0
            x0 = np.array([0., 0., 0.])
            for i in nodes:
                n = self.nodes[i]
                m += n.coords
                m_s += self.midline[n.midline_indx]
                x0 += self.midline_x[n.midline_indx]

            n_phi_list = [self.nodes[i].phi for i in nodes]
            if (np.max(n_phi_list) - np.min(n_phi_list))>np.pi:
                def correct_angle(x):
                    if x < np.pi:
                        return x+2*np.pi
                    else:
                        return x
                n_phi_list = list(map(correct_angle, n_phi_list))
            m_phi = sum(n_phi_list)

            e.midpoint = m / len(nodes)
            e.midpoint_phi = m_phi / len(nodes)

            if hole_instance(x0 / len(nodes), m_s / len(nodes), e):
                e.active = False

    def degenerate_crack(self, crack_instance : RadialCrack):

        crack_mid_idx, left_idx, right_idx = crack_instance.affected_idx(self.midline)

        left_factor = np.linspace(1./(len(left_idx)+1), 1. - 1./(len(left_idx)+1), len(left_idx)).tolist()
        right_factor = np.linspace(1./(len(right_idx)+1), 1. - 1./(len(right_idx)+1), len(right_idx)).tolist()

        crack_idxs = []
        for n in self.nodes:
            if n.midline_indx == crack_mid_idx:
                if crack_instance.is_in_crack(n):
                    crack_idxs.append(n.global_id)

        from copy import deepcopy

        new_nnodes = self.nnodes-1
        for i in crack_idxs:
            new_nnodes += 1
            degen_node = deepcopy(self.nodes[i])
            degen_node.global_id = new_nnodes
            self.nodes.append(degen_node)

        for el_idx, e in enumerate(self.elements):
            if (e.midline_indx == crack_mid_idx):
                for crack_idx, id in enumerate(crack_idxs):
                    if id in e.list_of_nodes:
                        self.elements[el_idx].list_of_nodes = list(map(lambda x: self.nnodes + crack_idx if x == id else x, self.elements[el_idx].list_of_nodes))

        for n in self.nodes:
            if n.midline_indx in (left_idx+right_idx) or n.midline_indx == crack_mid_idx:
                if crack_instance.is_in_wedge(n):
                    if n.midline_indx == crack_mid_idx:
                        zfactor = 1.
                    elif n.midline_indx in left_idx:
                        zfactor = left_factor[left_idx.index(n.midline_indx)]
                    else:
                        zfactor = -1 * right_factor[right_idx.index(n.midline_indx)]
                    
                    dz, dr = crack_instance(n, zfactor, self.nnodes)

                    idx = n.midline_indx 
                    x0 = self.midline_x[idx]
                    x = n.coords

                    n.coords += dz * self.midline_triad[crack_mid_idx][0] + dr*((x - x0) / np.linalg.norm(x - x0))

        self.nnodes = len(self.nodes)

    def partition_pipe(self, nparts, size_overlap):

        self.nparts = nparts
        self.size_overlap = size_overlap

        # Build Initial Partition
        elem_midline = len(self.midline) - 1
        midline_elements = np.arange(elem_midline)
        self.parts = np.array_split(midline_elements, nparts)

        # Produce Overlapping Domain Decomposition
        e_indx = []
        n_indx = []
        data = []
        M = []
        for i, e in enumerate(midline_elements):
            for n in [e, e+1]:
                e_indx.append(i)
                n_indx.append(n)
                data.append(1)
        M = coo_matrix((data, (n_indx, e_indx)), shape=(elem_midline + 1, elem_midline))

        self.p_e = []
        self.p_nodes = []
        for i, p in enumerate(self.parts): # For each partition
            elem = p.copy()
            for j in np.arange(size_overlap):  # For each layer of the overlap
                e = np.zeros((elem_midline,))
                e[elem] = 1.0
                n = M @ e
                nodes = np.where(n > 0)
                n[nodes] = 1
                e = M.transpose() @ n
                elem = np.where(e > 0)

            self.p_e.append(elem[0])
            self.p_nodes.append(np.where(M @ e > 0))

        ## Construct Partition of Unity Operator
        xi = np.zeros((elem_midline + 1,))
        for i, n in enumerate(self.p_nodes):
            for id in n:
                xi[id] += 1.0
        self.xi = 1. / xi

       
        self.Omg = [[] for _ in range(self.nparts)]
        for i, e in enumerate(self.elements):
            for j, elem_list in enumerate(self.p_e):            
                if e.midline_indx in elem_list:
                    self.Omg[j].append(e.global_id)

    
    def build_midline(self):

        self.midline_triad = [ self.tri ] * len(self.midline)
        
        self.midline_x = [np.array([0.0, 0.0, 0.0])] * len(self.midline)

        self.section_starts = list(self.section_ends)
        self.section_starts.insert(0, 0.0) # At the start
        self.section_starts.pop() # remove the section end

        for i, id in enumerate(self.midline_to_section) : # For each of the midline points

            x0 = self.x_sec_ends[id]
            t = self.triads[id] # This is the triad at the start of the section

            s = self.midline[i] - self.section_starts[id] # Distane into this section (arclength)

            if self.section_list[id]['type'].lower() == 'straight':
                self.midline_x[i] = self.get_coords_straight(s, x0, t[0])
                self.midline_triad[i] = t.copy() # For a straight pipe this doesn't change
            elif self.section_list[id]['type'].lower() == 'straight_new':
                self.midline_x[i] = self.get_coords_straight(s, x0, t[0])
                self.midline_triad[i] = t.copy() # For a straight pipe this doesn't change
            elif self.section_list[id]['type'].lower() == 'bend':

                roc = self.section_list[id]['param']['radius']
                axis, v = self.which_rotation_axis(self.section_list[id]['param']['axis'], id)
                angle = s / roc * np.sign(self.section_list[id]['param']['angle'])

                # Need a point to rotate about
                p = x0 + np.sign(angle) * roc * v
                self.midline_x[i] = rotate_point_about_point(x0, axis, angle, p)
                self.midline_triad[i] = rotate_triad(t, axis, angle)

            elif self.section_list[id]['type'].lower() == 'bend_new':

                roc = self.section_list[id]['param']['radius']
                axis, v = self.which_rotation_axis_new(self.section_list[id]['param']['dir1'],self.section_list[id]['param']['dir2'], id)
                angle = s / roc

                # Need a point to rotate about
                p = x0 + roc * v
                self.midline_x[i] = rotate_point_about_point(x0, axis, angle, p)
                self.midline_triad[i] = rotate_triad(t, axis, angle)

            else:
                raise Exception("Unknown section type " + self.section_list[i]['type'])

    def transform_mesh(self):

        """
            This function transforms the straight mesh

            It takes no inputs or returns nothing. It updates each update each "node" with the index of the section.
        """
        for i, n in enumerate(self.nodes): # For each node.
            idx = n.midline_indx # This is the index of the midline
            dr = self.midline_triad[idx][1]
            dn = self.midline_triad[idx][2]
            n.coords = self.midline_x[idx] + n.v[1] * dr + n.v[2] * dn 
            
    def mideline_to_node(self,
                         node : Node):
        
        return node.coords - np.array([self.mesh.s[node.midline_indx], 0.0, 0.0])

    def get_coords_straight(self, 
                            s : float,
                            x0 : np.array, 
                            ds : np.array):
        return x0 + s * ds
    
    def export(self, 
               filename: str = "foo.vtk",
               point_data : dict = {},
               cell_data : dict = {},
               ):
        
        points = [] 
        for n in self.nodes:
            points.append(n.coords)
            
        connectivity = []
        for e in self.elements:
            if(e.active):
                connectivity.append(e.list_of_nodes)

        # Get the element type right for meshio
        elem_type = self.elem_type
    
        if(elem_type == "tri"):
            elem_type = "triangle"
        nodel = len(self.elements[0].list_of_nodes)

        if(nodel > 4 and elem_type != "hex"):
            elem_type += str(nodel)
        elif (elem_type == "hex"):
            elem_type = "hexahedron"

        cells = [
            (elem_type, connectivity),
        ]

        # Alternative with the same options
        meshio.write_points_cells(filename, points, cells, point_data=point_data, cell_data = cell_data)

