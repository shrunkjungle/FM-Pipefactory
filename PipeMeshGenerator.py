import numpy as np
import pipefactory as pf

########## Pipe Parameter Class ###############

class PipeParam:
    def __init__(self,
                 outer_radius, 
                 thickness,
                 element_size,
                 element_around_circum, 
                 elements_through_thickness,
                 initial_direction = [1.0,0.0,0.0],
                 origin = [0.,0.,0.]):
        
        self.outer_radius = outer_radius
        self.thickness=thickness
        self.element_size=element_size
        self.element_around_circum=element_around_circum
        self.elements_through_thickness=elements_through_thickness

        self.section_list = [] # for pipefactory - still with bend_new and straight_new
        self.mesh_sections = [] # identical list but with renamed type
        self.current_dir = initial_direction
        self.origin = origin

    def add_straight(self, length):

        s_dict = {'length':length,
                  'dir': np.array(self.current_dir),
                  'type': 'Straight_new'}
        
        ss_dict = {'length':length,
                  'direction': self.current_dir,
                  'type': 'Straight'}

        self.section_list.append(s_dict)
        self.mesh_sections.append(ss_dict)

    def add_bend(self, radius, direction):

        b_dict = {'type': 'Bend_new',
                  'param': {'dir1' : np.array(self.current_dir),
                            'dir2' : np.array(direction),
                            'radius': radius}
                            }
        
        bb_dict = {'type': 'Bend',
                   'direction_begin' : self.current_dir,
                   'direction_end' : direction,
                   'radius': radius}

        self.section_list.append(b_dict)
        self.mesh_sections.append(bb_dict)

        self.current_dir = direction

    def save_to_json(self, name, midline):

        self.pipe_parameters = {
            'Outer Radius': self.outer_radius,
            'Thickness': self.thickness,
            'Element Size': self.element_size,
            'Elements through Thickness': self.elements_through_thickness,
            'Elements around Circumference': self.element_around_circum,
            'Pipe Mesh Origin': self.origin
        }

        data_to_save = {
            'Pipe Parameters': self.pipe_parameters,
            'Mesh Sections': self.mesh_sections,
            'Defects': None,
            'Midline': midline
        }

        from json import dump

        with open(f'{name}.json', 'w') as file:
            dump(data_to_save, file, indent=4)

########### Input Parameters #############

name = "foo"

mesh_info = PipeParam(outer_radius = 0.0365, 
                      thickness = 0.01, 
                      element_size = 0.01,
                      element_around_circum = 48, 
                      elements_through_thickness = 2)

mesh_info.add_straight(2.0)
# mesh_info.add_bend(0.3,[1.0,1.0,0.0])



mesh = pf.Pipe(outer_radius = mesh_info.outer_radius, 
               thickness = mesh_info.thickness, 
               section_list=mesh_info.section_list, 
               elem_type=("hex", False), 
               element_size = mesh_info.element_size,
               element_around_circum = mesh_info.element_around_circum, 
               elements_through_thickness = mesh_info.elements_through_thickness)
            #    mesh_refinement=pf.AxialRefinement(0.5,0.0025, pf.Ramp(0.1,0.3)))

mesh.degenerate_crack(pf.RadialCrack(s0=0.5,phi0=67.6,phi1=112.4,crack_width=0.0,crack_depth=0.01,outer_radius=0.0365,thickness = 0.01))

mesh.export(f'{name}.xdmf')
mesh_info.save_to_json(f'{name}', mesh.midline.tolist())