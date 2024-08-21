import pipefactory as pf
from json import load
import numpy as np

def PartitionROM(input_json):
    """
    Partition mesh into separate bends in order to do Reduced Order FEM.
    
    Parameters
    ----------
    input_json : str
        Name of the pipe to partition
    """
    with open(f"{input_json}.json", 'r') as file:
        data = load(file)

    param_dict = data["Pipe Parameters"]

    radius = param_dict["Outer Radius"]
    thickness = param_dict["Thickness"]
    el_len =  param_dict["Element Size"]
    el_thruthick = param_dict["Elements through Thickness"]
    el_circum = param_dict["Elements around Circumference"]

    bend_count = 0
    current_dir = [1.0,0.0,0.0]

    for dict in data["Mesh Sections"]:
        if dict["type"] == "Bend":
            mesh_info = PipeParam(outer_radius = radius, 
                                thickness = thickness, 
                                element_size = el_len,
                                element_around_circum = el_circum, 
                                elements_through_thickness = el_thruthick,
                                initial_direction = current_dir)

            mesh_info.add_straight(1.6)
            mesh_info.add_bend(dict["radius"], dict["direction_end"])
            mesh_info.add_straight(1.6)

            mesh = pf.Pipe(outer_radius = mesh_info.outer_radius, 
                        thickness = mesh_info.thickness, 
                        section_list=mesh_info.section_list, 
                        elem_type=("hex", False), 
                        element_size = mesh_info.element_size,
                        element_around_circum = mesh_info.element_around_circum, 
                        elements_through_thickness = mesh_info.elements_through_thickness)

            mesh.export(f'{input_json}_bend{bend_count}.xdmf')
            mesh_info.save_to_json(f'{input_json}_bend{bend_count}', mesh.midline.tolist())

            current_dir = dict["direction_end"]
            bend_count += 1

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

        direction = np.array(direction)/np.linalg.norm(np.array(direction))

        angle = np.arccos(np.dot(np.array(self.current_dir),direction))

        b_dict = {'type': 'Bend_new',
                  'param': {'dir1' : np.array(self.current_dir),
                            'dir2' : np.array(direction),
                            'radius': radius,
                            'angle' : angle}
                            }
        
        bb_dict = {'type': 'Bend',
                   'direction_begin' : self.current_dir,
                   'direction_end' : direction.tolist(),
                   'radius': radius,
                   'length' : radius*angle}

        self.section_list.append(b_dict)
        self.mesh_sections.append(bb_dict)

        self.current_dir = direction.tolist()

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