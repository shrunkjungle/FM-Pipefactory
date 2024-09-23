# from pipefactory import Pipe
from .Pipe import Pipe
from json import load
import numpy as np

def PartitionROM(input_json : str, limb_len : float = 1.6):
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

    mesh_info = PipeParam(outer_radius = radius, 
                        thickness = thickness, 
                        element_size = el_len,
                        element_around_circum = el_circum, 
                        elements_through_thickness = el_thruthick,
                        initial_direction = current_dir)

    mesh_info.add_straight(2*limb_len)

    mesh = Pipe(outer_radius = mesh_info.outer_radius, 
                thickness = mesh_info.thickness, 
                section_list=mesh_info.section_list, 
                elem_type=("hex", False), 
                element_size = mesh_info.element_size,
                element_around_circum = mesh_info.element_around_circum, 
                elements_through_thickness = mesh_info.elements_through_thickness)

    mesh.export(f'{input_json}_straight.xdmf')
    mesh_info.save_to_json(f'{input_json}_straight', mesh.midline.tolist())

    for dict in data["Mesh Sections"]:
        if dict["type"] == "Bend":
            mesh_info = PipeParam(outer_radius = radius, 
                                thickness = thickness, 
                                element_size = el_len,
                                element_around_circum = el_circum, 
                                elements_through_thickness = el_thruthick,
                                initial_direction = current_dir)

            mesh_info.add_straight(limb_len)
            mesh_info.add_bend(dict["radius"], dict["direction_end"])
            mesh_info.add_straight(limb_len)

            mesh = Pipe(outer_radius = mesh_info.outer_radius, 
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
                 outer_radius : float, 
                 thickness : float,
                 element_size : float,
                 element_around_circum : int, 
                 elements_through_thickness : int,
                 initial_direction : list = [1.0,0.0,0.0],
                 origin : list = [0.,0.,0.]):
        
        self.outer_radius = outer_radius
        self.thickness=thickness
        self.element_size=element_size
        self.element_around_circum=element_around_circum
        self.elements_through_thickness=elements_through_thickness

        self.section_list = [] # for pipefactory - still with bend_new and straight_new
        self.mesh_sections = [] # identical list but with renamed type
        self.current_dir = initial_direction
        self.origin = origin

    def add_straight(self, length : float):

        s_dict = {'length':length,
                  'dir': np.array(self.current_dir),
                  'type': 'Straight_new'}
        
        ss_dict = {'length':length,
                  'direction': self.current_dir,
                  'type': 'Straight'}

        self.section_list.append(s_dict)
        self.mesh_sections.append(ss_dict)

    def add_bend(self, radius : float, direction : list | np.ndarray):

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

    def save_to_json(self, name : str, midline : np.ndarray | None):

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

    def load_json(self, name : str):

        from json import load
        with open(f"{name}.json", 'r') as file:
            data = load(file)

        params = data["Pipe Parameters"]

        self.outer_radius = params['Outer Radius']
        self.thickness=params['Thickness']
        self.element_size=params['Element Size']
        self.element_around_circum=params['Elements around Circumference']
        self.elements_through_thickness=params['Elements through Thickness']

        self.section_list = []
        for dic in data["Mesh Sections"]:
            if dic["type"] == "Straight":
                self.section_list.append({'length':dic['length'],
                                          'dir': np.array(dic['direction']),
                                          'type': 'Straight_new'})
            else:
                self.section_list.append({'type': 'Bend_new',
                                          'param': {'dir1' : np.array(dic['direction_begin']),
                                                    'dir2' : np.array(dic['direction_end']),
                                                    'radius': dic['radius'],
                                                    'angle' : dic['length']/dic['radius']}
                                                    })
        self.mesh_sections = data
