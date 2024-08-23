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


name = f"al056004"

mesh_info = PipeParam(outer_radius = 0.01220, 
                        thickness = 0.00643, 
                        element_size = 0.005,
                        element_around_circum = 60, 
                        elements_through_thickness = 4)

mesh_info.add_straight(5.004)
#mesh_info.add_bend(0.2,[0.,1.,0.])

# sample_positions = [0.02, 0.50, 1.00, 1.50, 2.20, 2.90, 3.60, 4.30, 4.90]
# temps = [100.747, 66.518, 53.437, 43.54, 34.922, 29.784, 26.264, 24.532, 23.477]
# temp_func = lambda x : np.interp(x, sample_positions, temps)
# therm_opt = {"distribution": temp_func, "lin_exp_coeff": 24 * 1e-5, "ambient_temp": 20}

mesh = pf.Pipe(outer_radius = mesh_info.outer_radius, 
            thickness = mesh_info.thickness, 
            section_list=mesh_info.section_list, 
            elem_type=("hex", False), 
            element_size = mesh_info.element_size,
            element_around_circum = mesh_info.element_around_circum, 
            elements_through_thickness = mesh_info.elements_through_thickness,
            mesh_refinement=pf.AxialRefinement(1.50,0.001, pf.Ramp(0.1,0.3)),
            #thermal_expansion_opt=therm_opt
            )

#mesh.degenerate_crack(pf.RadialCrack(s0=3.5,phi0=0,phi1=361,crack_width=0.002,crack_depth=0.004, smoothing_dist=0.0,outer_radius=mesh_info.outer_radius,thickness = mesh_info.thickness, el_thru_thick=mesh_info.elements_through_thickness))

def slit_profile(ds_slit, z):

    thickness = 0.0063
    slit_depth = 0.004
    depth_ratio = slit_depth/thickness
    #z ranges from -1 to 1 through the thickness
    if z < 1 - 2*depth_ratio:
        return False
    else:
        return True

mesh.remove_elements( pf.Radial_Slit(s0=1.50, phi0=0, phi1=360, slit_width=0.002,outer_radius=mesh_info.outer_radius,thickness = mesh_info.thickness, partial = True, profile=slit_profile))

mesh.export(f'FM-FEM/input_xdmf/{name}.xdmf')
#mesh.export_walls(f'{name}_walls.xdmf')
mesh_info.save_to_json(f'FM-FEM/input_json/{name}', mesh.midline.tolist())