import numpy as np
import pipefactory as pf

########### Input Parameters #############

name = "foo"

mesh_info = pf.PipeParam(outer_radius = 0.0366, 
                      thickness = 0.00305, 
                      element_size = 0.01,
                      element_around_circum = 48, 
                      elements_through_thickness = 3)

mesh_info.add_straight(1.6)
mesh_info.add_bend(1.0, [1.,0.,-1.])
mesh_info.add_straight(1.6)


mesh = pf.Pipe(outer_radius = mesh_info.outer_radius, 
               thickness = mesh_info.thickness, 
               section_list=mesh_info.section_list, 
               elem_type=("hex", False), 
               element_size = mesh_info.element_size,
               element_around_circum = mesh_info.element_around_circum, 
               elements_through_thickness = mesh_info.elements_through_thickness)

mesh.export(f'{name}.xdmf')
# mesh_info.save_to_json(f'{name}', mesh.midline.tolist())
