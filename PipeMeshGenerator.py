import numpy as np
from pipefactory import PipeParam, Pipe, RadialCrack, PartitionROM

########### Input Parameters #############

name = f"foo"

mesh_info = PipeParam(outer_radius = 0.0365, 
                      thickness = 0.00305, 
                      element_size = 0.01,
                      element_around_circum = 48, 
                      elements_through_thickness = 3)

mesh_info.add_straight(1.6)

mesh = Pipe(outer_radius = mesh_info.outer_radius, 
            thickness = mesh_info.thickness, 
            section_list=mesh_info.section_list, 
            elem_type=("hex", False), 
            element_size = mesh_info.element_size,
            element_around_circum = mesh_info.element_around_circum, 
            elements_through_thickness = mesh_info.elements_through_thickness,
            #mesh_refinement=pf.AxialRefinement(1.50,0.001, pf.Ramp(0.1,0.3)),
            #thermal_expansion_opt=therm_opt,
            )

# mesh.degenerate_crack(RadialCrack(0.4,0.,np.pi,0., 0.01,0.))
mesh.export(f'{name}.xdmf')
mesh_info.save_to_json(f'{name}', mesh.midline.tolist())

