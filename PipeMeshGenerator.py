import numpy as np
import pipefactory as pf

########### Input Parameters #############

name = f"bar"

mesh_info = pf.PipeParam(outer_radius = 0.0365, 
                      thickness = 0.01, 
                      element_size = 0.01,
                      element_around_circum = 32, 
                      elements_through_thickness = 3,
                      )

mesh_info.add_straight(0.1)
mesh_info.add_bend(1.0, [1.,1.,1.])
mesh_info.add_straight(0.1)
mesh_info.add_bend(0.2, [-1.0,0.6,0.8])

# name = f"foo"

# mesh_info = pf.PipeParam(outer_radius = 0.0365, 
#                       thickness = 0.00305, 
#                       element_size = 0.01,
#                       element_around_circum = 48, 
#                       elements_through_thickness = 3,
#                       initial_direction=[np.cos(np.deg2rad(8)), 0., -np.sin(np.deg2rad(8))],
#                       origin=[0.,0.,0.])

# mesh_info.add_straight(1.015)
# mesh_info.add_bend(1.0, [np.cos(np.deg2rad(14)), 0., -np.sin(np.deg2rad(14))])
# mesh_info.add_straight(0.095)
# mesh_info.add_bend(1.0, [1.,0.,0.])
# mesh_info.add_straight(0.802)
# mesh_info.add_bend(1.0, [np.cos(np.deg2rad(18)), 0., -np.sin(np.deg2rad(18))])
# mesh_info.add_straight(0.110)
# mesh_info.add_bend(1.0, [1.,0.,0.])
# mesh_info.add_straight(1.192)

# mesh_info.add_bend(1.0, [np.cos(np.deg2rad(45)), -np.sin(np.deg2rad(45)),0.])
# mesh_info.add_straight(0.539)
# mesh_info.add_bend(1.0, [1.,0.,0.])
# mesh_info.add_straight(2.023)
# mesh_info.add_bend(1.0, [np.cos(np.deg2rad(23)), -np.sin(np.deg2rad(23)),0.])
# mesh_info.add_straight(0.258)
# mesh_info.add_bend(1.0, [1.,0.,0.])
# mesh_info.add_straight(1.239)

mesh = pf.Pipe(outer_radius = mesh_info.outer_radius, 
            thickness = mesh_info.thickness, 
            section_list=mesh_info.section_list, 
            elem_type=("hex", False), 
            element_size = mesh_info.element_size,
            element_around_circum = mesh_info.element_around_circum, 
            elements_through_thickness = mesh_info.elements_through_thickness,
            origin= mesh_info.origin,
            #mesh_refinement=pf.AxialRefinement(1.50,0.001, pf.Ramp(0.1,0.3)),
            #thermal_expansion_opt=therm_opt,
            )

mesh.degenerate_crack2(pf.AxialCrack(0.1, 0.0,0.05, 0.003,0.005, 0.03))
mesh.degenerate_crack(pf.RadialCrack(0.2,6.0,1.0, 0.004, 0.01, 0.03))
mesh.add_defect_displacement(pf.Dimple(0.3,0.0,0.02, 0.04))
mesh.add_elements(pf.Cuboid(0.6, 6.0, 1.0, 0.05, 0.02))

mesh.quality_check()

mesh.export(f'{name}.xdmf', save_point_data=True)
# mesh_info.save_to_json(f'{name}', mesh.midline.tolist())

# PartitionROM("ITER_M6")

