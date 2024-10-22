import numpy as np
from pipefactory import PipeParam, Pipe, AxialRefinement, Ramp, Hole, Dimple, Weld, RadialCrack, AxialCrack, Cuboid, PartitionROM
import json
import sys
import os
import shutil
from scipy.stats import qmc

def generate_parameter_combinations(parameters, sample_num, ele_around_circum):
        """
        Generate all combinations of parameters based on lower and upper
        bounds for each, and the total number of samples
        """
        defects = parameters.pop('defects')
        defect_types = []
        for defect in defects.items():
            defect_type = defect[0].split('_')[0]
            defect_types.append(defect_type)
            for param, values in defect[1].items():
                if f'{defect_type}_{param}' not in parameters:
                    if 'phi' in param:
                        parameters[f'{defect_type}_{param}'] = list(map(np.deg2rad, list(map(float, values))))
                        if parameters[f'{defect_type}_{param}'][0] > parameters[f'{defect_type}_{param}'][1]:
                            parameters[f'{defect_type}_{param}'][1] += 2*np.pi
                    elif 's0' not in param:
                        parameters[f'{defect_type}_{param}'] = [value / 1000 for value in list(map(float, values))]
                    else:
                        parameters[f'{defect_type}_{param}'] = values
                else:
                    added = False
                    i = 1
                    while added == False:
                        if f'{defect_type}_{param}_{i}' not in parameters:
                            if 'phi' in param:
                                parameters[f'{defect_type}_{param}_{i}'] = list(map(np.deg2rad, list(map(float, values))))
                                if parameters[f'{defect_type}_{param}'][0] > parameters[f'{defect_type}_{param}'][1]:
                                    parameters[f'{defect_type}_{param}'][1] += 2*np.pi
                            elif 's0' not in param:
                                parameters[f'{defect_type}_{param}_{i}'] = [value / 1000 for value in list(map(float, values))]
                            else:
                                parameters[f'{defect_type}_{param}_{i}'] = values
                            added = True
                        else:
                            continue
        dimensionality = len(parameters.items())
        constant_params = {}
        for param in parameters.items():
            #values = param[1].split(',')
            values = param[1]
            if values[0] == values[1]:
                constant_params[param[0]] = values[0]
        for consts in constant_params:
            parameters.pop(consts)
            dimensionality -= 1
        delta = 360 / float(ele_around_circum)
        if len(parameters) != 0:
            thresholds_low = [float(value[0]) if 'RadialCrack_phi' not in param and 'AttachedCuboid_phi' not in param else 0 if 'phi_span' not in param else -0.5 for param, value in parameters.items()]
            thresholds_high = [float(value[1]) if 'RadialCrack_phi' not in param and 'AttachedCuboid_phi' not in param else int((value[1] - value[0]) / np.deg2rad(delta/2)) if 'phi_span' not in param else (int((value[1] - value[0]) / np.deg2rad(delta)) + 0.5) for param, value in parameters.items()]
            sampler = qmc.LatinHypercube(d=dimensionality)
            params_sample = sampler.random(sample_num)
            params_sample_scaled = qmc.scale(params_sample, thresholds_low, thresholds_high)

            combinations_list = [list(zip(list(parameters.keys()), sample)) for sample in params_sample_scaled]
            combinations_dict_list = [dict(combinations_list[i]) for i in range(len(combinations_list))]

            for const, const_value in constant_params.items():
                for sample in combinations_dict_list:
                    sample[const] = float(const_value)

            i = 0
            for param, value in parameters.items():
                if ('RadialCrack_phi' in param or 'AttachedCuboid_phi' in param) and 'phi_span' not in param:
                    for sample in combinations_dict_list:
                        sample[param] = value[0] + (round(sample[param]) * np.deg2rad(delta/2))
                if 'phi' in param and 'phi_span' not in param:
                    for sample in combinations_dict_list:
                        sample[param] = np.mod(sample[param], 2*np.pi)
                i += 1

            for sample in combinations_dict_list:
                sample["Material Properties"] = {"density": sample.pop("density"), "ym": sample.pop("ym"), "poisson": sample.pop("poisson")}
                sample["Defects"] = {}
                hole_params = ["s0", "phi0", "radius"]
                dimple_params = ["s0", "phi0", "depth", "radius"]
                weld_params = ["s0", "Aout", "Ain", "ell_out", "ell_in"]
                radialcrack_params = ["s0", "phi0", "phi_span", "width", "depth", "smoothing_dist"]
                axialcrack_params = ["s0", "phi", "length", "width", "depth", "smoothing_dist"]
                attachedcuboid_params = ["s0", "phi0", "phi_span", "length", "height"]
                for defect in defect_types:
                    if f'{defect}' not in sample["Defects"]:
                        sample["Defects"][f'{defect}'] = {}
                        if defect == "Hole":
                            for param in hole_params:
                                sample["Defects"][f'{defect}'][param] = sample.pop(f"{defect}_{param}")
                        elif defect == "Dimple":
                            for param in dimple_params:
                                sample["Defects"][f'{defect}'][param] = sample.pop(f"{defect}_{param}")
                        elif defect == "Weld":
                            for param in weld_params:
                                sample["Defects"][f'{defect}'][param] = sample.pop(f"{defect}_{param}")
                        elif defect == "RadialCrack":
                            for param in radialcrack_params:
                                sample["Defects"][f'{defect}'][param] = sample.pop(f"{defect}_{param}")
                            sample["Defects"][f'{defect}']["phi"] = sample["Defects"][f'{defect}'].pop("phi0")

                            dphi_sample = sample["Defects"][f'{defect}'].pop("phi_span")
                            if f"{defect}_phi_span" in constant_params:
                                if np.mod(dphi_sample/np.deg2rad(delta), 2) == 0:
                                    sample["Defects"][f'{defect}']['phi'] = round(sample["Defects"][f'{defect}']['phi']/np.deg2rad(delta)) * np.deg2rad(delta)
                                elif np.mod(dphi_sample/np.deg2rad(delta), 2) == 1:
                                    sample["Defects"][f'{defect}']['phi'] = min(range(1, 2*ele_around_circum, 2), key=lambda x:abs(x-sample["Defects"][f'{defect}']['phi']/np.deg2rad(delta/2))) * np.deg2rad(delta/2)
                                    if sample["Defects"][f'{defect}']['phi'] > parameters[f'{defect}_phi0'][1]:
                                        sample["Defects"][f'{defect}']['phi'] -= 2*np.deg2rad(delta)
                                    elif sample["Defects"][f'{defect}']['phi'] < parameters[f'{defect}_phi0'][0]:
                                        sample["Defects"][f'{defect}']['phi'] += 2*np.deg2rad(delta)
                                else:
                                    raise ValueError("Some rounding has gone wrong")
                            else:
                                dphi_low = parameters[f"{defect}_phi_span"][0]
                                odds = range(1, int((parameters[f"{defect}_phi_span"][1] - parameters[f"{defect}_phi_span"][0]) / np.deg2rad(delta)), 2)
                                evens = range(0, int((parameters[f"{defect}_phi_span"][1] - parameters[f"{defect}_phi_span"][0]) / np.deg2rad(delta)), 2)
                                if np.absolute(np.mod(sample["Defects"][f'{defect}']["phi"], np.deg2rad(delta)) - np.deg2rad(delta/2)) > np.deg2rad(delta/4):
                                    if np.mod(np.ceil(dphi_sample), 2) == 0:
                                        if dphi_sample > int((parameters[f"{defect}_phi_span"][1] - parameters[f"{defect}_phi_span"][0]) / np.deg2rad(delta)):
                                            dphi_sample = np.random.choice(evens)
                                        else:
                                            dphi_sample = np.ceil(dphi_sample)
                                    else:
                                        dphi_sample = np.floor(dphi_sample)
                                else:
                                    if np.mod(np.ceil(dphi_sample), 2) != 0:
                                        if dphi_sample > int((parameters[f"{defect}_phi_span"][1] - parameters[f"{defect}_phi_span"][0]) / np.deg2rad(delta)):
                                            dphi_sample = np.random.choice(odds)
                                        else:
                                            dphi_sample = np.ceil(dphi_sample)
                                    else:
                                        dphi_sample = np.floor(dphi_sample)
                                dphi_sample = np.mod(dphi_low + (dphi_sample * np.deg2rad(delta)), 2*np.pi)
                            sample["Defects"][f'{defect}']["dphi"] = dphi_sample
                        elif defect == "AxialCrack":
                            for param in axialcrack_params:
                                sample["Defects"][f'{defect}'][param] = sample.pop(f"{defect}_{param}")
                        elif defect == "AttachedCuboid":
                            for param in attachedcuboid_params:
                                sample["Defects"][f'{defect}'][param] = sample.pop(f"{defect}_{param}")
                            sample["Defects"][f'{defect}']["phi"] = sample["Defects"][f'{defect}'].pop("phi0")
                            
                            dphi_sample = sample["Defects"][f'{defect}'].pop("phi_span")
                            if f"{defect}_phi_span" in constant_params:
                                if np.mod(dphi_sample/np.deg2rad(delta), 2) == 0:
                                    sample["Defects"][f'{defect}']['phi'] = round(sample["Defects"][f'{defect}']['phi']/np.deg2rad(delta)) * np.deg2rad(delta)
                                elif np.mod(dphi_sample/np.deg2rad(delta), 2) == 1:
                                    sample["Defects"][f'{defect}']['phi'] = min(range(1, 2*ele_around_circum, 2), key=lambda x:abs(x-sample["Defects"][f'{defect}']['phi']/np.deg2rad(delta/2))) * np.deg2rad(delta/2)
                                    if sample["Defects"][f'{defect}']['phi'] > parameters[f'{defect}_phi0'][1]:
                                        sample["Defects"][f'{defect}']['phi'] -= 2*np.deg2rad(delta)
                                    elif sample["Defects"][f'{defect}']['phi'] < parameters[f'{defect}_phi0'][0]:
                                        sample["Defects"][f'{defect}']['phi'] += 2*np.deg2rad(delta)
                                else:
                                    raise ValueError("Some rounding has gone wrong")
                            else:
                                dphi_low = parameters[f"{defect}_phi_span"][0]
                                odds = range(1, int((parameters[f"{defect}_phi_span"][1] - parameters[f"{defect}_phi_span"][0]) / np.deg2rad(delta)), 2)
                                evens = range(0, int((parameters[f"{defect}_phi_span"][1] - parameters[f"{defect}_phi_span"][0]) / np.deg2rad(delta)), 2)
                                if np.absolute(np.mod(sample["Defects"][f'{defect}']["phi"], np.deg2rad(delta)) - np.deg2rad(delta/2)) > np.deg2rad(delta/4):
                                    if np.mod(np.ceil(dphi_sample), 2) == 0:
                                        if dphi_sample > int((parameters[f"{defect}_phi_span"][1] - parameters[f"{defect}_phi_span"][0]) / np.deg2rad(delta)):
                                            dphi_sample = np.random.choice(evens)
                                        else:
                                            dphi_sample = np.ceil(dphi_sample)
                                    else:
                                        dphi_sample = np.floor(dphi_sample)
                                else:
                                    if np.mod(np.ceil(dphi_sample), 2) != 0:
                                        if dphi_sample > int((parameters[f"{defect}_phi_span"][1] - parameters[f"{defect}_phi_span"][0]) / np.deg2rad(delta)):
                                            dphi_sample = np.random.choice(odds)
                                        else:
                                            dphi_sample = np.ceil(dphi_sample)
                                    else:
                                        dphi_sample = np.floor(dphi_sample)
                                dphi_sample = np.mod(dphi_low + (dphi_sample * np.deg2rad(delta)), 2*np.pi)
                            sample["Defects"][f'{defect}']["dphi"] = dphi_sample
                    else:
                        added = False
                        i = 1
                        while added == False:
                            if f'{defect}_{i}' not in sample["Defects"]:
                                sample["Defects"][f'{defect}_{i}'] = {}
                                if defect == "Hole":
                                    for param in hole_params:
                                        sample["Defects"][f'{defect}_{i}'][param] = sample.pop(f"{defect}_{param}_{i}")
                                elif defect == "Dimple":
                                    for param in dimple_params:
                                        sample["Defects"][f'{defect}_{i}'][param] = sample.pop(f"{defect}_{param}_{i}")
                                elif defect == "Weld":
                                    for param in dimple_params:
                                        sample["Defects"][f'{defect}_{i}'][param] = sample.pop(f"{defect}_{param}_{i}")
                                elif defect == "RadialCrack":
                                    for param in radialcrack_params:
                                        sample["Defects"][f'{defect}_{i}'][param] = sample.pop(f"{defect}_{param}_{i}")
                                    sample["Defects"][f'{defect}_{i}']["phi"] = sample["Defects"][f'{defect}_{i}'].pop("phi0")
                                    
                                    dphi_sample = sample["Defects"][f'{defect}_{i}'].pop("phi_span")
                                    if f'{defect}_phi_span_{i}' in constant_params:
                                        if np.mod(dphi_sample/np.deg2rad(delta), 2) == 0:
                                            sample["Defects"][f'{defect}_{i}']['phi'] = round(sample["Defects"][f'{defect}_{i}']['phi']/np.deg2rad(delta)) * np.deg2rad(delta)
                                        elif np.mod(dphi_sample/np.deg2rad(delta), 2) == 1:
                                            sample["Defects"][f'{defect}_{i}']['phi'] = min(range(1, 2*ele_around_circum, 2), key=lambda x:abs(x-sample["Defects"][f'{defect}_{i}']['phi']/np.deg2rad(delta/2))) * np.deg2rad(delta/2)
                                            if sample["Defects"][f'{defect}_{i}']['phi'] > parameters[f'{defect}_phi0_{i}'][1]:
                                                sample["Defects"][f'{defect}_{i}']['phi'] -= 2*np.deg2rad(delta)
                                            elif sample["Defects"][f'{defect}_{i}']['phi'] < parameters[f'{defect}_phi0_{i}'][0]:
                                                sample["Defects"][f'{defect}_{i}']['phi'] += 2*np.deg2rad(delta)
                                        else:
                                            raise ValueError("Some rounding has gone wrong")
                                    else:
                                        dphi_low = parameters[f'{defect}_phi_span_{i}'][0]
                                        odds = range(1, int((parameters[f'{defect}_phi_span_{i}'][1] - parameters[f'{defect}_phi_span_{i}'][0]) / np.deg2rad(delta)), 2)
                                        evens = range(0, int((parameters[f'{defect}_phi_span_{i}'][1] - parameters[f'{defect}_phi_span_{i}'][0]) / np.deg2rad(delta)), 2)
                                        if np.absolute(np.mod(sample["Defects"][f'{defect}_{i}']["phi"], np.deg2rad(delta)) - np.deg2rad(delta/2)) > np.deg2rad(delta/4):
                                            if np.mod(np.ceil(dphi_sample), 2) == 0:
                                                if dphi_sample > int((parameters[f'{defect}_phi_span_{i}'][1] - parameters[f'{defect}_phi_span_{i}'][0]) / np.deg2rad(delta)):
                                                    dphi_sample = np.random.choice(evens)
                                                else:
                                                    dphi_sample = np.ceil(dphi_sample)
                                            else:
                                                dphi_sample = np.floor(dphi_sample)
                                        else:
                                            if np.mod(np.ceil(dphi_sample), 2) != 0:
                                                if dphi_sample > int((parameters[f'{defect}_phi_span_{i}'][1] - parameters[f'{defect}_phi_span_{i}'][0]) / np.deg2rad(delta)):
                                                    dphi_sample = np.random.choice(odds)
                                                else:
                                                    dphi_sample = np.ceil(dphi_sample)
                                            else:
                                                dphi_sample = np.floor(dphi_sample)
                                        dphi_sample = np.mod(dphi_low + (dphi_sample * np.deg2rad(delta)), 2*np.pi)
                                    sample["Defects"][f'{defect}_{i}']["dphi"] = dphi_sample
                                elif defect == "AxialCrack":
                                    for param in axialcrack_params:
                                        sample["Defects"][f'{defect}_{i}'][param] = sample.pop(f"{defect}_{param}_{i}")
                                elif defect == "AttachedCuboid":
                                    for param in attachedcuboid_params:
                                        sample["Defects"][f'{defect}_{i}'][param] = sample.pop(f"{defect}_{param}_{i}")
                                    sample["Defects"][f'{defect}_{i}']["phi"] = sample["Defects"][f'{defect}_{i}'].pop("phi0")
                                    
                                    dphi_sample = sample["Defects"][f'{defect}_{i}'].pop("phi_span")
                                    if f'{defect}_phi_span_{i}' in constant_params:
                                        if np.mod(dphi_sample/np.deg2rad(delta), 2) == 0:
                                            sample["Defects"][f'{defect}_{i}']['phi'] = round(sample["Defects"][f'{defect}_{i}']['phi']/np.deg2rad(delta)) * np.deg2rad(delta)
                                        elif np.mod(dphi_sample/np.deg2rad(delta), 2) == 1:
                                            sample["Defects"][f'{defect}_{i}']['phi'] = min(range(1, 2*ele_around_circum, 2), key=lambda x:abs(x-sample["Defects"][f'{defect}_{i}']['phi']/np.deg2rad(delta/2))) * np.deg2rad(delta/2)
                                            if sample["Defects"][f'{defect}_{i}']['phi'] > parameters[f'{defect}_phi0_{i}'][1]:
                                                sample["Defects"][f'{defect}_{i}']['phi'] -= 2*np.deg2rad(delta)
                                            elif sample["Defects"][f'{defect}_{i}']['phi'] < parameters[f'{defect}_phi0_{i}'][0]:
                                                sample["Defects"][f'{defect}_{i}']['phi'] += 2*np.deg2rad(delta)
                                        else:
                                            raise ValueError("Some rounding has gone wrong")
                                    else:
                                        dphi_low = parameters[f'{defect}_phi_span_{i}'][0]
                                        odds = range(1, int((parameters[f'{defect}_phi_span_{i}'][1] - parameters[f'{defect}_phi_span_{i}'][0]) / np.deg2rad(delta)), 2)
                                        evens = range(0, int((parameters[f'{defect}_phi_span_{i}'][1] - parameters[f'{defect}_phi_span_{i}'][0]) / np.deg2rad(delta)), 2)
                                        if np.absolute(np.mod(sample["Defects"][f'{defect}_{i}']["phi"], np.deg2rad(delta)) - np.deg2rad(delta/2)) > np.deg2rad(delta/4):
                                            if np.mod(np.ceil(dphi_sample), 2) == 0:
                                                if dphi_sample > int((parameters[f'{defect}_phi_span_{i}'][1] - parameters[f'{defect}_phi_span_{i}'][0]) / np.deg2rad(delta)):
                                                    dphi_sample = np.random.choice(evens)
                                                else:
                                                    dphi_sample = np.ceil(dphi_sample)
                                            else:
                                                dphi_sample = np.floor(dphi_sample)
                                        else:
                                            if np.mod(np.ceil(dphi_sample), 2) != 0:
                                                if dphi_sample > int((parameters[f'{defect}_phi_span_{i}'][1] - parameters[f'{defect}_phi_span_{i}'][0]) / np.deg2rad(delta)):
                                                    dphi_sample = np.random.choice(odds)
                                                else:
                                                    dphi_sample = np.ceil(dphi_sample)
                                            else:
                                                dphi_sample = np.floor(dphi_sample)
                                        dphi_sample = np.mod(dphi_low + (dphi_sample * np.deg2rad(delta)), 2*np.pi)
                                    sample["Defects"][f'{defect}_{i}']["dphi"] = dphi_sample
                                added = True
                            else:
                                i += 1
                                continue
        else:
            combinations_list = [[tuple()]] * sample_num
            combinations_dict_list = [{}] * sample_num
    
        return combinations_dict_list

name = sys.argv[1]
mesh_complexity = sys.argv[2]

if mesh_complexity == 'Basic':
    mesh_info = PipeParam(#outer_radius = 0.0365,
                          outer_radius = float(sys.argv[4]) / 1000,
                          #thickness = 0.00305,
                          thickness = float(sys.argv[5]) / 1000,
                          element_size = 0.1,
                          element_around_circum = 32,
                          elements_through_thickness = 1,
                          #initial_direction=[np.cos(np.deg2rad(8)), 0., -np.sin(np.deg2rad(8))],
                          initial_direction = [0.0, 0.0, 1.0],
                          origin=[0.,0.,0.])
elif mesh_complexity == 'Detailed':
    mesh_info = PipeParam(#outer_radius = 0.0365,
                          outer_radius = float(sys.argv[4]) / 1000,
                          #thickness = 0.00305,
                          thickness = float(sys.argv[5]) / 1000,
                          element_size = 0.01,
                          element_around_circum = 48, 
                          elements_through_thickness = 3,
                          #initial_direction=[np.cos(np.deg2rad(8)), 0., -np.sin(np.deg2rad(8))],
                          initial_direction = [0.0, 0.0, 1.0],
                          origin=[0.,0.,0.])
elif mesh_complexity == 'Sampling':
    params = json.loads(sys.argv[3])
    ele_size = float(params.pop("ele_size"))
    ele_thru_thick = int(params.pop("ele_thru_thick"))
    ele_around_circum = int(params.pop("ele_around_circum"))
    radius = float(params.pop("radius"))
    wall_thickness = float(params.pop("thickness"))
    sample_num = int(params.pop("sample_num"))
    sections = params.pop("pipe_geom")
    mesh_info = PipeParam(outer_radius = radius + (wall_thickness / 2),
                          thickness = wall_thickness,
                          element_size = ele_size,
                          element_around_circum = ele_around_circum,
                          elements_through_thickness = ele_thru_thick,
                          initial_direction = [0.0, 0.0, 1.0],
                          origin = [0.,0.,0.])
    for section in sections.values():
        if 'Straight' in section:
            mesh_info.add_straight(float(section['Straight'][0]))
        elif 'Bend' in section:
            mesh_info.add_bend(float(section['Bend'][0]), list(map(float, section['Bend'][2])))
        else:
            raise ValueError("Section must be of type 'Straight' or 'Bend'")
    combinations = generate_parameter_combinations(params, sample_num, ele_around_circum)
elif mesh_complexity != 'Existing':
    raise ValueError("Mesh complexity must be either 'Detailed', 'Basic' or 'Existing'")

    
if mesh_complexity != 'Existing' and mesh_complexity != 'Sampling':
    sections = json.loads(sys.argv[3])
    for section in sections.values():
        if 'Straight' in section:
            mesh_info.add_straight(float(section['Straight'][0]))
        elif 'Bend' in section:
            mesh_info.add_bend(float(section['Bend'][0]), list(map(float, section['Bend'][2])))
        else:
            raise ValueError("Section must be of type 'Straight' or 'Bend'")

    mesh = Pipe(outer_radius = mesh_info.outer_radius, 
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

    if mesh_complexity == 'Basic':
        if os.path.exists(f'{name}.json'):
            os.remove(f'{name}.json')
        if os.path.exists(f'{name}.xdmf'):
            os.remove(f'{name}.xdmf')
        if os.path.exists(f'{name}.h5'):
            os.remove(f'{name}.h5')
        mesh.export(f'{name}.xdmf')
        mesh_info.save_to_json(f'{name}')

    elif mesh_complexity == 'Detailed':
        if not os.path.exists('FM_FEM/fmfem/ReducedOrder/input_xdmf'):
            os.mkdir('FM_FEM/fmfem/ReducedOrder/input_xdmf')
        if not os.path.exists('FM_FEM/fmfem/ReducedOrder/input_json'):
            os.mkdir('FM_FEM/fmfem/ReducedOrder/input_json')
        mesh.export(f'FM_FEM/fmfem/ReducedOrder/input_xdmf/{name}.xdmf')
        mesh_info.save_to_json(f'FM_FEM/fmfem/ReducedOrder/input_json/{name}')
        PartitionROM(f"{name}")

elif mesh_complexity == 'Sampling':
    os.mkdir(f'{name}_meshsamples')
    i = 1
    for sample in combinations:
        # Sequentially add the defects to the pipe mesh
        for defect, params in sample["Defects"].items():
            # Sequentially add mesh refinement CURRENTLY WORKS ONLY WITH AT MOST ONE NON-RADIAL-CRACK DEFECT (not needed for radial cracks)
            if "RadialCrack" not in defect:
                mesh_refinement = AxialRefinement(params["s0"], ele_size/5.0, Ramp(0.1,0.3))
                break
            else:
                mesh_refinement = None
            
        mesh = Pipe(outer_radius = mesh_info.outer_radius, 
                thickness = mesh_info.thickness, 
                section_list=mesh_info.section_list, 
                elem_type=("hex", False), 
                element_size = mesh_info.element_size,
                element_around_circum = mesh_info.element_around_circum, 
                elements_through_thickness = mesh_info.elements_through_thickness,
                origin = mesh_info.origin,
                mesh_refinement = mesh_refinement
                #mesh_refinement=pf.AxialRefinement(1.50,0.001, pf.Ramp(0.1,0.3)),
                #thermal_expansion_opt=therm_opt,
                )

        for defect, params in sample["Defects"].items():
            # Add defect
            if "Hole" in defect:
                mesh.remove_elements(Hole(float(params["s0"]), float(params["phi0"]), float(params["radius"])))
            if "Dimple" in defect:
                mesh.add_defect_displacement(Dimple(float(params["s0"]), float(params["phi0"]), float(params["depth"]), float(params["radius"])))
            if "Weld" in defect:
                mesh.remove_elements(Weld(float(params["s0"]), float(params["Aout"]), float(params["Ain"]), float(params["ell_out"]), float(params["ell_in"])))
            if "RadialCrack" in defect:
                mesh.degenerate_crack(RadialCrack(float(params["s0"]), float(params["phi"]), float(params["dphi"]), float(params["width"]), float(params["depth"]),float(params["smoothing_dist"])))
            if "AxialCrack" in defect:
                mesh.degenerate_crack2(AxialCrack(float(params["s0"]), float(params["phi"]), float(params["length"]), float(params["width"]), float(params["depth"]), float(params["smoothing_dist"])))
            if "AttachedCuboid" in defect:
                mesh.add_elements(Cuboid(float(params["s0"]), float(params["phi"]), float(params["dphi"]), float(params["length"]), float(params["height"])))

        mesh.export(f'{name}_meshsamples/{name}_{i:0{7}d}.xdmf')
        mesh_info.save_to_json(f'{name}_meshsamples/{name}_{i:0{7}d}', materials=sample["Material Properties"], defects=list(sample["Defects"].items()))
        print(i)
        i += 1
else:
    if not os.path.exists('FM_FEM/fmfem/ReducedOrder/input_xdmf'):
        os.mkdir('FM_FEM/fmfem/ReducedOrder/input_xdmf')
    if not os.path.exists('FM_FEM/fmfem/ReducedOrder/input_json'):
        os.mkdir('FM_FEM/fmfem/ReducedOrder/input_json')
    shutil.move(f'tar_extract/{name}.json', 'FM_FEM/fmfem/ReducedOrder/input_json/')
    shutil.move(f'tar_extract/{name}.xdmf', 'FM_FEM/fmfem/ReducedOrder/input_xdmf/')
    shutil.move(f'tar_extract/{name}.h5', 'FM_FEM/fmfem/ReducedOrder/input_xdmf/')
    PartitionROM(f"{name}")
# PartitionROM("ITER_M6")

