from meshio import Mesh
import numpy as np

from neatmesh._reader import MeshReader3D
from neatmesh._analyzer import Analyzer3D

from prettytable import PrettyTable

class Analyzer:

    def __init__(self, meshio_mesh : Mesh):

        self.mesh = MeshReader3D(meshio_mesh)
        self.analyzer = Analyzer3D(self.mesh)
        self.analyze()

    def analyze(self):
        print("Collecting cell types..")
        self.analyzer.count_cell_types()

        print("Analyzing faces...")
        self.analyzer.analyze_faces()

        print("Analyzing cells...")
        self.analyzer.analyze_cells()

        print("Checking non-orthogonality...")
        self.analyzer.analyze_non_ortho()

        print("Checking adjacent cells volume ratio...")
        self.analyzer.analyze_adjacents_volume_ratio()

    def get_stats_from_array(self, array):

        arr_max = np.nanmax(array)
        arr_min = np.nanmin(array)
        arr_mean = np.nanmean(array)
        arr_std = np.nanstd(array)

        return [arr_max, arr_min, arr_mean, arr_std]

    def __call__(self):

        quality_metrics_dict = {
            "Face Area": {"array": self.analyzer.face_areas},
            "Face Aspect Ratio": {"array": self.analyzer.face_aspect_ratios},
            "Cell Volume": {"array": self.analyzer.cells_volumes},
            "Non-Orthogonality": {"array": self.analyzer.non_ortho},
            "Adjacent Cells Volume Ratio": {"array": self.analyzer.adj_ratio},
        }
        
        t = PrettyTable(["Stat.", "Max.", "Min.", "Mean", "Std."])

        for k, v in quality_metrics_dict.items():
            t.add_row([k] + self.get_stats_from_array(v["array"]))

        print(t)




        


