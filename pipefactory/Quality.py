from meshio import Mesh
import numpy as np
import pyvista as pv

from neatmesh._reader import MeshReader3D
from neatmesh._analyzer import Analyzer3D

from prettytable import PrettyTable

class Analyzer2:

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


class Analyzer:

    def __init__(self, points, cells):

        pv_cells = [x for sublist in cells[0][1] for x in ([8]+sublist)]

        grid = pv.UnstructuredGrid(pv_cells,[pv.CellType.HEXAHEDRON]*len(cells[0][1]), points)

        self.analyze(grid)

    def analyze(self, grid):
        self.skew = grid.compute_cell_quality(quality_measure="skew")
        self.dim = grid.compute_cell_quality(quality_measure="dimension")
        self.vol = grid.compute_cell_quality(quality_measure="volume")
        self.diag = grid.compute_cell_quality(quality_measure="diagonal")
        

    def get_stats_from_array(self, array):

        arr_max = np.nanmax(array)
        arr_min = np.nanmin(array)
        arr_mean = np.nanmean(array)
        arr_std = np.nanstd(array)

        return [arr_max, arr_min, arr_mean, arr_std]

    def __call__(self):

        quality_metrics_dict = {
            "volume": self.vol, "skew": self.skew, "dimension": self.dim, "diagonal": self.diag
        }
        
        acceptable_range = {
            "volume": None, "skew": [0,0.5], "dimension": None, "diagonal": [0.65, 1]
        }
        
        t = PrettyTable(["Stat.", "Max.", "Min.", "Mean", "Std.", "Status"])

        for k, v in quality_metrics_dict.items():
            listofstat = self.get_stats_from_array(pv.convert_array(v.active_scalars))
            status = [" "]
            if acceptable_range[k] is not None:
                colour = '\033[92m' if (acceptable_range[k][1]>=listofstat[0]) and (acceptable_range[k][0]<=listofstat[1]) else '\033[91m'
                status = [colour + f"{acceptable_range[k]}" + "\033[0m"]

            t.add_row([k] + listofstat + status)

        print(t)



        


